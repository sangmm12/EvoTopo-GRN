# EvoTopo-GRN — simulation and method comparison

Code accompanying the Methods section *Computer Simulation and Performance Metrics*.
It generates synthetic expression data indexed by an expression index (EI), reconstructs
the regulatory network with EvoTopo-GRN, and benchmarks it against seven alternatives.

```
simulation/   data generator (equivalent R and Python implementations)
methods/      EvoTopo-GRN implementation
benchmark/    8-method comparison over SNR 1–20
results/      completed runs (15 replicates per cell)
```

## Requirements

```bash
pip install numpy scipy pandas scikit-learn matplotlib pysindy xgboost
```

`simulation/sim_power_EI_all.R` depends only on base R — no third-party packages.
`pysindy` and `xgboost` are needed only to run the competing methods in `benchmark/`;
the simulator and EvoTopo-GRN itself require neither.

Tested with scikit-learn 1.8.0, SciPy 1.17.1, NumPy 2.4.4, PySINDy (commit `c4421fc`)
and XGBoost 3.5.0.

## 1. Generate simulated data

**R** — single file, prints a full self-check and writes a diagnostic figure:

```bash
Rscript simulation/sim_power_EI_all.R
```

**Python**:

```python
from sim_power_EI import sim_power_EI

d = sim_power_EI(p=20, n=100, snr=10.0, seed=101)
d["expr"]      # 20 x 100 counts for the focal genes
d["expr_all"]  # 500 x 100 full matrix — EI is computed from this
d["EI"]        # expression index, strictly increasing
d["adj_main"]  # ground-truth network (undirected, zero diagonal)
d["mu_true"]   # noise-free power-law trajectories
```

Each gene's log-count follows a power law in the expression index,
`log10(count + 1) = a·EI^b`, and the regulatory network is encoded in the
covariance of the residuals rather than in the dynamics. Counts are drawn by
multinomial sampling with library size `round(10^EI)`, so that the expression
index recomputed from the simulated data reproduces the design values exactly.

Key arguments (Table 1 of the Methods):

| Argument | Default | Meaning |
|---|---|---|
| `p` / `n` / `n_bg` | 20 / 100 / 480 | focal genes / samples / background genes |
| `EI_range` | `(6.95, 7.5)` | expression-index range, matched to the observed data |
| `net_sd` | `0.10` | **network signal strength — fixed, independent of SNR** |
| `net_sparsity` | `0.12` | network sparsity; mean degree 2.4 |
| `snr` | `1 … 20` | controls measurement noise only, `σ_g = sqrt(Var(μ_g)/snr)` |
| `net_type` | `"ggm"` | `"ggm"` undirected Gaussian graphical model; `"indeg1"` exactly one regulator per gene |
| `noise_sd` | `None` | legacy: scales the network with SNR (reverses the MCC trend) — for ablation only |

Setting `snr` governs the independent measurement error and leaves `net_sd` untouched.
This separation matters: the structured residual is the only carrier of network
information, so scaling it by SNR would make SNR an inverse measure of network strength.

## 2. Run EvoTopo-GRN

```python
from evotopo_grn import evotopo_grn, calculate_metrics

adj, fit = evotopo_grn(d["expr"].values)        # adj[target, regulator]
print(calculate_metrics(adj, d["adj_main"].values))
# {'TPR': ..., 'FPR': ..., 'MCC': ...}
```

Three steps, matching the Methods:

| Function | Step |
|---|---|
| `power_equation_fit(counts)` | fit `log10(count+1) = a·EI^b` per gene by nonlinear least squares |
| `detrend(fit)` | subtract the fitted power law, returning the residual matrix |
| `select_edges(resid, target)` | marginal-correlation screening, then adaptive lasso (ridge weights, one-standard-error rule) |

Running `python methods/evotopo_grn.py` applies the method to one simulated data set
and prints the metrics.

## 3. Run the full comparison

```bash
export MTL_SRC=/path/to/Compare_with_pySINDY_dynGenie3_nonlinearODE.py
cd benchmark
SNR_LO=1 SNR_HI=20 NREP=15 python bench.py
python mkplot.py
```

`bench.py` reuses the dynGENIE3, SINDy and nonlinear-ODE implementations from the
script named by `MTL_SRC`, exposing dynGENIE3's `tree_method` as a parameter so that
the random-forest and extremely-randomised-tree variants come from the same routine.

Results are appended row by row to `results/bench_results.csv`, so an interrupted run
loses nothing; re-running skips the `(snr, seed)` combinations already present.
`SNR_LO`, `SNR_HI` and `NREP` are read from the environment.

`results/bench_results.csv` already contains SNR 1–20 with 15 replicates each
(300 combinations). To extend to 50 replicates, simply set `NREP=50` and run again.

The eight methods compared:

| Method | Models | Needs a time axis |
|---|---|---|
| `EvoTopo-GRN` | covariation among samples, after detrending | no |
| `GENIE3` | covariation among samples (steady-state formulation) | no |
| `dynGENIE3_RF` / `_ET` | `dx/dt + αx` on the remaining genes; random forests / extra trees | yes |
| `SINDY_STLSQ` / `_SR3` / `_Lasso` | `dx/dt` on the library `{sin(x/2), x³}` | yes |
| `nonlinear_ODE` | `dx/dt + αx`, gradient-boosted trees | yes |

Methods returning continuous importance scores (GENIE3, both dynGENIE3 variants,
nonlinear ODE) are thresholded at 0.1, following the convention of the parent
comparison. All estimates are scored on the `p(p − 1)` off-diagonal positions.

## Notes

**The diagonal is excluded from scoring.** Some methods set the diagonal of the
estimated adjacency to unity by construction, while tree-based importance matrices
have zero diagonal by construction. Including it awards free true positives to the
former and imposes unavoidable false negatives on the latter — up to 0.2 in MCC, with
the sign depending on the method.

**SINDy is degenerate under its default threshold.** With `threshold = 0.1`, STLSQ and
SR3 retain essentially every candidate edge (TPR = FPR = 1, hence MCC = 0). Finite
differences over closely spaced samples amplify measurement error, leaving the
regression coefficients large and dense relative to the threshold. The same behaviour
appears when the identical code is applied to the parent study's own simulation, so it
is not specific to these data; meaningful SINDy results require the threshold to be
recalibrated per data set.

**Derivative-based methods cannot see this network by construction.** Because the
network lives in the residual covariance and not in the dynamics, `dY/dEI` is
determined entirely by the power-law trend. Near-zero MCC for SINDy, dynGENIE3 and the
nonlinear ODE model reflects a mismatch between their model assumptions and the
data-generating mechanism, not a deficiency of the methods.

**Reproducibility.** Every call is determined by a single `seed`. The R and Python
implementations use different random number streams, so a given seed does not
reproduce identical data across the two, although the statistical properties are
equivalent. Both carry built-in verification of the expression index, the power-law
R² of the noise-free trajectories, the zero-count rate, and the ground-truth degree
distribution.

`methods/mtode_py.py` is a Python port of the MTODE procedure of the parent study
(with ADSIHT replaced by group lasso with BIC selection). It is included for the
diagnostic comparison only and is not among the methods reported here.
