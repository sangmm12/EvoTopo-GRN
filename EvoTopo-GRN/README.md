# EvoTopo-GRN

Gene regulatory network inference by power-law detrending along the expression
index, followed by marginal screening and adaptive lasso selection. Includes a
simulator and a benchmark against seven published methods.

## Install

```bash
pip install -r requirements.txt
```

`dynGENIE3`, `SINDy` and `nonlinear_ODE` are executed from the reference
comparison script. Point `MTL_SRC` at it:

```bash
export MTL_SRC=/path/to/Compare_with_pySINDY_dynGenie3_nonlinearODE.py
```

The other five methods run without it.

## Infer a network

```python
import sys; sys.path.insert(0, "methods")
from evotopo_grn import evotopo_grn
import numpy as np

counts = ...                       # genes x samples, raw counts
EI = np.log10(counts.sum(0) + 1)   # from the FULL count matrix
adj, fit = evotopo_grn(counts, EI) # adj[target, regulator], binary
```

Pass `EI` explicitly whenever `counts` is a subset of the genes; the default
recomputes it from the subset and gives the wrong sample order.

## Simulate

```python
import sys; sys.path.insert(0, "simulation")
from simulate import simulate

d = simulate(p=20, n=100, snr=10, seed=1)
d["expr"]          # genes x samples counts
d["EI"]            # expression index
d["adj_main"]      # directed ground truth
d["adj_skeleton"]  # symmetrised ground truth
```

## Run the benchmark

```bash
cd benchmark

NS=100,200,500,1000 SNRS=1-20 REPS=10 INPUTS=detrend \
  python3 run_benchmark.py

python3 summarize.py benchmark_results.csv --rule topk --input detrend
```

Long grids: run in the background and poll. The driver skips cells already in
the CSV, so it resumes after any interruption.

```bash
nohup python3 run_benchmark.py > run.log 2>&1 &
```

## Environment variables

| Variable | Default | Effect |
|---|---|---|
| `NS` | `100,200,500,1000` | sample sizes; ranges as `100-200` |
| `SNRS` | `1-20` | signal-to-noise ratios |
| `REPS` | `10` | replicates per cell |
| `INPUTS` | `detrend` | `raw`, `detrend`, or both comma-separated |
| `METHODS` | all eight | comma-separated subset |
| `DYNGENIE3_VARIANTS` | `RF` | `RF,ET` |
| `SINDY_VARIANTS` | `STLSQ` | `STLSQ,SR3,Lasso` |
| `NTREES` | `500` | shared by GENIE3 and dynGENIE3 |
| `NGENE` / `NPAR` / `NETSD` | `20` / `1` / `0.10` | genes, max in-degree, network strength |
| `SEED_BASE` | `3000` | seed = SEED_BASE + replicate index |
| `BUDGET` | unlimited | seconds before a clean exit |
| `OUT` | `benchmark_results.csv` | output file |
| `MTL_SRC` | — | path to the reference comparison script |
| `REF_METHOD` | `EvoTopo-GRN` | reference for the paired tests |
| `LASSO_MAXITER` / `LASSO_SELECTION` | `20000` / `cyclic` | set to `3000` / `random` on near-collinear designs |

## Output

`summarize.py` writes `summary_overall.csv`, `summary_by_snr.csv`,
`paired_tests.csv`, `benchmark_summary.xlsx` and `benchmark_figure.png`.

Each result row carries a `rule` column with three thresholding conventions,
reported side by side:

- `native` — the method's own binary output (EvoTopo-GRN, SINDy)
- `cut0.1` — fixed importance threshold (GENIE3, dynGENIE3, nonlinear_ODE)
- `topk` — matched sparsity, k = the number of edges EvoTopo-GRN selected

Threshold-free scores (AUROC, AUPRC, EP@k) appear on every row.

## Files

```
simulation/simulate.py            simulator
methods/evotopo_grn.py            the method
methods/comparison_methods.py     all eight methods behind one interface
benchmark/scores.py               continuous scores, threshold-free metrics
benchmark/evaluate.py             confusion-matrix metrics, thresholding rules
benchmark/run_benchmark.py        grid driver
benchmark/summarize.py            tables, paired tests, figure
```

## Note on two comparison methods

`GRNBoost2` and `GRNVBEM` are equivalent reimplementations, not the original
packages: arboreto could not be installed and the original GRNVBEM is MATLAB
source. Results attributed to them should carry that caveat.
