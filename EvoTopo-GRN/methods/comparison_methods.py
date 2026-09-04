"""Unified interface to EvoTopo-GRN and seven comparison methods.

    EvoTopo-GRN     power-law detrending, screening, adaptive lasso
    GENIE3          steady-state random forests, x_j ~ x_{-j}
    dynGENIE3       dx/dt + alpha*x on the remaining genes, RF or ET
    SINDy           dx/dt on the library {sin(x/2), x**3}; STLSQ, SR3 or Lasso
    nonlinear_ODE   dx/dt + alpha*x, gradient-boosted trees
    TIGRESS         randomised LARS with stability selection
    GRNBoost2       per-target gradient boosting, gain importance
    GRNVBEM         lag-one autoregression with an ARD sparsity prior

Provenance
----------
GENIE3, dynGENIE3, SINDy and the nonlinear ODE model are executed from the
reference comparison script (path in MTL_SRC), unchanged except for the three
adjustments recorded below. TIGRESS follows Haury et al. (2012).
GRNBoost2 and GRNVBEM are EQUIVALENT REIMPLEMENTATIONS, not the original
packages: arboreto cannot be installed in this environment and the original
GRNVBEM is MATLAB source. Any result attributed to them carries that caveat.

Each method returns (V, E):
    V   p x p continuous score matrix, larger = more edge-like, or None
    E   p x p binary adjacency when the method has its own rule, or None
Rows are targets, columns are regulators; the diagonal is not evaluated.
"""
import os
import sys
import numpy as np
import warnings
warnings.filterwarnings("ignore")
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import lars_path, ARDRegression

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _d in ("methods", "simulation", "benchmark"):
    sys.path.insert(0, os.path.join(_ROOT, _d))

from evotopo_grn import power_equation_fit, detrend, select_edges   # noqa: E402
import scores as SC                                                 # noqa: E402

# ------------------------------------------------------------------ settings
NTREES = int(os.environ.get("NTREES", "500"))     # shared by GENIE3 and dynGENIE3
CUT = float(os.environ.get("CUT", "0.1"))
TIG_R = int(os.environ.get("TIGRESS_R", "300"))
TIG_L = int(os.environ.get("TIGRESS_L", "5"))
TIG_A = float(os.environ.get("TIGRESS_ALPHA", "0.2"))
DYN_VAR = os.environ.get("DYNGENIE3_VARIANTS", "RF").split(",")
SIN_VAR = os.environ.get("SINDY_VARIANTS", "STLSQ").split(",")

MTL_SRC = os.environ.get("MTL_SRC", "")

_UPSTREAM = {}


def _load_upstream():
    """Load the reference script lazily. Two edits are applied:
    (i) the module-level call to its own simulator is removed so that importing
        it does not run a simulation;
    (ii) dynGENIE3 exposes tree_method and its tree count is set to NTREES; the
        routine defaults to 1000, which would give the time-series variant twice
        the budget of steady-state GENIE3."""
    if _UPSTREAM:
        return _UPSTREAM
    if not MTL_SRC or not os.path.exists(MTL_SRC):
        raise FileNotFoundError(
            "Reference comparison script not found. Set MTL_SRC to its path. "
            "Without it, dynGENIE3, SINDy and nonlinear_ODE cannot be run; the "
            "remaining methods work unchanged.")
    raw = open(MTL_SRC).read().replace("X,X_obs,t,true_adj = sim20(snr)", "")
    ns = {"__name__": "upstream"}
    exec(compile(raw, MTL_SRC, "exec"), ns)
    i0 = raw.index("def run_dynGenie3(X,t):")
    i1 = raw.index("def run_nonliner_ODE(X,t):")
    body = (raw[i0:i1]
            .replace("def run_dynGenie3(X,t):", "def run_dynGenie3_tm(X,t,tree_method='RF'):")
            .replace("dynGENIE3(TS_data, time_points)",
                     f"dynGENIE3(TS_data, time_points, tree_method=tree_method, ntrees={NTREES})"))
    ns2 = dict(ns)
    exec(compile(body, "<dyngenie3-patched>", "exec"), ns2)
    _UPSTREAM.update(run_pySINDY=ns["run_pySINDY"],
                     run_nonliner_ODE=ns["run_nonliner_ODE"],
                     run_dynGenie3=ns2["run_dynGenie3_tm"])
    return _UPSTREAM


# ------------------------------------------------------- numerical safeguards
def _dedupe_time(X, t, rel=0.01):
    """Merge samples whose index values are closer than `rel` times the median
    interval. Finite differences otherwise divide by a near-zero interval."""
    t = np.asarray(t, float)
    order = np.argsort(t)
    t, X = t[order], np.asarray(X)[order]
    thr = rel * float(np.median(np.diff(t)))
    keep = [0]
    for i in range(1, len(t)):
        if t[i] - t[keep[-1]] >= thr:
            keep.append(i)
    keep = np.array(keep)
    return X[keep], t[keep], len(t) - len(keep)


def _robust(fn, X, t, *a):
    """Retry on de-duplicated samples if a time-series method fails. The catch
    must be broad: scikit-learn raises ValueError, XGBoost raises its own
    XGBoostError, which is not a ValueError subclass."""
    try:
        return fn(X, t, *a), 0
    except Exception:
        Xd, td, nd = _dedupe_time(X, t)
        return fn(Xd, td, *a), nd


def _dejitter_constant(D, rng=None):
    """Break ties in constant columns. A constant gene puts its minimum and
    maximum at the same sample, so the degradation-rate estimator divides by
    zero."""
    D = np.asarray(D, float).copy()
    sd = D.std(0, ddof=1)
    bad = np.where(sd < 1e-12)[0]
    if len(bad):
        rng = np.random.default_rng(0) if rng is None else rng
        scale = max(float(np.median(sd[sd > 0])) if np.any(sd > 0) else 1.0, 1e-9)
        D[:, bad] += rng.normal(0, 1e-6 * scale, (D.shape[0], len(bad)))
    return D


def _positivize(Y, mode):
    """dynGENIE3 takes logarithms of the per-gene extrema, so its input must be
    positive. On the raw scale log10(count+1) is zero exactly at zero counts, so
    a floor at 1e-3 binds only there. On detrended residuals, which take both
    signs, the series is shifted instead: clipping would collapse every negative
    residual onto one value."""
    Y = np.asarray(Y, float)
    if mode == "raw":
        return np.maximum(Y, 1e-3)
    m = float(Y.min())
    return Y + (1e-3 - m) if m < 1e-3 else Y


# ---------------------------------------------------------------- the methods
def m_evotopo(D, seed=0):
    p = D.shape[1]
    E = np.zeros((p, p), int)
    for c in range(p):
        try:
            E[c] = select_edges(D, c, seed=seed)
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            pass
    V = np.vstack([SC.evotopo_scores(D, c) for c in range(p)])
    return V, E


def m_genie3(D, seed=0):
    p = D.shape[1]
    V = np.zeros((p, p))
    for j in range(p):
        idx = [k for k in range(p) if k != j]
        m = RandomForestRegressor(n_estimators=NTREES, max_features="sqrt",
                                  random_state=seed, n_jobs=-1).fit(D[:, idx], D[:, j])
        fi = m.feature_importances_
        s = fi.sum()
        V[j, idx] = fi / s if s > 0 else fi
    return V, None


def m_dyngenie3(D, EI, mode, variant="RF"):
    up = _load_upstream()
    V, nd = _robust(up["run_dynGenie3"], _positivize(_dejitter_constant(D), mode), EI, variant)
    return np.asarray(V, float), None, nd


def m_sindy(D, EI, which=("STLSQ",)):
    """The reference routine returns (STLSQ, SR3, Lasso) binary matrices at once."""
    up = _load_upstream()
    (s0, s1, s2), nd = _robust(up["run_pySINDY"], D, EI)
    got = {"STLSQ": s0, "SR3": s1, "Lasso": s2}
    return {k: np.asarray(got[k], int) for k in which}, nd


def m_nonlinear_ode(D, EI):
    up = _load_upstream()
    V, nd = _robust(up["run_nonliner_ODE"], _dejitter_constant(D), EI)
    return np.asarray(V, float), None, nd


def m_tigress(D, seed=0, R=None, L=None, alpha=None):
    """Haury et al. (2012): R randomised LARS runs per target, each on a random
    half of the samples with predictor weights drawn from U(alpha, 1), the first
    L LARS steps, area scoring (a predictor in the active set at step t
    contributes 1/L per step)."""
    R = TIG_R if R is None else R
    L = TIG_L if L is None else L
    alpha = TIG_A if alpha is None else alpha
    rng = np.random.default_rng(seed)
    n, p = D.shape
    half = max(n // 2, 5)
    Z = (D - D.mean(0)) / np.maximum(D.std(0, ddof=1), 1e-12)
    V = np.zeros((p, p))
    for j in range(p):
        idx = [k for k in range(p) if k != j]
        X0, y0 = Z[:, idx], Z[:, j]
        f = np.zeros(len(idx))
        for _ in range(R):
            r = rng.choice(n, half, replace=False)
            w = rng.uniform(alpha, 1.0, len(idx))
            try:
                _, _, c = lars_path(X0[r] * w[None, :], y0[r], method="lar", max_iter=L)
            except Exception:
                continue
            st = min(L, c.shape[1] - 1)
            for t in range(1, st + 1):
                f[np.abs(c[:, t]) > 0] += 1.0 / st
        V[j, idx] = f / R
    return V, None


def m_grnboost2(D, seed=0):
    """Equivalent reimplementation (not arboreto): per-target gradient boosting
    with gain importance, in the configuration GRNBoost2 specifies."""
    from xgboost import XGBRegressor
    p = D.shape[1]
    V = np.zeros((p, p))
    for j in range(p):
        idx = [k for k in range(p) if k != j]
        m = XGBRegressor(n_estimators=500, learning_rate=0.01, max_depth=3,
                         subsample=0.9, colsample_bytree=0.1,
                         importance_type="gain", n_jobs=-1,
                         random_state=seed, verbosity=0).fit(D[:, idx], D[:, j])
        fi = m.feature_importances_
        s = fi.sum()
        V[j, idx] = fi / s if s > 0 else fi
    return V, None


def m_grnvbem(D):
    """Equivalent reimplementation (not the MATLAB original): x_j(t+1) ~ sum_k
    A_jk x_k(t) fitted by automatic relevance determination; the score is the
    coefficient magnitude. Samples must already be in index order."""
    p = D.shape[1]
    Xl, Xn = D[:-1], D[1:]
    V = np.zeros((p, p))
    for j in range(p):
        idx = [k for k in range(p) if k != j]
        m = ARDRegression(max_iter=300).fit(Xl[:, idx], Xn[:, j])
        V[j, idx] = np.abs(m.coef_)
    return V, None


# -------------------------------------------------------------- orchestration
def fit_all(data, input_mode="detrend", seed=0, methods=None):
    """Run the requested methods on one data set.

    input_mode  'raw'      every method receives log10(count + 1)
                'detrend'  every method receives the power-law residuals

    Returns {method_name: dict(V=, E=, dropped=)}, variant suffixes included.
    """
    EI = data["EI"]
    if input_mode == "raw":
        D = np.asarray(data["Y"], float)
    elif input_mode == "detrend":
        D = detrend(power_equation_fit(data["expr"], EI))
    else:
        raise ValueError(f"unknown input_mode {input_mode!r}")
    order = np.argsort(EI, kind="stable")
    D, EIs = D[order], np.asarray(EI, float)[order]

    want = set(methods or ["EvoTopo-GRN", "GENIE3", "dynGENIE3", "SINDy",
                           "nonlinear_ODE", "TIGRESS", "GRNBoost2", "GRNVBEM"])
    out = {}
    if "EvoTopo-GRN" in want:
        V, E = m_evotopo(D, seed=seed)
        out["EvoTopo-GRN"] = dict(V=V, E=E, dropped=0)
    if "GENIE3" in want:
        V, E = m_genie3(D, seed=seed)
        out["GENIE3"] = dict(V=V, E=E, dropped=0)
    if "dynGENIE3" in want:
        for v in DYN_VAR:
            V, E, nd = m_dyngenie3(D, EIs, input_mode, v)
            out[f"dynGENIE3_{v}"] = dict(V=V, E=None, dropped=nd)
    if "SINDy" in want:
        got, nd = m_sindy(D, EIs, tuple(SIN_VAR))
        for k, E in got.items():
            out[f"SINDy_{k}"] = dict(V=E.astype(float), E=E, dropped=nd)
    if "nonlinear_ODE" in want:
        V, E, nd = m_nonlinear_ode(D, EIs)
        out["nonlinear_ODE"] = dict(V=V, E=None, dropped=nd)
    if "TIGRESS" in want:
        V, E = m_tigress(D, seed=seed)
        out["TIGRESS"] = dict(V=V, E=None, dropped=0)
    if "GRNBoost2" in want:
        V, E = m_grnboost2(D, seed=seed)
        out["GRNBoost2"] = dict(V=V, E=None, dropped=0)
    if "GRNVBEM" in want:
        V, E = m_grnvbem(D)
        out["GRNVBEM"] = dict(V=V, E=None, dropped=0)
    return out
