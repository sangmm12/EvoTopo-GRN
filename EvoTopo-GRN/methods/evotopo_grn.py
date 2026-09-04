"""EvoTopo-GRN: power-law detrending along the expression index, followed by
marginal screening and adaptive lasso selection.

Pipeline
--------
1. power_equation_fit  fit log10(count + 1) = a_j * EI ** b_j per gene
2. detrend             subtract the fitted trend, giving a residual matrix
3. select_edges        for each target: screen to floor((m) / log m) candidates
                       by absolute marginal correlation, obtain ridge weights by
                       10-fold CV at the one-standard-error rule, then fit an
                       adaptive lasso with penalty factors 1 / |beta_ridge| and
                       take its support at the one-standard-error rule.

Environment variables
---------------------
LASSO_MAXITER    coordinate-descent iteration cap (default 20000)
LASSO_SELECTION  'cyclic' (default) or 'random'; use 'random' on near-collinear
                 designs, where cyclic descent can stall.
"""
import os
import numpy as np
from scipy.optimize import curve_fit
from sklearn.linear_model import Lasso, Ridge
from sklearn.model_selection import KFold

_MAXIT = int(os.environ.get("LASSO_MAXITER", "20000"))
_SEL = os.environ.get("LASSO_SELECTION", "cyclic")


# --------------------------------------------------------------- 1. detrending
def _fit_trend(X, y):
    """Fit y = a * X ** b with three fallbacks: nonlinear least squares, the
    closed-form log-log fit, and a constant. Returns (fitted values, (a, b))."""
    pos = y > 0
    if pos.sum() < 5:
        return np.full_like(X, float(y.mean())), (np.nan, np.nan)
    b0, la0 = np.polyfit(np.log(X[pos]), np.log(y[pos]), 1)
    cands = []
    try:
        (a, b), _ = curve_fit(lambda x, a, b: a * x ** b, X, y,
                              p0=[np.exp(la0), b0], maxfev=200000)
        cands.append((a, b))
    except Exception:
        pass
    cands.append((np.exp(la0), b0))
    for a, b in cands:
        if not np.isfinite([a, b]).all() or abs(b) > 50:
            continue
        tr = a * X ** b
        if np.isfinite(tr).all():
            return tr, (float(a), float(b))
    return np.full_like(X, float(y.mean())), (np.nan, np.nan)


def power_equation_fit(counts, EI=None):
    """counts: genes x samples. EI must be computed from the FULL count matrix.
    Passing a focal-gene subset without EI uses the wrong sample order and
    destroys every downstream step that depends on the ordering."""
    counts = np.asarray(counts, float)
    if EI is None:
        EI = np.log10(counts.sum(0) + 1.0)
    EI = np.asarray(EI, float)
    order = np.argsort(EI, kind="stable")
    dat, X = counts[:, order], EI[order]
    td = np.log10(dat + 1.0)
    trend = np.empty_like(td)
    par = np.empty((dat.shape[0], 2))
    n_fallback = 0
    for g in range(dat.shape[0]):
        tr, ab = _fit_trend(X, td[g])
        if np.allclose(tr, tr[0]):
            n_fallback += 1
        trend[g], par[g] = tr, ab
    return dict(counts=dat, Time=X, order=order, trans_data=td,
                trend=trend, power_par=par, n_fallback=n_fallback)


def detrend(fit):
    """Return the residual matrix as samples x genes, sorted by index."""
    return np.nan_to_num((fit["trans_data"] - fit["trend"]).T,
                         nan=0.0, posinf=0.0, neginf=0.0)


# --------------------------------------------------------------- 2. selection
def _lasso(a):
    return Lasso(alpha=a, max_iter=_MAXIT, selection=_SEL, random_state=0)


def _cv_coef(X, y, alphas, ridge=False, rule="1se", nfolds=10, seed=0):
    """Cross-validated coefficients at the one-standard-error rule (or minimum)."""
    kf = KFold(n_splits=nfolds, shuffle=True, random_state=seed)
    mse = np.zeros((len(alphas), nfolds))
    for f, (tr, te) in enumerate(kf.split(X)):
        for i, a in enumerate(alphas):
            m = Ridge(alpha=a) if ridge else _lasso(a)
            m.fit(X[tr], y[tr])
            mse[i, f] = np.mean((y[te] - m.predict(X[te])) ** 2)
    cvm, cvsd = mse.mean(1), mse.std(1, ddof=1) / np.sqrt(nfolds)
    i_min = int(np.argmin(cvm))
    if rule == "1se":
        cand = np.where(cvm <= cvm[i_min] + cvsd[i_min])[0]
        i_sel = int(cand[np.argmax(alphas[cand])]) if len(cand) else i_min
    else:
        i_sel = i_min
    m = Ridge(alpha=alphas[i_sel]) if ridge else _lasso(alphas[i_sel])
    m.fit(X, y)
    return m.coef_


def select_edges(D, col, reduction=True, rule="1se", seed=0):
    """Regulators of gene `col` from the residual matrix D (samples x genes).
    Returns a binary row of length p."""
    p = D.shape[1]
    y = D[:, col]
    others = np.array([k for k in range(p) if k != col])
    out = np.zeros(p, int)
    if y.std() < 1e-12:
        return out
    x = D[:, others]
    if reduction:
        vec = np.nan_to_num(np.abs([np.corrcoef(x[:, c], y)[0, 1]
                                    for c in range(x.shape[1])]))
        nk = max(int(x.shape[1] / np.log(x.shape[1])), 1)
        sel = np.argsort(-vec)[:nk]
        x, others = x[:, sel], others[sel]
    xs = (x - x.mean(0)) / np.maximum(x.std(0, ddof=1), 1e-12)
    ys = y - y.mean()
    rc = _cv_coef(xs, ys, np.logspace(3, -2, 40), ridge=True, rule="1se", seed=seed)
    w = np.maximum(np.abs(rc), 1e-8)
    xa = xs * w[None, :]
    amax = np.max(np.abs(xa.T @ ys)) / len(ys)
    bc = _cv_coef(xa, ys, amax * np.logspace(0, -3, 40), rule=rule, seed=seed)
    out[others[np.abs(bc) > 1e-10]] = 1
    return out


def evotopo_grn(counts, EI=None, reduction=True, seed=0):
    """End-to-end: counts (genes x samples) -> binary adjacency [target, regulator]."""
    fit = power_equation_fit(counts, EI)
    D = detrend(fit)
    p = D.shape[1]
    adj = np.zeros((p, p), int)
    for c in range(p):
        try:
            adj[c] = select_edges(D, c, reduction=reduction, seed=seed)
        except (np.linalg.LinAlgError, ValueError, FloatingPointError):
            pass
    return adj, fit
