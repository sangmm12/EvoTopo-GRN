"""Continuous scores for EvoTopo-GRN and threshold-free metrics."""
import numpy as np
import warnings
warnings.filterwarnings("ignore")
from sklearn.linear_model import lasso_path, Ridge
from sklearn.metrics import roc_auc_score, average_precision_score


def evotopo_scores(D, col, reduction=True, n_alphas=100):
    """Score each candidate regulator by the penalty value at which it first
    enters the adaptive lasso path (larger = stronger). Candidates removed by
    screening keep a negligible score derived from their marginal correlation,
    so the ranking stays complete for AUROC/AUPRC."""
    p = D.shape[1]
    y = D[:, col]
    others = np.array([k for k in range(p) if k != col])
    x = D[:, others]
    out = np.zeros(p)
    if y.std() < 1e-12:
        return out
    if reduction:
        vec = np.nan_to_num(np.abs([np.corrcoef(x[:, c], y)[0, 1]
                                    for c in range(x.shape[1])]))
        nk = max(int(x.shape[1] / np.log(x.shape[1])), 1)
        sel = np.argsort(-vec)[:nk]
        out[others] = vec * 1e-6
        x, others = x[:, sel], others[sel]
    xs = (x - x.mean(0)) / np.maximum(x.std(0, ddof=1), 1e-12)
    ys = y - y.mean()
    rc = Ridge(alpha=1.0).fit(xs, ys).coef_
    w = np.maximum(np.abs(rc), 1e-8)
    xa = xs * w[None, :]
    alphas, coefs, _ = lasso_path(xa, ys, n_alphas=n_alphas)
    for j in range(xa.shape[1]):
        nz = np.where(np.abs(coefs[j]) > 1e-12)[0]
        if len(nz):
            out[others[j]] = float(alphas[nz[0]])
    return out


def threshold_free(S, A):
    """AUROC, AUPRC and precision among the top-k scores, k = number of true
    edges. Evaluated on the off-diagonal ordered pairs only."""
    p = A.shape[0]
    o = ~np.eye(p, dtype=bool)
    s = np.asarray(S, float)[o]
    a = np.asarray(A, int)[o]
    k = int(a.sum())
    if k == 0 or k == len(a):
        return dict(AUROC=np.nan, AUPRC=np.nan, EPk=np.nan, AUPRC_base=np.nan)
    s = np.nan_to_num(s, nan=0.0, posinf=0.0, neginf=0.0)
    top = np.argsort(-s)[:k]
    return dict(AUROC=float(roc_auc_score(a, s)),
                AUPRC=float(average_precision_score(a, s)),
                EPk=float(a[top].mean()),
                AUPRC_base=k / len(a))


def shd(E, A):
    """Structural Hamming distance, counting a reversal once rather than twice."""
    E = np.asarray(E, int).copy()
    A = np.asarray(A, int).copy()
    np.fill_diagonal(E, 0)
    np.fill_diagonal(A, 0)
    rev = int(((E == 1) & (A == 0) & (A.T == 1)).sum())
    fp = int(((E == 1) & (A == 0)).sum())
    fn = int(((E == 0) & (A == 1)).sum())
    return dict(SHD=fp + fn - rev, reversed_edges=rev,
                extra_edges=fp - rev, missed_edges=fn - rev)
