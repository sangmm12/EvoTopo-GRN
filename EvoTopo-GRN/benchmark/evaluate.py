"""Turn method output into comparable metric rows.

Three thresholding rules are reported side by side.

  native   the method's own binary output. Defined for EvoTopo-GRN (adaptive
           lasso path truncated at the one-standard-error rule) and SINDy
           (non-zero coefficients).
  cut0.1   a fixed importance threshold of 0.1, the convention of the reference
           comparison study. Applied only to GENIE3, dynGENIE3 and the nonlinear
           ODE model, whose scores are normalised variable importances.
  topk     matched sparsity: the k highest-scoring off-diagonal entries, where k
           is the number of edges EvoTopo-GRN selected on the same data set. All
           methods are then compared at equal graph density.

AUROC, AUPRC and EP@k are threshold-free and repeated on every row.
"""
import numpy as np
import scores as SC

UPSTREAM_CUT = ("GENIE3", "dynGENIE3", "nonlinear_ODE")   # prefix match


def confusion(E, A):
    p = A.shape[0]
    o = ~np.eye(p, dtype=bool)
    e = np.asarray(E, int)[o]
    t = np.asarray(A, int)[o]
    TP = int(((e == 1) & (t == 1)).sum())
    FP = int(((e == 1) & (t == 0)).sum())
    TN = int(((e == 0) & (t == 0)).sum())
    FN = int(((e == 0) & (t == 1)).sum())
    return TP, FP, TN, FN


def metrics(E, A, S=None):
    """Directed metrics against A; skeleton MCC against S when supplied.
    Precision is undefined when no edge is declared and is recorded as NaN
    rather than zero."""
    TP, FP, TN, FN = confusion(E, A)
    den = np.sqrt(float(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    TPR = TP / (TP + FN) if TP + FN else np.nan
    PPV = TP / (TP + FP) if TP + FP else np.nan
    F1 = (2 * PPV * TPR / (PPV + TPR)) if (PPV and TPR and PPV + TPR > 0) else np.nan
    r = dict(TP=TP, FP=FP, TN=TN, FN=FN,
             TPR=TPR, FPR=FP / (FP + TN) if FP + TN else np.nan,
             Precision=PPV, FDR=(1 - PPV) if PPV == PPV else np.nan, F1=F1,
             Specificity=TN / (TN + FP) if TN + FP else np.nan,
             Accuracy=(TP + TN) / (TP + FP + TN + FN),
             MCC=(TP * TN - FP * FN) / den if den else 0.0,
             edges=int(np.asarray(E, int).sum() - np.trace(np.asarray(E, int))))
    if S is not None:
        Ei = np.asarray(E, int)
        r["skel_MCC"] = metrics(((Ei + Ei.T) > 0).astype(int), S)["MCC"]
    return r


def topk(V, k):
    """Keep the k highest off-diagonal scores."""
    V = np.asarray(V, float).copy()
    np.fill_diagonal(V, -np.inf)
    p = V.shape[0]
    E = np.zeros((p, p), int)
    if k > 0:
        for t in np.argsort(-V, axis=None)[:int(k)]:
            E[t // p, t % p] = 1
    return E


def rows_for(fitted, A, S, base, matched_k=None):
    """fitted: output of comparison_methods.fit_all. base: fields written on
    every row. Returns one row per (method, thresholding rule)."""
    if matched_k is None:
        e = fitted.get("EvoTopo-GRN", {}).get("E")
        matched_k = int(e.sum() - np.trace(e)) if e is not None else int(A.sum())
    out = []
    for name, r in fitted.items():
        V, E, nd = r["V"], r["E"], r.get("dropped", 0)
        tf = SC.threshold_free(V, A) if V is not None else dict(
            AUROC=np.nan, AUPRC=np.nan, EPk=np.nan, AUPRC_base=np.nan)
        cand = {}
        if E is not None:
            cand["native"] = np.asarray(E, int)
        if V is not None and name.startswith(UPSTREAM_CUT):
            cand["cut0.1"] = np.where(np.asarray(V, float) >= 0.1, 1, 0)
        if V is not None:
            cand["topk"] = topk(V, matched_k)
        for rule, Ei in cand.items():
            out.append(dict(**base, method=name, rule=rule, matched_k=matched_k,
                            time_dedup=nd, **metrics(Ei, A, S), **tf))
    return out
