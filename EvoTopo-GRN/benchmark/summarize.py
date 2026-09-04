"""Summarise a results CSV into tables, paired tests and a figure.

    python3 summarize.py [results.csv] [--rule topk] [--input detrend]

Writes next to the input file:
    summary_overall.csv     mean and sd per (input, rule, n, method)
    summary_by_snr.csv      the same, expanded by SNR
    paired_tests.csv        paired Wilcoxon against REF_METHOD
    benchmark_summary.xlsx
    benchmark_figure.png

Methods within a cell analyse the same data set, so comparisons are paired on
(SNR, replicate); an unpaired test would ignore the shared variation between
data sets.
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

METS = ["MCC", "skel_MCC", "TPR", "FPR", "Precision", "FDR", "F1",
        "Specificity", "Accuracy", "AUROC", "AUPRC", "EPk", "edges"]
ORDER = ["EvoTopo-GRN", "TIGRESS", "GRNBoost2", "GRNVBEM", "GENIE3",
         "dynGENIE3_RF", "dynGENIE3_ET", "SINDy_STLSQ", "SINDy_SR3",
         "SINDy_Lasso", "nonlinear_ODE"]
REF = os.environ.get("REF_METHOD", "EvoTopo-GRN")


def _agg(df, keys):
    mets = [m for m in METS if m in df.columns]
    g = df.groupby(keys, observed=True)[mets].agg(["mean", "std"])
    g.columns = [f"{a}_{'sd' if b == 'std' else 'mean'}" for a, b in g.columns]
    g = g.reset_index()
    g["method"] = pd.Categorical(g["method"], [m for m in ORDER if m in set(df.method)],
                                 ordered=True)
    return g.sort_values(keys).round(4)


def paired(df, metric="MCC"):
    out = []
    for (inp, rule, n), sub in df.groupby(["input", "rule", "n"], observed=True):
        if REF not in set(sub.method):
            continue
        a = sub[sub.method == REF].set_index(["snr", "seed"])[metric]
        for m in sub.method.unique():
            if m == REF:
                continue
            b = sub[sub.method == m].set_index(["snr", "seed"])[metric]
            x, y = a.align(b, join="inner")
            ok = x.notna() & y.notna()
            x, y = x[ok], y[ok]
            if len(x) < 10 or np.allclose(x, y):
                continue
            out.append(dict(input=inp, rule=rule, n=n, metric=metric, ref=REF,
                            method=m, n_pairs=len(x),
                            ref_mean=round(x.mean(), 4), mean=round(y.mean(), 4),
                            diff=round(x.mean() - y.mean(), 4),
                            ref_win_rate=round(float((x > y).mean()), 3),
                            p=stats.wilcoxon(x, y).pvalue))
    return pd.DataFrame(out)


def figure(df, path, rule, inp):
    sub = df[(df.rule == rule) & (df.input == inp)]
    if sub.empty:
        return None
    ns = sorted(sub.n.unique())
    ms = [m for m in ORDER if m in set(sub.method)]
    cmap = plt.get_cmap("tab10")
    col = {m: cmap(i % 10) for i, m in enumerate(ms)}
    col[REF] = "#c0392b"
    fig, ax = plt.subplots(2, len(ns), figsize=(4.3 * len(ns), 7.5),
                           squeeze=False, sharex=True)
    for c, n in enumerate(ns):
        for r, met in enumerate(["MCC", "skel_MCC"]):
            a = ax[r][c]
            for m in ms:
                s = sub[(sub.n == n) & (sub.method == m)].groupby("snr")[met]
                mu, sd = s.mean(), s.std()
                a.plot(mu.index, mu.values, "-o", ms=3, lw=1.6, color=col[m], label=m)
                a.fill_between(mu.index, mu - sd, mu + sd, color=col[m], alpha=.10, lw=0)
            a.grid(alpha=.25)
            a.set_ylim(-0.2, 1.0)
            if r == 0:
                a.set_title(f"n = {n}")
            if c == 0:
                a.set_ylabel("Directed MCC" if met == "MCC" else "Skeleton MCC")
            if r == 1:
                a.set_xlabel("SNR")
    ax[0][0].legend(fontsize=8, loc="upper left", ncol=2, framealpha=.9)
    fig.suptitle(f"input = {inp}   thresholding = {rule}")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    src = args[0] if args else "benchmark_results.csv"
    opt = {}
    for i, a in enumerate(sys.argv[1:], 1):
        if a.startswith("--") and i < len(sys.argv) - 1:
            opt[a.lstrip("-")] = sys.argv[i + 1]
    out_dir = os.path.dirname(os.path.abspath(src))
    df = pd.read_csv(src).drop_duplicates(
        subset=["input", "n", "snr", "seed", "method", "rule"], keep="last")

    ov = _agg(df, ["input", "rule", "n", "method"])
    bs = _agg(df, ["input", "rule", "n", "snr", "method"])
    pt = pd.concat([paired(df, "MCC"), paired(df, "skel_MCC")], ignore_index=True)
    ov.to_csv(os.path.join(out_dir, "summary_overall.csv"), index=False)
    bs.to_csv(os.path.join(out_dir, "summary_by_snr.csv"), index=False)
    pt.to_csv(os.path.join(out_dir, "paired_tests.csv"), index=False)

    rule = opt.get("rule", "topk")
    inp = opt.get("input", df.input.iloc[0])
    figure(df, os.path.join(out_dir, "benchmark_figure.png"), rule, inp)

    with pd.ExcelWriter(os.path.join(out_dir, "benchmark_summary.xlsx"),
                        engine="openpyxl") as w:
        ov.to_excel(w, sheet_name="overall", index=False)
        bs.to_excel(w, sheet_name="by_snr", index=False)
        if not pt.empty:
            pt.to_excel(w, sheet_name="paired_tests", index=False)
        df.to_excel(w, sheet_name="raw", index=False)

    cells = df.groupby(["input", "n", "snr", "seed"], observed=True).ngroups
    print(f"{cells} cells, {len(df)} rows, methods {sorted(df.method.unique())}")
    show = ov[(ov.rule == rule) & (ov.input == inp)]
    print(show[["n", "method", "MCC_mean", "MCC_sd", "skel_MCC_mean", "AUROC_mean"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
