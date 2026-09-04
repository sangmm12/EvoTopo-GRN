"""Grid driver: input protocol x sample size x SNR x replicate, appended to CSV.

Three properties matter for long runs.

  * Each finished cell is appended immediately, so an interrupted run keeps
    everything already computed.
  * On startup the existing CSV is read and completed (input, n, snr, seed)
    cells are skipped, so re-running is idempotent and a run split into
    segments gives the same result as a single run.
  * The column layout must not change once a run has started; adding a column
    mid-run leaves earlier rows short and misaligns the file. Rewrite the whole
    file with the new column back-filled instead.

The seed depends on the replicate index only, so a given replicate is an
independent data set at each sample size rather than a subsample of a common one.
"""
import os
import sys
import time
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "simulation"))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "methods"))

from simulate import simulate                 # noqa: E402
import comparison_methods as M                # noqa: E402
import evaluate as EV                         # noqa: E402


def _ints(s):
    out = []
    for part in str(s).split(","):
        part = part.strip()
        if "-" in part and not part.startswith("-"):
            a, b = part.split("-")
            out += list(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out


NS = _ints(os.environ.get("NS", "100,200,500,1000"))
SNRS = _ints(os.environ.get("SNRS", "1-20"))
REPS = int(os.environ.get("REPS", "10"))
INPUTS = [s.strip() for s in os.environ.get("INPUTS", "detrend").split(",") if s.strip()]
METHODS = [s.strip() for s in os.environ.get(
    "METHODS", "EvoTopo-GRN,GENIE3,dynGENIE3,SINDy,nonlinear_ODE,"
               "TIGRESS,GRNBoost2,GRNVBEM").split(",") if s.strip()]
OUT = os.environ.get("OUT", "benchmark_results.csv")
SEED_BASE = int(os.environ.get("SEED_BASE", "3000"))
BUDGET = float(os.environ.get("BUDGET", "1e9"))
NGENE = int(os.environ.get("NGENE", "20"))
NPAR = int(os.environ.get("NPAR", "1"))
NETSD = float(os.environ.get("NETSD", "0.10"))

KEY = ["input", "n", "snr", "seed"]


def done_keys(path):
    if not os.path.exists(path):
        return set()
    df = pd.read_csv(path)
    if not set(KEY) <= set(df.columns):
        raise SystemExit(f"{path} has an incompatible column layout; use a new OUT file.")
    return set(map(tuple, df[KEY].astype(str).values))


def one_cell(input_mode, n, snr, seed):
    d = simulate(p=NGENE, n=n, snr=snr, seed=seed, net_sd=NETSD, n_parents=NPAR)
    fitted = M.fit_all(d, input_mode, seed=seed, methods=METHODS)
    base = dict(input=input_mode, n=n, snr=snr, seed=seed,
                true_edges=int(d["adj_main"].sum()))
    return EV.rows_for(fitted, d["adj_main"], d["adj_skeleton"], base)


def main():
    todo = [(im, n, s, SEED_BASE + r)
            for n in NS for s in SNRS for r in range(REPS) for im in INPUTS]
    have = done_keys(OUT)
    todo = [c for c in todo if tuple(map(str, c)) not in have]
    total = len(NS) * len(SNRS) * REPS * len(INPUTS)
    print(f"grid {total} cells, {total - len(todo)} done, {len(todo)} to run", flush=True)

    t0 = time.time()
    k = 0
    for cell in todo:
        if time.time() - t0 > BUDGET:
            print(f"budget spent: {k} cells this run, {total - len(todo) + k}/{total}", flush=True)
            return
        try:
            rows = one_cell(*cell)
        except Exception as e:
            print(f"  skipped {cell}: {type(e).__name__}: {e}", flush=True)
            continue
        pd.DataFrame(rows).to_csv(OUT, mode="a", header=not os.path.exists(OUT), index=False)
        k += 1
        if k % 10 == 0:
            print(f"  {k}/{len(todo)}  {(time.time() - t0) / k:.1f}s per cell", flush=True)
    print(f"ALL DONE: {k} cells this run", flush=True)


if __name__ == "__main__":
    main()
