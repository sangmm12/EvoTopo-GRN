"""Simulator for EvoTopo-GRN benchmarking.

Generative model
----------------
1. Expression index    EI_i = 6.95 + 0.55 * u_(i),  u_(i) sorted uniforms.
                       Design library size N_i = round(10 ** EI_i).
2. Allometric trend    mu_j(EI) = a_j * EI ** b_j on the log10(count + 1) scale,
                       with (a_j, b_j) obtained by inverting the endpoint values
                       mu_lo ~ U(0.75, 2.70) and mu_hi = min(mu_lo + s * 0.55, 3.10),
                       s ~ U(0.35, 1.05).
3. Network             A random permutation fixes a topological order; each
                       non-root gene draws at most `n_parents` parents from the
                       genes preceding it, with weights W_jk = +/- U(0.4, 0.8).
                       Residuals follow eps = W eps + u, i.e. eps = (I - W)^-1 u.
                       Columns are standardised afterwards, which removes
                       varsortability so that orientation cannot be read off the
                       marginal variances.
4. Noise               snr controls independent measurement noise only; the
                       network signal is fixed at net_sd and does not scale with snr.
5. Counts              Focal intensities lam = max(10 ** mu - 1, 1e-9) * 10 ** eps.
                       n_bg background genes carry the library size; the whole
                       intensity vector is drawn as Multinomial(N_i, lam / sum lam),
                       so the realised index reproduces the design values.

Ground truth is the support of W (rows = targets, columns = regulators).
"""
import numpy as np
import pandas as pd


def power_from_endpoints(mu_lo, mu_hi, x_lo, x_hi):
    """Invert mu(x) = a * x ** b from its values at two endpoints."""
    b = np.log(mu_hi / mu_lo) / np.log(x_hi / x_lo)
    return mu_lo / x_lo ** b, b


def simulate(p=20, n=100, n_bg=480,
             EI_range=(6.95, 7.5),
             mu_lo_range=(0.75, 2.7),
             slope_range=(0.35, 1.05),
             mu_cap=3.1,
             snr=10.0,
             net_sd=0.10,
             w_range=(0.4, 0.8),
             n_parents=1,
             seed=0):
    """Return one simulated data set.

    Parameters
    ----------
    p, n        focal genes, samples
    n_bg        background genes carrying the library size
    snr         signal-to-noise ratio; controls measurement noise only
    net_sd      dispersion of the structured residual (the network signal)
    w_range     magnitude of the non-zero edge weights, signs equiprobable
    n_parents   maximum in-degree of the directed acyclic graph

    Returns
    -------
    dict with keys
        expr         (p, n)      focal gene counts
        expr_all     (p+n_bg, n) full count matrix
        Y            (n, p)      log10(expr + 1)
        EI           (n,)        expression index recomputed from expr_all
        W            (p, p)      true edge weights
        adj_main     (p, p)      directed ground truth, rows = targets
        adj_skeleton (p, p)      symmetrised ground truth
    """
    rng = np.random.default_rng(seed)
    x_lo, x_hi = EI_range

    # 1. Ordered, unevenly spaced index grid
    EI = x_lo + (x_hi - x_lo) * np.sort(rng.uniform(0, 1, n))
    N = np.round(10 ** EI).astype(np.int64)

    # 2. Per-gene allometric trend on the log10(count + 1) scale
    mu_lo = rng.uniform(*mu_lo_range, p)
    slope = rng.uniform(*slope_range, p)
    mu_hi = np.minimum(mu_lo + slope * (x_hi - x_lo), mu_cap)
    a0, b0 = power_from_endpoints(mu_lo, mu_hi, x_lo, x_hi)
    mu = a0[None, :] * EI[:, None] ** b0[None, :]

    # 3. Directed acyclic structural model on the residuals: eps = W eps + u
    perm = rng.permutation(p)
    W = np.zeros((p, p))
    for r in range(1, p):
        j = perm[r]
        for k in rng.choice(perm[:r], min(int(n_parents), r), replace=False):
            W[j, k] = rng.uniform(*w_range) * rng.choice([1, -1])
    A_inv = np.linalg.inv(np.eye(p) - W)
    u_sd = rng.uniform(0.6, 1.6, p)
    eps_net = (rng.normal(0, 1, (n, p)) * u_sd[None, :]) @ A_inv.T
    eps_net = (eps_net - eps_net.mean(0)) / eps_net.std(0, ddof=1)

    eps_net = eps_net * float(net_sd)
    sigma_g = np.sqrt(mu.var(0, ddof=1) / snr) if snr else np.zeros(p)
    eps = eps_net + rng.normal(0, 1, (n, p)) * sigma_g[None, :]
    sd_g = np.sqrt(net_sd ** 2 + sigma_g ** 2)

    # 4. Background genes bring the library total up to 10 ** EI
    bg_lo = np.maximum(rng.normal(3.4, 0.9, n_bg), 0.05)
    bg_hi = np.maximum(bg_lo + rng.uniform(0.7, 1.1, n_bg) * (x_hi - x_lo), 0.06)
    a_bg, b_bg = power_from_endpoints(bg_lo, bg_hi, x_lo, x_hi)
    lam_bg = np.maximum(10 ** (a_bg[None, :] * EI[:, None] ** b_bg[None, :]) - 1.0, 1e-9)
    lam_focal = np.maximum(10 ** mu - 1.0, 1e-9)
    lam_bg *= np.median((10 ** EI - lam_focal.sum(1)) / lam_bg.sum(1))

    # 5. Multinomial sampling at the design library size
    lam_obs = np.hstack([lam_focal * 10 ** eps,
                         lam_bg * 10 ** rng.normal(0, float(np.median(sd_g)), lam_bg.shape)])
    P = lam_obs / lam_obs.sum(1, keepdims=True)
    counts = np.vstack([rng.multinomial(N[i], P[i]) for i in range(n)])

    genes = [f"G{i + 1}" for i in range(p)]
    samples = [f"S{i + 1}" for i in range(n)]
    expr_all = pd.DataFrame(counts.T, index=genes + [f"BG{i + 1}" for i in range(n_bg)],
                            columns=samples)
    expr = expr_all.loc[genes]

    adj = (np.abs(W) > 1e-12).astype(int)
    np.fill_diagonal(adj, 0)
    return dict(
        p=p, n=n, snr=snr, seed=seed,
        expr=expr.values.astype(float), expr_all=expr_all,
        Y=np.log10(expr.values.astype(float) + 1.0).T,
        EI=np.log10(expr_all.values.sum(0) + 1.0),
        W=pd.DataFrame(W, index=genes, columns=genes),
        adj_main=adj,
        adj_skeleton=((adj + adj.T) > 0).astype(int),
        genes=genes, samples=samples)


if __name__ == "__main__":
    d = simulate(n=100, snr=10, seed=1)
    EI = d["EI"]
    print("EI range      : %.4f - %.4f" % (EI.min(), EI.max()))
    print("EI monotone   :", bool(np.all(np.diff(EI) > 0)))
    print("counts        :", d["expr"].shape, " background:", d["expr_all"].shape)
    print("directed edges:", int(d["adj_main"].sum()))
