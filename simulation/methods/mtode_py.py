"""mtode_py.py — snr_1.R 里 MTODE 的 Python 移植（用于在没有 R 的环境里核算 MCC）

逐段对应 R 版：
  预平滑(smoothing spline) -> 多项式基 poly_basis_1d -> 梯形积分 cumtrapz
  -> 每个 target 分组正交化(QR) -> 块对角设计阵 -> 组稀疏回归 -> 还原 beta
  -> 组贡献 f_group -> adj_est

唯一替换：ADSIHT 换成 **group lasso + BIC 选 lambda**。
理由：ADSIHT 是 R 包，本环境装不上；snr_1.R 里注释掉的备选正是
      cv.sparsegl / cv.grpreg，也就是 group lasso。
因为每组做过 QR 正交化，组内设计阵是标准正交，块坐标下降有闭式解。
"""
import numpy as np
from scipy.interpolate import UnivariateSpline, make_smoothing_spline
from scipy.integrate import cumulative_trapezoid


def calculate_metrics(predicted, true):
    """与 snr_1.R 的 calculate_metrics 完全一致"""
    pr = np.asarray(predicted).ravel().astype(int)
    tr = np.asarray(true).ravel().astype(int)
    TP = np.sum((pr == 1) & (tr == 1)); FP = np.sum((pr == 1) & (tr == 0))
    TN = np.sum((pr == 0) & (tr == 0)); FN = np.sum((pr == 0) & (tr == 1))
    TPR = TP / (TP + FN) if (TP + FN) else 0.0
    FPR = FP / (FP + TN) if (FP + TN) else 0.0
    den = np.sqrt(float(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    MCC = (TP * TN - FP * FN) / den if den else 0.0
    return dict(TPR=TPR, FPR=FPR, MCC=MCC, TP=TP, FP=FP, TN=TN, FN=FN)


def poly_basis_1d(x, degree=3):
    return np.column_stack([x**d for d in range(1, degree + 1)])


def _group_lasso_orth(Q, y, groups, lam, w, beta=None, tol=1e-8, maxit=300):
    """组内正交（Q_g'Q_g = I）时的 group lasso 块坐标下降；支持热启动"""
    beta = np.zeros(Q.shape[1]) if beta is None else beta.copy()
    r = y - Q @ beta
    for _ in range(maxit):
        delta = 0.0
        for g, cols in groups.items():
            bg = beta[cols]
            z = Q[:, cols].T @ r + bg              # 组内正交 → 直接加回
            nz = np.linalg.norm(z)
            new = z * max(0.0, 1.0 - lam * w[g] / nz) if nz > 0 else np.zeros_like(z)
            d = new - bg
            if np.any(d):
                r -= Q[:, cols] @ d
                beta[cols] = new
                delta = max(delta, np.max(np.abs(d)))
        if delta < tol:
            break
    return beta


def _group_iht(Q, y, groups, ic_mult=1.0, ebic_gamma=0.0, maxit=50):
    """迭代组硬阈值（IHT）+ IC 选组数 —— 对应 ADSIHT 的行为。
    关键：每步用当前残差的梯度重新筛组并refit，而不是只做一次边际筛选。
    组间并不正交，一次性边际筛选在共线设计下会选错。"""
    n = len(y); G = list(groups.keys())
    gl = {g: groups[g] for g in G}
    kmax = min(len(G), max(1, n // 6))
    best = (np.inf, np.zeros(Q.shape[1]))
    prev_keep = None
    for k in range(1, kmax + 1):
        b = np.zeros(Q.shape[1])
        keep = None
        for _ in range(maxit):
            u = b + Q.T @ (y - Q @ b)                      # 梯度步（步长 1）
            norms = {g: np.linalg.norm(u[gl[g]]) for g in G}
            new_keep = sorted(sorted(G, key=lambda g: -norms[g])[:k])
            cols = np.concatenate([gl[g] for g in new_keep])
            if len(cols) >= n - 1:
                break
            bk, *_ = np.linalg.lstsq(Q[:, cols], y, rcond=None)
            b = np.zeros(Q.shape[1]); b[cols] = bk          # 硬阈值 + debias
            if new_keep == keep:
                break
            keep = new_keep
        if keep is None:
            break
        cols = np.concatenate([gl[g] for g in keep])
        rss = np.sum((y - Q[:, cols] @ b[cols])**2)
        pen = ic_mult * len(cols) * np.log(n)
        if ebic_gamma > 0:
            from math import lgamma
            GG = len(G)
            pen += 2.0 * ebic_gamma * (lgamma(GG + 1) - lgamma(k + 1) - lgamma(GG - k + 1))
        val = n * np.log(max(rss, 1e-12) / n) + pen
        if val < best[0]:
            best = (val, b.copy())
    return best[1]


def MTODE(Y, times, M=3, n_expand=50, smooth_s=None, verbose=False,
          selector="grplasso", ic_mult=1.0, ebic_gamma=0.0):
    Y = np.asarray(Y, float)
    n, p = Y.shape
    t0, rng_t = times.min(), times.max() - times.min()
    tn = (times - t0) / rng_t
    t_new = np.linspace(tn.min(), tn.max(), n_expand)
    t_restored = t_new * rng_t + t0

    # ---- 1. 预平滑（对应 R 的 smooth.spline，默认 GCV 定平滑量） ----
    x_smooth = np.zeros((n, p)); x_smooth2 = np.zeros((n_expand, p))
    for j in range(p):
        if smooth_s == "none":
            sp = UnivariateSpline(tn, Y[:, j], s=0, k=3)     # 插值，不平滑
        elif smooth_s is None:
            sp = make_smoothing_spline(tn, Y[:, j])          # lam=None → GCV
        else:
            sp = UnivariateSpline(tn, Y[:, j], s=smooth_s)
        x_smooth[:, j] = sp(tn); x_smooth2[:, j] = sp(t_new)

    # ---- 2. 多项式基 + 梯形积分 ----
    phi = np.hstack([poly_basis_1d(x_smooth[:, j], M) for j in range(p)])
    phi_int = np.column_stack([cumulative_trapezoid(phi[:, c], tn, initial=0)
                               for c in range(phi.shape[1])])
    phi2 = np.hstack([poly_basis_1d(x_smooth2[:, j], M) for j in range(p)])
    phi_int2 = np.column_stack([cumulative_trapezoid(phi2[:, c], t_new, initial=0)
                                for c in range(phi2.shape[1])])

    # ---- 3. 分组正交化（对应 construct_group；各 target 结构相同，只算一次） ----
    X_local = np.column_stack([tn, phi_int])          # n x (1 + p*M)
    gid = np.concatenate([[0], np.repeat(np.arange(p), M)])   # 时间列并入第 1 组
    groups, Qs, Rs, mXs, srcs = {}, [], [], [], []
    off = 0
    for g in range(p):
        cols = np.where(gid == g)[0]
        Xg = X_local[:, cols]
        mXg = Xg.mean(0); Xc = Xg - mXg
        Qg, Rg = np.linalg.qr(Xc)
        groups[g] = np.arange(off, off + Qg.shape[1]); off += Qg.shape[1]
        Qs.append(Qg); Rs.append(Rg); mXs.append(mXg); srcs.append(cols)
    Q_all = np.hstack(Qs)
    w = {g: np.sqrt(len(groups[g])) for g in groups}

    # ---- 4. 标准化 + 组稀疏回归（ADSIHT → group lasso + BIC） ----
    mY, sY = Y.mean(0), Y.std(0, ddof=1)
    Ys = (Y - mY) / sY

    lam_max = max(np.linalg.norm(Q_all[:, groups[g]].T @ Ys[:, j]) / w[g]
                  for j in range(p) for g in groups)
    lams = lam_max * np.logspace(0, -3, 40)          # 由大到小，热启动

    beta_Q = np.zeros((Q_all.shape[1], p))
    if selector == "iht":
        for j in range(p):
            beta_Q[:, j] = _group_iht(Q_all, Ys[:, j], groups,
                                      ic_mult=ic_mult, ebic_gamma=ebic_gamma)
    else:
        for j in range(p):
            best = (np.inf, None); warm = None
            for lam in lams:
                b = _group_lasso_orth(Q_all, Ys[:, j], groups, lam, w, beta=warm)
                warm = b
                rss = np.sum((Ys[:, j] - Q_all @ b)**2)
                df = np.sum(np.abs(b) > 1e-10)
                bic = n * np.log(max(rss, 1e-12) / n) + df * np.log(n)
                if bic < best[0]:
                    best = (bic, b)
            beta_Q[:, j] = best[1]

    # ---- 5. 还原原始尺度的 beta（对应 recover_beta） ----
    ncol_local = X_local.shape[1]
    beta_all = np.zeros((ncol_local + 1, p))          # 第 0 行 = 截距
    for j in range(p):
        bX = np.zeros(ncol_local); mX = np.zeros(ncol_local)
        for g in groups:
            bg = np.linalg.solve(Rs[g], beta_Q[groups[g], j])
            bX[srcs[g]] = bg; mX[srcs[g]] = mXs[g]
        beta_all[1:, j] = sY[j] * bX
        beta_all[0, j] = mY[j] - sY[j] * float(mX @ bX)

    # ---- 6. 组贡献 f_group → adj_est ----
    adj_est = np.zeros((p, p), int)
    f_group_est = []
    for jj in range(p):
        betam = beta_all[2:, jj]                      # 去掉截距和时间项 → p*M
        est = np.zeros((n_expand, p))
        for g in range(p):
            cols = np.arange(g * M, (g + 1) * M)
            est[:, g] = phi_int2[:, cols] @ betam[cols]
        est[:, jj] += beta_all[1, jj] * t_new + beta_all[0, jj]
        f_group_est.append(est)
        adj_est[jj, np.where(est.sum(0) != 0)[0]] = 1

    return dict(beta_all=beta_all, adj_est=adj_est, f_group_est=f_group_est,
                times=times, times_restored=t_restored, x_smooth=x_smooth)
