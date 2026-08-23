"""sim_power_EI.py — 20 基因 x 50 样本模拟数据（幂律版，对齐 Figure1）

信号：log10(count + 1) = mu_g(EI) + eps_g ,  mu_g(x) = a_g * x^(b_g)
      → log–log 图上单调上升的曲线，散点紧贴曲线，与 Figure1 一致

网络：编码在**残差**里。理由：ode_1.R 的 qdODE_parallel 第一步是
      get_interaction(result$original_data, ...) —— 在原始计数矩阵上跑
      adaptive LASSO 做变量选择，网络来自样本间共变，而不是曲线形状。
      （若把网络编码进曲线形状，20 条幂律曲线两两相关 > 0.998，
        设计矩阵条件数 ~1e17，任何随机支撑都能拟合，真值不可辨识。）
      残差用稀疏精度矩阵 Omega 的高斯图模型：eps ~ N(0, Omega^{-1})，
      真值 = support(Omega)，正是邻域选择 LASSO 的目标。

计数：libsize = round(10^EI) 的 multinomial 抽样
      → colSums 精确等于设计值 → EI 可复算、严格单调
"""
import numpy as np, pandas as pd


def power_from_endpoints(mu_lo, mu_hi, x_lo, x_hi):
    """由端点值反解 mu(x) = a * x^b 的 (a, b)"""
    b = np.log(mu_hi / mu_lo) / np.log(x_hi / x_lo)
    return mu_lo / x_lo**b, b


def make_precision(p, sparsity, w_range, rng):
    """稀疏、对角占优（因而正定）的精度矩阵"""
    Om = np.zeros((p, p))
    for i in range(p):
        for j in range(i + 1, p):
            if rng.uniform() < sparsity:
                Om[i, j] = Om[j, i] = rng.uniform(*w_range) * rng.choice([1, -1])
    np.fill_diagonal(Om, np.abs(Om).sum(1) + 0.1)
    return Om


def sim_power_EI(p=20, n=50, n_bg=480,
                 EI_range=(6.95, 7.5),      # 对齐 Figure1 横轴 10^6.95 ~ 10^7.5
                 mu_lo_range=(0.75, 2.7),   # EI 下端的 log10(count+1)
                 slope_range=(0.35, 1.05),  # d mu / d EI
                 mu_cap=3.1,                # log10(count+1) 上限
                 snr=None,                  # 信噪比：控制**独立测量噪声**（网络信号不受影响）
                 net_sd=0.10,               # 结构化残差 sd —— 网络信号，固定，不随 snr 变
                 noise_sd=None,             # 兼容旧接口：给了它就退回“网络随 snr 缩放”的老行为
                 net_sparsity=0.12,         # 残差网络稀疏度（net_type='ggm' 时用）
                 w_range=(0.4, 0.8),        # 边强度
                 net_type="ggm",            # 'ggm' = 无向高斯图模型；'indeg1' = 每基因恰好 1 个调控者
                 seed=42):
    rng = np.random.default_rng(seed)
    x_lo, x_hi = EI_range

    # ---------- 1. EI 网格：有序、不等距（像真实样本） ----------
    u = np.sort(rng.uniform(0, 1, n))
    EI = x_lo + (x_hi - x_lo) * u
    N = np.round(10**EI).astype(np.int64)                  # 设计文库大小

    # ---------- 2. 焦点基因的幂律曲线 mu_j(EI) = a_j * EI^b_j ----------
    mu_lo = rng.uniform(*mu_lo_range, p)
    slope = rng.uniform(*slope_range, p)
    mu_hi = np.minimum(mu_lo + slope * (x_hi - x_lo), mu_cap)
    a0, b0 = power_from_endpoints(mu_lo, mu_hi, x_lo, x_hi)
    mu = a0[None, :] * EI[:, None]**b0[None, :]            # n x p

    # ---------- 3. 残差：拆成“网络信号”和“测量噪声”两部分 ----------
    #   eps_net  ~ N(0, Omega^{-1})，sd 固定为 net_sd  → 网络信号，**不随 snr 变**
    #   eps_meas ~ N(0, sigma_g^2)  独立同分布        → 测量噪声，sigma_g = sqrt(var(mu_g)/snr)
    # 这样 snr 越大 = 测量噪声越小 = 网络越容易被看见，MCC 才会随 snr 上升。
    # 老接口 noise_sd 保留：给了它就退回“网络随 snr 一起缩小”的旧行为（趋势会是反的）。
    if net_type == "indeg1":
        # 有向 SEM：每个基因恰好 1 个父节点（根节点除外），eps = (I-W)^{-1} u
        # in-degree <= 1 => 没有 v-structure => 道德图 == 真骨架，邻域选择的目标就是真值
        perm = rng.permutation(p)
        W = np.zeros((p, p))
        for r in range(1, p):
            j = perm[r]; k = perm[rng.integers(0, r)]          # 父节点必须序在前 => 无环
            W[j, k] = rng.uniform(*w_range) * rng.choice([1, -1])
        Omega = W                                              # 复用返回位（真值取 support）
        A_inv = np.linalg.inv(np.eye(p) - W)
        eps_net = rng.normal(0, 1, (n, p)) @ A_inv.T
    else:
        Omega = make_precision(p, net_sparsity, w_range, rng)
        L = np.linalg.cholesky(np.linalg.inv(Omega))
        eps_net = rng.normal(0, 1, (n, p)) @ L.T
    eps_net = (eps_net - eps_net.mean(0)) / eps_net.std(0, ddof=1)

    if noise_sd is not None:                      # 旧行为
        sd_g = (np.sqrt(mu.var(0, ddof=1) / snr) if snr is not None
                else np.full(p, float(noise_sd)))
        eps = eps_net * sd_g[None, :]
        sigma_g = np.zeros(p)
    else:                                          # 新行为（默认）
        eps_net = eps_net * float(net_sd)
        sigma_g = (np.sqrt(mu.var(0, ddof=1) / snr) if snr is not None else np.zeros(p))
        eps = eps_net + rng.normal(0, 1, (n, p)) * sigma_g[None, :]
        sd_g = np.sqrt(net_sd**2 + sigma_g**2)

    # ---------- 4. 背景基因：把文库总量撑到 10^7 量级 ----------
    bg_lo = np.maximum(rng.normal(3.4, 0.9, n_bg), 0.05)
    bg_hi = np.maximum(bg_lo + rng.uniform(0.7, 1.1, n_bg) * (x_hi - x_lo), 0.06)
    a_bg, b_bg = power_from_endpoints(bg_lo, bg_hi, x_lo, x_hi)
    lam_bg = np.maximum(10**(a_bg[None, :] * EI[:, None]**b_bg[None, :]) - 1.0, 1e-9)
    lam_focal = np.maximum(10**mu - 1.0, 1e-9)
    # 标定背景总量，使 sum(lambda) ≈ 10^EI（归一化偏移 c(EI) ≈ 0）
    lam_bg *= np.median((10**EI - lam_focal.sum(1)) / lam_bg.sum(1))

    # ---------- 5. 比例 + multinomial 抽样 ----------
    lam_obs = np.hstack([lam_focal * 10**eps,
                         lam_bg * 10**rng.normal(0, float(np.median(sd_g)), lam_bg.shape)])
    P = lam_obs / lam_obs.sum(1, keepdims=True)
    counts = np.vstack([rng.multinomial(N[i], P[i]) for i in range(n)])

    genes = [f"G{i+1}" for i in range(p)]
    samples = [f"S{i+1}" for i in range(n)]
    expr_all = pd.DataFrame(counts.T, index=genes + [f"BG{i+1}" for i in range(n_bg)],
                            columns=samples)
    expr = expr_all.loc[genes]
    EI_real = np.log10(expr_all.values.sum(0) + 1.0)

    # ---------- 6. 真值曲线与幂律参数 ----------
    c_off = EI - np.log10(np.hstack([lam_focal, lam_bg]).sum(1))   # 归一化共同偏移
    mu_true = mu + c_off[:, None]
    lx = np.log(EI)
    a_true, b_true, r2 = np.zeros(p), np.zeros(p), np.zeros(p)
    for j in range(p):
        s, i0 = np.polyfit(lx, np.log(mu_true[:, j]), 1)
        a_true[j], b_true[j] = np.exp(i0), s
        fit = a_true[j] * EI**b_true[j]
        r2[j] = 1 - ((mu_true[:, j] - fit)**2).sum() / ((mu_true[:, j] - mu_true[:, j].mean())**2).sum()

    adj = (np.abs(Omega) > 1e-12).astype(int)
    np.fill_diagonal(adj, 0)
    if net_type == "indeg1":
        adj = ((adj + adj.T) > 0).astype(int)     # 观测数据无法定向，真值按无向骨架评价

    return dict(
        p=p, n=n, EI=EI_real, times=EI_real, N=N, noise_sd=sd_g, snr=snr, c_off=c_off,
        net_sd_used=net_sd, meas_sd=sigma_g, eps_net_only=None,
        expr=expr, expr_all=expr_all,
        Y=pd.DataFrame(np.log10(expr.values + 1.0).T, columns=genes),
        mu_true=pd.DataFrame(mu_true, columns=genes),
        eps_true=pd.DataFrame(eps, columns=genes),
        power_par=pd.DataFrame({"a": a_true, "b": b_true, "R2_powerlaw": r2}, index=genes),
        Omega=pd.DataFrame(Omega, index=genes, columns=genes),
        adj_main=pd.DataFrame(adj, index=genes, columns=genes),
        genes=genes, samples=samples)


if __name__ == "__main__":
    d = sim_power_EI()
    e, EI = d["expr"].values, d["EI"]
    print("EI 范围        : %.4f ~ %.4f  (10^%.2f ~ 10^%.2f)" % (EI.min(), EI.max(), EI.min(), EI.max()))
    print("EI 严格单调    :", bool(np.all(np.diff(EI) > 0)))
    print("EI 可复算      :", np.allclose(EI, np.log10(d["expr_all"].values.sum(0) + 1)))
    print("归一化偏移跨度 : %.4f" % (d["c_off"].max() - d["c_off"].min()))
    print("焦点基因计数   : %d ~ %d" % (e.min(), e.max()))
    print("幂指数 b       : %.2f ~ %.2f" % (d["power_par"].b.min(), d["power_par"].b.max()))
    print("真曲线幂律 R^2 : 最小 %.6f" % d["power_par"].R2_powerlaw.min())
    print("残差网络边数   : %d 条（无向） / %d" % (d["adj_main"].values.sum() // 2, 20 * 19 // 2))

    d["expr"].to_csv("pw_expr_20g_50s.csv")
    d["expr_all"].to_csv("pw_expr_all_500g_50s.csv")
    pd.DataFrame({"sample": d["samples"], "libsize": d["expr_all"].values.sum(0), "EI": EI}
                 ).to_csv("pw_EI.csv", index=False)
    d["Y"].round(6).to_csv("pw_Y_50s_20g.csv", index=False)
    d["mu_true"].round(6).to_csv("pw_mu_true_50s_20g.csv", index=False)
    d["eps_true"].round(6).to_csv("pw_eps_true_50s_20g.csv", index=False)
    d["power_par"].round(8).to_csv("pw_power_par.csv")
    d["Omega"].round(6).to_csv("pw_Omega_true.csv")
    d["adj_main"].to_csv("pw_adj_true.csv")
