"""
evotopo_grn.py — EvoTopo-GRN: 沿表达指数(EI)的基因调控网络重构

三步：
  1. power_equation_fit  每个基因用 nls 拟合 log10(count+1) = a * EI^b
  2. detrend             扣掉拟合的幂律趋势，得到残差矩阵
  3. select_edges        对每个 target 做自适应 LASSO 变量选择

扣趋势是关键：原始计数上 20 个回归量被共同的幂律趋势主导、两两相关 > 0.99，
LASSO 等于在一堆几乎相同的预测变量里挑；扣掉之后剩下的结构化残差才是网络信号。

依赖：numpy, scipy, scikit-learn
"""
import numpy as np
from scipy.optimize import curve_fit
from sklearn.linear_model import Lasso, Ridge
from sklearn.model_selection import KFold


# ---------------------------------------------------------------- 第 1 步
def power_equation_fit(counts, n_grid=20):
    """counts: 基因 x 样本 的计数矩阵

    返回 dict:
      original_data  按文库大小升序重排后的计数矩阵
      Time           每个样本的表达指数 EI = log10(colSums + 1)
      power_par      每个基因的 (a, b)，满足 log10(count+1) ≈ a * EI^b
      power_fit      在 n_grid 个等距 EI 点上的插值曲线
      keep           成功拟合的基因下标
    """
    order = np.argsort(counts.sum(0), kind="stable")     # 按文库大小升序
    counts = counts[:, order]
    EI = np.log10(counts.sum(0) + 1.0)                   # ★ 表达指数
    logc = np.log10(counts + 1.0)

    keep, pars = [], []
    for g in range(counts.shape[0]):
        y = logc[g]
        try:                                             # 先 log-lm 取初值，再 nls
            b0, la0 = np.polyfit(np.log(EI), np.log(np.maximum(y, 1e-9)), 1)
            (a, b), _ = curve_fit(lambda x, a, b: a * x**b, EI, y,
                                  p0=[np.exp(la0), b0], maxfev=200000)
            keep.append(g); pars.append((a, b))
        except Exception:
            pass                                         # 拟合失败的基因剔除
    keep = np.asarray(keep); pars = np.asarray(pars)

    grid = np.linspace(EI.min(), EI.max(), n_grid)
    return dict(original_data=counts[keep], Time=EI, power_par=pars,
                power_fit=pars[:, 0:1] * grid[None, :]**pars[:, 1:2],
                grid=grid, keep=keep, sample_order=order)


# ---------------------------------------------------------------- 第 2 步
def detrend(fit):
    """扣掉幂律趋势，返回 样本 x 基因 的残差矩阵"""
    a, b = fit["power_par"][:, 0:1], fit["power_par"][:, 1:2]
    trend = a * fit["Time"][None, :]**b
    return (np.log10(fit["original_data"] + 1.0) - trend).T


# ---------------------------------------------------------------- 第 3 步
def _cv_coef(X, y, alphas, ridge=False, rule="1se", nfolds=10, seed=0):
    """沿正则化路径做 K 折交叉验证，按 lambda.min 或 lambda.1se 取系数"""
    kf = KFold(n_splits=nfolds, shuffle=True, random_state=seed)
    mse = np.zeros((len(alphas), nfolds))
    for f, (tr, te) in enumerate(kf.split(X)):
        for i, a in enumerate(alphas):
            m = Ridge(alpha=a) if ridge else Lasso(alpha=a, max_iter=20000)
            m.fit(X[tr], y[tr])
            mse[i, f] = np.mean((y[te] - m.predict(X[te]))**2)
    cvm = mse.mean(1); cvse = mse.std(1, ddof=1) / np.sqrt(nfolds)
    i_min = int(np.argmin(cvm))
    if rule == "1se":                                    # glmnet 的 1se 规则
        cand = np.where(cvm <= cvm[i_min] + cvse[i_min])[0]
        i_sel = int(cand[np.argmax(alphas[cand])]) if len(cand) else i_min
    else:
        i_sel = i_min
    m = Ridge(alpha=alphas[i_sel]) if ridge else Lasso(alpha=alphas[i_sel], max_iter=20000)
    m.fit(X, y)
    return m.coef_


def select_edges(D, col, reduction=True, rule="1se", seed=0):
    """对 target `col` 做自适应 LASSO。D: 样本 x 基因。返回长度 p 的 0/1 入边指示"""
    p = D.shape[1]
    y = D[:, col]
    others = np.array([k for k in range(p) if k != col])
    x = D[:, others]

    if reduction:                                        # 按边际相关预筛到 (p-1)/log(p-1) 个
        with np.errstate(invalid="ignore"):
            r = np.abs([np.corrcoef(x[:, c], y)[0, 1] for c in range(x.shape[1])])
        r = np.nan_to_num(r)
        nk = max(int(x.shape[1] / np.log(x.shape[1])), 1)
        sel = np.argsort(-r)[:nk]
        x, others = x[:, sel], others[sel]

    xs = (x - x.mean(0)) / np.maximum(x.std(0, ddof=1), 1e-12)
    ys = y - y.mean()
    out = np.zeros(p, int)
    if ys.std() < 1e-12:
        return out

    # ridge 给自适应权重
    rc = _cv_coef(xs, ys, np.logspace(3, -2, 40), ridge=True, rule="1se", seed=seed)
    w = np.maximum(np.abs(rc), 1e-8)
    # 变量替换 x_k -> x_k * w_k，等价于 penalty.factor = 1/w
    xa = xs * w[None, :]
    amax = np.max(np.abs(xa.T @ ys)) / len(ys)
    bc = _cv_coef(xa, ys, amax * np.logspace(0, -3, 40), rule=rule, seed=seed)

    out[others[np.abs(bc) > 1e-10]] = 1
    return out


# ---------------------------------------------------------------- 入口
def evotopo_grn(counts, n_grid=20, reduction=True, rule="1se", seed=0):
    """counts: 基因 x 样本 计数矩阵 → (p x p 邻接矩阵, fit)
    adj[target, regulator] = 1"""
    fit = power_equation_fit(counts, n_grid=n_grid)
    resid = detrend(fit)
    p = resid.shape[1]
    adj = np.vstack([select_edges(resid, c, reduction, rule, seed) for c in range(p)])
    return adj, fit


# ---------------------------------------------------------------- 评价
def calculate_metrics(pred, true):
    """在非对角位置上算 TPR / FPR / MCC"""
    pred = np.asarray(pred, int); true = np.asarray(true, int)
    if pred.ndim == 2:
        o = ~np.eye(pred.shape[0], dtype=bool)
        pred, true = pred[o], true[o]
    TP = int(((pred == 1) & (true == 1)).sum()); FP = int(((pred == 1) & (true == 0)).sum())
    TN = int(((pred == 0) & (true == 0)).sum()); FN = int(((pred == 0) & (true == 1)).sum())
    den = np.sqrt(float(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    return dict(TP=TP, FP=FP, TN=TN, FN=FN,
                TPR=TP / (TP + FN) if TP + FN else 0.0,
                FPR=FP / (FP + TN) if FP + TN else 0.0,
                MCC=(TP * TN - FP * FN) / den if den else 0.0)


if __name__ == "__main__":
    from sim_power_EI import sim_power_EI
    d = sim_power_EI(p=20, n=100, snr=10.0, seed=101)
    adj, fit = evotopo_grn(d["expr"].values)
    m = calculate_metrics(adj, d["adj_main"].values)
    print(f"EvoTopo-GRN  TPR={m['TPR']:.3f}  FPR={m['FPR']:.3f}  MCC={m['MCC']:+.3f}")
