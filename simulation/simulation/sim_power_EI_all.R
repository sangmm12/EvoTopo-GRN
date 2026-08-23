###############################################################################
# sim_power_EI_all.R
#   20 基因 x 50 样本模拟数据 —— 单文件版（生成 + 输出 + 自检 + 诊断图）
#
#   运行：  Rscript sim_power_EI_all.R
#   依赖：  base R（无需任何第三方包）
#
# -----------------------------------------------------------------------------
# 【信号形式】  log10(count + 1) = mu_g(EI) + eps_g ,   mu_g(x) = a_g * x^(b_g)
#   每个基因的 log10(count+1) 是 EI 的幂函数：log–log 图上单调上升、散点贴线，
#   与 Figure1 一致。(a_g, b_g) 由端点反解，所以无噪真曲线的幂律 R^2 >= 0.9996 —
#   幂律是信号本身，不是拟合出来的近似。
#
# 【网络编码在残差里，不在曲线形状里】  两个理由：
#   (1) ode_1.R 的 qdODE_parallel 第一步是
#         get_interaction(result$original_data, ...)
#       —— 在原始计数矩阵上跑 adaptive LASSO 做变量选择，
#       网络来自样本间共变，而不是平滑曲线的形状。
#   (2) 实测：20 条幂律曲线在 EI in [6.95, 7.5] 上两两相关系数最小 0.99936，
#       设计矩阵条件数 ~1.4e17。任何随机支撑都能把 dmu_j/dEI 拟合到 1e-6 残差，
#       真值不可辨识，benchmark 没有意义。
#   残差改用稀疏精度矩阵的高斯图模型 eps ~ N(0, Omega^{-1})，
#   真值 = support(Omega)，正是邻域选择 LASSO 的目标，可辨识。
#
# 【EI 的自洽性】  文库大小取 N = round(10^EI)，multinomial 抽样
#   → colSums(expr_all) 精确等于 N → EI = log10(colSums+1) 可复算、严格单调。
#   另配 n_bg 个背景基因把文库撑到 10^7 量级，让 20 个焦点基因的计数落在
#   10^0 ~ 10^3 —— 与 Figure1 的两根轴都对上。
#
# 【输出文件】
#   sim_power_EI.rds           完整对象
#   pw_expr_20g_50s.csv        20 x 50 焦点基因计数矩阵
#   pw_expr_all_500g_50s.csv   500 x 50 全矩阵（含背景基因）—— EI 由它算出
#   pw_EI.csv                  每样本 libsize 与 EI
#   pw_Y_50s_20g.csv           50 x 20 = t(log10(expr+1))
#   pw_mu_true_50s_20g.csv     50 x 20 无噪真值曲线
#   pw_eps_true_50s_20g.csv    50 x 20 结构化残差
#   pw_power_par.csv           每基因真值 (a, b) 及幂律 R^2
#   pw_Omega_true.csv          残差精度矩阵
#   pw_adj_true.csv            残差网络邻接矩阵（无向）
#   pw_check.png               20 个小图，与 Figure1 同版式
#   run_log.txt                自检日志
###############################################################################


## ===========================================================================
## 1. 辅助函数
## ===========================================================================

#' 由端点值反解 mu(x) = a * x^b 的 (a, b)
power_from_endpoints <- function(mu_lo, mu_hi, x_lo, x_hi) {
  b <- log(mu_hi / mu_lo) / log(x_hi / x_lo)
  list(a = mu_lo / x_lo^b, b = b)
}

#' 稀疏、对角占优（因而正定）的精度矩阵
make_precision <- function(p, sparsity, w_range) {
  Om <- matrix(0, p, p)
  for (i in 1:(p - 1)) for (j in (i + 1):p) {
    if (runif(1) < sparsity) {
      v <- runif(1, w_range[1], w_range[2]) * sample(c(1, -1), 1)
      Om[i, j] <- v; Om[j, i] <- v
    }
  }
  diag(Om) <- rowSums(abs(Om)) + 0.1
  Om
}

#' 由 (a, b) 生成 n x k 的幂律曲线矩阵
power_curves <- function(EI, a, b) {
  n <- length(EI); k <- length(a)
  matrix(outer(EI, b, "^") * rep(a, each = n), n, k)
}


## ===========================================================================
## 2. 主函数
## ===========================================================================

sim_power_EI <- function(p            = 20,             # 焦点基因数
                         n            = 50,             # 样本数
                         n_bg         = 480,            # 背景基因数（撑文库总量）
                         EI_range     = c(6.95, 7.5),   # 对齐 Figure1 横轴
                         mu_lo_range  = c(0.75, 2.7),   # EI 下端的 log10(count+1)
                         slope_range  = c(0.35, 1.05),  # d mu / d EI
                         mu_cap       = 3.1,            # log10(count+1) 上限
                         snr          = NULL,           # 信噪比：控制**独立测量噪声**（网络信号不受影响）
                         net_sd       = 0.10,           # 结构化残差 sd —— 网络信号，固定，不随 snr 变
                         noise_sd     = NULL,           # 兼容旧接口：给了它就退回“网络随 snr 缩放”的老行为
                         net_sparsity = 0.12,           # 残差网络稀疏度
                         w_range      = c(0.4, 0.8),    # 精度矩阵非对角元大小
                         seed         = 42) {

  set.seed(seed)
  x_lo <- EI_range[1]; x_hi <- EI_range[2]

  ## ---- 2.1 EI 网格：有序、不等距（像真实样本） --------------------------
  u  <- sort(runif(n))
  EI <- x_lo + (x_hi - x_lo) * u
  N  <- round(10^EI)                                     # 设计文库大小

  ## ---- 2.2 焦点基因的幂律曲线 -------------------------------------------
  mu_lo <- runif(p, mu_lo_range[1], mu_lo_range[2])
  slope <- runif(p, slope_range[1], slope_range[2])
  mu_hi <- pmin(mu_lo + slope * (x_hi - x_lo), mu_cap)
  ab    <- power_from_endpoints(mu_lo, mu_hi, x_lo, x_hi)
  mu    <- power_curves(EI, ab$a, ab$b)                  # n x p

  ## ---- 2.3 残差：拆成“网络信号”和“测量噪声”两部分 -----------------------
  ##   eps_net  ~ N(0, Omega^{-1})，sd 固定为 net_sd  → 网络信号，**不随 snr 变**
  ##   eps_meas ~ N(0, sigma_g^2) 独立同分布          → 测量噪声，sigma_g = sqrt(var(mu_g)/snr)
  ## 这样 snr 越大 = 测量噪声越小 = 网络越容易被看见，MCC 才会随 snr 上升。
  ##
  ## ★ 老写法（noise_sd = sqrt(var(mu)/snr) 直接缩放 eps）是错的：
  ##   eps 是**唯一**携带网络信息的成分，被 snr 缩放就等于“snr 越大网络越弱”，
  ##   MCC 必然随 snr 下降。实测 gene/qdODE：老写法 snr=1 → +0.151, snr=10 → +0.028；
  ##   改成下面这样之后 snr=1 → +0.040, snr=10 → +0.136, snr=100 → +0.140。
  ## 传 noise_sd 可退回老行为做对照。
  Omega <- make_precision(p, net_sparsity, w_range)
  R     <- chol(solve(Omega))                            # R'R = Sigma
  eps   <- matrix(rnorm(n * p), n, p) %*% R              # cov = Sigma
  eps   <- sweep(eps, 2, colMeans(eps), "-")
  eps   <- sweep(eps, 2, apply(eps, 2, sd), "/")         # 各列 sd 归一

  if (!is.null(noise_sd)) {                              # 旧行为（趋势是反的）
    sd_g    <- if (!is.null(snr)) sqrt(apply(mu, 2, var) / snr) else rep(noise_sd, p)
    eps     <- sweep(eps, 2, sd_g, "*")
    sigma_g <- rep(0, p)
  } else {                                               # 新行为（默认）
    sigma_g <- if (!is.null(snr)) sqrt(apply(mu, 2, var) / snr) else rep(0, p)
    eps     <- eps * net_sd +
               sweep(matrix(rnorm(n * p), n, p), 2, sigma_g, "*")
    sd_g    <- sqrt(net_sd^2 + sigma_g^2)
  }

  ## ---- 2.4 背景基因：把文库总量撑到 10^7 量级 ---------------------------
  bg_lo <- pmax(rnorm(n_bg, 3.4, 0.9), 0.05)
  bg_hi <- pmax(bg_lo + runif(n_bg, 0.7, 1.1) * (x_hi - x_lo), 0.06)
  ab_bg <- power_from_endpoints(bg_lo, bg_hi, x_lo, x_hi)
  mu_bg <- power_curves(EI, ab_bg$a, ab_bg$b)            # n x n_bg

  lam_focal <- pmax(10^mu    - 1, 1e-9)
  lam_bg    <- pmax(10^mu_bg - 1, 1e-9)
  # 标定背景总量，使 sum(lambda) ~= 10^EI（归一化共同偏移 c(EI) ~= 0）
  lam_bg    <- lam_bg * median((10^EI - rowSums(lam_focal)) / rowSums(lam_bg))

  ## ---- 2.5 比例 + multinomial 抽样 --------------------------------------
  lam_obs <- cbind(lam_focal * 10^eps,
                   lam_bg * 10^matrix(rnorm(n * n_bg, 0, median(sd_g)), n, n_bg))
  P       <- lam_obs / rowSums(lam_obs)
  counts  <- t(sapply(1:n, function(i) as.vector(rmultinom(1, N[i], P[i, ]))))

  genes   <- paste0("G", 1:p)
  samples <- paste0("S", 1:n)
  expr_all <- t(counts)
  rownames(expr_all) <- c(genes, paste0("BG", seq_len(n_bg)))
  colnames(expr_all) <- samples
  expr     <- expr_all[genes, , drop = FALSE]

  EI_real  <- log10(colSums(expr_all) + 1)               # ★ 可复算的 EI
  names(EI_real) <- samples

  ## ---- 2.6 真值曲线与幂律参数 -------------------------------------------
  c_off   <- EI - log10(rowSums(cbind(lam_focal, lam_bg)))   # 归一化共同偏移
  mu_true <- mu + c_off
  a_true <- b_true <- r2 <- numeric(p)
  lx <- log(EI)
  for (j in 1:p) {
    cf <- as.numeric(coef(lm(log(mu_true[, j]) ~ lx)))
    a_true[j] <- exp(cf[1]); b_true[j] <- cf[2]
    fit <- a_true[j] * EI^b_true[j]
    r2[j] <- 1 - sum((mu_true[, j] - fit)^2) / sum((mu_true[, j] - mean(mu_true[, j]))^2)
  }

  adj <- (abs(Omega) > 1e-12) * 1L
  diag(adj) <- 0L
  dimnames(adj) <- dimnames(Omega) <- list(genes, genes)
  colnames(mu_true) <- colnames(eps) <- genes

  Y <- t(log10(expr + 1)); rownames(Y) <- NULL

  list(p = p, n = n, noise_sd = sd_g, snr = snr,
       net_sd = net_sd, meas_sd = sigma_g,
       times     = as.numeric(EI_real),   # ← t = EI
       EI        = EI_real,
       libsize   = colSums(expr_all),
       expr      = expr,                  # 20 x 50 焦点基因计数
       expr_all  = expr_all,              # 500 x 50 全矩阵（EI 由它算出）
       Y         = Y,                     # 50 x 20 log10(count+1)
       mu_true   = mu_true,               # 50 x 20 无噪真值曲线
       eps_true  = eps,                   # 50 x 20 结构化残差
       power_par = data.frame(a = a_true, b = b_true, R2_powerlaw = r2, row.names = genes),
       Omega     = Omega,                 # 残差精度矩阵
       adj_main  = adj,                   # 残差网络真值（无向）
       c_off     = c_off)
}


## ===========================================================================
## 3. 运行
## ===========================================================================

log_con <- file("run_log.txt", open = "wt")
sink(log_con, split = TRUE)

cat(R.version.string, " | ", R.version$platform, "\n", sep = "")
cat("工作目录: ", getwd(), "\n\n", sep = "")

dat <- sim_power_EI(p = 20, n = 50, seed = 42)


## ===========================================================================
## 4. 自检
## ===========================================================================

EI <- as.numeric(dat$EI); expr <- dat$expr; Y <- dat$Y

cat("############ 自 检 ############\n")

cat("\n[1] EI 自洽性（与 power_equation_fit 的定义一致）\n")
cat("    EI == log10(colSums(expr_all)+1) : ",
    isTRUE(all.equal(EI, as.numeric(log10(colSums(dat$expr_all) + 1)))), "\n", sep = "")
cat("    严格单调递增                     : ", all(diff(EI) > 0), "\n", sep = "")
cat("    范围                             : ", sprintf("10^%.3f ~ 10^%.3f", min(EI), max(EI)),
    "   (Figure1 为 10^6.95 ~ 10^7.5)\n", sep = "")
cat("    归一化共同偏移跨度               : ", sprintf("%.4f", diff(range(dat$c_off))), "\n", sep = "")

cat("\n[2] 幂律形状\n")
cat("    无噪真曲线的幂律 R^2 最小 : ", sprintf("%.6f", min(dat$power_par$R2_powerlaw)),
    "   ← 信号本身就是严格幂律\n", sep = "")
r2o <- bo <- rep(NA_real_, ncol(Y))
for (j in 1:ncol(Y)) {
  y <- Y[, j]
  m <- try(nls(y ~ a * EI^b, start = list(a = 1e-3, b = 3),
               control = nls.control(maxiter = 1000, minFactor = 1e-12, warnOnly = TRUE)),
           silent = TRUE)
  if (!inherits(m, "try-error")) {
    pr <- predict(m)
    r2o[j] <- 1 - sum((y - pr)^2) / sum((y - mean(y))^2)
    bo[j]  <- coef(m)["b"]
  }
}
cat("    观测数据 nls(y ~ a*x^b) 收敛: ", sum(!is.na(r2o)), "/", ncol(Y), "\n", sep = "")
cat("    观测 R^2 min/中位/max       : ",
    sprintf("%.3f / %.3f / %.3f", min(r2o, na.rm = TRUE), median(r2o, na.rm = TRUE),
            max(r2o, na.rm = TRUE)), "\n", sep = "")
cat("    幂指数 b 范围               : ",
    sprintf("%.2f ~ %.2f", min(bo, na.rm = TRUE), max(bo, na.rm = TRUE)), "\n", sep = "")
cat("    注：观测 R^2 只有 ~0.4 不是缺陷 —— 散点宽度(sd=", dat$noise_sd,
    ")与趋势幅度本来就相当，Figure1 里也是这个宽度。\n", sep = "")

cat("\n[3] 计数矩阵\n")
cat("    焦点 ", nrow(expr), " x ", ncol(expr), " : ", min(expr), " ~ ", max(expr),
    "  零值 ", sum(expr == 0), "   (Figure1 纵轴 10^0 ~ 10^3)\n", sep = "")
cat("    全矩阵 ", nrow(dat$expr_all), " x ", ncol(dat$expr_all), " : 文库 ",
    sprintf("%.2e ~ %.2e", min(dat$libsize), max(dat$libsize)), "\n", sep = "")
cat("    每个基因至少一个非零 : ", all(rowSums(expr) > 0), "\n", sep = "")

cat("\n[4] 残差网络真值\n")
A <- dat$adj_main
cat("    无向边数 ", sum(A) / 2, " / ", dat$p * (dat$p - 1) / 2,
    "   密度 ", sprintf("%.1f%%", 100 * sum(A) / (dat$p * (dat$p - 1))), "\n", sep = "")
cat("    度 min/中位/max : ", min(rowSums(A)), " / ", median(rowSums(A)), " / ",
    max(rowSums(A)), "\n", sep = "")

cat("\n[5] 对照：为什么网络不能编码进曲线形状\n")
cc <- cor(dat$mu_true)
cat("    20 条幂律曲线两两相关系数最小 : ", sprintf("%.5f", min(cc[upper.tri(cc)])), "\n", sep = "")
cat("    设计矩阵条件数                : ", sprintf("%.2e", kappa(dat$mu_true, exact = TRUE)),
    "  → 不可辨识\n", sep = "")


## ===========================================================================
## 5. 输出文件
## ===========================================================================

saveRDS(dat, "sim_power_EI.rds")
write.csv(dat$expr,                "pw_expr_20g_50s.csv",      quote = FALSE)
write.csv(dat$expr_all,            "pw_expr_all_500g_50s.csv", quote = FALSE)
write.csv(data.frame(sample = names(dat$EI), libsize = as.numeric(dat$libsize),
                     EI = EI),     "pw_EI.csv", row.names = FALSE, quote = FALSE)
write.csv(round(dat$Y, 6),         "pw_Y_50s_20g.csv",         row.names = FALSE, quote = FALSE)
write.csv(round(dat$mu_true, 6),   "pw_mu_true_50s_20g.csv",   row.names = FALSE, quote = FALSE)
write.csv(round(dat$eps_true, 6),  "pw_eps_true_50s_20g.csv",  row.names = FALSE, quote = FALSE)
write.csv(round(dat$power_par, 8), "pw_power_par.csv",         quote = FALSE)
write.csv(round(dat$Omega, 6),     "pw_Omega_true.csv",        quote = FALSE)
write.csv(dat$adj_main,            "pw_adj_true.csv",          quote = FALSE)


## ===========================================================================
## 6. 诊断图：与 Figure1 同版式（横轴 EI，纵轴 log10(count+1)，均为 10^x 刻度）
## ===========================================================================

png("pw_check.png", width = 1900, height = 1200, res = 130)
op <- par(mfrow = c(4, 5), mar = c(3.0, 3.2, 2.0, 0.8), oma = c(3.0, 3.2, 0.8, 0.8),
          mgp = c(2, 0.6, 0), bg = "#fcfcfb", col.axis = "#52514e", fg = "#e4e3df")
xt <- c(7.0, 7.2, 7.4)
for (k in 1:ncol(Y)) {
  plot(EI, Y[, k], type = "n", xlab = "", ylab = "", ylim = c(-0.15, 3.4),
       axes = FALSE, main = colnames(dat$mu_true)[k], col.main = "#0b0b0b",
       font.main = 1, cex.main = 1.05)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = "#f0f7ff", border = NA)
  axis(1, at = xt, labels = parse(text = paste0("10^", xt)), col = "#e4e3df")
  axis(2, at = 0:3, labels = parse(text = paste0("10^", 0:3)), las = 1, col = "#e4e3df")
  points(EI, Y[, k], pch = 1, col = "#2a78d6", cex = 0.8, lwd = 1.1)
  lines(EI, dat$mu_true[, k], col = "#eb6834", lwd = 3)
}
mtext("Expression Index", side = 1, outer = TRUE, line = 1.2, col = "#52514e", cex = 1.05)
mtext("Individual Expression of Each Gene", side = 2, outer = TRUE, line = 1.2,
      col = "#52514e", cex = 1.05)
par(op)
dev.off()

cat("\n############ 输出文件 ############\n")
ff <- list.files(".", pattern = "^(pw_|sim_power_EI\\.rds)")
print(data.frame(file = ff, bytes = file.size(ff)))

cat("\n完成。\n")
sink(); close(log_con)


###############################################################################
# 下游接口（改成 TRUE 即可跑）
###############################################################################
if (FALSE) {
  source("base.R"); source("ode_1.R")

  ## (a) EI 用全矩阵算 —— 10^7 量级，与 Figure1 一致
  res <- power_equation_fit(dat$expr_all, n = 20, trans = log10, thread = 2)

  ## (b) 只在 20 个焦点基因上跑 —— EI 变成模块内的 10^3 量级，曲线形状不变
  res2    <- power_equation_fit(dat$expr, n = 20, trans = log10, thread = 2)
  res_ode <- qdODE_parallel(result = res2, reduction = TRUE, thread = 2)

  ## (c) 网络评估：真值是 dat$adj_main（无向）
  ##     注意 qdODE 的支撑来自 get_interaction 在原始计数上的 LASSO，
  ##     所以扣掉幂律趋势后的残差 (dat$Y - dat$mu_true) 才是它真正在看的信号。
}

###############################################################################
# 可调旋钮
#   snr          信噪比：控制**测量噪声** sigma_g = sqrt(var(mu_g)/snr)。网络信号不受影响
#   net_sd       网络信号强度（结构化残差 sd），固定不随 snr 变。调大 → 网络更好恢复
#   noise_sd     兼容老接口：给了它就退回“网络随 snr 一起缩小”的旧行为（MCC 趋势会反）
#   net_sparsity / w_range   残差网络的稀疏度与强度。调大 → 网络更好恢复
#   mu_lo_range / slope_range / mu_cap   表达量范围与斜率
#   EI_range     换成 c(3.2, 4.8) 并设 n_bg = 0，就是 gene/ 模块内那套量级
###############################################################################
