## =====================================================================
## HACK + MED 事件发生时间序列: 聚集程度诊断指标与图形
##
## 数据: PRC1_HACK_MED_daily_interval.csv 的 230 个事件日(同一天只算一次)
## 观测窗: [0, L), L = 4410 天 (= 630 个 7 天窗口, 与上一轮一致)
##
## 两个零模型(置换检验, B = 999, 均为整数天无放回抽样, 保持"每天最多 1 事件"):
##   H0 齐次随机     : 事件日在窗内均匀随机  -> 检验"是否非均匀"
##   H1 保留趋势随机 : 按高斯核平滑强度(带宽 180 天)加权抽样
##                     -> 检验"是否超出长期趋势的短程聚集"
## 只有在 H1 下仍显著, 才能把聚集归因于短程相依(自激发), 而非事件率上升。
##
## 点 Source 直接出图; 换图改顶部 FIG
## =====================================================================

DIR <- "E:/paper/rdata/hackexp"
f   <- file.path(DIR, "PRC1_HACK_MED_daily_interval.csv")
d   <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
d$Date <- as.Date(d$Date)
Date0  <- min(d$Date)
x      <- sort(as.numeric(d$Date - Date0))     # 0..4404, 整数天
N      <- length(x)
L      <- 4410L
it     <- diff(x)                              # 229 个事件间隔
stopifnot(N == 230L, max(x) < L)

## ------------------------------------------------------------ 设定
DELTA <- c(7, 10, 14, 21, 30, 49, 70, 105, 147, 210, 315, 441, 630, 882)
stopifnot(all(L %% DELTA == 0))                # 全部整除 L, 分箱无残窗
SMAX  <- 365L
sgrid <- 1:SMAX
SG    <- c(1:30, seq(35, 200, by = 5), 250, 300, 400, 500)
B     <- 999L
H     <- 180                                   # 平滑强度带宽(天)
set.seed(20261108)

## ------------------------------------------------------- 指标函数
fano_curve <- function(xx) {
  sapply(DELTA, function(D) {
    cc <- tabulate(floor(xx / D) + 1L, nbins = L / D)
    var(cc) / mean(cc)
  })
}
inter_stats <- function(xx) {
  ii <- diff(sort(xx)); cv <- sd(ii) / mean(ii)
  c(cv = cv, burst = (cv - 1) / (cv + 1),
    p_le7 = mean(ii <= 7), p_le30 = mean(ii <= 30))
}
K_curve <- function(xx) {                      # 一维 Ripley K(s), 有序对
  pd <- as.integer(as.vector(dist(xx)))       # 事件时间均为整数天
  2 * cumsum(tabulate(pd, nbins = SMAX)) * L / (N * N)
}
surv_curve <- function(xx) {
  ii <- diff(sort(xx))
  sapply(SG, function(t) mean(ii > t))
}
ks_stat <- function(xx) max(abs((1:N) / N - sort(xx) / L))
KL_excess <- function(k) sum(k[1:60] - 2 * (1:60))

fo <- fano_curve(x); io <- inter_stats(x); ko <- K_curve(x)
so <- surv_curve(x); kso <- ks_stat(x); klo <- KL_excess(ko)

## ------------------------------------- 保留趋势的零模型: 平滑强度
M   <- outer(0:(L - 1L), x, "-") / H
lam <- rowSums(dnorm(M))
lam <- lam / mean(lam) + 0.02                 # 加背景, 避免权重为 0

sim_hom   <- function() sort(sample.int(L, N) - 1L)
sim_trend <- function() {                     # Gumbel-top-k = 带权无放回抽样
  key <- log(pmax(runif(L), 1e-12)) / lam
  sort(order(key, decreasing = TRUE)[1:N] - 1L)
}

## ------------------------------------------------------- 置换检验
FH <- FT <- matrix(0, B, length(DELTA))
IH <- IT <- matrix(0, B, 4)
KH <- KT <- matrix(0, B, SMAX)
SH <- ST <- matrix(0, B, length(SG))
KSH <- KST <- KLH <- KLT <- numeric(B)

for (b in seq_len(B)) {
  xh <- sim_hom(); xt <- sim_trend()
  FH[b, ] <- fano_curve(xh);   FT[b, ] <- fano_curve(xt)
  IH[b, ] <- inter_stats(xh);  IT[b, ] <- inter_stats(xt)
  kh <- K_curve(xh); kt <- K_curve(xt)
  KH[b, ] <- kh; KT[b, ] <- kt
  KLH[b] <- KL_excess(kh); KLT[b] <- KL_excess(kt)
  SH[b, ] <- surv_curve(xh);  ST[b, ] <- surv_curve(xt)
  KSH[b] <- ks_stat(xh);      KST[b] <- ks_stat(xt)
}

p_up  <- function(obs, null) (1 + sum(null >= obs)) / (length(null) + 1)
qband <- function(M) apply(M, 2, quantile, probs = c(0.025, 0.975))

## ------------------------------------------------------- 指标汇总表
make_row <- function(name, obs, nh, nt) {
  data.frame(metric = name,
             observed = as.numeric(obs),
             hom_mean = mean(nh), hom_lo = as.numeric(quantile(nh, .025)),
             hom_hi = as.numeric(quantile(nh, .975)), hom_sd = sd(nh),
             p_hom = p_up(obs, nh),
             trend_mean = mean(nt), trend_lo = as.numeric(quantile(nt, .025)),
             trend_hi = as.numeric(quantile(nt, .975)), trend_sd = sd(nt),
             p_trend = p_up(obs, nt))
}
k30 <- which(DELTA == 30); k105 <- which(DELTA == 105)
MET <- rbind(
  make_row("CV of intervals",        io[1], IH[, 1], IT[, 1]),
  make_row("Burstiness B",           io[2], IH[, 2], IT[, 2]),
  make_row("P(interval <= 7d)",      io[3], IH[, 3], IT[, 3]),
  make_row("P(interval <= 30d)",     io[4], IH[, 4], IT[, 4]),
  make_row("Fano (Delta = 7d)",      fo[1], FH[, 1], FT[, 1]),
  make_row("Fano (Delta = 30d)",     fo[k30], FH[, k30], FT[, k30]),
  make_row("Fano (Delta = 105d)",    fo[k105], FH[, k105], FT[, k105]),
  make_row("K(s)-2s area, s<=60d",   klo, KLH, KLT),
  make_row("KS deviation from uniform", kso, KSH, KST)
)
rownames(MET) <- NULL

FANOTAB <- data.frame(
  Delta = DELTA, observed = fo,
  hom_mean = colMeans(FH), hom_lo = qband(FH)[1, ], hom_hi = qband(FH)[2, ],
  p_hom = sapply(seq_along(DELTA), function(j) p_up(fo[j], FH[, j])),
  trend_mean = colMeans(FT), trend_lo = qband(FT)[1, ], trend_hi = qband(FT)[2, ],
  p_trend = sapply(seq_along(DELTA), function(j) p_up(fo[j], FT[, j]))
)
rownames(FANOTAB) <- NULL

cat(sprintf("事件日 N = %d, 观测窗 %d 天 (%s ~ %s), 置换 B = %d\n",
            N, L, Date0, Date0 + L - 1, B))
cat(sprintf("间隔: 均值 %.2f 天, 中位 %.1f 天, CV %.3f, Burstiness %.3f, P(<=7d) %.3f\n",
            mean(it), median(it), io[1], io[2], io[3]))
cat("\n指标 (p_hom / p_trend 为单侧置换 p 值, 越小越偏离零模型):\n")
print(MET[, c("metric", "observed", "hom_mean", "trend_mean", "p_hom", "p_trend")],
      row.names = FALSE, digits = 4)

## ------------------------------------------------------------- 配色
COL_OBS   <- "#2C7FB8"
COL_HOM   <- "grey55"
COL_TREND <- "#E39B2D"
COL_REF   <- "#D94841"

open_fig <- function(width = 10.5, height = 5.6) {
  cur <- dev.cur(); dl <- dev.list()
  if (!is.null(dl) && names(cur) %in% c("png", "jpeg", "tiff", "pdf", "svg"))
    dev.off(cur)
  if (is.null(dev.list())) dev.new(width = width, height = height)
  par(mfrow = c(1, 1), mar = c(4.3, 4.8, 2.3, 1.1),
      cex.lab = 1.12, cex.axis = 0.92, cex.main = 1.08, bg = "white")
}

## 日期刻度(数值天坐标 -> 日期标签), 避免 Date 对象在 segments 中的类型问题
date_axis <- function() {
  at <- seq(0, L, by = 730)
  axis(1, at = at, labels = format(Date0 + at, "%Y-%m"))
}

## ---------- 图 1: 事件发生时刻(地毯图) + 平滑强度(显示长期趋势)
draw_events <- function() {
  open_fig(width = 11, height = 4.8)
  rate <- lam / mean(lam) * (N / L) * 365
  ymax <- max(rate) * 1.18
  plot(0:(L - 1L), rate, type = "n", xaxt = "n", ylim = c(0, ymax),
       xlab = "Date", ylab = "Smoothed event rate (events / year)",
       main = "HACK + MED: event occurrence times and long-run rate trend")
  date_axis()
  segments(x, 0, x, ymax * 0.10, col = "grey30", lwd = 0.8)
  lines(0:(L - 1L), rate, col = COL_TREND, lwd = 2.2)
  abline(h = (N / L) * 365, col = COL_REF, lty = 2, lwd = 1.5)
  legend("topleft", bty = "n", cex = 0.88,
         legend = c("event day (rug)", "smoothed rate (h = 180 d)", "homogeneous mean"),
         col = c("grey30", COL_TREND, COL_REF), lty = c(1, 1, 2), lwd = c(1, 2.2, 1.5))
  box(); invisible(NULL)
}

## ---------- 图 2: 累积计数曲线 vs 齐次期望
draw_cum <- function() {
  open_fig(width = 9.6, height = 5.4)
  plot(c(0, L), c(0, N), type = "n", xaxt = "n",
       xlab = "Date", ylab = "Cumulative event days N(t)",
       main = "Cumulative count vs homogeneous Poisson expectation")
  date_axis()
  lines(c(0, x, L), c(0, 1:N, N), type = "s", col = COL_OBS, lwd = 2)
  lines(c(0, L), c(0, N), col = COL_REF, lty = 2, lwd = 1.8)
  legend("topleft", bty = "n", cex = 0.9,
         legend = c("observed N(t)", "homogeneous expectation"),
         col = c(COL_OBS, COL_REF), lty = c(1, 2), lwd = c(2, 1.8))
  text(L * 0.98, N * 0.06, pos = 2, cex = 0.9,
       labels = sprintf("KS D = %.3f, p(hom) = %.3f", kso, p_up(kso, KSH)))
  box(); invisible(NULL)
}

## ---------- 图 3: 多尺度 Fano 因子(核心)
draw_fano <- function() {
  open_fig()
  eh <- qband(FH); et <- qband(FT)
  yr <- range(c(fo, as.vector(eh), as.vector(et), 1))
  plot(DELTA, fo, type = "n", log = "x", xaxt = "n",
       ylim = c(max(0, yr[1] * 0.94), yr[2] * 1.08),
       xlab = expression("window length "~Delta~" (days, log scale)"),
       ylab = "Fano factor  Var(N) / E(N)",
       main = "Multi-scale clustering: Fano factor of window counts")
  axis(1, at = DELTA, labels = as.character(DELTA), cex.axis = 0.85)
  polygon(c(DELTA, rev(DELTA)), c(eh[1, ], rev(eh[2, ])),
          col = adjustcolor(COL_HOM, 0.22), border = NA)
  polygon(c(DELTA, rev(DELTA)), c(et[1, ], rev(et[2, ])),
          col = adjustcolor(COL_TREND, 0.26), border = NA)
  abline(h = 1, col = COL_REF, lty = 2, lwd = 1.5)
  lines(DELTA, fo, type = "b", pch = 16, col = COL_OBS, lwd = 2.2)
  legend("topleft", bty = "n", cex = 0.86,
         legend = c("observed", "95% band: homogeneous null",
                    "95% band: trend-preserving null", "Poisson (Fano = 1)"),
         col = c(COL_OBS, adjustcolor(COL_HOM, 0.6),
                 adjustcolor(COL_TREND, 0.85), COL_REF),
         lty = c(1, 1, 1, 2), lwd = c(2.2, 8, 8, 1.5), pch = c(16, NA, NA, NA))
  box(); invisible(NULL)
}

## ---------- 图 4: 间隔生存函数 vs 指数(Poisson)基准
draw_inter <- function() {
  open_fig()
  eh <- qband(SH); et <- qband(ST)
  lo <- 1 / (length(it) + 1)
  plot(SG, pmax(so, lo), type = "n", log = "y", ylim = c(lo, 1.05),
       xlab = "Inter-event interval  t  (days)",
       ylab = "P(interval > t)   [log scale]",
       main = "Inter-event interval survival: short gaps and long tail")
  polygon(c(SG, rev(SG)), c(pmax(eh[1, ], lo), rev(pmax(eh[2, ], lo))),
          col = adjustcolor(COL_HOM, 0.22), border = NA)
  polygon(c(SG, rev(SG)), c(pmax(et[1, ], lo), rev(pmax(et[2, ], lo))),
          col = adjustcolor(COL_TREND, 0.26), border = NA)
  lines(SG, exp(-SG / mean(it)), col = COL_REF, lty = 2, lwd = 1.8)
  lines(SG, pmax(so, lo), type = "b", pch = 16, cex = 0.65,
        col = COL_OBS, lwd = 1.8)
  legend("topright", bty = "n", cex = 0.84,
         legend = c("observed", "95% band: homogeneous", "95% band: trend-preserving",
                    "exponential (Poisson)"),
         col = c(COL_OBS, adjustcolor(COL_HOM, 0.6),
                 adjustcolor(COL_TREND, 0.85), COL_REF),
         lty = c(1, 1, 1, 2), lwd = c(1.8, 8, 8, 1.8), pch = c(16, NA, NA, NA))
  box(); invisible(NULL)
}

## ---------- 图 5: 成对时间差(一维 Ripley K)超出量
draw_K <- function() {
  open_fig()
  KHd <- sweep(KH, 2, 2 * sgrid, "-")
  KTd <- sweep(KT, 2, 2 * sgrid, "-")
  eo  <- ko - 2 * sgrid
  eh  <- qband(KHd); et <- qband(KTd)
  yr  <- range(c(eo, as.vector(eh), as.vector(et), 0))
  plot(sgrid, eo, type = "n", ylim = yr * c(1.05, 1.05),
       xlab = "Lag  s  (days)", ylab = "K(s) - 2s",
       main = "Pairwise time-difference excess (1D Ripley K)")
  polygon(c(sgrid, rev(sgrid)), c(eh[1, ], rev(eh[2, ])),
          col = adjustcolor(COL_HOM, 0.22), border = NA)
  polygon(c(sgrid, rev(sgrid)), c(et[1, ], rev(et[2, ])),
          col = adjustcolor(COL_TREND, 0.26), border = NA)
  abline(h = 0, col = COL_REF, lty = 2, lwd = 1.5)
  lines(sgrid, eo, col = COL_OBS, lwd = 1.8)
  legend("topleft", bty = "n", cex = 0.86,
         legend = c("observed", "95% band: homogeneous",
                    "95% band: trend-preserving", "Poisson (K(s)=2s)"),
         col = c(COL_OBS, adjustcolor(COL_HOM, 0.6),
                 adjustcolor(COL_TREND, 0.85), COL_REF),
         lty = c(1, 1, 1, 2), lwd = c(1.8, 8, 8, 1.5))
  box(); invisible(NULL)
}

## ---------- 图 6: 指标汇总(相对保留趋势零模型的标准化偏离)
draw_summary <- function() {
  open_fig(width = 9.6, height = 5.8)
  par(mar = c(4.3, 9.6, 2.3, 1.1))
  z  <- (MET$observed - MET$trend_mean) / MET$trend_sd
  nn <- nrow(MET); yy <- nn:1
  xl <- range(c(-2.4, z * 1.18))
  plot(xl, c(0.5, nn + 0.5), type = "n", yaxt = "n", ylab = "",
       xlab = "Standardized deviation from trend-preserving null (z)",
       main = "Clustering diagnostics: observed vs null models")
  segments(-1.96, yy, 1.96, yy, col = "grey80", lwd = 9)
  abline(v = 0, col = COL_REF, lty = 2, lwd = 1.4)
  points(z, yy, pch = 16, col = COL_OBS, cex = 1.35)
  text(z, yy, labels = sprintf("  %.1f", z), pos = 4, cex = 0.8)
  axis(2, at = yy, labels = MET$metric, las = 1, cex.axis = 0.82)
  legend("bottomright", bty = "n", cex = 0.82,
         legend = c("observed z", "null 95% range (|z| < 1.96)"),
         col = c(COL_OBS, "grey70"), pch = c(16, NA), lty = c(0, 1),
         lwd = c(1, 6))
  box(); invisible(NULL)
}

if (!exists("FIG")) FIG <- "fano"
FIGS <- list(events = draw_events, cum = draw_cum, fano = draw_fano,
             inter = draw_inter, K = draw_K, summary = draw_summary,
             none = function() invisible(NULL))
FIGS[[FIG]]()
