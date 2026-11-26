## =====================================================================
## 损失标记 Z 与事件发生频率 / 间隔序列的相关性: 设计与置换检验
##
## 数据: PRC1_HACK_MED_event_marks.csv (230 个事件日)
##       Z_j = 9.59 + 0.57 ln(Total.Records)  (单调变换, 秩相关等价于 ln 记录数)
##       gaps[j] = t_{j+1} - t_j  (事件 j 之后的间隔, 229 个)
##
## 待检假设 H0: 给定事件时间与长期趋势, 损失标记与后续发生过程独立。
## 主统计量: rho_s(Z_j, log gaps[j])  (预设方向: 损失越大, 后续间隔越短 -> rho < 0)
##
## 混杂控制: 块置换(把 Z 按事件顺序切成连续块, 块内重排)保留年际趋势,
##           只破坏"哪个损失配哪个间隔"; 另用全局置换/循环移位做稳健性。
## 补充: 偏相关(控制时间 t 与前一间隔)、lag 1..10 交叉相关、未来 30/90/180 天窗口计数。
## 点 Source 直接出图; 换图改顶部 FIG
## =====================================================================

DIR <- "E:/paper/rdata/hackexp"
d  <- read.csv(file.path(DIR, "PRC1_HACK_MED_event_marks.csv"),
               stringsAsFactors = FALSE, check.names = FALSE)
d$Date <- as.Date(d$Date)
etimes <- as.numeric(d$Date - min(d$Date))     # 0..4404, 230 个事件日
Z      <- d$Z
NE     <- length(etimes)
gaps   <- diff(etimes)                          # 229 个后续间隔
z      <- Z[-NE]                                # 与 gaps 对齐的 229 个标记
lg     <- log(gaps)
L      <- 4410L
KMAX   <- 10L
WS     <- c(30, 90, 180)
stopifnot(NE == 230L, all(gaps > 0))

rho_of <- function(a, b) suppressWarnings(cor(a, b, method = "spearman"))

## ------------------------------------------- lag k: Z_j vs 第 k 个后续间隔
Y_LAG <- lapply(1:KMAX, function(k) lg[k:(NE - 1L)])
Z_LAG <- lapply(1:KMAX, function(k) Z[1:(NE - k)])

## ------------------------------------------- 未来 W 天窗口内的事件日数
## 只保留窗口完整落在观测期内的事件, 避免右端截断造成的假低计数
WINS <- lapply(WS, function(W) {
  ok  <- which(etimes + W <= max(etimes))
  cnt <- vapply(ok, function(j) sum(etimes > etimes[j] & etimes <= etimes[j] + W),
                numeric(1))
  list(ok = ok, n = cnt)
})
names(WINS) <- paste0("W", WS)

## ------------------------------------------- 偏相关: 控制 t_j 与前一间隔
jj  <- 2:(NE - 1L)
X   <- cbind(1, splines::ns(etimes[jj], df = 4), log(gaps[jj - 1L]))
qrx <- qr(X)
ry  <- qr.resid(qrx, lg[jj])
partial_rho <- function(zfull) rho_of(qr.resid(qrx, zfull[jj]), ry)

obs_main <- rho_of(z, lg)
obs_lag  <- sapply(1:KMAX, function(k) rho_of(Z_LAG[[k]], Y_LAG[[k]]))
obs_win  <- sapply(seq_along(WS), function(i) rho_of(Z[WINS[[i]]$ok], WINS[[i]]$n))
obs_part <- partial_rho(Z)

## ------------------------------------------------------- 置换零模型
block_perm <- function(v, m) {
  bid <- (seq_along(v) - 1L) %/% m
  unname(unlist(lapply(split(v, bid), sample), use.names = FALSE))
}
shift_perm <- function(v) {
  k <- sample.int(length(v) - 1L, 1L)
  v[c((k + 1L):length(v), 1L:k)]
}

p_lower <- function(obs, null) (1 + sum(null <= obs)) / (length(null) + 1)
p_upper <- function(obs, null) (1 + sum(null >= obs)) / (length(null) + 1)
p_two   <- function(obs, null) {
  m0 <- mean(null)
  (1 + sum(abs(null - m0) >= abs(obs - m0))) / (length(null) + 1)
}
qband <- function(v) as.numeric(quantile(v, c(0.025, 0.975)))

B <- 1999L
set.seed(20261108)
NC <- 1L + KMAX + length(WS) + 1L
NM <- matrix(0, B, NC)
for (b in seq_len(B)) {
  zp <- block_perm(Z, 20)                      # 主零模型: 块长 20 个事件(约 1 年)
  NM[b, 1L] <- rho_of(zp[-NE], lg)
  for (k in 1:KMAX) NM[b, 1L + k] <- rho_of(zp[1:(NE - k)], Y_LAG[[k]])
  for (i in seq_along(WS))
    NM[b, 1L + KMAX + i] <- rho_of(zp[WINS[[i]]$ok], WINS[[i]]$n)
  NM[b, NC] <- partial_rho(zp)
}

null_main <- NM[, 1L]
null_part <- NM[, NC]

## ------------------------------------------- 稳健性(仅主统计量)
rob <- function(gen, B2 = 999L) replicate(B2, rho_of(gen(Z)[-NE], lg))
R10 <- rob(function(v) block_perm(v, 10))
R40 <- rob(function(v) block_perm(v, 40))
RG  <- rob(sample)
RS  <- rob(shift_perm)

## ------------------------------------------------------- 指标汇总表
make_test <- function(name, obs, null, dir) {
  q <- qband(null)
  data.frame(test = name, observed = as.numeric(obs), n = NA_integer_,
             null_mean = mean(null), null_lo = q[1], null_hi = q[2],
             p_directional = if (dir == "low") p_lower(obs, null) else p_upper(obs, null),
             p_two_sided = p_two(obs, null))
}
MET <- rbind(
  make_test("rho(Z, log next gap) [MAIN]", obs_main, null_main, "low"),
  make_test("partial rho | t, prev gap", obs_part, null_part, "low"),
  do.call(rbind, lapply(1:KMAX, function(k)
    make_test(sprintf("lag %d: rho(Z_j, log gap_{j+%d})", k, k - 1),
              obs_lag[k], NM[, 1L + k], "two"))),
  do.call(rbind, lapply(seq_along(WS), function(i)
    make_test(sprintf("N_%dd vs Z", WS[i]), obs_win[i], NM[, 1L + KMAX + i], "up")))
)
MET$n <- c(NE - 1L, length(jj),
           sapply(1:KMAX, function(k) NE - k),
           sapply(WINS, function(w) length(w$ok)))
rownames(MET) <- NULL
## lag/窗口用双尾 p 覆盖方向列
idx_lag <- 3:(2 + KMAX)
MET$p_directional[idx_lag] <- MET$p_two_sided[idx_lag]

ROB <- do.call(rbind, lapply(list(
  list("block m=20 (main)", null_main), list("block m=10", R10),
  list("block m=40", R40), list("global shuffle", RG),
  list("circular shift", RS)), function(o) {
    q <- qband(o[[2]])
    data.frame(null_model = o[[1]], observed = obs_main,
               null_mean = mean(o[[2]]), null_lo = q[1], null_hi = q[2],
               p_lower = p_lower(obs_main, o[[2]]),
               p_two = p_two(obs_main, o[[2]]))
  }))
rownames(ROB) <- NULL

cat(sprintf("事件日 %d, 后续间隔 %d; 主检验 rho_s = %.4f (块置换 B=%d)\n",
            NE, length(gaps), obs_main, B))
cat(sprintf("  H0 均值 %.4f, 95%% 区间 [%.4f, %.4f], 单尾 p = %.4f, 双尾 p = %.4f\n",
            mean(null_main), qband(null_main)[1], qband(null_main)[2],
            p_lower(obs_main, null_main), p_two(obs_main, null_main)))
cat(sprintf("偏相关(控制 t 与前一间隔) = %.4f, p = %.4f\n",
            obs_part, p_lower(obs_part, null_part)))
cat("\n各检验:\n")
print(MET, row.names = FALSE, digits = 4)
cat("\n稳健性:\n")
print(ROB, row.names = FALSE, digits = 4)

## ------------------------------------------------------------- 配色
COL     <- "#2C7FB8"
COL_REF <- "#D94841"
COL_ALT <- "#E39B2D"

open_fig <- function(width = 10.5, height = 5.6) {
  cur <- dev.cur(); dl <- dev.list()
  if (!is.null(dl) && names(cur) %in% c("png", "jpeg", "tiff", "pdf", "svg"))
    dev.off(cur)
  if (is.null(dev.list())) dev.new(width = width, height = height)
  par(mfrow = c(1, 1), mar = c(4.3, 4.8, 2.3, 1.1),
      cex.lab = 1.08, cex.axis = 0.92, cex.main = 1.02, bg = "white")
}

grp3 <- cut(z, breaks = quantile(z, c(0, 1/3, 2/3, 1)),
            include.lowest = TRUE, labels = c("Low Z", "Mid Z", "High Z"))
kw_p <- kruskal.test(lg ~ grp3)$p.value

## ---------- 图 1: 损失 vs 后续间隔散点
draw_scatter <- function() {
  open_fig()
  plot(z, lg, pch = 16, cex = 0.8, col = adjustcolor(COL, 0.75),
       xlab = expression("Transformed loss  " * Z[j]),
       ylab = "log(next inter-event gap, days)",
       main = "Loss mark vs subsequent waiting time")
  o <- order(z); fit <- loess(lg ~ z, span = 0.6)
  lines(z[o], fitted(fit)[o], col = COL_REF, lwd = 2)
  legend("topright", bty = "n", cex = 0.82,
         legend = c(sprintf("Spearman rho = %.3f", obs_main),
                    "H1: rho < 0 (large loss -> shorter gap)",
                    sprintf("p(one-sided, H1) = %.3f", p_lower(obs_main, null_main)),
                    sprintf("p(two-sided) = %.3f", p_two(obs_main, null_main)),
                    "loess fit"),
         col = c("white", "white", "white", "white", COL_REF),
         lty = c(0, 0, 0, 0, 1), lwd = c(1, 1, 1, 1, 2))
  box(); invisible(NULL)
}

## ---------- 图 2: Z 三分位组的后续间隔箱线图
draw_groups <- function() {
  open_fig(width = 8.6, height = 5.6)
  boxplot(lg ~ grp3, col = c("#9EC5E8", "#5B9BD5", COL), border = "grey35",
          xlab = "Transformed loss group (tertiles)",
          ylab = "log(next inter-event gap, days)",
          main = "Waiting time after low / mid / high loss events")
  med <- tapply(lg, grp3, median)
  text(1:3, med, labels = sprintf("median %.1f d", exp(med)), pos = 3, cex = 0.85)
  legend("topright", bty = "n", cex = 0.88,
         legend = sprintf("Kruskal-Wallis p = %.3f", kw_p))
  box(); invisible(NULL)
}

## ---------- 图 3: 高/低损失组的间隔生存曲线
draw_surv <- function() {
  open_fig()
  sgrid <- 1:200
  hi <- gaps[grp3 == "High Z"]; lo <- gaps[grp3 == "Low Z"]
  S_hi <- sapply(sgrid, function(s) mean(hi > s))
  S_lo <- sapply(sgrid, function(s) mean(lo > s))
  plot(sgrid, S_lo, type = "s", col = "#9EC5E8", lwd = 2.2,
       ylim = c(0, 1), xlab = "Waiting time s (days)",
       ylab = "P(next gap > s)",
       main = "Survival of waiting time by loss group")
  lines(sgrid, S_hi, type = "s", col = COL, lwd = 2.2)
  abline(h = 0.5, col = "grey60", lty = 3)
  legend("topright", bty = "n", cex = 0.9,
         legend = c(sprintf("Low Z  (n=%d, median %.0f d)", length(lo), median(lo)),
                    sprintf("High Z (n=%d, median %.0f d)", length(hi), median(hi))),
         col = c("#9EC5E8", COL), lty = 1, lwd = 2.2)
  box(); invisible(NULL)
}

## ---------- 图 4: 滞后交叉相关 + 块置换 95% 带
draw_lags <- function() {
  open_fig()
  band <- t(apply(NM[, 1L + 1:KMAX, drop = FALSE], 2, qband))
  yr <- range(c(obs_lag, as.vector(band), 0))
  pad <- 0.12 * diff(yr)
  plot(1:KMAX, obs_lag, type = "n", ylim = c(yr[1] - pad, yr[2] + pad),
       xlab = "Lag k  (k=1 is the immediately following gap)",
       ylab = "Spearman rho  (Z_j vs k-th subsequent gap)",
       main = "Cross-correlation between loss mark and later gaps")
  polygon(c(1:KMAX, KMAX:1), c(band[, 1], rev(band[, 2])),
          col = adjustcolor(COL_ALT, 0.30), border = NA)
  abline(h = 0, col = COL_REF, lty = 2, lwd = 1.4)
  lines(1:KMAX, obs_lag, type = "b", pch = 16, col = COL, lwd = 2)
  sig <- MET$p_two_sided[idx_lag] < 0.05
  points((1:KMAX)[sig], obs_lag[sig], pch = 21, cex = 1.8,
         bg = COL_REF, col = COL_REF)
  legend("topright", bty = "n", cex = 0.86,
         legend = c("observed rho", "95% block-permutation band",
                    "p(two-sided) < 0.05 (uncorrected)"),
         col = c(COL, adjustcolor(COL_ALT, 0.8), COL_REF),
         lty = c(1, 1, 0), lwd = c(2, 8, 1), pch = c(16, NA, 21))
  box(); invisible(NULL)
}

## ---------- 图 5: 主检验的置换零分布
draw_null <- function() {
  open_fig(width = 8.8, height = 5.6)
  h <- hist(null_main, breaks = 45, plot = FALSE)
  plot(h, col = "grey78", border = "white", freq = FALSE,
       xlab = "Spearman rho under block-permutation H0",
       ylab = "Density",
       main = "Null distribution of the main correlation statistic",
       xlim = range(c(null_main, obs_main)) * c(1.05, 1.05))
  q <- qband(null_main)
  abline(v = q, col = COL_ALT, lty = 2, lwd = 1.6)
  abline(v = obs_main, col = COL_REF, lwd = 2.2)
  legend("topright", bty = "n", cex = 0.86,
         legend = c(sprintf("observed rho = %.3f", obs_main),
                    sprintf("H0 95%% band [%.3f, %.3f]", q[1], q[2]),
                    sprintf("p(one-sided) = %.3f", p_lower(obs_main, null_main))),
         col = c(COL_REF, COL_ALT, "white"), lty = c(1, 2, 0), lwd = c(2.2, 1.6, 1))
  box(); invisible(NULL)
}

## ---------- 图 6: 未来窗口事件率 vs 损失分位组
draw_windows <- function() {
  open_fig()
  grp5 <- cut(Z, breaks = quantile(Z, probs = seq(0, 1, length.out = 6)),
              include.lowest = TRUE, labels = FALSE)
  xc <- seq_along(levels(factor(grp5)))
  rates <- sapply(seq_along(WS), function(i) {
    w <- WINS[[i]]
    tapply(w$n / WS[i], factor(grp5[w$ok], levels = 1:5), mean)
  })
  yr <- range(rates)
  plot(xc, rates[, 1], type = "n", ylim = c(0, yr[2] * 1.18),
       xaxt = "n", xlab = "Transformed loss Z (quintile group, low -> high)",
       ylab = "Subsequent event rate (event days / day)",
       main = "Future event frequency by loss mark")
  axis(1, at = xc, labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))
  cols <- c("#9EC5E8", COL, COL_ALT)
  for (i in seq_along(WS)) lines(xc, rates[, i], type = "b", pch = 16,
                                 col = cols[i], lwd = 2)
  legend("topleft", bty = "n", cex = 0.85,
         title = "p_dir: one-sided p for H1 rho > 0",
         legend = sprintf("W = %d d   rho = %.3f   p_dir = %.3f",
                          WS, obs_win, MET$p_directional[2 + KMAX + seq_along(WS)]),
         col = cols, lty = 1, lwd = 2, pch = 16)
  box(); invisible(NULL)
}

if (!exists("FIG")) FIG <- "scatter"
FIGS <- list(scatter = draw_scatter, groups = draw_groups, surv = draw_surv,
             lags = draw_lags, null = draw_null, windows = draw_windows,
             none = function() invisible(NULL))
FIGS[[FIG]]()
