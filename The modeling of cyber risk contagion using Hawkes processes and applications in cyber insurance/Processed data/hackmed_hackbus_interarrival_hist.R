## =====================================================================
## 仿论文 Figure 4 —— HACK-MED / HACK-BUS 事件到达间隔(interarrival time)直方图
##
## 数据: hackmed_daily.csv / hackbus_daily.csv
##       (2010-01-01 ~ 2018-12-30 窗口, 同一天合并为 1 个事件日, 记录数累加)
##       日序列 inter 列首行 0 是起点, 真实间隔 = inter[inter > 0]
##
## 样式: 沿用原论文 hist_med.pdf / hist_bus.pdf 的 ggplot2 做法
##       (steelblue 柱 + 红色粗虚线均值 + 绿色点线中位数 + 彩色数值标注)
## 修正: 原作图代码 (prc0823.R) 用 binwidth=5 却把横轴标成
##       "Interarrival times (7 days)"; 本脚本统一为 7 天分箱, 横轴明确写 days。
##
## 点 Source 直接出图; 只改顶部 FIG:
##   "both" 双面板全范围(最贴近原图) | "zoom" 双面板放大 0~70 天(推荐阅读)
##   "bus" / "med" 单独一张 | "exp" 叠加 i.i.d. 指数间隔的期望箱计数
## 本脚本不写任何文件; 导出 PNG/CSV 用 hackmed_hackbus_interarrival_export.R
## =====================================================================
rm(list=ls())
DIR <- "E:/paper/rdata/hackexp"
if (!exists("FIG", inherits = TRUE)) FIG <- "both"

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

BINW    <- 7L
ZOOM_TO <- 70L
COL_BAR  <- "steelblue"
COL_MEAN <- "#D94841"
COL_MED  <- "#238B45"
COL_EXP  <- "#6A51A3"

## --------------------------------------------------------------- 1. 读数据
read_gaps <- function(file, label) {
  d <- read.csv(file.path(DIR, file), stringsAsFactors = FALSE,
                check.names = FALSE)
  stopifnot(all(c("Date", "inter", "Total.Records") %in% names(d)))
  d$Date <- as.Date(d$Date)
  stopifnot(!anyNA(d$inter), d$inter[1] == 0, all(d$inter >= 0),
            !anyDuplicated(d$Date),
            all(as.integer(diff(d$Date)) == d$inter[-1]))
  g <- as.integer(d$inter[d$inter > 0])
  list(label = label, gaps = g, daily = d,
       stats = list(n = length(g), mean = mean(g), median = median(g),
                    min = min(g), max = max(g), sd = sd(g),
                    cv = sd(g) / mean(g),
                    p_le7 = mean(g <= 7), p_le30 = mean(g <= 30)))
}

## --------------------------------------------------------------- 2. 7 天分箱
## 右闭区间 (0,7], (7,14], ... ; 与 "间隔 <= 7 天" 的口径一致
make_bins <- function(x, binw = BINW) {
  g  <- x$gaps
  hi <- binw * ceiling(max(g) / binw)
  br <- seq(0, hi, by = binw)
  idx <- findInterval(g, br, left.open = TRUE)
  idx <- pmin(pmax(idx, 1L), length(br) - 1L)
  cnt <- tabulate(idx, nbins = length(br) - 1L)
  lo <- head(br, -1); hi_v <- tail(br, -1)
  rate <- 1 / mean(g)
  b <- data.frame(dataset = x$label, bin_lo = lo, bin_hi = hi_v,
                  bin_mid = (lo + hi_v) / 2,
                  count = as.integer(cnt),
                  exp_count = length(g) * (pexp(hi_v, rate) - pexp(lo, rate)))
  stopifnot(sum(b$count) == length(g),
            b$count[1] == sum(g <= binw),      # 首箱 = (0, 7]
            all(b$count >= 0))
  x$bins <- b
  x$rate <- rate
  x
}

D <- list(
  bus = make_bins(read_gaps("hackbus_daily.csv", "HACK-BUS")),
  med = make_bins(read_gaps("hackmed_daily.csv", "HACK-MED"))
)

## --------------------------------------------------------------- 3. 单面板图
draw_one <- function(x, panel = "", exp_ref = FALSE, zoom = FALSE) {
  b <- x$bins; s <- x$stats
  xmax <- if (zoom) ZOOM_TO else max(b$bin_hi)
  ymax <- max(c(b$count, if (exp_ref) b$exp_count else 0)) * 1.22
  ntail <- sum(x$gaps > ZOOM_TO)

  p <- ggplot(b, aes(x = bin_mid, y = count)) +
    geom_col(width = BINW, fill = COL_BAR, colour = "white",
             alpha = 0.85, linewidth = 0.35) +
    geom_vline(xintercept = s$mean,   linetype = "dashed",
               colour = COL_MEAN, linewidth = 1.1) +
    geom_vline(xintercept = s$median, linetype = "dotted",
               colour = COL_MED,  linewidth = 1.3) +
    annotate("text", x = xmax * 0.985, y = ymax * 0.91, hjust = 1,
             size = 4.2, colour = COL_MEAN,
             label = sprintf("Mean = %.1f days", s$mean)) +
    annotate("text", x = xmax * 0.985, y = ymax * 0.81, hjust = 1,
             size = 4.2, colour = COL_MED,
             label = sprintf("Median = %g days", s$median)) +
    annotate("text", x = xmax * 0.985, y = ymax * 0.71, hjust = 1,
             size = 3.3, colour = "grey35",
             label = sprintf("", s$n, BINW)) +
    labs(x = "Interarrival times (days)", y = "Breaches number",
         cex.axis=1.5,   title = panel) +
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    scale_x_continuous(breaks = pretty(c(0, xmax), n = 6),
                       expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          axis.title  = element_text(size = 14),
          axis.text   = element_text(size = 11),
          plot.title  = element_text(size = 13, face = "bold"),
          plot.margin = margin(6, 12, 6, 6))

  if (exp_ref) {
    p <- p +
      ## 阶梯线: 每个 7 天箱内的期望计数为常数, 与柱的口径完全一致
      geom_step(data = b, aes(x = bin_lo, y = exp_count),
                inherit.aes = FALSE, colour = COL_EXP,
                linewidth = 1.15, linetype = "dotdash", direction = "hv") +
      annotate("text", x = xmax * 0.985, y = ymax * 0.61, hjust = 1,
               size = 3.3, colour = COL_EXP,
               label = sprintf("Exponential gaps: rate = %.4f /day", x$rate)) +
      annotate("text", x = xmax * 0.985, y = ymax * 0.53, hjust = 1,
               size = 3.1, colour = "grey35",
               label = "Dotdash step: expected bin counts (i.i.d. exponential)")
  }
  if (zoom) {
    ## 用视窗裁剪而非 scale limits: 不丢弃数据, 也不产生 "removed rows" 警告
    p <- p +
      coord_cartesian(xlim = c(0, ZOOM_TO), expand = FALSE) +
      annotate("text", x = xmax * 0.985, y = ymax * 0.43, hjust = 1,
               size = 3.1, colour = "grey35",
               label = sprintf("%d gaps > %d d outside the displayed range",
                               ntail, ZOOM_TO))
  }
  p
}

## --------------------------------------------------------------- 4. 出图
FIGS <- list(
  both = function() grid.arrange(
    draw_one(D$bus, "(a) HACK-BUS"), draw_one(D$med, "(b) HACK-MED"), ncol = 2),
  zoom = function() grid.arrange(
    draw_one(D$bus, "(a) HACK-BUS", zoom = TRUE),
    draw_one(D$med, "(b) HACK-MED", zoom = TRUE), ncol = 2),
  bus  = function() print(draw_one(D$bus, "HACK-BUS")),
  med  = function() print(draw_one(D$med, "HACK-MED")),
  exp  = function() grid.arrange(
    draw_one(D$bus, "(a) HACK-BUS", exp_ref = TRUE, zoom = TRUE),
    draw_one(D$med, "(b) HACK-MED", exp_ref = TRUE, zoom = TRUE), ncol = 2)
)

if (FIG %in% names(FIGS)) FIGS[[FIG]]()
