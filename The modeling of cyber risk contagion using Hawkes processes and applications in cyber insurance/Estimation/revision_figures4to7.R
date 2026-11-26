## ---------------------------------------------------------------
## 修订版 Figures 4--7（送审稿应用部分）统一重绘 / 归集
##
## 送审稿原图定义（paper_text.txt 1123--1334 行）:
##   Figure 4  间隔直方图 (business (a) / medical (b))
##   Figure 5  间隔对指数分布的 QQ 图 (business (a) / medical (b))
##   Figure 6  固定箱长 tau=10、滞后 delta=1 的计数 ACF (business (a) / medical (b))
##   Figure 7  M1/M2 累计损失期望曲线互比  ->  修订为 held-out 实测 vs 预测
##
## 唯一数据源: hackmed_daily.csv / hackbus_daily.csv
##   PRC1.csv 2010-01-01~2018-12-30 + HACK + Total.Records>0
##   -> MED / BSO+BSF+BSR -> 同日合并累加 -> Z = 9.59 + 0.57 log Y
## 口径: 时间单位=天; 真实间隔 = inter[inter > 0] (首行 0 为起点)
##       ACF 窗口按日历 2010-01-01 ~ 2018-12-30 (日历口径, 见说明文档)
##
## 运行: Rscript --vanilla revision_figures4to7.R
## 输出: E:/paper/rdata/hackexp/hawkesplots/revised_figures4to7/
## ---------------------------------------------------------------
rm(list=ls())
DIR <- "E:/paper/rdata/hackexp"
OUT <- file.path(DIR, "hawkesplots", "revised_figures4to7")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(ggplot2))

FILES <- c("HACK-BUS" = "hackbus_daily.csv", "HACK-MED" = "hackmed_daily.csv")
COL   <- c("HACK-BUS" = "#2C7FB8", "HACK-MED" = "#D95F02")
ORIGIN <- as.Date("2010-01-01")
END    <- as.Date("2018-12-30")
TAU    <- 10L          # 箱长, 天
LAG    <- 30L          # 最大滞后, 箱
CI_Z   <- 1.96

## 坐标轴标题字号倍率 (等价于 base R 的 cex.lab = 1.5)
AXIS_CEX <- 1.5

## 统一主题：轴标题 = 基础字号 12 x 1.5 = 18 (用户明确要求)，去次要网格
## 刻度标签保持 10.5，避免与放大后的轴标题争夺空间 / 日期刻度拥挤
TH <- theme_minimal(base_size = 12) +
  theme(axis.title = element_text(size = 12 * AXIS_CEX),
        axis.text  = element_text(size = 10.5),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        panel.grid.minor = element_blank())

SER <- lapply(names(FILES), function(nm) {
  x <- read.csv(file.path(DIR, FILES[[nm]]), stringsAsFactors = FALSE)
  x$Date <- as.Date(x$Date)
  x <- x[order(x$Date), ]
  list(label = nm, Date = x$Date, gap = x$inter[x$inter > 0],
       Y = x$Total.Records, Z = x$Z)
})
names(SER) <- names(FILES)

save_both <- function(p, stem, w = 9.2, h = 4.8) {
  ggsave(file.path(OUT, paste0(stem, ".png")), p,
         width = w, height = h, dpi = 180, bg = "white")
  ggsave(file.path(OUT, paste0(stem, ".pdf")), p,
         width = w, height = h, device = "pdf", bg = "white")
  invisible(NULL)
}

cat("[0] 数据概览\n")
for (nm in names(SER)) {
  s <- SER[[nm]]
  cat(sprintf("  %-9s 事件日=%3d 间隔=%3d 首/末=%s/%s  均值=%.3f 中位=%g 最大=%g\n",
              nm, length(s$Date), length(s$gap),
              as.character(min(s$Date)), as.character(max(s$Date)),
              mean(s$gap), median(s$gap), max(s$gap)))
}

## ---------------------------------------------------------------
## Figure 4: 7 天箱间隔直方图 (右闭 (0,7],(7,14],...)
## ---------------------------------------------------------------
cat("\n[1] Figure 4 间隔直方图\n")

DH <- do.call(rbind, lapply(names(SER), function(nm)
  data.frame(dataset = nm, gap = SER[[nm]]$gap)))
DH$dataset <- factor(DH$dataset, levels = names(FILES))

ST <- do.call(rbind, lapply(names(SER), function(nm) {
  d <- SER[[nm]]$gap
  data.frame(dataset = nm,
             mean_d = mean(d), median_d = median(d),
             n_gaps = length(d),
             omitted_beyond_84 = sum(d > 84),
             label = sprintf("Mean = %.1f days\nMedian = %g days",
                             mean(d), median(d)))
}))
ST$dataset <- factor(ST$dataset, levels = names(FILES))

BRK <- seq(0, 7 * ceiling(max(DH$gap) / 7), by = 7)
p4 <- ggplot(DH, aes(x = gap)) +
  geom_histogram(aes(fill = dataset), breaks = BRK, closed = "right",
                 colour = "white", linewidth = 0.3, alpha = 0.9) +
  geom_vline(data = ST, aes(xintercept = mean_d),
             colour = "#B2182B", linetype = "dashed", linewidth = 0.9) +
  geom_vline(data = ST, aes(xintercept = median_d),
             colour = "#1A9850", linetype = "dotted", linewidth = 1.0) +
  geom_label(data = ST,
             aes(x = Inf, y = Inf, label = label),
             hjust = 1.05, vjust = 1.15, size = 3.5,
             fill = "white", alpha = 0.9, colour = "grey25",
             lineheight = 1.05) +
  facet_wrap(~ dataset, scales = "free_y") +
  scale_fill_manual(values = COL) +
  labs(x = "Interarrival times (days)", y = "Number of breaches",
       subtitle = "Bins are seven days wide and right-closed: (0,7], (7,14], ...") +
  TH + theme(legend.position = "none")
save_both(p4, "Figure4_interarrival_full")

## 正文推荐版：横轴放大到 0--84 天，尾部省略数写进 subtitle
p4z <- p4 +
  coord_cartesian(xlim = c(0, 84)) +
  labs(subtitle = sprintf(paste0("Seven-day bins; %d BUS and %d MED gaps beyond ",
                                 "84 days are omitted from the display only"),
                          ST$omitted_beyond_84[ST$dataset == "HACK-BUS"],
                          ST$omitted_beyond_84[ST$dataset == "HACK-MED"]))
save_both(p4z, "Figure4_interarrival_zoom")

## ---------------------------------------------------------------
## Figure 5: 间隔对指数分布的 QQ 图
## ---------------------------------------------------------------
cat("[2] Figure 5 指数 QQ 图\n")

QQ <- do.call(rbind, lapply(names(SER), function(nm) {
  d <- sort(SER[[nm]]$gap); n <- length(d)
  th <- qexp(ppoints(n), rate = 1 / mean(d))
  data.frame(dataset = nm, theoretical = th, sample = d)
}))
QQ$dataset <- factor(QQ$dataset, levels = names(FILES))

KS <- do.call(rbind, lapply(names(SER), function(nm) {
  d <- SER[[nm]]$gap
  kt <- suppressWarnings(ks.test(d, "pexp", rate = 1 / mean(d)))
  data.frame(dataset = nm, D = as.numeric(kt$statistic),
             label = sprintf("KS D = %.3f", as.numeric(kt$statistic)))
}))
KS$dataset <- factor(KS$dataset, levels = names(FILES))

p5 <- ggplot(QQ, aes(x = theoretical, y = sample)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey40", linewidth = 0.9) +
  geom_point(aes(colour = dataset), size = 1.5, alpha = 0.75) +
  geom_label(data = KS, aes(x = 0, y = Inf, label = label),
             hjust = -0.05, vjust = 1.15, size = 3.5,
             fill = "white", alpha = 0.9, colour = "grey25") +
  facet_wrap(~ dataset, scales = "free") +
  scale_colour_manual(values = COL) +
  labs(x = "Theoretical exponential quantiles (days)",
       y = "Sample interarrival quantiles (days)",
       subtitle = "") +
  TH + theme(legend.position = "none")
save_both(p5, "Figure5_qq_exponential")

## ---------------------------------------------------------------
## Figure 6: 计数 ACF (tau = 10 days, delta = 1 bin) + 按年去趋势对照
## ---------------------------------------------------------------
cat("[3] Figure 6 计数 ACF\n")

EDGES <- seq(ORIGIN, END + 1, by = paste0(TAU, " days"))
NBIN  <- length(EDGES) - 1L
YR    <- format(EDGES[-length(EDGES)], "%Y")

AC <- do.call(rbind, lapply(names(SER), function(nm) {
  id <- findInterval(as.numeric(SER[[nm]]$Date), as.numeric(EDGES))
  cnt <- tabulate(id, nbins = NBIN)
  yr_mean <- tapply(cnt, YR, mean)
  det <- cnt - as.numeric(yr_mean[YR])
  ar <- as.numeric(acf(cnt, lag.max = LAG, plot = FALSE)$acf)[-1L]
  ad <- as.numeric(acf(det, lag.max = LAG, plot = FALSE)$acf)[-1L]
  rbind(data.frame(dataset = nm, series = "raw counts",
                   lag = seq_len(LAG), acf = ar),
        data.frame(dataset = nm, series = "year-detrended",
                   lag = seq_len(LAG), acf = ad))
}))
AC$dataset <- factor(AC$dataset, levels = names(FILES))
AC$series  <- factor(AC$series, levels = c("raw counts", "year-detrended"))

CI <- CI_Z / sqrt(NBIN)
p6 <- ggplot(AC, aes(x = lag, y = acf, colour = series, linetype = series)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.5) +
  geom_hline(yintercept = c(-CI, CI), linetype = "dotted",
             colour = "grey45", linewidth = 0.7) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 1.6) +
  facet_wrap(~ dataset) +
  scale_colour_manual(values = c("raw counts" = "#2C7FB8",
                                 "year-detrended" = "#7A7A7A")) +
  scale_linetype_manual(values = c("raw counts" = "solid",
                                   "year-detrended" = "dashed")) +
  scale_x_continuous(breaks = seq(0, LAG, by = 5)) +
  labs(x = expression("Lag (bins of " * tau * " = 10 days; " * delta * " = 1 bin)"),
       y = "Autocorrelation of breach counts") +
  TH
save_both(p6, "Figure6_acf_counts", w = 9.6, h = 4.9)

## ---------------------------------------------------------------
## Figure 7: 归集已验证的 held-out 图（模型 vs 实测）
## ---------------------------------------------------------------
cat("[4] Figure 7 归集 held-out 图\n")
SRC <- file.path(DIR, "hawkesplots")
cp <- function(from, to) {
  s <- file.path(SRC, from); d <- file.path(OUT, to)
  if (file.exists(s)) { file.copy(s, d, overwrite = TRUE); TRUE } else FALSE
}
ok7 <- c(
  cp("fig_holdout_BUS.pdf", "Figure7_heldout_BUS.pdf"),
  cp("fig_holdout_BUS.png", "Figure7_heldout_BUS.png"),
  cp("fig_holdout_MED.pdf", "Figure7_heldout_MED.pdf"),
  cp("fig_holdout_MED.png", "Figure7_heldout_MED.png"),
  cp("fig_holdout_facets.png", "Figure7_heldout_facets.png"),
  cp("fig_holdout_rolling.png", "Figure7_heldout_rolling.png"))
cat("  已复制 Figure 7 文件数:", sum(ok7), "/", length(ok7), "\n")

## ---------------------------------------------------------------
## 统计汇总（供图注与回复信核对）
## ---------------------------------------------------------------
LAG1 <- do.call(rbind, lapply(names(SER), function(nm) {
  a <- AC[AC$dataset == nm, ]
  data.frame(dataset = nm,
             acf_lag1_raw = a$acf[a$lag == 1 & a$series == "raw counts"],
             acf_lag1_detrended = a$acf[a$lag == 1 & a$series == "year-detrended"])
}))
SUM <- merge(merge(ST[, c("dataset", "n_gaps", "mean_d", "median_d",
                          "omitted_beyond_84")],
                   KS[, c("dataset", "D")], by = "dataset"),
             LAG1, by = "dataset")
SUM$n_event_days <- vapply(as.character(SUM$dataset),
                           function(nm) length(SER[[nm]]$Date), integer(1))
SUM$ci_acf <- CI
write.csv(SUM, file.path(OUT, "figures4to6_stats.csv"), row.names = FALSE)
cat("\n-- 统计汇总 --\n")
print(SUM, row.names = FALSE, digits = 4)

cat("\n-- 输出文件 --\n")
print(list.files(OUT, pattern = "\\.(png|pdf|csv)$"))
cat("\nDONE\n")
