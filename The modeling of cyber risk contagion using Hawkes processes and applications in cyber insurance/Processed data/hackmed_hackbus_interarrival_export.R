## =====================================================================
## 导出: 仿论文 Figure 4 的间隔直方图 PNG + 分箱/统计表
## 运行: Rscript hackmed_hackbus_interarrival_export.R
## 主脚本只负责出图; 这里设 FIG="none" 后再 source, 避免重复弹图
## =====================================================================

DIR <- "E:/paper/rdata/hackexp"
OUT <- file.path(DIR, "hawkesplots")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

FIG <- "none"
source("E:/paper/rdata/hackexp")

## --------------------------------------------------------------- PNG
ggsave(file.path(OUT, "fig_interarrival_hist_both.png"),
       plot = arrangeGrob(grobs = list(
         draw_one(D$bus, "(a) HACK-BUS"),
         draw_one(D$med, "(b) HACK-MED")), ncol = 2),
       width = 13.2, height = 5.4, dpi = 150, bg = "white")

ggsave(file.path(OUT, "fig_interarrival_hist_zoom.png"),
       plot = arrangeGrob(grobs = list(
         draw_one(D$bus, "(a) HACK-BUS", zoom = TRUE),
         draw_one(D$med, "(b) HACK-MED", zoom = TRUE)), ncol = 2),
       width = 13.2, height = 5.4, dpi = 150, bg = "white")

ggsave(file.path(OUT, "fig_interarrival_hist_bus.png"),
       plot = draw_one(D$bus, "HACK-BUS"),
       width = 7.4, height = 5.2, dpi = 150, bg = "white")

ggsave(file.path(OUT, "fig_interarrival_hist_med.png"),
       plot = draw_one(D$med, "HACK-MED"),
       width = 7.4, height = 5.2, dpi = 150, bg = "white")

ggsave(file.path(OUT, "fig_interarrival_hist_exponential.png"),
       plot = arrangeGrob(grobs = list(
         draw_one(D$bus, "(a) HACK-BUS", exp_ref = TRUE, zoom = TRUE),
         draw_one(D$med, "(b) HACK-MED", exp_ref = TRUE, zoom = TRUE)), ncol = 2),
       width = 13.2, height = 5.4, dpi = 150, bg = "white")

## --------------------------------------------------------------- CSV
bins <- do.call(rbind, lapply(D, function(x) x$bins))
row.names(bins) <- NULL
write.csv(bins, file.path(DIR, "hackmed_hackbus_interarrival_bins.csv"),
          row.names = FALSE)

summ <- do.call(rbind, lapply(D, function(x) {
  s <- x$stats
  ## 间隔以整数天记录, 存在大量并列, KS 的渐近 p 值只作参考
  ks <- suppressWarnings(stats::ks.test(x$gaps, "pexp", rate = x$rate))
  data.frame(dataset = x$label,
             n_event_days = nrow(x$daily), n_gaps = s$n,
             min_d = s$min, median_d = s$median, mean_d = round(s$mean, 4),
             max_d = s$max, sd_d = round(s$sd, 4), CV = round(s$cv, 4),
             P_le7d = round(s$p_le7, 4), P_le30d = round(s$p_le30, 4),
             exp_rate_per_day = round(x$rate, 6),
             bin1_observed = x$bins$count[1],
             bin1_expected_exponential = round(x$bins$exp_count[1], 2),
             bin1_obs_exp_ratio = round(x$bins$count[1] / x$bins$exp_count[1], 3),
             ks_D = round(as.numeric(ks$statistic), 4),
             ks_p_approx = signif(ks$p.value, 4))
}))
write.csv(summ, file.path(DIR, "hackmed_hackbus_interarrival_stats.csv"),
          row.names = FALSE)

## --------------------------------------------------------------- 核对打印
cat("=== 仿 Figure 4 间隔直方图 ===\n")
print(summ[, c("dataset", "n_gaps", "mean_d", "median_d", "max_d",
               "P_le7d", "bin1_observed", "bin1_expected_exponential",
               "bin1_obs_exp_ratio", "ks_D", "ks_p_approx")])
cat("\n首箱 (0,7] 与 i.i.d. 指数间隔的期望对比: 观测/期望 > 1 即短间隔超额。\n")
for (nm in names(D)) {
  x <- D[[nm]]
  cat(sprintf("%s: 7 天箱 %d 个, 最大间隔 %d 天, 超过 %d 天的间隔 %d 个\n",
              x$label, nrow(x$bins), max(x$gaps), ZOOM_TO,
              sum(x$gaps > ZOOM_TO)))
}
cat("\n已导出 5 张 PNG 到 ", OUT, " 及 2 张 CSV 到 ", DIR, "\n", sep = "")
