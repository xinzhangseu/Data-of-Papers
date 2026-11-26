## =====================================================================
##  Export MED / BUS arrival series and transformed-loss series
##  for the training and held-out windows, as CSV files.
##
##  数据: hackmed_daily.csv / hackbus_daily.csv (同日合并, Y 累加,
##        Z = 9.59 + 0.57*ln(Y)), 与全文 PRC1.csv 2010-2018 口径一致.
##  切分: 与 referee_heldout_recompute.R 的 make_split() 完全一致
##        origin   = 各序列首个事件日 (MED 2010-03-02, BUS 2010-01-01)
##        t_global = Date - origin (模型内部时间, 单位: 天)
##        train    = t_global <  2016-01-01 - origin
##        test     = 2016-01-01 - origin <= t_global <= 2018-12-30 - origin
##
##  输出 (全部在 DIR):
##    MED_train_arrivals.csv / MED_test_arrivals.csv   内部到达(时刻/间隔)序列
##    BUS_train_arrivals.csv / BUS_test_arrivals.csv
##    MED_train_loss.csv     / MED_test_loss.csv       变换后损失序列 Z
##    BUS_train_loss.csv     / BUS_test_loss.csv
##    series_metadata.csv                            每段窗口与汇总, 便于复核
##
##  两类序列按 idx / Date 一一对齐, 每行 = 一个事件日。
## =====================================================================

DIR   <- "E:/paper/rdata/hackexp/era_hawkes"
CUT   <- as.Date("2016-01-01")
END   <- as.Date("2018-12-30")
FILES <- c(MED = "hackmed_daily.csv", BUS = "hackbus_daily.csv")
LABEL <- c(MED = "HACK-MED", BUS = "HACK-BUS")

meta <- list()
for (nm in names(FILES)) {
  d <- read.csv(file.path(DIR, FILES[[nm]]),
                stringsAsFactors = FALSE, check.names = FALSE)
  d$Date <- as.Date(d$Date)
  d <- d[order(d$Date), ]
  stopifnot(!anyDuplicated(d$Date), all(d$Z > 0))

  origin <- min(d$Date)
  T_tr   <- as.numeric(CUT - origin)
  T_end  <- as.numeric(END - origin)

  d$t_global   <- as.numeric(d$Date - origin)
  d$inter_prev <- c(NA_integer_, as.integer(diff(d$Date)))
  d$Z_check    <- 9.59 + 0.57 * d$ln_Total.Records
  stopifnot(max(abs(d$Z_check - d$Z)) < 1e-8)   # 变换口径校验

  seg <- list(train = d[d$t_global <  T_tr, ],
              test  = d[d$t_global >= T_tr & d$t_global <= T_end, ])
  seg_start <- list(train = origin, test = CUT)
  seg_t0    <- c(train = 0, test = T_tr)
  seg_t1    <- c(train = T_tr, test = T_end)

  for (sp in c("train", "test")) {
    s <- seg[[sp]]
    n <- nrow(s)
    idx <- seq_len(n)

    ## --- 内部到达序列: 全局内部时刻 + 段内时刻 + 两种间隔 ---
    arr <- data.frame(
      idx           = idx,
      Date          = s$Date,
      t_global      = s$t_global,                             # 距序列首个事件日
      t_segment     = as.numeric(s$Date - seg_start[[sp]]),  # 距本段起点
      inter_prev    = s$inter_prev,                           # 与全局上一事件日的间隔(天)
      inter_within  = c(NA_integer_, as.integer(diff(s$Date))), # 段内相邻间隔
      N_cum         = idx,                                    # 段内累计事件数
      stringsAsFactors = FALSE)

    ## --- 变换后损失序列: 原始记录数 -> ln -> Z -> 累计损失 ---
    loss <- data.frame(
      idx             = idx,
      Date            = s$Date,
      n_events        = s$n_events,            # 同日合并的明细条数
      Total.Records   = s$Total.Records,       # 同日累加后的原始记录数
      ln_Total.Records= s$ln_Total.Records,    # ln(Y)
      Z               = s$Z,                   # 变换后损失 = 9.59 + 0.57 ln(Y)
      dZ              = c(NA_real_, diff(s$Z)), # 相邻事件日 Z 的变化量(首行 NA)
      L_cum           = cumsum(s$Z),           # 段内累计损失
      stringsAsFactors = FALSE)

    fa <- file.path(DIR, sprintf("%s_%s_arrivals.csv", nm, sp))
    fl <- file.path(DIR, sprintf("%s_%s_loss.csv",      nm, sp))
    write.csv(arr,  fa, row.names = FALSE)
    write.csv(loss, fl, row.names = FALSE)

    meta[[length(meta) + 1L]] <- data.frame(
      dataset        = LABEL[[nm]],
      split          = sp,
      origin         = as.character(origin),
      window_start   = as.character(seg_start[[sp]]),
      window_end     = as.character(if (sp == "train") CUT - 1L else END),
      first_event_date = as.character(min(s$Date)),
      last_event_date  = as.character(max(s$Date)),
      T_start_global = seg_t0[[sp]],
      T_end_global   = seg_t1[[sp]],
      n_events       = n,
      sum_Z          = sum(s$Z),
      mean_Z         = mean(s$Z),
      mean_inter_prev   = mean(s$inter_prev,   na.rm = TRUE),
      mean_inter_within = mean(arr$inter_within, na.rm = TRUE),
      max_inter_prev    = max(s$inter_prev,    na.rm = TRUE),
      Z_min          = min(s$Z),
      Z_max          = max(s$Z),
      stringsAsFactors = FALSE)
  }
}
META <- do.call(rbind, meta)
write.csv(META, file.path(DIR, "series_metadata.csv"), row.names = FALSE)

## ---------------- 与既有留出段结果对拍 ----------------------------------
SUM <- read.csv(file.path(DIR, "table34_heldout_summary.csv"),
                stringsAsFactors = FALSE)
chk <- META[META$split == "test", c("dataset", "n_events", "sum_Z")]
ref <- do.call(rbind, lapply(unique(chk$dataset), function(ds) {
  r <- SUM[SUM$dataset == ds, ][1L, ]
  data.frame(dataset = ds, n_events = r$N_obs, sum_Z = r$L_obs_end)
}))
m <- merge(chk, ref, by = "dataset", suffixes = c("_new", "_ref"))
cat("=== 测试段对拍 (与 table34_heldout_summary.csv) ===\n")
print(m, row.names = FALSE, digits = 8)
cat("最大绝对偏差: n =", max(abs(m$n_events_new - m$n_events_ref)),
    ", sum_Z =", sprintf("%.3g", max(abs(m$sum_Z_new - m$sum_Z_ref))), "\n")

cat("\n=== 四段窗口汇总 ===\n")
print(META[, c("dataset", "split", "origin", "window_start", "window_end",
               "T_start_global", "T_end_global", "n_events", "sum_Z", "mean_Z")],
      row.names = FALSE, digits = 6)
cat("\n已导出 8 个序列 CSV + series_metadata.csv 到:", DIR, "\n")
