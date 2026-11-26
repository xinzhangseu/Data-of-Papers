## =====================================================================
## hackmed.csv / hackbus.csv -> 同日合并(Total.Records 累加)
##                           -> ln(累加值) + 论文变换 Z = 9.59 + 0.57*ln(Y)
## 纯 base R; 点 Source 直接运行。原始明细文件只读取、不覆写。
##
## 顺序不可颠倒: 先累加再取对数, ln(sum_j Y_j) != sum_j ln(Y_j)。
## =====================================================================

DIR <- "E:/paper/rdata/hackexp"

process <- function(fname, outname, label) {
  f <- file.path(DIR, fname)
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  d$Date.Made.Public <- as.Date(d$Date.Made.Public)
  d$Total.Records    <- suppressWarnings(as.numeric(d$Total.Records))
  stopifnot(all(!is.na(d$Total.Records) & d$Total.Records > 0))

  ## ---- 同日整合: 记录数求和, 并记下当日合并了几条明细
  daily <- aggregate(
    cbind(Total.Records = d$Total.Records, n_events = rep(1L, nrow(d))),
    by = list(Date = d$Date.Made.Public), FUN = sum)
  daily <- daily[order(daily$Date), ]
  row.names(daily) <- NULL

  ## ---- 变换: 对累加后的 Total.Records 取自然对数
  daily$ln_Total.Records <- log(daily$Total.Records)
  daily$Z                <- 9.59 + 0.57 * log(daily$Total.Records)
  daily$inter            <- c(0L, as.integer(diff(daily$Date)))  # 首行 0 表示起点

  ## ---- 校验
  stopifnot(nrow(daily) == length(unique(d$Date.Made.Public)),
            isTRUE(all.equal(sum(daily$Total.Records), sum(d$Total.Records))),
            sum(daily$n_events) == nrow(d),
            all(daily$Total.Records > 0),
            all(diff(daily$Date) >= 1),
            sum(daily$inter) == as.integer(max(d$Date.Made.Public) -
                                           min(d$Date.Made.Public)))

  write.csv(daily, file.path(DIR, outname), row.names = FALSE)

  ## ---- 核对打印
  dup  <- sum(daily$n_events > 1)
  imin <- which.min(daily$Z); imax <- which.max(daily$Z)
  cat(sprintf("\n===== %s (%s) =====\n", label, fname))
  cat(sprintf("明细 %d 条 -> 事件日 %d 个 (合并 %d 条; 有 %d 天重复, 单日最多 %d 条)\n",
              nrow(d), nrow(daily), nrow(d) - nrow(daily), dup,
              max(daily$n_events)))
  cat(sprintf("日期 %s ~ %s | 累加前后记录数合计均为 %s\n",
              min(daily$Date), max(daily$Date),
              format(sum(daily$Total.Records), scientific = FALSE)))
  cat(sprintf("真实间隔 %d 个, 累计 %d 天 (等于起止日期跨度); 间隔为 1 天的有 %d 个\n",
              nrow(daily) - 1L, sum(daily$inter), sum(daily$inter == 1)))
  cat(sprintf("Z 范围 %.4f ~ %.4f | 均值 %.4f | 中位 %.4f\n",
              min(daily$Z), max(daily$Z), mean(daily$Z), median(daily$Z)))
  cat(sprintf("Z 最小: %s, Y=%s ; Z 最大: %s, Y=%s\n",
              daily$Date[imin], format(daily$Total.Records[imin], scientific = FALSE),
              daily$Date[imax], format(daily$Total.Records[imax], scientific = FALSE)))
  cat(sprintf("已导出: %s\n", file.path(DIR, outname)))
  invisible(daily)
}

hackmed_daily <- process("hackmed.csv", "hackmed_daily.csv", "hackmed")
hackbus_daily <- process("hackbus.csv", "hackbus_daily.csv", "hackbus")

## 结果对象: hackmed_daily / hackbus_daily
## 列: Date, Total.Records(累加), n_events, ln_Total.Records, Z, inter
