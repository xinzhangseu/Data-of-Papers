## =====================================================================
## HACK + MED 数据集 -> 以天为单位的间隔序列
## 规则: 同一天的多条记录只算一个事件日, Total.Records 累加
## 输出列与既有 prc_daily.csv 一致: Date, Total.Records, lossnew, inter
##       (inter 首行 = 0, 之后 = 相邻事件日的天数差)
##       lossnew = 9.59 + 0.57 * ln(Total.Records)
## 纯 base R; 点 Source 直接运行
## =====================================================================

DIR <- "E:/paper/rdata/hackexp"
f   <- file.path(DIR, "PRC1_HACK_MED_positive.csv")

d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
d$Date.Made.Public <- as.Date(d$Date.Made.Public)
d$Total.Records    <- suppressWarnings(as.numeric(d$Total.Records))
stopifnot(all(d$Type.of.breach == "HACK"),
          all(d$Type.of.organization == "MED"),
          all(!is.na(d$Total.Records) & d$Total.Records > 0))

## --------------------------------------------- 按日聚合(同一天只算一次)
## 用非 formula 形式, 避免列名含点/首列无名时 aggregate 报错
daily <- aggregate(d$Total.Records,
                   by = list(Date = d$Date.Made.Public), FUN = sum)
names(daily) <- c("Date", "Total.Records")
daily <- daily[order(daily$Date), ]
cnt <- table(d$Date.Made.Public)
daily$n_events <- as.integer(cnt[as.character(daily$Date)])

## ------------------------------------------------- 损失变换与日间隔(天)
daily$lossnew <- 9.59 + 0.57 * log(daily$Total.Records)
daily$inter   <- c(0L, as.integer(diff(daily$Date)))   # 首行 0 表示起点

## --------------------------------------------------------------- 校验
stopifnot(nrow(daily) == length(unique(d$Date.Made.Public)),
          isTRUE(all.equal(sum(daily$Total.Records), sum(d$Total.Records))),
          sum(daily$n_events) == nrow(d),
          all(daily$inter[-1] > 0))                    # 去重后不应有 0 间隔

cat(sprintf("明细 %d 条 -> 事件日 %d 个 (合并了 %d 条同日记录)\n",
            nrow(d), nrow(daily), nrow(d) - nrow(daily)))
cat(sprintf("  有重复的日期 %d 天, 单日最多 %d 条\n",
            sum(daily$n_events > 1), max(daily$n_events)))
cat(sprintf("  记录数合计(前后一致): %s\n",
            format(sum(daily$Total.Records), scientific = FALSE)))
cat(sprintf("  日期范围: %s ~ %s\n", min(daily$Date), max(daily$Date)))
it <- daily$inter[-1]                                  # 真实的 229 个间隔
cat(sprintf("间隔(天, 共 %d 个): 最小 %d | 中位 %.1f | 均值 %.2f | 最大 %d | 间隔=1 的 %d 个\n",
            length(it), min(it), median(it), mean(it), max(it), sum(it == 1)))
cat("\n前 8 个事件日:\n")
print(head(daily[, c("Date", "Total.Records", "lossnew", "inter", "n_events")], 8),
      row.names = FALSE)

## --------------------------------------------------------------- 导出
out <- file.path(DIR, "PRC1_HACK_MED_daily_interval.csv")
write.csv(daily[, c("Date", "Total.Records", "lossnew", "inter", "n_events")],
          out, row.names = FALSE)
cat(sprintf("\n已导出 %d 个事件日到: %s\n", nrow(daily), out))

## 结果对象: d(273 条明细) / daily(230 个事件日) / it(229 个间隔向量)
