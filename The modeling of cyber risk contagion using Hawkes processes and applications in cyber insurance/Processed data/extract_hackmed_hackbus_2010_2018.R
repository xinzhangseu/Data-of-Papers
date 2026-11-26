## =====================================================================
## PRC1.csv -> 2010-01-01 ~ 2018-12-30, Total.Records > 0, HACK
##          -> hackmed (MED) / hackbus (BSO, BSF, BSR)
## 纯 base R; 点 Source 直接运行
##
## 列名口径: MED / BSO / BSF / BSR 属于 Type.of.organization(机构类型);
##           Type.of.breach 只有 HACK/PORT/DISC/PHYS/INSD/STAT/CARD。
## 日期口径: 闭区间 [2010-01-01, 2018-12-30], 单位 = 天(Date.Made.Public)。
## =====================================================================

DIR  <- "E:/paper/rdata/hackexp"
f    <- file.path(DIR, "PRC1.csv")
from <- as.Date("2010-01-01")
to   <- as.Date("2018-12-30")

d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
if (names(d)[1] == "") names(d)[1] <- "row_id"       # 首列是无名行号
d$Total.Records    <- suppressWarnings(as.numeric(d$Total.Records))  # 含 3e+05 写法
d$Date.Made.Public <- as.Date(d$Date.Made.Public)

## ---------------------------------------------- 条件 1: 三个条件一次筛完
keep <- with(d, Type.of.breach == "HACK" &
                 !is.na(Total.Records) & Total.Records > 0 &
                 !is.na(Date.Made.Public) &
                 Date.Made.Public >= from & Date.Made.Public <= to)
hack <- d[which(keep), ]
hack <- hack[order(hack$Date.Made.Public, hack$row_id), ]
row.names(hack) <- NULL

## ---------------------------------------------- 条件 2: 按机构类型拆两个集
hackmed <- hack[hack$Type.of.organization == "MED", ]
hackbus <- hack[hack$Type.of.organization %in% c("BSO", "BSF", "BSR"), ]
row.names(hackmed) <- NULL
row.names(hackbus) <- NULL

## -------------------------------------------------------------------- 校验
stopifnot(all(hack$Type.of.breach == "HACK"),
          all(hack$Total.Records > 0),
          all(hack$Date.Made.Public >= from & hack$Date.Made.Public <= to),
          all(hackmed$Type.of.organization == "MED"),
          all(hackbus$Type.of.organization %in% c("BSO", "BSF", "BSR")),
          length(intersect(hackmed$row_id, hackbus$row_id)) == 0,
          nrow(hackmed) + nrow(hackbus) <= nrow(hack))

## -------------------------------------------------------------------- 导出
write.csv(hack,    file.path(DIR, "hack_2010_2018.csv"), row.names = FALSE)
write.csv(hackmed, file.path(DIR, "hackmed.csv"),        row.names = FALSE)
write.csv(hackbus, file.path(DIR, "hackbus.csv"),        row.names = FALSE)

## -------------------------------------------------------------------- 核对
info <- function(x, label) {
  cat(sprintf("\n%s: %d 条 (占条件1 %.1f%%)\n",
              label, nrow(x), 100 * nrow(x) / nrow(hack)))
  cat(sprintf("  日期 %s ~ %s | 记录数合计 %s | 最小正值 %s\n",
              min(x$Date.Made.Public), max(x$Date.Made.Public),
              format(sum(x$Total.Records), scientific = FALSE),
              min(x$Total.Records)))
  cat("  机构类型: ")
  print(table(x$Type.of.organization))
}
cat(sprintf("原始 %d 条 -> 条件1 (2010-01-01~2018-12-30 & HACK & >0): %d 条\n",
            nrow(d), nrow(hack)))
cat(sprintf("条件1 日期范围: %s ~ %s | 记录数合计 %s\n",
            min(hack$Date.Made.Public), max(hack$Date.Made.Public),
            format(sum(hack$Total.Records), scientific = FALSE)))
cat(sprintf("边界检查: 2010-01-01 当天 %d 条; 数据集中晚于 2018-12-30 的记录 %d 条\n",
            sum(d$Date.Made.Public == from, na.rm = TRUE),
            sum(d$Date.Made.Public > to,   na.rm = TRUE)))
info(hackmed, "hackmed (HACK + MED)")
info(hackbus, "hackbus (HACK + BSO/BSF/BSR)")
cat(sprintf("\n[校验] hackmed %d + hackbus %d + 其余机构 %d = 条件1 %d\n",
            nrow(hackmed), nrow(hackbus),
            nrow(hack) - nrow(hackmed) - nrow(hackbus), nrow(hack)))
cat("\n结果对象: hack(条件1) / hackmed / hackbus\n")
cat(sprintf("导出: %s, %s, %s\n",
            file.path(DIR, "hack_2010_2018.csv"),
            file.path(DIR, "hackmed.csv"),
            file.path(DIR, "hackbus.csv")))
