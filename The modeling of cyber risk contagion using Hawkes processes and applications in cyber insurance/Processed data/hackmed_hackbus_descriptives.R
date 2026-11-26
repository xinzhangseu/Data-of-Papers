## =====================================================================
## hackmed / hackbus 的统计描述
## 输入: hackmed.csv, hackbus.csv (明细); hackmed_daily.csv, hackbus_daily.csv (日序列)
## 输出: hackmed_hackbus_descriptives.csv  (主要变量描述统计长表)
##       hackmed_hackbus_descriptives_extra.csv (概览/间隔/年度/机构分组)
## 纯 base R; 点 Source 直接运行。所有输入只读、不覆写。
## =====================================================================

DIR <- "E:/paper/rdata/hackexp"
raw <- list(hackmed = "hackmed.csv",      hackbus = "hackbus.csv")
day <- list(hackmed = "hackmed_daily.csv", hackbus = "hackbus_daily.csv")

stats1 <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  m <- mean(x); s <- sd(x)
  c(n = length(x), min = min(x), q1 = unname(quantile(x, .25)),
    median = median(x), mean = m, q3 = unname(quantile(x, .75)),
    max = max(x), sd = s, cv = s / m,
    skewness = mean(((x - m) / s)^3), sum = sum(x))
}

main_rows <- list(); extra_rows <- list()
add_extra <- function(ds, group, item, value)
  extra_rows[[length(extra_rows) + 1L]] <<-
    data.frame(dataset = ds, group = group, item = item, value = as.character(value))

for (nm in names(day)) {
  r <- read.csv(file.path(DIR, raw[[nm]]), stringsAsFactors = FALSE, check.names = FALSE)
  d <- read.csv(file.path(DIR, day[[nm]]),  stringsAsFactors = FALSE, check.names = FALSE)
  r$Date.Made.Public <- as.Date(r$Date.Made.Public)
  d$Date <- as.Date(d$Date)
  d$Total.Records <- suppressWarnings(as.numeric(d$Total.Records))
  r$Total.Records <- suppressWarnings(as.numeric(r$Total.Records))

  ## ---------------- 主要变量描述统计 (inter 去掉首行起点 0)
  vars <- list(
    "Total.Records (daily sum)"   = d$Total.Records,
    "ln_Total.Records"            = d$ln_Total.Records,
    "Z (loss mark)"               = d$Z,
    "inter (days, excl. origin)"  = d$inter[-1L],
    "n_events (records per day)"  = d$n_events,
    "Total.Records (raw record)"  = r$Total.Records
  )
  for (v in names(vars)) {
    s <- stats1(vars[[v]])
    main_rows[[length(main_rows) + 1L]] <-
      cbind(data.frame(dataset = nm, variable = v), as.data.frame(as.list(s)))
  }

  ## ---------------- 概览
  add_extra(nm, "overview", "raw records", nrow(r))
  add_extra(nm, "overview", "distinct event days", nrow(d))
  add_extra(nm, "overview", "merged records", nrow(r) - nrow(d))
  add_extra(nm, "overview", "dates with >1 record", sum(d$n_events > 1))
  add_extra(nm, "overview", "max records on one day", max(d$n_events))
  add_extra(nm, "overview", "first date", as.character(min(d$Date)))
  add_extra(nm, "overview", "last date", as.character(max(d$Date)))
  add_extra(nm, "overview", "span (days)", as.integer(max(d$Date) - min(d$Date)))
  add_extra(nm, "overview", "sum of daily records",
            format(sum(d$Total.Records), scientific = FALSE))
  add_extra(nm, "overview", "sum equals raw sum",
            isTRUE(all.equal(sum(d$Total.Records), sum(r$Total.Records))))

  ## ---------------- 间隔 / 聚集的描述性指标
  g <- d$inter[-1L]
  cv <- sd(g) / mean(g)
  x <- as.numeric(d$Date - min(d$Date))
  nwin <- floor(max(x) / 7) + 1L
  cnt <- hist(x, breaks = seq(0, 7 * nwin, by = 7), right = FALSE, plot = FALSE)$counts
  add_extra(nm, "intervals", "n gaps", length(g))
  add_extra(nm, "intervals", "mean gap (days)", round(mean(g), 4))
  add_extra(nm, "intervals", "median gap (days)", median(g))
  add_extra(nm, "intervals", "sd / CV", sprintf("%.4f / %.4f", sd(g), cv))
  add_extra(nm, "intervals", "burstiness B=(CV-1)/(CV+1)",
            sprintf("%.4f", (cv - 1) / (cv + 1)))
  add_extra(nm, "intervals", "gaps = 1 day", sum(g == 1))
  add_extra(nm, "intervals", "P(gap <= 7 days)",
            sprintf("%.4f", mean(g <= 7)))
  add_extra(nm, "intervals", "Fano (7-day windows)",
            sprintf("%.4f", var(cnt) / mean(cnt)))
  add_extra(nm, "intervals", "non-empty 7-day windows", sum(cnt > 0))
  add_extra(nm, "intervals", "max event days in one window", max(cnt))

  ## ---------------- 分年度事件日数
  yr <- table(format(d$Date, "%Y"))
  for (k in names(yr)) add_extra(nm, "by year", k, as.integer(yr[k]))

  ## ---------------- 机构类型构成 (明细条数 / 事件日数)
  org <- table(r$Type.of.organization)
  for (k in names(org)) {
    sub <- r[r$Type.of.organization == k, ]
    add_extra(nm, "organization", sprintf("%s: records / event days", k),
              sprintf("%d / %d", as.integer(org[k]),
                      length(unique(sub$Date.Made.Public))))
  }
}

main  <- do.call(rbind, main_rows)
extra <- do.call(rbind, extra_rows)
num <- sapply(main, is.numeric)
main[num] <- lapply(main[num], function(z) round(z, 4))

write.csv(main,  file.path(DIR, "hackmed_hackbus_descriptives.csv"),       row.names = FALSE)
write.csv(extra, file.path(DIR, "hackmed_hackbus_descriptives_extra.csv"), row.names = FALSE)

cat("=== 主要变量描述统计 ===\n")
print(main, row.names = FALSE)
cat("\n=== 概览 / 间隔 / 分组 ===\n")
print(extra, row.names = FALSE)
cat("\n结果对象: main / extra\n")
