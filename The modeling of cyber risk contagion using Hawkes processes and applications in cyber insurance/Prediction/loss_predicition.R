rm(list=ls())
suppressPackageStartupMessages(library(ggplot2))

## ------------------------------- 0. 数据目录与开关 -----------------------
DATA_HOME  <- "E:/paper/rdata/hackexp/era_hawkes"

script_dir <- function() {
  for (fr in sys.frames()) {
    of <- fr$ofile
    if (!is.null(of)) return(dirname(normalizePath(of, mustWork = FALSE)))
  }
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f))
    return(dirname(normalizePath(sub("^--file=", "", f[1L]), mustWork = FALSE)))
  getwd()
}
DIR <- if (dir.exists(DATA_HOME)) DATA_HOME else script_dir()

WHICH      <- "ALL"       # "MED" | "BUS" | "TERMINAL" | "ALL"
BAND_SPEC  <- "M1_tv"
AXIS_CEX   <- 1.5         # 轴标题放大倍数 (与 Figures 4-7 一致)
LEGEND_CEX <-1.4
NEW_WINDOW <- FALSE       # TRUE: 每张图单独一个窗口; FALSE: 依次显示, 可翻页
CUT        <- as.Date("2016-01-01")
END        <- as.Date("2018-12-30")

f_cur <- file.path(DIR, "table34_heldout_curves.csv")
f_sum <- file.path(DIR, "table34_heldout_summary.csv")
if (!all(file.exists(c(f_cur, f_sum))))
  stop("缺少结果 CSV。请把两个 table34_heldout_*.csv 放到: ", DIR,
       "\n或先运行同目录下的 referee_table34_heldout_predict.R")

CUR <- read.csv(f_cur, stringsAsFactors = FALSE); CUR$Date <- as.Date(CUR$Date)
SUM <- read.csv(f_sum, stringsAsFactors = FALSE)

## ------------------------------- 1. 规格、颜色、线型 ---------------------
SPECS <- c("P", "Pa", "M1", "M2", "M1_tv", "M2_tv")
SPECLAB <- c(P     = "Poisson homogeneous (5.3)",
             Pa    = "Poisson, annual (5.4), last-rate",
             M1    = "Hawkes gamma=0, const mu",
             M2    = "Hawkes gamma=1, const mu",
             M1_tv = "Hawkes gamma=0, annual mu(t)",
             M2_tv = "Hawkes gamma=1, annual mu(t)")
COLS <- c(P = "#7A7A7A", Pa = "#B0B0B0", M1 = "#2C7FB8",
          M2 = "#7FB2D9", M1_tv = "#D94841", M2_tv = "#EFA9A3")
LTS  <- c(P = "dotted", Pa = "dashed", M1 = "solid",
          M2 = "dashed", M1_tv = "solid", M2_tv = "dashed")
SHORT <- c(P = "P", Pa = "Pa",
           M1 = "M1\nconst, k=0", M2 = "M2\nconst, k=1",
           M1_tv = "M1_tv\nannual, k=0", M2_tv = "M2_tv\nannual, k=1")
DSOF <- c(MED = "HACK-MED", BUS = "HACK-BUS")

CUR$spec <- factor(CUR$spec, levels = SPECS)
SUM$spec <- factor(SUM$spec, levels = SPECS)

## ------------------------------- 2. 实测累计损失阶梯线 -------------------
obs_step <- function(ds) {
  fn <- file.path(DIR, if (ds == "HACK-MED") "hackmed_daily.csv"
                  else "hackbus_daily.csv")
  if (!file.exists(fn)) return(NULL)
  d <- read.csv(fn, stringsAsFactors = FALSE)
  d$Date <- as.Date(d$Date)
  d <- d[d$Date >= CUT & d$Date <= END, ]
  d <- d[order(d$Date), ]
  data.frame(Date = c(CUT, d$Date, END),
             L = c(0, cumsum(d$Z), sum(d$Z)))
}

obs_grid <- function(ds) {
  cc <- CUR[CUR$dataset == ds & CUR$spec == "P", ]
  data.frame(Date = cc$Date, L = cc$L_obs)
}

## ------------------------------- 3. 图 A: 累计损失轨迹 -------------------
plot_loss_curve <- function(nm) {
  ds <- DSOF[[nm]]
  cc <- CUR[CUR$dataset == ds, ]
  ob <- obs_step(ds); if (is.null(ob)) ob <- obs_grid(ds)
  band <- cc[cc$spec == BAND_SPEC, ]
  s <- SUM[SUM$dataset == ds, ]
  Lobs <- s$L_obs_end[1L]
  
  cat(sprintf("[%s] 阶梯线终点 L = %.3f, summary 中 L_obs_end = %.3f, 差 = %.2g\n",
              ds, max(ob$L), Lobs, max(ob$L) - Lobs))
  
  ggplot() +
    geom_ribbon(data = band, aes(x = Date, ymin = L_lo, ymax = L_hi),
                fill = "#D94841", alpha = 0.12, colour = NA) +
    geom_line(data = cc, aes(x = Date, y = L_pred, colour = spec,
                             linetype = spec), linewidth = 1.05) +
    geom_step(data = ob, aes(x = Date, y = L), direction = "hv",
              colour = "black", linewidth = 1.45) +
    scale_colour_manual(values = COLS, labels = unname(SPECLAB[SPECS]),
                        name = NULL) +
    scale_linetype_manual(values = LTS, labels = unname(SPECLAB[SPECS]),
                          name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    guides(colour = guide_legend(nrow = 2L, byrow = TRUE)) +
    labs(x = " Date",
         y = "Cumulative loss") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(size = 10 * LEGEND_CEX),
          axis.title = element_text(size = 13 * AXIS_CEX),
          axis.text = element_text(size = 11),
          plot.margin = margin(6, 16, 6, 30, unit = "pt"),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 9.5, colour = "grey30"),
          panel.grid.minor = element_blank())
}

## ------------------------------- 4. 图 B: 终点损失与 90% 区间 -----------
plot_loss_terminal <- function() {
  obs <- unique(SUM[, c("dataset", "L_obs_end")])
  d <- SUM
  lev <- c("HACK-MED", "HACK-BUS")     # MED 在左, 与正文/表 3、表 4 顺序一致
  d$dataset <- factor(d$dataset, levels = lev)
  obs$dataset <- factor(obs$dataset, levels = lev)
  d$lab <- sprintf("%+.0f%%", 100 * (d$L_pred_end - d$L_obs_end) / d$L_obs_end)
  
  ggplot(d, aes(x = spec, y = L_pred_end, colour = spec)) +
    geom_hline(data = obs, aes(yintercept = L_obs_end),
               colour = "grey25", linetype = "dashed", linewidth = 0.9) +
    geom_errorbar(aes(ymin = L_lo_end, ymax = L_hi_end),
                  width = 0.16, linewidth = 0.85) +
    geom_point(size = 3.1) +
    geom_text(aes(y = L_hi_end, label = lab), vjust = -0.5,
              size = 3.1, colour = "grey20") +
    facet_wrap(~ dataset, ncol = 2L, scales = "free_y") +
    scale_colour_manual(values = COLS, guide = "none") +
    scale_x_discrete(labels = unname(SHORT[SPECS])) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
    labs(x = "Specification (estimated on 2010-2015 only)",
         y = "Terminal loss L(T)  (sum of Z)",
         title = "",
         subtitle = "Points = mean of 1,000 simulated paths; bars = 5th/95th percentiles; dashed line = realised; labels = relative bias.") +
    theme_minimal(base_size = 13) +
    theme(axis.title = element_text(size = 13 * AXIS_CEX),
          axis.text.x = element_text(size = 9.5, lineheight = 0.95),
          axis.text.y = element_text(size = 10.5),
          plot.margin = margin(6, 16, 6, 34, unit = "pt"),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 9.5, colour = "grey30"),
          strip.text = element_text(size = 12, face = "bold"),
          panel.grid.minor = element_blank())
}

## ------------------------------- 5. 直接显示 (不写文件) ------------------
PLOTS <- list()
if (WHICH %in% c("MED", "ALL"))       PLOTS[["MED"]]      <- plot_loss_curve("MED")
if (WHICH %in% c("BUS", "ALL"))       PLOTS[["BUS"]]      <- plot_loss_curve("BUS")
if (WHICH %in% c("TERMINAL", "ALL"))  PLOTS[["terminal"]] <- plot_loss_terminal()

if (interactive()) {
  for (nm in names(PLOTS)) {
    if (NEW_WINDOW) dev.new()
    print(PLOTS[[nm]])
  }
} else {
  cat("非交互环境 (Rscript): 已生成", length(PLOTS),
      "个图对象, 未写文件; 在 RStudio 中 Source 本文件即可直接显示图形。\n")
}

## ------------------------------- 6. 控制台核对表 -------------------------

invisible(PLOTS)
