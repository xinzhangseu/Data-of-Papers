## =====================================================================
##  Tables 3 / 4 : training-window MLE for the six candidate specifications
##  Self-contained, base R only.  Source this file -> estimates + CSV outputs.
##
##  Data (same directory): hackmed_daily.csv, hackbus_daily.csv
##    PRC1.csv, 2010-01-01..2018-12-30, HACK, Total.Records > 0, MED / BSO+BSF+BSR,
##    same-day records summed, Z = 9.59 + 0.57 log(Total.Records).
##  Time unit = day; t_i = Date_i - first event day of that series.
##  Training window: t < CUT (2016-01-01); no held-out event enters any fit.
##
##  Specifications (lambda(t) = mu(t) + alpha * sum_{T_j<t} w_j exp(-beta(t-T_j))):
##    P      mu(t)=mu (homogeneous Poisson, eq. 5.3)         alpha = 0
##    Pa     mu(t)=mu_k, k = calendar-year segment (eq. 5.4) alpha = 0
##    M1     kappa = 0, w_j = 1,     mu(t)=mu                (lambda^(1))
##    M2     kappa = 1, w_j = Z_j,   mu(t)=mu                (lambda^(2))
##    M1_tv  kappa = 0, annual piecewise-constant mu(t)
##    M2_tv  kappa = 1, annual piecewise-constant mu(t)
##
##  Parametrisation (kept identical to referee_timevarying_background.R so that
##  the numbers reproduce Tables 3/4 exactly):
##    mu_k   = exp(gamma_k),  beta = exp(theta_b),
##    eta    = 0.999 * plogis(theta_e),                 (subcriticality, eta<1)
##    alpha  = eta * beta / E[w],   w = (Z/Zbar)^kappa if kappa is supplied,
##                                   w = Z (kappa=1) or 1 (kappa=0).
##    => eta = alpha*E[w]/beta is the branching ratio, NOT "stationarity".
##  logL = sum_i log lambda(T_i) - int_0^T mu(s) ds
##         - (alpha/beta) * sum_j w_j (1 - exp(-beta(T-T_j))).
##  Annual backgrounds are integrated analytically (sum_k mu_k * L_k), so the
##  compensator identity Lambda(T_tr) = N_tr can be checked at machine precision.
##  Optimisation: 6 starting values, BFGS + Nelder-Mead refinement, best value kept.
##
##  Outputs (all written to DIR):
##    table34_model_estimates.csv   parameter estimates, SE, logL/AIC/BIC/LR
##    table34_annual_mu.csv         annual background rates mu_k per segment
##    table34_estimate_check.csv    reproduction vs existing result tables
##    table34_kappa_profile.csv     optional profile likelihood over kappa>=0
##  No figures are written.  Set CUT / STUDY_END to change the estimation window.
## =====================================================================

## ------------------------------------------------------------------ 0. 设置
DIR <- "E:/paper/rdata/hackexp/era_hawkes"
if (!dir.exists(DIR)) {
  of <- tryCatch(normalizePath(sys.frame(1L)$ofile), error = function(e) NA_character_)
  DIR <- if (!is.na(of)) dirname(of) else getwd()
}
OUT <- DIR

CUT                <- as.Date("2016-01-01")   # 训练段右端点(不含)
STUDY_END          <- as.Date("2018-12-30")   # 研究终点, 只用于切分记录
NSTART             <- 6L                       # 多起点数(与表 3/4 一致)
COMPUTE_SE         <- TRUE                      # 观测信息阵 + delta method 标准误
RUN_KAPPA_PROFILE  <- Sys.getenv("T34_KAPPA", "0") == "1"  # 剖面似然; T34_KAPPA=1 开启
KAPPA_GRID         <- seq(0, 2, by = 0.25)

FILES <- c(MED = "hackmed_daily.csv", BUS = "hackbus_daily.csv")
LABEL <- c(MED = "HACK-MED",          BUS = "HACK-BUS")
SPECS <- c("P", "Pa", "M1", "M2", "M1_tv", "M2_tv")

## ------------------------------------------------------------- 1. 数据与切分
load_full <- function(file, label) {
  d <- read.csv(file.path(DIR, file), stringsAsFactors = FALSE)
  d$Date <- as.Date(d$Date)
  stopifnot(!anyDuplicated(d$Date), all(d$Z > 0))
  origin <- min(d$Date)
  list(label = label, file = file, origin = origin,
       t = as.numeric(d$Date - origin), Z = d$Z, N = nrow(d),
       Date = d$Date, year = as.integer(format(d$Date, "%Y")),
       T_event = max(as.numeric(d$Date - origin)),
       T_end = as.numeric(STUDY_END - origin))
}

make_split <- function(full, cut) {
  T_tr <- as.numeric(as.Date(cut) - full$origin)
  itr  <- which(full$t < T_tr)
  ite  <- which(full$t >= T_tr & full$t <= full$T_end)
  xtr  <- list(label = full$label, t = full$t[itr], Z = full$Z[itr],
               N = length(itr), T = T_tr,
               Date = full$origin + full$t[itr], year = full$year[itr])
  list(label = full$label, origin = full$origin, xtr = xtr,
       T_tr = T_tr, T_end = full$T_end, cut = as.Date(cut),
       n_tr = length(itr), n_te = length(ite), te_Z = full$Z[ite])
}

## ------------------------------------------------------------- 2. 背景设计
year_bounds <- function(x) {
  yrs <- min(x$year):max(x$year)
  b <- as.numeric(as.Date(paste0(yrs, "-01-01")) - min(x$Date))
  sort(unique(c(0, b[b > 0 & b < x$T], x$T)))
}

make_bg <- function(x, type = c("const", "piece9")) {
  type <- match.arg(type)
  if (type == "const")
    return(list(B = matrix(1, x$N, 1L), L_k = x$T, K = 1L,
                type = type, bounds = c(0, x$T)))
  bounds <- year_bounds(x)
  K <- length(bounds) - 1L
  k <- findInterval(x$t, bounds, left.open = TRUE)
  k <- pmin(pmax(k, 1L), K)
  B <- matrix(0, x$N, K)
  B[cbind(seq_len(x$N), k)] <- 1
  list(B = B, L_k = diff(bounds), K = K, type = type, bounds = bounds)
}

## ------------------------------------------------------------- 3. 似然与拟合
mark_weights <- function(x, mark = FALSE, kappa = NULL) {
  if (!is.null(kappa)) {
    w <- (x$Z / mean(x$Z))^kappa
    return(list(w = w, mw = mean(w)))
  }
  if (mark) return(list(w = x$Z, mw = mean(x$Z)))
  list(w = rep(1, x$N), mw = 1)
}

make_ll <- function(x, bg, mark = FALSE, kappa = NULL) {
  t <- x$t; T <- x$T; N <- x$N; B <- bg$B; K <- bg$K
  mw0 <- mark_weights(x, mark, kappa)
  w <- mw0$w; mw <- mw0$mw
  function(p) {
    gam  <- p[seq_len(K)]
    beta <- exp(p[K + 1L])
    eta  <- 0.999 * plogis(p[K + 2L])
    if (!is.finite(beta) || beta > 1e4 || beta < 1e-6) return(1e10)
    mu_i <- as.vector(exp(B %*% gam))
    if (any(!is.finite(mu_i)) || any(mu_i <= 0)) return(1e10)
    bg_int <- sum(exp(gam) * bg$L_k)
    if (!is.finite(bg_int)) return(1e10)
    alpha <- eta * beta / mw
    R <- numeric(N)
    if (N > 1L) for (i in 2:N)
      R[i] <- exp(-beta * (t[i] - t[i - 1L])) * (R[i - 1L] + w[i - 1L])
    lam <- mu_i + alpha * R
    if (any(lam <= 0) || any(!is.finite(lam))) return(1e10)
    LL <- sum(log(lam)) - bg_int -
          (alpha / beta) * sum(w * (1 - exp(-beta * (T - t))))
    if (!is.finite(LL)) 1e10 else -LL
  }
}

summarise_fit <- function(x, bg, p, mark = FALSE, kappa = NULL) {
  K <- bg$K; t <- x$t; T <- x$T; N <- x$N
  mw0 <- mark_weights(x, mark, kappa)
  w <- mw0$w; mw <- mw0$mw
  gam  <- p[seq_len(K)]
  beta <- exp(p[K + 1L])
  eta  <- 0.999 * plogis(p[K + 2L])
  mu_k <- exp(gam)
  bg_int <- sum(mu_k * bg$L_k)
  alpha  <- eta * beta / mw
  R <- numeric(N)
  if (N > 1L) for (i in 2:N)
    R[i] <- exp(-beta * (t[i] - t[i - 1L])) * (R[i - 1L] + w[i - 1L])
  lam  <- as.vector(exp(bg$B %*% gam)) + alpha * R
  Lam  <- bg_int + (alpha / beta) * sum(w * (1 - exp(-beta * (T - t))))
  logL <- sum(log(lam)) - Lam
  npar <- K + 2L
  list(mu_k = mu_k, mu_last = tail(mu_k, 1L),
       mu_const = if (K == 1L) mu_k[1L] else NA_real_,
       mu_avg = bg_int / T, alpha = alpha, beta = beta, eta = eta,
       logL = logL, AIC = -2 * logL + 2 * npar,
       BIC = -2 * logL + log(N) * npar, npar = npar,
       Lambda = Lam, Lam_over_N = Lam / N, resid = abs(Lam - N) / N,
       bg_int = bg_int, memory_days = 1 / beta, par = p)
}

fit_hawkes <- function(x, bg, mark = FALSE, kappa = NULL,
                       nstart = NSTART, refine = TRUE) {
  fn <- make_ll(x, bg, mark, kappa)
  K  <- bg$K
  mu0 <- x$N / x$T
  ## 低/高/中 eta 与不同 beta 起点交错, 覆盖近临界脊
  hv <- c(0.05, 0.85, 0.35, 0.05, 0.85, 0.35)
  bv <- c(0.20, 0.05, 0.50, 0.80, 0.30, 0.10)
  starts <- lapply(seq_len(min(nstart, 6L)), function(j)
    c(rep(log(mu0), K), log(bv[j]), qlogis(hv[j] / 0.999)))
  best <- NULL; bestv <- Inf
  for (p0 in starts) {
    o <- try(optim(p0, fn, method = "BFGS",
                   control = list(maxit = 800L, reltol = 1e-10)), silent = TRUE)
    if (!inherits(o, "try-error") && is.finite(o$value) && o$value < bestv) {
      bestv <- o$value; best <- o$par
    }
    if (refine) {
      o2 <- try(optim(p0, fn, method = "Nelder-Mead",
                      control = list(maxit = 1500L, reltol = 1e-10)), silent = TRUE)
      if (!inherits(o2, "try-error") && is.finite(o2$value) && o2$value < bestv) {
        bestv <- o2$value; best <- o2$par
      }
    }
  }
  s <- summarise_fit(x, bg, best, mark, kappa)
  s$fn <- fn; s$value <- bestv; s$bg <- bg
  s$mark <- mark; s$kappa <- kappa
  s
}

## 同背景的纯 Poisson MLE: 常数背景 mu=N/T; 分段背景 mu_k=n_k/L_k
fit_bg_only <- function(x, bg) {
  K <- bg$K
  if (K == 1L) {
    mu <- x$N / x$T
    logL <- x$N * log(mu) - mu * x$T
    return(list(mu_k = mu, n_k = x$N, L_k = x$T, logL = logL, npar = 1L,
                AIC = -2 * logL + 2, BIC = -2 * logL + log(x$N),
                Lambda = mu * x$T))
  }
  kk <- pmin(pmax(findInterval(x$t, bg$bounds, left.open = TRUE), 1L), K)
  n_k <- tabulate(kk, nbins = K); L_k <- diff(bg$bounds)
  mu_k <- (n_k + 1e-9) / L_k
  logL <- sum(n_k * log(mu_k)) - sum(mu_k * L_k)
  list(mu_k = mu_k, n_k = n_k, L_k = L_k, logL = logL, npar = K,
       AIC = -2 * logL + 2 * K, BIC = -2 * logL + log(x$N) * K,
       Lambda = sum(mu_k * L_k))
}

## 观测信息阵 + delta method: SE(mu_last, alpha, beta, eta)
delta_se <- function(fit, x) {
  p <- fit$par; K <- fit$bg$K
  mw <- mark_weights(x, fit$mark, fit$kappa)$mw
  trans <- function(q) {
    gam <- q[seq_len(K)]
    beta <- exp(q[K + 1L]); eta <- 0.999 * plogis(q[K + 2L])
    c(mu_last = exp(gam[K]), alpha = eta * beta / mw,
      beta = beta, eta = eta)
  }
  H <- tryCatch(optimHess(p, fit$fn), error = function(e) NULL)
  V <- if (!is.null(H))
    tryCatch(solve(H), error = function(e) matrix(NA_real_, length(p), length(p)))
    else matrix(NA_real_, length(p), length(p))
  h <- 1e-5
  J <- sapply(seq_along(p), function(j) {
    q1 <- p; q0 <- p; q1[j] <- q1[j] + h; q0[j] <- q0[j] - h
    (trans(q1) - trans(q0)) / (2 * h)
  })
  v <- tryCatch(diag(J %*% V %*% t(J)), error = function(e) rep(NA_real_, 4L))
  setNames(sqrt(pmax(v, 0)), names(trans(p)))
}

## --------------------------------------------------------------- 4. 主流程
SER <- lapply(names(FILES), function(nm) load_full(FILES[[nm]], LABEL[[nm]]))
names(SER) <- names(FILES)
SP  <- lapply(SER, make_split, cut = CUT)

ROWS <- list(); ANNUAL <- list(); FITS <- list(); PROF <- list()

for (nm in names(SP)) {
  sp <- SP[[nm]]; x <- sp$xtr

  for (md in SPECS) {
    if (md %in% c("P", "Pa")) {
      bg <- make_bg(x, if (md == "P") "const" else "piece9")
      b0 <- fit_bg_only(x, bg)
      r <- list(mu_k = b0$mu_k, mu_last = tail(b0$mu_k, 1L),
                mu_const = if (bg$K == 1L) b0$mu_k[1L] else NA_real_,
                mu_avg = sum(b0$mu_k * b0$L_k) / x$T,
                alpha = 0, beta = NA_real_, eta = 0,
                memory_days = NA_real_, logL = b0$logL,
                AIC = b0$AIC, BIC = b0$BIC, npar = b0$npar,
                Lambda = b0$Lambda, Lam_over_N = b0$Lambda / x$N,
                resid = abs(b0$Lambda - x$N) / x$N,
                bg_int = b0$Lambda, logL_bg = b0$logL, AIC_bg = b0$AIC,
                LR_vs_bg_only = 0)
      se <- c(mu_last = if (bg$K == 1L) sqrt(x$N) / x$T else
                          sqrt(tail(b0$n_k, 1L)) / tail(b0$L_k, 1L),
              alpha = 0, beta = NA_real_, eta = 0)
      if (md == "Pa") {
        ANNUAL[[length(ANNUAL) + 1L]] <- data.frame(
          dataset = sp$label, model = md, segment = seq_len(bg$K),
          start_date = as.character(sp$origin + bg$bounds[-length(bg$bounds)]),
          end_date   = as.character(sp$origin + bg$bounds[-1L] - 1L),
          L_k = diff(bg$bounds), n_k = b0$n_k,
          mu_k = b0$mu_k, mu_k_poisson = b0$mu_k,
          stringsAsFactors = FALSE)
      }
    } else {
      ty <- if (grepl("_tv", md)) "piece9" else "const"
      bg <- make_bg(x, ty)
      b0 <- fit_bg_only(x, bg)
      f  <- fit_hawkes(x, bg, mark = grepl("M2", md))
      FITS[[paste(nm, md)]] <- f
      r <- list(mu_k = f$mu_k, mu_last = f$mu_last,
                mu_const = f$mu_const, mu_avg = f$mu_avg,
                alpha = f$alpha, beta = f$beta, eta = f$eta,
                memory_days = f$memory_days, logL = f$logL,
                AIC = f$AIC, BIC = f$BIC, npar = f$npar,
                Lambda = f$Lambda, Lam_over_N = f$Lam_over_N,
                resid = f$resid, bg_int = f$bg_int,
                logL_bg = b0$logL, AIC_bg = b0$AIC,
                LR_vs_bg_only = 2 * (f$logL - b0$logL))
      se <- if (COMPUTE_SE) delta_se(f, x) else
        c(mu_last = NA_real_, alpha = NA_real_, beta = NA_real_, eta = NA_real_)
      ## 按年背景明细(仅 M1_tv / M2_tv 与分段 Poisson)
      if (ty == "piece9") {
        ANNUAL[[length(ANNUAL) + 1L]] <- data.frame(
          dataset = sp$label, model = md, segment = seq_len(bg$K),
          start_date = as.character(sp$origin + bg$bounds[-length(bg$bounds)]),
          end_date   = as.character(sp$origin + bg$bounds[-1L] - 1L),
          L_k = diff(bg$bounds),
          n_k = if (!is.null(b0$n_k)) b0$n_k else NA_integer_,
          mu_k = f$mu_k, mu_k_poisson = b0$mu_k,
          stringsAsFactors = FALSE)
      }
    }

    ## 表 3/4 中展示的 mu: P/M1/M2 用 mu_hat, 按年背景用窗口均值 mu_avg
    mu_table <- if (md %in% c("P", "M1", "M2")) r$mu_const else
                if (grepl("_tv", md)) r$mu_avg else NA_real_
    ROWS[[length(ROWS) + 1L]] <- data.frame(
      dataset = sp$label, model = md,
      kappa = if (grepl("M2", md)) 1 else 0,
      background = if (md %in% c("P", "M1", "M2")) "const" else "annual",
      K_bg = if (md %in% c("P", "M1", "M2")) 1L else length(r$mu_k),
      cut = as.character(CUT), n_train = x$N, T_train = x$T,
      mu_const = r$mu_const, mu_last = r$mu_last,
      mu_avg = r$mu_avg, mu_table = mu_table,
      alpha = r$alpha, beta = r$beta, eta = r$eta,
      memory_days = r$memory_days,
      se_mu = unname(se["mu_last"]), se_alpha = unname(se["alpha"]),
      se_beta = unname(se["beta"]), se_eta = unname(se["eta"]),
      logL = r$logL, npar = r$npar, AIC = r$AIC, BIC = r$BIC,
      logL_bg = r$logL_bg, AIC_bg = r$AIC_bg,
      LR_vs_bg_only = r$LR_vs_bg_only,
      Lambda = r$Lambda, Lambda_over_N = r$Lam_over_N,
      compensator_resid = r$resid,
      stringsAsFactors = FALSE)
  }
}

EST <- do.call(rbind, ROWS)
## dAIC(kappa=1 - kappa=0): 同一背景下两个标记规格的拟合差
## Pa 仅为 bootstrap/预测用的分段 Poisson DGP, 论文表中其 AIC 记为 "-",
## 此处仍计算并输出, 但不参与规格选择.
EST$dAIC_kappa1_vs_0 <- NA_real_
for (nm in names(SP)) for (ty in c("const", "annual")) {
  a <- if (ty == "const") c("M1", "M2") else c("M1_tv", "M2_tv")
  i1 <- which(EST$dataset == LABEL[[nm]] & EST$model == a[1L])
  i2 <- which(EST$dataset == LABEL[[nm]] & EST$model == a[2L])
  EST$dAIC_kappa1_vs_0[c(i1, i2)] <- EST$AIC[i2] - EST$AIC[i1]
}
if (length(ANNUAL)) ANNTAB <- do.call(rbind, ANNUAL)

## 可选: kappa>=0 的剖面似然(常数背景; 用于 §5.2 的单边 score / LR 诊断)
## 注意: 自由估计落入边界 kappa=0 时只能写 absence of evidence;
##       截尾/负值不能解释为「负的损失依赖」。
if (RUN_KAPPA_PROFILE) {
  cat("[kappa] 剖面似然 (常数背景)...\n")
  for (nm in names(SP)) {
    x <- SP[[nm]]$xtr; bg <- make_bg(x, "const")
    for (kap in KAPPA_GRID) {
      f <- fit_hawkes(x, bg, mark = FALSE, kappa = kap,
                      nstart = 2L, refine = FALSE)
      PROF[[length(PROF) + 1L]] <- data.frame(
        dataset = LABEL[[nm]], kappa = kap, alpha = f$alpha,
        beta = f$beta, eta = f$eta, mu = f$mu_const,
        logL = f$logL, AIC = f$AIC, stringsAsFactors = FALSE)
    }
  }
  PT <- do.call(rbind, PROF)
  for (ds in unique(PT$dataset)) {
    s <- PT[PT$dataset == ds, ]
    l0 <- s$logL[s$kappa == 0]
    cat(sprintf("  %s: kappa_hat = %.2f, LR_+ = %.3f\n", ds,
                s$kappa[which.max(s$logL)], 2 * (max(s$logL) - l0)))
  }
  write.csv(PT, file.path(OUT, "table34_kappa_profile.csv"), row.names = FALSE)
}

## --------------------------------------------------------------- 5. 数值对拍
CK <- list()
ck_row <- function(lab, md, q, ref, val)
  data.frame(dataset = lab, model = md, quantity = q, reference = ref,
             reproduced = val, abs_diff = abs(ref - val),
             stringsAsFactors = FALSE)

ref_file <- file.path(DIR, "heldout_training_fits.csv")
if (file.exists(ref_file)) {
  REF <- read.csv(ref_file, stringsAsFactors = FALSE)
  REF <- REF[REF$is_primary & REF$cut == as.character(CUT), ]
  for (i in seq_len(nrow(REF))) {
    e <- EST[EST$dataset == REF$dataset[i] & EST$model == REF$model[i], ]
    if (!nrow(e)) next
    for (q in c("mu", "alpha", "beta", "eta", "AIC_train", "logL_train")) {
      rv <- switch(q, mu = e$mu_last, alpha = e$alpha, beta = e$beta,
                   eta = e$eta, AIC_train = e$AIC, logL_train = e$logL)
      if (is.na(REF[[q]][i]) || is.na(rv)) next
      CK[[length(CK) + 1L]] <- ck_row(REF$dataset[i], REF$model[i], q, REF[[q]][i], rv)
    }
  }
}
s52_file <- file.path(DIR, "s52_background_training.csv")
if (file.exists(s52_file)) {
  S52 <- read.csv(s52_file, stringsAsFactors = FALSE)
  for (nm in names(SP)) for (md in c("M1_tv", "M2_tv")) {
    e <- EST[EST$dataset == LABEL[[nm]] & EST$model == md, ]
    s <- S52[S52$dataset == LABEL[[nm]] &
             S52$background == "annual" &
             S52$model == if (md == "M1_tv") "M1 (lambda1)" else "M2 (lambda2)", ]
    if (nrow(e) && nrow(s)) CK[[length(CK) + 1L]] <- ck_row(LABEL[[nm]], md, "s52_mu_avg",
                                   s$mu_avg[1L], e$mu_avg)
  }
}
## 截图表 3/4 的展示值(仅用于核对, 不参与拟合)
SCREEN <- data.frame(
  dataset = rep(c("HACK-MED", "HACK-BUS"), each = 5L),
  model   = rep(c("P", "M1", "M2", "M1_tv", "M2_tv"), 2L),
  mu_shown = c(0.05772, 0.01417, 0.01310, 0.03986, 0.04074,
               0.09995, 0.04166, 0.04020, 0.07119, 0.07406),
  alpha = c(NA, 0.017611, 0.001291, 0.004112, 0.000280,
            NA, 0.013024, 0.000837, 0.009584, 0.000530),
  beta  = c(NA, 0.023037, 0.021591, 0.013049, 0.012223,
            NA, 0.021919, 0.020384, 0.032902, 0.030002),
  eta   = c(0, 0.7645, 0.7845, 0.3151, 0.3006,
            0, 0.5942, 0.6105, 0.2913, 0.2627),
  AIC   = c(949.6, 857.2, 862.6, 844.3, 845.4,
            1448.7, 1420.8, 1423.5, 1425.5, 1426.5),
  stringsAsFactors = FALSE)
for (i in seq_len(nrow(SCREEN))) {
  e <- EST[EST$dataset == SCREEN$dataset[i] & EST$model == SCREEN$model[i], ]
  if (!nrow(e)) next
  vals <- list(mu_shown = e$mu_table, alpha = e$alpha, beta = e$beta,
               eta = e$eta, AIC = e$AIC)
  for (q in names(vals))
    if (!is.na(SCREEN[[q]][i]) && !is.na(vals[[q]]))
      CK[[length(CK) + 1L]] <- ck_row(SCREEN$dataset[i], SCREEN$model[i], paste0("screen_", q),
             SCREEN[[q]][i], vals[[q]])
}
CK <- if (length(CK)) do.call(rbind, CK) else
  data.frame(dataset = character(), model = character(), quantity = character(),
             reference = numeric(), reproduced = numeric(), abs_diff = numeric())

if (nrow(CK)) {
  exact <- CK[!grepl("^screen_", CK$quantity), ]
  cat(sprintf("[check] 与既有结果表逐项对拍: %d 项, max |diff| = %.3g\n",
              nrow(exact), max(exact$abs_diff, na.rm = TRUE)))
  cat(sprintf("[check] 与截图展示值(含四舍五入)对拍: %d 项, max |diff| = %.3g\n",
              sum(grepl("^screen_", CK$quantity)),
              max(CK$abs_diff[grepl("^screen_", CK$quantity)], na.rm = TRUE)))
}

## --------------------------------------------------------------- 6. 输出
write.csv(EST, file.path(OUT, "table34_model_estimates.csv"), row.names = FALSE)
if (exists("ANNTAB"))
  write.csv(ANNTAB, file.path(OUT, "table34_annual_mu.csv"), row.names = FALSE)
write.csv(CK,  file.path(OUT, "table34_estimate_check.csv"), row.names = FALSE)

cat("\n=== Training-window MLE (cut = ", as.character(CUT), ") ===\n", sep = "")
print(EST[, c("dataset", "model", "mu_table", "alpha", "beta", "eta",
              "AIC", "logL", "LR_vs_bg_only")],
      row.names = FALSE, digits = 5)
cat("\n=== Standard errors (observed information, delta method) ===\n")
print(EST[, c("dataset", "model", "se_mu", "se_alpha", "se_beta", "se_eta")],
      row.names = FALSE, digits = 4)
cat("\nCompensator check, Hawkes fits: max |Lambda(T_tr)/N_tr - 1| = ",
    sprintf("%.3g", max(EST$compensator_resid[EST$alpha > 0], na.rm = TRUE)), "\n")
cat("Outputs: table34_model_estimates.csv, table34_annual_mu.csv, ",
    "table34_estimate_check.csv", if (RUN_KAPPA_PROFILE)
      ", table34_kappa_profile.csv" else "", "\n", sep = "")

invisible(list(estimates = EST, annual = if (exists("ANNTAB")) ANNTAB else NULL,
               check = CK, splits = SP))
