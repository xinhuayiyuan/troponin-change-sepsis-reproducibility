# =============================================================================
# Remaining quantitative items
#   V1 optimism-corrected Harrell C (bootstrap)
#   V2 calibration curve at 24 days (apparent + bootstrap optimism-corrected)
#   V3 NRI / IDI at 24 days (naive, uncensored subset only; documented)
#   B1 measurement-intensity IPW (multinomial propensity for the measurement count)
#   B2 tipping-point: confounding strength needed to shift the estimate to the MCID
#   B3 regression-dilution correction across a grid of exposure reliabilities
#      (a literature intra-individual CV is NOT assumed here: the reliability is
#       varied explicitly and the required CV is reported instead)
#   D1 dose-response curve data for the rate (dual-spline model, reference r = 0)
#
# Input : data/landmark_analysis_dataset.csv
# Output: output/lock/sensitivity/*.csv
# Contents:
#   - V2 includes the bootstrap optimism-corrected calibration slope
#   - V3 NRI/IDI (naive, uncensored subset, no IPCW)
#   - B2 (HR > RR for non-rare outcomes => conservative E-value)
# =============================================================================
suppressPackageStartupMessages({
  library(survival); library(splines); library(pec); library(riskRegression)
})
set.seed(20261005)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT  <- file.path(BASE, 'output', 'lock', 'sensitivity')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(BASE, 'data', 'landmark_analysis_dataset.csv'), stringsAsFactors = FALSE)
for (bc in c('alive_at_72', 'in_hosp_at_72', 'landmark')) d[[bc]] <- as.logical(d[[bc]])
d <- d[which(d$landmark), ]
d$time_days <- d$time_landmark / 24
d$female <- as.integer(d$gender == 'F')
d$log_baseline <- log(d$first_tn)
q <- quantile(d$r, c(0.01, 0.99), na.rm = TRUE)
d$r_w <- pmin(pmax(d$r, q[1]), q[2])

cov_adm <- c('age', 'female', 'sofa_24h', 'invasive_vent', 'norepinephrine',
             'creatinine_max', 'lactate_max')
cov_lm  <- c('sofa72', 'vaso72')
need <- c('time_days', 'death28', 'r_w', 'n_meas', 'log_baseline', cov_adm, cov_lm)
dd <- d[complete.cases(d[, need]), ]
cat('analysis n =', nrow(dd), ' deaths =', sum(dd$death28), '\n\n')

f0 <- Surv(time_days, death28) ~ n_meas + log_baseline + age + female +
  sofa_24h + invasive_vent + norepinephrine + creatinine_max + lactate_max + sofa72 + vaso72
f1 <- update(f0, . ~ . + r_w)
m0 <- coxph(f0, data = dd, x = TRUE, y = TRUE)
m1 <- coxph(f1, data = dd, x = TRUE, y = TRUE)

# ---- V1: optimism-corrected Harrell C ---------------------------------------
B <- 200
c_app <- concordance(m1)$concordance
opt <- numeric(B)
for (b in 1:B) {
  idx <- sample(nrow(dd), replace = TRUE)
  mb <- coxph(f1, data = dd[idx, ], x = TRUE, y = TRUE)
  c_in   <- concordance(mb)$concordance
  c_test <- concordance(mb, newdata = dd)$concordance
  opt[b] <- c_in - c_test
}
c_corr <- c_app - mean(opt, na.rm = TRUE)
c_opt_se <- sd(opt, na.rm = TRUE) / sqrt(B)
cat(sprintf('V1 Harrell C apparent = %.4f | optimism = %.4f (MC SE %.4f) | corrected = %.4f (B = %d)\n',
            c_app, mean(opt, na.rm = TRUE), c_opt_se, c_corr, B))
write.csv(data.frame(model = '+rate (full covariate set)', B = B,
                     apparent_C = c_app, optimism = mean(opt, na.rm = TRUE),
                     optimism_mc_se = c_opt_se, corrected_C = c_corr),
          file.path(OUT, 'val_optimism_cindex.csv'), row.names = FALSE)

# ---- V2: calibration curve at 24 days (apparent + optimism-corrected) -------
TAU <- 24
p1v <- as.numeric(predictSurvProb(m1, newdata = dd, times = TAU))
risk1 <- 1 - p1v
brk <- unique(quantile(risk1, probs = seq(0, 1, 0.1), na.rm = TRUE))
grp <- cut(risk1, breaks = brk, include.lowest = TRUE)
apparent <- do.call(rbind, tapply(seq_len(nrow(dd)), grp, function(i) {
  s <- survfit(Surv(time_days, death28) ~ 1, data = dd[i, ])
  pred_mean <- mean(risk1[i])
  obs <- 1 - summary(s, times = TAU)$surv
  data.frame(decile = as.character(grp[i][1]), pred = pred_mean, obs = obs, n = length(i))
}, simplify = FALSE))
write.csv(apparent, file.path(OUT, 'val_calibration_curve.csv'), row.names = FALSE)
cat('\nV2 apparent calibration curve (deciles):\n'); print(apparent, row.names = FALSE)

# bootstrap optimism-corrected calibration slope
B2 <- 200
cal_slope <- function(y, lp) {
  mm <- glm(y ~ lp, family = binomial)
  unname(coef(mm)[2])
}
lp1 <- log(risk1 / (1 - risk1))
slope_app <- cal_slope(dd$death28, lp1)
slope_opt <- numeric(B2)
for(b in 1:B2) {
  idx <- sample(nrow(dd), replace = TRUE)
  mb <- coxph(f1, data = dd[idx, ], x = TRUE, y = TRUE)
  risk_b <- 1 - as.numeric(predictSurvProb(mb, newdata = dd[idx, ], times = TAU))
  lp_b   <- log(risk_b / (1 - risk_b))
  s_app  <- cal_slope(dd$death28[idx], lp_b)
  risk_t <- 1 - as.numeric(predictSurvProb(mb, newdata = dd, times = TAU))
  lp_t   <- log(risk_t / (1 - risk_t))
  s_test <- cal_slope(dd$death28, lp_t)
  slope_opt[b] <- s_app - s_test
}
opt_mean <- mean(slope_opt, na.rm = TRUE)
opt_se <- sd(slope_opt, na.rm = TRUE) / sqrt(B2)
slope_corr <- slope_app - opt_mean
# Monte Carlo uncertainty of the AVERAGED estimate (normal approximation);
# not the spread of individual replicate estimates.
ci_corr <- slope_corr + c(-1, 1) * 1.96 * opt_se
cat(sprintf('V2 optimism-corrected calibration slope: apparent=%.4f | optimism=%.4f (MC SE %.4f) | corrected=%.4f (95%% MC interval %.4f-%.4f, B=%d)\n',
            slope_app, opt_mean, opt_se, slope_corr, ci_corr[1], ci_corr[2], B2))
write.csv(data.frame(model = 'reference+rate', B = B2, apparent_slope = slope_app,
                     optimism = opt_mean, optimism_mc_se = opt_se,
                     corrected_slope = slope_corr,
                     corrected_mc_lo = ci_corr[1], corrected_mc_hi = ci_corr[2]),
          file.path(OUT, 'val_calibration_slope_optimism.csv'), row.names = FALSE)

# ---- V3: NRI / IDI at 24 days (naive, uncensored subset only) ---------------
# Known 24-day status: either died within 24d or followed >=24d (censored or alive)
known24 <- dd[dd$time_days >= 24 | dd$death28 == 1, ]
y24 <- known24$death28
p0_24 <- 1 - as.numeric(predictSurvProb(m0, newdata = known24, times = 24))
p1_24 <- 1 - as.numeric(predictSurvProb(m1, newdata = known24, times = 24))
# continuous NRI
nri_event   <- mean((p1_24 > p0_24)[y24 == 1]) - mean((p1_24 < p0_24)[y24 == 1])
nri_nonevent<- mean((p1_24 < p0_24)[y24 == 0]) - mean((p1_24 > p0_24)[y24 == 0])
nri_cont    <- nri_event + nri_nonevent
# IDI
idi <- (mean(p1_24[y24 == 1]) - mean(p0_24[y24 == 1])) -
       (mean(p1_24[y24 == 0]) - mean(p0_24[y24 == 0]))
cat(sprintf('\nV3 NRI/IDI (naive, 24-day known status): n_known=%d, events=%d | continuous NRI=%.4f | IDI=%.4e\n',
            nrow(known24), sum(y24), nri_cont, idi))
write.csv(data.frame(time_days = 24, n_known = nrow(known24), events = sum(y24),
                     continuous_NRI = nri_cont, IDI = idi,
                     note = 'naive (uncensored subset only; no IPCW)'),
          file.path(OUT, 'val_nri_idi.csv'), row.names = FALSE)

# ---- B1: measurement-intensity IPW ------------------------------------------
# The number of cTnT measurements is itself a clinical decision, so patients are
# weighted by the inverse probability of their observed measurement count
# (2, 3, 4, 5+) given baseline covariates, and the main model is refitted.
if (requireNamespace('nnet', quietly = TRUE)) {
  library(nnet)
  dd$kcat <- factor(pmin(dd$n_meas, 5))
  pw <- multinom(kcat ~ age + female + sofa_24h + log_baseline + creatinine_max,
                 data = dd, trace = FALSE)
  pr <- predict(pw, type = 'probs')
  dd$w_ipw <- 1 / pr[cbind(seq_len(nrow(dd)), as.integer(dd$kcat))]
  dd$w_ipw <- pmin(dd$w_ipw, quantile(dd$w_ipw, 0.99))     # stabilise the tail
  m_ipw <- coxph(f1, data = dd, weights = w_ipw, x = TRUE, y = TRUE)
  cf_ipw <- coef(m_ipw)['r_w']; se_ipw <- sqrt(vcov(m_ipw)['r_w', 'r_w'])
  cat(sprintf('
B1 measurement-intensity IPW: n = %d | HR(r_w) = %.4f (%.4f-%.4f), P = %.4f | weight range %.3f-%.3f
',
              m_ipw$n, exp(cf_ipw), exp(cf_ipw - 1.96 * se_ipw), exp(cf_ipw + 1.96 * se_ipw),
              2 * pnorm(-abs(cf_ipw / se_ipw)), min(dd$w_ipw), max(dd$w_ipw)))
  write.csv(data.frame(analysis = 'B1 measurement-intensity IPW',
                       n = m_ipw$n, events = m_ipw$nevent, term = 'r_w',
                       HR = exp(cf_ipw), lo = exp(cf_ipw - 1.96 * se_ipw),
                       hi = exp(cf_ipw + 1.96 * se_ipw),
                       p = 2 * pnorm(-abs(cf_ipw / se_ipw)),
                       w_min = min(dd$w_ipw), w_max = max(dd$w_ipw)),
            file.path(OUT, 'bias_ipw_measurement.csv'), row.names = FALSE)
}

# ---- B2: tipping-point analysis (confounding needed to reach the MCID) ------
hr_to_rr <- function(x) (1 - 0.5^sqrt(x)) / (1 - 0.5^sqrt(1 / x))
evalue_rr <- function(rr) { if (rr < 1) rr <- 1 / rr; rr + sqrt(rr * (rr - 1)) }
cf <- coef(m1)['r_w']; se <- sqrt(vcov(m1)['r_w', 'r_w'])
hr <- exp(cf)
sd_r <- sd(dd$r_w, na.rm = TRUE)
tp <- data.frame(
  target = c('observed HR per 1 SD', 'preset MCID HR per 1 SD = 1.20',
             'observed HR per 1 unit'),
  HR = c(hr^sd_r, 1.20, hr),
  evalue_point = c(evalue_rr(hr_to_rr(hr^sd_r)), evalue_rr(hr_to_rr(1.20)),
                   evalue_rr(hr_to_rr(hr))))
cat('\nB2 tipping-point (E-values):\n'); print(tp, row.names = FALSE)
write.csv(tp, file.path(OUT, 'bias_tipping_point.csv'), row.names = FALSE)
# NOTE: outcome prevalence is 26% (non-rare). For non-rare outcomes HR > RR;
# treating HR as RR yields a slightly conservative (larger) E-value.

# ---- B3: regression-dilution correction over a reliability grid -------------
lam <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60)
rd <- data.frame(reliability = lam,
                 HR_corrected = exp(cf / lam),
                 lo_corrected = exp((cf - 1.96 * se) / lam),
                 hi_corrected = exp((cf + 1.96 * se) / lam))
rd$HR_corrected_per_SD <- rd$HR_corrected^sd_r
cat('\nB3 regression-dilution grid (per 1 unit r):\n'); print(rd, row.names = FALSE)
write.csv(rd, file.path(OUT, 'bias_regression_dilution.csv'), row.names = FALSE)

# ---- D1: dose-response curve of the rate (dual-spline model) ----------------
ms <- coxph(as.formula(paste('Surv(time_days, death28) ~ ns(r_w, 4) + ns(log_baseline, 4) +',
                             'r_w:log_baseline + n_meas +',
                             paste(c(cov_adm, cov_lm), collapse = ' + '))),
            data = dd, x = TRUE, y = TRUE)
grid <- seq(quantile(dd$r_w, 0.025), quantile(dd$r_w, 0.975), length.out = 60)
nd0 <- dd[rep(1, length(grid)), ]
nd0$r_w <- 0
nd1 <- nd0; nd1$r_w <- grid
X <- model.matrix(delete.response(terms(ms)), data = rbind(nd0[1, ], nd1),
                  xlev = ms$xlevels)
b <- coef(ms); V <- vcov(ms)
X <- X[, names(b), drop = FALSE]
L <- X[-1, , drop = FALSE] - matrix(X[1, ], nrow = nrow(X) - 1, ncol = ncol(X),
                                    byrow = TRUE)
fit <- as.numeric(L %*% b)
se <- sqrt(rowSums((L %*% V) * L))
curve <- data.frame(r_w = grid,
                    logHR = fit,
                    HR = exp(fit),
                    HR_lo = exp(fit - 1.96 * se),
                    HR_hi = exp(fit + 1.96 * se))
write.csv(curve, file.path(OUT, 'dose_response_curve.csv'), row.names = FALSE)
cat(sprintf('\nD1 dose-response curve: %d grid points over r_w %.2f to %.2f (reference r = 0); HR range %.3f-%.3f\n',
            nrow(curve), min(grid), max(grid), min(curve$HR), max(curve$HR)))
# The anova() row-name regex was not robust (it printed NA).
# Recompute the joint LR for the rate block explicitly (spline + interaction)
# so the console value matches the locked landmark_rate_jointLR.csv.
f_norate <- as.formula(paste('Surv(time_days, death28) ~ ns(log_baseline, 4) + n_meas +',
                             paste(c(cov_adm, cov_lm), collapse = ' + ')))
m_norate <- coxph(f_norate, data = dd, x = TRUE, y = TRUE)
lr_joint <- 2 * (ms$loglik[2] - m_norate$loglik[2])
df_joint <- length(coef(ms)) - length(coef(m_norate))
p_joint <- pchisq(lr_joint, df_joint, lower.tail = FALSE)
cat(sprintf('   rate-block joint LR (spline + interaction) = %.4f, df = %d, P = %.4f\n',
            lr_joint, df_joint, p_joint))

cat('\nDONE. outputs in', OUT, '\n')
