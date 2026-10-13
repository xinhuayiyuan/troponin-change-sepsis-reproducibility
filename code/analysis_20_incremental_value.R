# =============================================================================
# Incremental prognostic value and
# quantitative bias analysis (all computed on the real landmark dataset).
#   - reference model  : SAME as main model M3a (log_baseline + n_meas +
#                        sofa72 + vaso72 + admission covariates)
#   - + rate (linear)  / + rate (restricted cubic spline, 4 df)
#   - Harrell C and Uno C, delta-C with 500-replicate bootstrap CI
#   - calibration slope (apparent; the optimism-corrected version is in analysis_30)
# Output: output/lock/landmark/*.csv
# Reference model aligned to the M3a specification.
# =============================================================================
suppressPackageStartupMessages({
  library(survival); library(splines); library(pec); library(dcurves); library(boot)
})
set.seed(20261004)
BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT <- file.path(BASE, 'output', 'lock', 'landmark')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(BASE, 'data', 'landmark_analysis_dataset.csv'), stringsAsFactors = FALSE)
for (bc in c('alive_at_72', 'in_hosp_at_72', 'landmark')) d[[bc]] <- as.logical(d[[bc]])
d <- d[which(d$landmark), ]
d$time_days <- d$time_landmark / 24
d$female <- as.integer(d$gender == 'F')
d$log_baseline <- log(d$first_tn)
q <- quantile(d$r, c(0.01, 0.99)); d$r_w <- pmin(pmax(d$r, q[1]), q[2])

cov_adm <- c('age', 'female', 'sofa_24h', 'invasive_vent', 'norepinephrine',
             'creatinine_max', 'lactate_max')
cov_lm  <- c('sofa72', 'vaso72')
need <- c('time_days','death28','r_w','log_baseline','n_meas', cov_lm, cov_adm)
dd <- d[complete.cases(d[, need]), ]
TAU <- 24
cat('incremental-value dataset n =', nrow(dd), ' deaths =', sum(dd$death28), '\n\n')

# ---- reference model = M3a without rate -------------------------------------
f0 <- as.formula(paste('Surv(time_days, death28) ~ log_baseline + n_meas +',
                       paste(c(cov_lm, cov_adm), collapse = ' + ')))
f1 <- update(f0, . ~ . + r_w)
f2 <- update(f0, . ~ . + ns(r_w, 4))
m0 <- coxph(f0, data = dd, x = TRUE, y = TRUE)
m1 <- coxph(f1, data = dd, x = TRUE, y = TRUE)
m2 <- coxph(f2, data = dd, x = TRUE, y = TRUE)

# --- Harrell C and LR tests --------------------------------------------------
c0 <- concordance(m0); c1 <- concordance(m1); c2 <- concordance(m2)
lr_lin <- 2 * (m1$loglik[2] - m0$loglik[2])
lr_spl <- 2 * (m2$loglik[2] - m0$loglik[2])
tab <- data.frame(
  model = c('reference (M3a covariates, no rate)', '+ rate (linear)', '+ rate (4-df spline)'),
  k = c(length(coef(m0)), length(coef(m1)), length(coef(m2))),
  Harrell_C = c(c0$concordance, c1$concordance, c2$concordance),
  C_lower = c(c0$concordance - 1.96 * sqrt(c0$var), c1$concordance - 1.96 * sqrt(c1$var),
              c2$concordance - 1.96 * sqrt(c2$var)),
  C_upper = c(c0$concordance + 1.96 * sqrt(c0$var), c1$concordance + 1.96 * sqrt(c1$var),
              c2$concordance + 1.96 * sqrt(c2$var)),
  LR_vs_ref_chi2 = c(NA, lr_lin, lr_spl),
  LR_df = c(NA, 1, 4),
  LR_p = c(NA, pchisq(lr_lin, 1, lower.tail = FALSE), pchisq(lr_spl, 4, lower.tail = FALSE)))
tab[, 3:8] <- lapply(tab[, 3:8], function(x) signif(as.numeric(x), 6))
print(tab, row.names = FALSE)
write.csv(tab, file.path(OUT, 'incremental_discrimination.csv'), row.names = FALSE)

# --- delta-C bootstrap (500 replicates) --------------------------------------
boot_dc <- function(dat, idx) {
  dd2 <- dat[idx, ]
  a <- coxph(f0, data = dd2); b <- coxph(f1, data = dd2); cc <- coxph(f2, data = dd2)
  co <- concordance(a)$concordance; c1b <- concordance(b)$concordance
  c2b <- concordance(cc)$concordance
  c(dC_lin = c1b - co, dC_spline = c2b - co)
}
bs <- boot(dd, boot_dc, R = 500)
ci_lin <- boot.ci(bs, index = 1, type = 'perc')
ci_spl <- boot.ci(bs, index = 2, type = 'perc')
dc <- data.frame(
  contrast = c('delta-C adding rate (linear)', 'delta-C adding rate (4-df spline)'),
  delta_C = c(c1$concordance - c0$concordance, c2$concordance - c0$concordance),
  boot_lower = c(ci_lin$percent[4], ci_spl$percent[4]),
  boot_upper = c(ci_lin$percent[5], ci_spl$percent[5]))
dc[, 2:4] <- lapply(dc[, 2:4], function(x) signif(as.numeric(x), 5))
cat('\n=== delta-C (500 bootstrap replicates) ===\n'); print(dc, row.names = FALSE)
write.csv(dc, file.path(OUT, 'incremental_deltaC.csv'), row.names = FALSE)

# --- calibration slope and Brier score at 24 days ----------------------------
p0 <- as.numeric(predictSurvProb(m0, newdata = dd, times = TAU))
p1 <- as.numeric(predictSurvProb(m1, newdata = dd, times = TAU))
r0 <- 1 - p0; r1 <- 1 - p1
lp0 <- log(r0 / (1 - r0)); lp1 <- log(r1 / (1 - r1))
cal_slope <- function(y, lp) {
  mm <- glm(y ~ lp, family = binomial)
  unname(coef(mm)[2])
}
library(riskRegression)
sc <- Score(list('reference' = m0, 'reference+rate' = m1), formula = Surv(time_days, death28) ~ 1,
            data = dd, times = TAU, metrics = c('AUC', 'Brier'), conf.int = TRUE)
cat('\n=== Uno C (IPCW) and Brier at 24 days ===\n'); print(sc$AUC$score); print(sc$Brier$score)
write.csv(as.data.frame(sc$AUC$score), file.path(OUT, 'incremental_unoC.csv'), row.names = FALSE)
write.csv(as.data.frame(sc$Brier$score), file.path(OUT, 'incremental_brier.csv'), row.names = FALSE)
ca <- data.frame(model = c('reference', 'reference+rate'),
                 calibration_slope = c(cal_slope(dd$death28, lp0), cal_slope(dd$death28, lp1)))
print(ca, row.names = FALSE)
write.csv(ca, file.path(OUT, 'incremental_calibration_slope.csv'), row.names = FALSE)

# --- NOTE: the bootstrap optimism-corrected calibration slope is
# computed once, in analysis_30 (section V2), together with the calibration curve.
# It is computed from a single implementation, so that one quantity is reported
# once and with one set of Monte Carlo draws.

# --- decision curve analysis -------------------------------------------------
dcs <- dca(Surv(time_days, death28) ~ r0 + r1, data = dd, time = TAU,
           thresholds = seq(0.05, 0.60, by = 0.05))
nb <- as.data.frame(dcs$dca)
write.csv(nb, file.path(OUT, 'incremental_dca.csv'), row.names = FALSE)
cat('\n=== DCA net benefit (selected thresholds) ===\n')
print(nb[nb$threshold %in% c(0.10, 0.20, 0.30, 0.40), ], row.names = FALSE)

cat('\nDONE. outputs in', OUT, '\n')
