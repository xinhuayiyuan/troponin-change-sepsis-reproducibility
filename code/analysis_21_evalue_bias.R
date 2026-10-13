# =============================================================================
# E-value (quantitative bias analysis) and IPCW discrimination.
# Output: output/lock/landmark/*.csv
# =============================================================================
suppressPackageStartupMessages({
  library(survival); library(splines); library(riskRegression); library(EValue)
})
BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT <- file.path(BASE, 'output', 'lock', 'landmark')

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

f0 <- as.formula(paste('Surv(time_days, death28) ~ log_baseline + n_meas +',
                       paste(c(cov_lm, cov_adm), collapse = ' + ')))
f1 <- update(f0, . ~ . + r_w)
m1 <- coxph(f1, data = dd, x = TRUE, y = TRUE)
m0 <- coxph(f0, data = dd, x = TRUE, y = TRUE)

# --- E-value for the fully adjusted rate HR -------------------------------
# Implemented directly from the published formula (VanderWeele & Ding,
# Ann Intern Med 2017) so the calculation is fully transparent:
#   HR -> RR  : RR = (1 - 0.5^sqrt(HR)) / (1 - 0.5^sqrt(1/HR))
#   E-value   : E = RR + sqrt(RR * (RR - 1))
cf <- coef(m1)['r_w']; se <- sqrt(diag(vcov(m1)))['r_w']
hr <- exp(cf); lo <- exp(cf - 1.96 * se); hi <- exp(cf + 1.96 * se)
hr_to_rr <- function(x) (1 - 0.5^sqrt(x)) / (1 - 0.5^sqrt(1 / x))
evalue_rr <- function(rr) {
  if (rr < 1) rr <- 1 / rr
  rr + sqrt(rr * (rr - 1))
}
e_point <- evalue_rr(hr_to_rr(hr))
e_ci <- if (lo > 1 || hi < 1) evalue_rr(hr_to_rr(if (lo > 1) lo else hi)) else 1
cat('\n=== E-value (adjusted rate HR) ===\n')
cat(sprintf('HR for r (per 1 unit) = %.4f (%.4f-%.4f)\n', hr, lo, hi))
cat(sprintf('RR equivalent        = %.4f\n', hr_to_rr(hr)))
cat(sprintf('E-value (point)      = %.4f\n', e_point))
cat(sprintf('E-value (CI bound)   = %.4f\n', e_ci))
write.csv(data.frame(term = 'rate r_w (adjusted, ref+rate model)',
                     HR = signif(hr, 5), lower = signif(lo, 5), upper = signif(hi, 5),
                     RR_equivalent = signif(hr_to_rr(hr), 5),
                     evalue_point = signif(e_point, 5), evalue_ci = signif(e_ci, 5)),
          file.path(OUT, 'evalue_rate.csv'), row.names = FALSE)

# --- corrected IPCW discrimination at tau = 24 d (avoid boundary censoring) --
TAU <- 24
sc <- Score(list('reference' = m0, 'reference+rate' = m1),
            formula = Surv(time_days, death28) ~ 1, data = dd, times = TAU,
            metrics = c('AUC', 'Brier'), conf.int = TRUE)
cat(sprintf('\n=== IPCW Uno C and Brier at %d days ===\n', TAU))
print(sc$AUC$score); print(sc$Brier$score)
# write.csv(as.data.frame(sc$AUC$score), file.path(OUT, 'incremental_unoC.csv'), row.names = FALSE)
# write.csv(as.data.frame(sc$Brier$score), file.path(OUT, 'incremental_brier.csv'), row.names = FALSE)

# --- NRI / IDI at 24 days (riskRegression) ----------------------------------
nri <- try(Score(list('reference' = m0, 'reference+rate' = m1),
                 formula = Surv(time_days, death28) ~ 1, data = dd, times = TAU,
                 metrics = c('AUC'), conf.int = FALSE), silent = TRUE)
cat('\nDone.\n')
