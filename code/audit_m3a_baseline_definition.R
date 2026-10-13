# =============================================================================
# AUDIT: which "baseline" does the primary landmark model
# M3a actually adjust for?
#   * the model code (analysis_19_landmark_models.R):
#       d$log_baseline <- log(d$first_tn)
#   * the manuscript calls this covariate "log(baseline troponin + 0.001)"
#   * the supplement calls it "log(first troponin)"
# This script quantifies the difference and refits M3a under each definition.
# Read-only: writes only to output/.
# =============================================================================
suppressPackageStartupMessages({ library(survival) })

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT  <- file.path(BASE, 'output')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, 'audit_m3a_baseline.log'))

d <- read.csv(file.path(BASE, 'data', 'landmark_analysis_dataset.csv'),
              stringsAsFactors = FALSE)
for (bc in c('alive_at_72', 'in_hosp_at_72', 'landmark')) d[[bc]] <- as.logical(d[[bc]])
d <- d[which(d$landmark), ]
d$time_days <- d$time_landmark / 24
d$female <- as.integer(d$gender == 'F')

cat('landmark cohort n =', nrow(d), ' deaths =', sum(d$death28), '\n\n')

# --- how often does first_tn differ from the pre-specified window-min baseline?
dif <- abs(d$first_tn - d$baseline_troponin)
cat('first_tn vs baseline_troponin (0-6 h window minimum):\n')
cat('  rows differing      :', sum(dif > 1e-9), '/', nrow(d), '\n')
cat('  rows differing >10% :', sum(dif > 0.1 * pmax(d$baseline_troponin, 1e-9)), '\n')
cat('  median |difference| :', signif(median(dif), 4), '\n')
cat('  max |difference|    :', signif(max(dif), 4), '\n\n')

q <- quantile(d$r, c(0.01, 0.99), na.rm = TRUE)
d$r_w <- pmin(pmax(d$r, q[1]), q[2])

cov_adm <- c('age', 'female', 'sofa_24h', 'invasive_vent', 'norepinephrine',
             'creatinine_max', 'lactate_max')
cov_lm  <- c('sofa72', 'vaso72')

fit_m3a <- function(dd, base_var, label, offset = 0) {
  dd <- dd[complete.cases(dd[, c('r_w', base_var, 'n_meas', cov_adm, cov_lm)]), ]
  dd$log_base_used <- log(dd[[base_var]] + offset)
  f <- as.formula(paste('Surv(time_days, death28) ~ r_w + n_meas + log_base_used +',
                        paste(c(cov_adm, cov_lm), collapse = ' + ')))
  m <- coxph(f, data = dd, x = TRUE, y = TRUE)
  cf <- coef(m); se <- sqrt(diag(vcov(m)))
  i <- which(names(cf) == 'r_w'); j <- which(names(cf) == 'log_base_used')
  cat(sprintf('%-46s n=%d events=%d\n', label, nrow(dd), sum(dd$death28)))
  cat(sprintf('   r_w        HR %.4f (%.4f-%.4f)\n', exp(cf[i]),
              exp(cf[i] - 1.96 * se[i]), exp(cf[i] + 1.96 * se[i])))
  cat(sprintf('   log(baseline) HR %.4f (%.4f-%.4f)\n', exp(cf[j]),
              exp(cf[j] - 1.96 * se[j]), exp(cf[j] + 1.96 * se[j])))
  invisible(list(n = nrow(dd), cf = cf, se = se))
}

cat('=== M3a under alternative baseline definitions ===\n')
a <- fit_m3a(d, 'first_tn',        'A) locked: log(first_tn)')
b <- fit_m3a(d, 'first_tn',        'B) log(first_tn + 0.001)', 0.001)
c <- fit_m3a(d, 'baseline_troponin','C) log(baseline_troponin + 0.001)', 0.001)

cat('\nDelta in r_w log-HR (C - A):',
    signif(c$cf[which(names(c$cf) == 'r_w')] - a$cf[which(names(a$cf) == 'r_w')], 4), '\n')
cat('\nR version:', R.version.string, '\n')
sink()
cat('log written to', file.path(OUT, 'audit_m3a_baseline.log'), '\n')
