# =============================================================================
# AUDIT: minimum detectable effect at 80% power for the prespecified
# continuous-rate (per-SD) primary analysis, from the fitted landmark M3a model.
# Read-only; writes only to output/.
# =============================================================================
suppressPackageStartupMessages({ library(survival) })
BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT  <- file.path(BASE, 'output')
d <- read.csv(file.path(BASE, 'data', 'landmark_analysis_dataset.csv'), stringsAsFactors = FALSE)
for (bc in c('alive_at_72','in_hosp_at_72','landmark')) d[[bc]] <- as.logical(d[[bc]])
d <- d[which(d$landmark), ]
d$time_days <- d$time_landmark / 24
d$female <- as.integer(d$gender == 'F')
d$log_baseline <- log(d$first_tn)
q <- quantile(d$r, c(0.01, 0.99), na.rm = TRUE)
d$r_w <- pmin(pmax(d$r, q[1]), q[2])
cov_adm <- c('age','female','sofa_24h','invasive_vent','norepinephrine','creatinine_max','lactate_max')
cov_lm  <- c('sofa72','vaso72')
dd <- d[complete.cases(d[, c('r_w','log_baseline','n_meas', cov_adm, cov_lm)]), ]
f <- as.formula(paste('Surv(time_days, death28) ~ r_w + n_meas + log_baseline +',
                      paste(c(cov_adm, cov_lm), collapse = ' + ')))
m <- coxph(f, data = dd, x = TRUE, y = TRUE)
se_unit <- sqrt(diag(vcov(m)))['r_w']
sd_rw   <- sd(dd$r_w)
se_sd   <- se_unit * sd_rw
z       <- qnorm(0.975) + qnorm(0.80)     # 2.8016 for 80% power, 2-sided 5%
cat('M3a sample n =', nrow(dd), ' deaths =', sum(dd$death28), '\n')
cat('SD(r_w) in M3a sample                 :', round(sd_rw, 4), '\n')
cat('SE(log HR) per 1 unit of r_w          :', round(se_unit, 6), '\n')
cat('SE(log HR) per 1 SD of r_w            :', round(se_sd, 6), '\n')
cat('95% CI per SD (should be 0.92-1.09)   :',
    round(exp(coef(m)['r_w']*sd_rw - 1.96*se_sd), 3), '-',
    round(exp(coef(m)['r_w']*sd_rw + 1.96*se_sd), 3), '\n')
cat('Minimum detectable HR (80% power, z=2.8016) per 1 SD :', round(exp(z*se_sd), 3), '\n')
cat('Minimum detectable HR (80% power) per 1 unit          :', round(exp(z*se_unit), 3), '\n')
cat('Power to detect MCID 1.20 per SD                      :',
    round(1 - pnorm(qnorm(0.975) - log(1.20)/se_sd), 4), '\n')
