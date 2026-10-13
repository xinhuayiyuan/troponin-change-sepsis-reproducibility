# =============================================================================
# Reproducibility cross-check: log(baseline cTnT) association with mortality
# under the landmark framing, and the four-group contrast.
# Output: output/lock/landmark/legacy_recheck_landmark.csv
# =============================================================================
suppressPackageStartupMessages({ library(survival) })
BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT <- file.path(BASE, 'output', 'lock', 'landmark')

d <- read.csv(file.path(BASE, 'data', 'landmark_analysis_dataset.csv'), stringsAsFactors = FALSE)
for (bc in c('alive_at_72', 'in_hosp_at_72', 'landmark')) d[[bc]] <- as.logical(d[[bc]])
d <- d[which(d$landmark), ]
d$time_days <- d$time_landmark / 24
d$female <- as.integer(d$gender == 'F')
d$log_baseline <- log(d$first_tn)
d$grp <- factor(d$trop_dynamic_group,
                levels = c('Stable-low', 'Falling', 'Stable-high', 'Rising'))

covs <- c('age', 'female', 'sofa_24h', 'invasive_vent', 'norepinephrine',
          'creatinine_max', 'myocardial_infarct', 'congestive_heart_failure', 'lactate_max')

get <- function(m, nm) {
  cf <- coef(m); se <- sqrt(diag(vcov(m)))
  i <- which(names(cf) == nm)
  if (!length(i)) return(c(NA, NA, NA))
  c(exp(cf[i]), exp(cf[i] - 1.96 * se[i]), exp(cf[i] + 1.96 * se[i]))
}

# A. legacy Model 3 specification, landmark time origin, full complete cases
cc <- d[complete.cases(d[, c('grp', covs, 'log_baseline', 'min_lvef_72h')]), ]
mA <- coxph(as.formula(paste('Surv(time_days, death28) ~ grp +',
                             paste(c(covs, 'log_baseline', 'min_lvef_72h'), collapse = ' + '))),
            data = cc, x = TRUE, y = TRUE)
cat('A. legacy Model 3 spec on landmark cohort: n=', nrow(cc), ' deaths=', sum(cc$death28), '\n')
print(summary(mA)$coefficients[c('grpRising', 'log_baseline', 'min_lvef_72h'), c('exp(coef)','Pr(>|z|)')])

# B. same but WITHOUT the echo/LVEF gate (all landmark complete cases w/o LVEF)
covs2 <- setdiff(covs, c('myocardial_infarct', 'congestive_heart_failure'))
cc2 <- d[complete.cases(d[, c('grp', covs2, 'log_baseline')]), ]
mB <- coxph(as.formula(paste('Surv(time_days, death28) ~ grp +',
                             paste(c(covs2, 'log_baseline'), collapse = ' + '))),
            data = cc2, x = TRUE, y = TRUE)
cat('\nB. no echo gate: n=', nrow(cc2), ' deaths=', sum(cc2$death28), '\n')
print(summary(mB)$coefficients[c('grpRising', 'log_baseline'), c('exp(coef)','Pr(>|z|)')])

# C. baseline alone in the landmark cohort (no group factor, no change) 
cc3 <- d[complete.cases(d[, c('log_baseline')]), ]
mC <- coxph(Surv(time_days, death28) ~ log_baseline, data = cc3, x = TRUE, y = TRUE)
cat('\nC. baseline alone: n=', nrow(cc3), ' deaths=', sum(cc3$death28), '\n')
print(summary(mC)$coefficients[c('log_baseline'), c('exp(coef)','Pr(>|z|)')])

rows <- rbind(
  data.frame(spec = 'A legacy M3 spec + LVEF (landmark)', n = nrow(cc), events = sum(cc$death28),
             term = 'Rising vs Stable-low', HR = get(mA, 'grpRising')[1], lo = get(mA, 'grpRising')[2], hi = get(mA, 'grpRising')[3]),
  data.frame(spec = 'A legacy M3 spec + LVEF (landmark)', n = nrow(cc), events = sum(cc$death28),
             term = 'log(baseline)', HR = get(mA, 'log_baseline')[1], lo = get(mA, 'log_baseline')[2], hi = get(mA, 'log_baseline')[3]),
  data.frame(spec = 'B no echo gate (landmark)', n = nrow(cc2), events = sum(cc2$death28),
             term = 'Rising vs Stable-low', HR = get(mB, 'grpRising')[1], lo = get(mB, 'grpRising')[2], hi = get(mB, 'grpRising')[3]),
  data.frame(spec = 'B no echo gate (landmark)', n = nrow(cc2), events = sum(cc2$death28),
             term = 'log(baseline)', HR = get(mB, 'log_baseline')[1], lo = get(mB, 'log_baseline')[2], hi = get(mB, 'log_baseline')[3]),
  data.frame(spec = 'C baseline alone (landmark)', n = nrow(cc3), events = sum(cc3$death28),
             term = 'log(baseline)', HR = get(mC, 'log_baseline')[1], lo = get(mC, 'log_baseline')[2], hi = get(mC, 'log_baseline')[3]))
rows[, 5:7] <- lapply(rows[, 5:7], function(x) signif(as.numeric(x), 5))
write.csv(rows, file.path(OUT, 'legacy_recheck_landmark.csv'), row.names = FALSE)
cat('\nDone.\n')
