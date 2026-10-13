# =============================================================================
# Sensitivity analyses that need the widened extraction
#   S1  alternative landmark times: 24 h and 48 h (rate measured over the same
#       window, always standardised to a 72-h fold-change so the three landmark
#       times are comparable; specification held constant = admission covariates
#       + n_meas + log(baseline), i.e. no 72-h state term)
#   S8  widened inclusion (>=1 cTnT in 0-72 h); the rate is missing for
#       single-measurement patients and is imputed in two documented ways
#   S12 negative control outcome: 30-day readmission with a NON-cardiac
#       principal diagnosis (logistic model; binary flag)
#
# Input : data/alt_cohort_vars.csv  (analysis_27; anchored on the original
#         cohort pipeline's own intermediate tables)
# Output: output/lock/sensitivity/*.csv
# =============================================================================
suppressPackageStartupMessages({ library(survival) })
set.seed(20261005)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT  <- file.path(BASE, 'output', 'lock', 'sensitivity')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

a <- read.csv(file.path(BASE, 'data', 'alt_cohort_vars.csv'), stringsAsFactors = FALSE)
a$female <- as.integer(a$gender == 'F')
cat('rows in the widened (>=1 cTnT) source:', nrow(a), '\n\n')

cov_adm <- c('age', 'female', 'sofa_24h', 'invasive_vent', 'norepinephrine',
             'creatinine_max', 'lactate_max')

grab <- function(m, term) {
  cf <- coef(m); se <- sqrt(diag(vcov(m)))
  i <- which(names(cf) == term)
  if (!length(i)) return(data.frame(term = term, HR = NA, lo = NA, hi = NA, p = NA))
  data.frame(term = term, HR = exp(cf[i]), lo = exp(cf[i] - 1.96 * se[i]),
             hi = exp(cf[i] + 1.96 * se[i]),
             p = summary(m)$coefficients[i, ncol(summary(m)$coefficients)])
}

# ---------------------------------------------------------------- S1 ---------
res1 <- list()
for (W in c(24, 48, 72)) {
  a[[paste0('alive', W)]]  <- (a$time_28d_true > W) | (a$mortality_28 == 0)
  a[[paste0('inhosp', W)]] <- a$hours_to_disch > W
  a[[paste0('lm', W)]]     <- a[[paste0('alive', W)]] & a[[paste0('inhosp', W)]]
  a[[paste0('tw', W)]]     <- pmax(a$time_28d_true - W, 0) / 24
  rv <- paste0('r', W)
  keep <- a[[paste0('lm', W)]] & a[[paste0('n', W)]] >= 2 &
    !is.na(a[[rv]]) & a[[paste0('first_v', W)]] > 0 & a[[paste0('dt_h', W)]] > 0
  dW <- a[keep, ]
  dW$log_baseline <- log(dW[[paste0('first_v', W)]])
  q <- quantile(dW[[rv]], c(0.01, 0.99), na.rm = TRUE)
  dW$r_w <- pmin(pmax(dW[[rv]], q[1]), q[2])
  sd_r <- sd(dW$r_w, na.rm = TRUE)
  names(dW)[names(dW) == sprintf('n%d', W)] <- 'n_meas'
  f <- as.formula(paste(sprintf('Surv(tw%d, death28) ~ r_w + n_meas + log_baseline +', W),
                        paste(cov_adm, collapse = ' + ')))
  m <- coxph(f, data = dW, x = TRUE, y = TRUE)
  g <- grab(m, 'r_w')
  res1[[as.character(W)]] <- data.frame(
    landmark_h = W, n = m$n, events = m$nevent, sd_r = sd_r,
    HR = g$HR, lo = g$lo, hi = g$hi, p = g$p,
    HR_per_SD = g$HR^sd_r, lo_per_SD = g$lo^sd_r, hi_per_SD = g$hi^sd_r,
    cohort = 'widened (>=1 cTnT in 0-72 h), landmark-eligible, >=2 measures within the window',
    specification = 'admission covariates + n_meas + log_baseline (no 72-h state term)',
    fu_days_from_landmark = round((672 - W) / 24, 2),
    note = 'NOT the main analysis: cohort and specification differ from the locked main model (n=2107)')
  cat(sprintf('S1 landmark %2d h: n = %4d, events = %3d | HR/unit = %.4f (%.4f-%.4f) | HR/SD = %.3f (%.3f-%.3f)\n',
              W, m$n, m$nevent, g$HR, g$lo, g$hi, g$HR^sd_r, g$lo^sd_r, g$hi^sd_r))
}
s1 <- do.call(rbind, res1)
write.csv(s1, file.path(OUT, 'sens_landmark_time.csv'), row.names = FALSE)

# ---------------------------------------------------------------- S8 ---------
d <- a[a$lm72 & a$n72 >= 1 & a$first_v72 > 0, ]
d$log_baseline <- log(d$first_v72)
q <- quantile(d$r72, c(0.01, 0.99), na.rm = TRUE)
d$r_w <- pmin(pmax(d$r72, q[1]), q[2])
d$time_days <- d$time_landmark / 24
d$r_missing <- as.integer(is.na(d$r72))
d$n_meas <- d$n72
cat(sprintf('\nS8 widened cohort: n = %d (>=2 measures %d, exactly 1 measure %d), deaths = %d\n',
            nrow(d), sum(d$n72 >= 2), sum(d$n72 == 1), sum(d$death28)))

f8  <- as.formula(paste('Surv(time_days, death28) ~ r_w + n_meas + log_baseline +',
                        paste(cov_adm, collapse = ' + ')))
f8i <- as.formula(paste('Surv(time_days, death28) ~ r_w + r_missing + n_meas + log_baseline +',
                        paste(cov_adm, collapse = ' + ')))
s8rows <- list()
d2 <- d[d$n72 >= 2 & !is.na(d$r_w), ]
m <- coxph(f8, data = d2, x = TRUE, y = TRUE); g <- grab(m, 'r_w')
s8rows[[1]] <- data.frame(approach = 'S8a >=2 measurements only (reference)',
                          n = m$n, events = m$nevent, g)
d0 <- d; d0$r_w[is.na(d0$r_w)] <- 0
m <- coxph(f8, data = d0, x = TRUE, y = TRUE); g <- grab(m, 'r_w')
s8rows[[2]] <- data.frame(approach = 'S8b widened, rate imputed as 0 (no observed change)',
                          n = m$n, events = m$nevent, g)
m <- coxph(f8i, data = d0, x = TRUE, y = TRUE); g <- grab(m, 'r_w')
gm <- grab(m, 'r_missing')
s8rows[[3]] <- data.frame(
  approach = sprintf('S8c widened + missing-rate indicator (indicator HR %.3f, P %.3f)', gm$HR, gm$p),
  n = m$n, events = m$nevent, g)
s8 <- do.call(rbind, s8rows)
sd8 <- sd(d2$r_w, na.rm = TRUE)
s8$HR_per_SD <- s8$HR^sd8; s8$lo_per_SD <- s8$lo^sd8; s8$hi_per_SD <- s8$hi^sd8
print(s8, row.names = FALSE)
write.csv(s8, file.path(OUT, 'sens_widened_inclusion.csv'), row.names = FALSE)

# ---------------------------------------------------------------- S12 --------
# Negative control: 30-day readmission with a NON-cardiac principal diagnosis.
# Report the fitted sample, not the supplied frame.
d0$negctrl <- d0$negctrl_readmit_noncardiac
gm <- glm(as.formula(paste('negctrl ~ r_w + n_meas + log_baseline +',
                           paste(cov_adm, collapse = ' + '))),
          data = d0, family = binomial)
co <- summary(gm)$coefficients
i <- which(rownames(co) == 'r_w')
orr <- exp(co[i, 1]); lo <- exp(co[i, 1] - 1.96 * co[i, 2]); hi <- exp(co[i, 1] + 1.96 * co[i, 2])
mf12 <- model.frame(gm)
n_fit12 <- nrow(mf12); ev_fit12 <- sum(mf12$negctrl)
cat(sprintf('\nS12 negative control (30-d non-cardiac readmission): fitted n = %d (supplied %d), events = %d (supplied %d) | OR/unit = %.4f (%.4f-%.4f), P = %.4f\n',
            n_fit12, nrow(d0), ev_fit12, sum(d0$negctrl, na.rm = TRUE), orr, lo, hi, co[i, 4]))
write.csv(data.frame(analysis = 'S12 negative control: 30-day non-cardiac readmission',
                     n = n_fit12, events = ev_fit12,
                     supplied_n = nrow(d0),
                     supplied_events = sum(d0$negctrl, na.rm = TRUE),
                     OR_per_unit = orr, lo = lo, hi = hi, p = co[i, 4]),
          file.path(OUT, 'sens_negative_control.csv'), row.names = FALSE)

cat('\nDONE. outputs in', OUT, '\n')
