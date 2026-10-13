# =============================================================================
# Landmark (72 h) dual-spline analysis:
#   - attenuation ladder  M1 -> M2 -> M3a -> M3b
#   - continuous rate r with restricted cubic splines on BOTH r and log(baseline)
#   - measurement-count covariate variants
#   - joint likelihood-ratio test on the rate term
#   - Schoenfeld PH diagnostics
#   - RMST difference
# Every number written out is produced by these models on the real dataset.
# Output: output/lock/landmark/*.csv
# =============================================================================
suppressPackageStartupMessages({
  library(survival); library(rms); library(Hmisc)
})
set.seed(20261004)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT <- file.path(BASE, 'output', 'lock', 'landmark')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(BASE, 'data', 'landmark_analysis_dataset.csv'),
              stringsAsFactors = FALSE)
# R's read.csv keeps "True"/"False" as character: coerce explicitly
for (bc in c('alive_at_72', 'in_hosp_at_72', 'landmark')) {
  d[[bc]] <- as.logical(d[[bc]])
}
d <- d[which(d$landmark), ]
d$time_days <- d$time_landmark / 24
d$female <- as.integer(d$gender == 'F')
d$log_baseline <- log(d$first_tn + 0.001)
ddist <- datadist(d); options(datadist = 'ddist')

# --- pre-specified stabilisation of the rate ---------------------------------
q <- quantile(d$r, c(0.01, 0.99), na.rm = TRUE)
d$r_w <- pmin(pmax(d$r, q[1]), q[2])          # winsorised at 1st/99th percentile
cat('rate winsor limits:', round(q, 3), '\n')
cat('landmark n =', nrow(d), ' deaths =', sum(d$death28), '\n\n')

cov_adm <- c('age', 'female', 'sofa_24h', 'invasive_vent', 'norepinephrine',
             'creatinine_max', 'lactate_max')
cov_lm <- c('sofa72', 'vaso72')

fit_and_report <- function(formula, data, label, terms_of_interest) {
  m <- try(cph(formula, data = data, x = TRUE, y = TRUE, surv = TRUE), silent = TRUE)
  if (inherits(m, 'try-error')) {
    cat('  [FAIL]', label, ':', attr(m, 'condition')$message, '\n')
    return(NULL)
  }
  cf <- coef(m); se <- sqrt(diag(vcov(m)))
  an <- anova(m)
  # --- Report the FITTED sample, not the supplied frame:
  # rms::cph applies na.delete by default, so nrow(data) overstates the analysed n.
  # model.frame(m) cannot be used here (the data argument is no longer in scope),
  # so the complete-case rows are computed explicitly from the variables used.
  vars_used <- all.vars(formula)
  cc <- complete.cases(data[, vars_used, drop = FALSE])
  used_n <- sum(cc)
  used_events <- sum(data$death28[cc])
  rows <- list()
  for (t in terms_of_interest) {
    i <- which(names(cf) == t)
    if (!length(i)) i <- grep(paste0('^', t), names(cf))[1]
    if (is.na(i) || !length(i)) next
    rn <- which(trimws(rownames(an)) == t)
    rows[[t]] <- data.frame(model = label, term = t, n = used_n, events = used_events,
                            supplied_n = nrow(data), supplied_events = sum(data$death28),
                            HR = exp(cf[i]),
                            lower = exp(cf[i] - 1.96 * se[i]),
                            upper = exp(cf[i] + 1.96 * se[i]),
                            p_term = if (length(rn)) an[rn, 'P'] else NA_real_)
  }
  do.call(rbind, rows)
}

res <- list()
res[['M1']]  <- fit_and_report(Surv(time_days, death28) ~ r_w, d, 'M1 crude r', 'r_w')
res[['M2']]  <- fit_and_report(as.formula(paste('Surv(time_days, death28) ~ r_w +',
                                                paste(cov_adm, collapse = ' + '))),
                               d, 'M2 +admission', 'r_w')
res[['M2n']] <- fit_and_report(as.formula(paste('Surv(time_days, death28) ~ r_w + n_meas +',
                                                paste(cov_adm, collapse = ' + '))),
                               d, 'M2 +admission +n_meas', 'r_w')
res[['M3a']] <- fit_and_report(as.formula(paste('Surv(time_days, death28) ~ r_w + n_meas +',
                                                'log_baseline +',
                                                paste(c(cov_adm, cov_lm), collapse = ' + '))),
                               d, 'M3a +landmark state (main)', c('r_w', 'log_baseline', 'n_meas'))
res[['M3b']] <- {
  dd <- d[!is.na(d$min_lvef_72h), ]
  fit_and_report(as.formula(paste('Surv(time_days, death28) ~ r_w + n_meas +',
                                  'log_baseline + min_lvef_72h +',
                                  paste(c(cov_adm, cov_lm), collapse = ' + '))),
                 dd, 'M3b +LVEF (echo subgroup)', c('r_w', 'log_baseline', 'n_meas'))
}
out <- do.call(rbind, res)
num <- c('HR', 'lower', 'upper', 'p_term')
out[num] <- lapply(out[num], function(x) signif(as.numeric(x), 6))
write.csv(out, file.path(OUT, 'landmark_attenuation_ladder.csv'), row.names = FALSE)
cat('=== attenuation ladder ===\n'); print(out, row.names = FALSE)

# --- dual spline model (rate + baseline both flexible) -----------------------
# use survival::coxph + splines::ns (restricted/natural cubic spline) to avoid
# rms Design interactions between two rcs() terms
suppressPackageStartupMessages(library(splines))
dd <- d[complete.cases(d[, c('r_w','log_baseline','n_meas', cov_adm, cov_lm)]), ]
cat('\ndual-spline complete-case n =', nrow(dd), ' deaths =', sum(dd$death28), '\n')

f_spline <- as.formula(paste('Surv(time_days, death28) ~ ns(r_w, 4) + ns(log_baseline, 4) +',
                             'r_w:log_baseline + n_meas +',
                             paste(c(cov_adm, cov_lm), collapse = ' + ')))
m_sp <- coxph(f_spline, data = dd, x = TRUE, y = TRUE)
cat('\n=== dual-spline model: LR tests per term ===\n')
an_sp <- drop1(m_sp, . ~ ., test = 'Chisq')
print(an_sp)
write.csv(data.frame(term = rownames(an_sp), as.data.frame(an_sp), row.names = NULL),
          file.path(OUT, 'landmark_dualspline_lrtests.csv'), row.names = FALSE)

f_norate <- as.formula(paste('Surv(time_days, death28) ~ ns(log_baseline, 4) + n_meas +',
                             paste(c(cov_adm, cov_lm), collapse = ' + ')))
m_norate <- coxph(f_norate, data = dd, x = TRUE, y = TRUE)
lr <- 2 * (m_sp$loglik[2] - m_norate$loglik[2])
df_lr <- length(coef(m_sp)) - length(coef(m_norate))
p_lr <- pchisq(lr, df_lr, lower.tail = FALSE)
cat(sprintf('\nJoint LR on rate block: chi2 = %.3f, df = %d, P = %.4f\n', lr, df_lr, p_lr))
write.csv(data.frame(test = 'joint LR on rate block (spline + interaction)',
                     chi2 = round(lr, 4), df = df_lr, p = signif(p_lr, 6)),
          file.path(OUT, 'landmark_rate_jointLR.csv'), row.names = FALSE)

# linear-terms version of the same joint test (plan-specified form)
f_lin <- as.formula(paste('Surv(time_days, death28) ~ r_w + log_baseline +',
                          'r_w:log_baseline + n_meas +',
                          paste(c(cov_adm, cov_lm), collapse = ' + ')))
m_lin <- coxph(f_lin, data = dd, x = TRUE, y = TRUE)
f_lin0 <- as.formula(paste('Surv(time_days, death28) ~ log_baseline + n_meas +',
                           paste(c(cov_adm, cov_lm), collapse = ' + ')))
m_lin0 <- coxph(f_lin0, data = dd, x = TRUE, y = TRUE)
lr2 <- 2 * (m_lin$loglik[2] - m_lin0$loglik[2])
p2 <- pchisq(lr2, 2, lower.tail = FALSE)
cat(sprintf('Joint LR on rate block (linear form, 2 df): chi2 = %.3f, P = %.4f\n', lr2, p2))
write.csv(data.frame(test = 'joint LR on rate block (linear form, 2 df)',
                     chi2 = round(lr2, 4), df = 2, p = signif(p2, 6)),
          file.path(OUT, 'landmark_rate_jointLR_linear.csv'), row.names = FALSE)

# --- PH diagnostics ----------------------------------------------------------
zph <- cox.zph(m_sp)
cat('\n=== Schoenfeld PH (dual-spline model) ===\n'); print(zph$table)
write.csv(data.frame(term = rownames(zph$table), as.data.frame(zph$table)),
          file.path(OUT, 'landmark_ph_schoenfeld.csv'), row.names = FALSE)

# --- sensitivity: require dt >= 6 h ------------------------------------------
d6 <- d[which(d$dt_hours >= 6), ]
d6 <- d6[complete.cases(d6[, c('r_w','log_baseline','n_meas', cov_adm, cov_lm)]), ]
m6 <- coxph(as.formula(paste('Surv(time_days, death28) ~ r_w + n_meas + log_baseline +',
                             paste(c(cov_adm, cov_lm), collapse = ' + '))),
            data = d6, x = TRUE, y = TRUE)
cf6 <- coef(m6); se6 <- sqrt(diag(vcov(m6)))
r6 <- which(names(cf6) == 'r_w')
a6 <- drop1(m6, . ~ ., test = 'Chisq'); p6 <- a6['r_w', 'Pr(>Chi)']
cat(sprintf('\nSensitivity dt>=6h: n=%d deaths=%d  HR(r) per 1 unit = %.3f (%.3f-%.3f), P=%.4f\n',
            nrow(d6), sum(d6$death28), exp(cf6[r6]), exp(cf6[r6] - 1.96 * se6[r6]),
            exp(cf6[r6] + 1.96 * se6[r6]), p6))
write.csv(data.frame(analysis = 'dt>=6h sensitivity', n = nrow(d6),
                     events = sum(d6$death28),
                     HR = signif(exp(cf6[r6]), 6),
                     lower = signif(exp(cf6[r6] - 1.96 * se6[r6]), 6),
                     upper = signif(exp(cf6[r6] + 1.96 * se6[r6]), 6),
                     p = signif(p6, 6)),
          file.path(OUT, 'landmark_dt6_sensitivity.csv'), row.names = FALSE)

# --- RMST difference (top vs bottom third of rate) ---------------------------
if (requireNamespace('survRM2', quietly = TRUE)) {
  library(survRM2)
  dd$r_tert <- cut(dd$r_w, quantile(dd$r_w, c(0, 2/3, 1)), labels = c('lowmid','high'),
                   include.lowest = TRUE)
  rm2 <- rmst2(dd$time_days, dd$death28, as.integer(dd$r_tert == 'high'), tau = 25)
  cat('\n=== RMST (high vs low/mid rate, tau = 25 days) ===\n'); print(rm2$unadjusted.result)
  write.csv(as.data.frame(rm2$unadjusted.result), file.path(OUT, 'landmark_rmst.csv'))
}
cat('\nDONE. outputs in', OUT, '\n')
