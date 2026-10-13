# =============================================================================
# Sensitivity analyses that need no new extraction
#   S2  exact/grace landmark (exclude deaths or discharges within 6 h after 72 h)
#       + Fisher exact test on the primary 2x2 contrast
#   S3  competing risk (Fine-Gray) at the landmark; discharge alive = competing
#   S4  complete case vs multiple imputation (m = 20) for the main model
#   S5  time-dependent (IPCW) AUC at 7/14/21/28 days, reference vs +rate
#   S6  excluding CKD or renal replacement therapy
#   S7  in-hospital death as the outcome instead of 28-day all-cause death
#   S9  echo subgroup (M3b) - reported here for a single comparison table
#   S10 LOD/2 sensitivity (values at/below the assay floor)
#   S11 Schoenfeld PH diagnostics for EVERY rung of the attenuation ladder
#
# Input : data/landmark_analysis_dataset.csv  (frozen derived dataset)
# Output: output/lock/sensitivity/*.csv
# Every number is produced by the models below; nothing is imputed by hand.
# =============================================================================
suppressPackageStartupMessages({
  library(survival); library(splines)
})
set.seed(20261005)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT  <- file.path(BASE, 'output', 'lock', 'sensitivity')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(BASE, 'data', 'landmark_analysis_dataset.csv'),
              stringsAsFactors = FALSE)
for (bc in c('alive_at_72', 'in_hosp_at_72', 'landmark')) d[[bc]] <- as.logical(d[[bc]])
d <- d[which(d$landmark), ]
d$time_days <- d$time_landmark / 24
d$female <- as.integer(d$gender == 'F')
d$log_baseline <- log(d$first_tn + 0.001)
qq <- quantile(d$r, c(0.01, 0.99), na.rm = TRUE)
d$r_w <- pmin(pmax(d$r, qq[1]), qq[2])
cat('landmark cohort n =', nrow(d), ' deaths =', sum(d$death28), '\n\n')

cov_adm <- c('age', 'female', 'sofa_24h', 'invasive_vent', 'norepinephrine',
             'creatinine_max', 'lactate_max')
cov_lm  <- c('sofa72', 'vaso72')
f_main  <- as.formula(paste('Surv(time_days, death28) ~ r_w + n_meas + log_baseline +',
                            paste(c(cov_adm, cov_lm), collapse = ' + ')))

grab <- function(m, term) {
  cf <- coef(m); se <- sqrt(diag(vcov(m)))
  i <- which(names(cf) == term)
  if (!length(i)) return(data.frame(term = term, HR = NA, lo = NA, hi = NA, p = NA))
  data.frame(term = term, HR = exp(cf[i]), lo = exp(cf[i] - 1.96 * se[i]),
             hi = exp(cf[i] + 1.96 * se[i]),
             p = summary(m)$coefficients[i, ncol(summary(m)$coefficients)])
}
rows <- list()
addrow <- function(analysis, note, m, term = 'r_w') {
  g <- grab(m, term)
  rows[[length(rows) + 1]] <<- data.frame(analysis = analysis, note = note,
                                          n = m$n, events = m$nevent, g)
}

# ---- S0: reference model (recomputed here for comparability) ----------------
m_ref <- coxph(f_main, data = d, x = TRUE, y = TRUE)
addrow('S0 primary model (landmark 72h)',
       paste('same specification as the locked main model; this table lists only the sensitivity',
             'analyses expressed as Cox hazard ratios - S1 (landmark timing), S3 (Fine-Gray sHR),',
             'S5 (IPCW AUC), S8 (widened inclusion) and S12 (negative-control OR) are in separate files'),
       m_ref)

# ---- S2: exact / grace landmark + Fisher exact test -------------------------
d$t_disch_d <- (d$hours_to_disch - 72) / 24
grace <- 6 / 24                                     # 6 hours
boundary <- ((d$death28 == 1) & (d$time_days < grace)) | (d$t_disch_d < grace)
d2 <- d[!boundary, ]
addrow('S2 exact landmark (6h grace window)',
       sprintf('excluded %d patients dying/discharged within 6h after the landmark', sum(boundary)),
       coxph(f_main, data = d2, x = TRUE, y = TRUE))
sub <- d[d$trop_dynamic_group %in% c('Rising', 'Stable-low'), ]
ft <- fisher.test(table(factor(sub$trop_dynamic_group, levels = c('Stable-low', 'Rising')),
                        factor(sub$death28, levels = c(0, 1))))
cat(sprintf('S2 Fisher exact (Rising vs Stable-low, crude): OR = %.3f (%.3f-%.3f), P = %.4f  [2x2]\n',
            ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value))
write.csv(data.frame(test = 'Fisher exact, Rising vs Stable-low, crude 2x2',
                     n = nrow(sub), deaths = sum(sub$death28),
                     OR = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2],
                     p = ft$p.value),
          file.path(OUT, 'sens_fisher_exact.csv'), row.names = FALSE)

# ---- S3: competing risk (Fine-Gray) at the landmark -------------------------
# Event time for deaths must be the actual time of death
# (time_landmark), not the discharge time.  Discharge-alive remains the competing
# event timed at hospital discharge (capped at 600 h = 25 days from landmark).
if (requireNamespace('cmprsk', quietly = TRUE)) {
  library(cmprsk)
  d$st <- ifelse(d$death28 == 1, 1,
                 ifelse((d$hours_to_disch - 72) < 600, 2, 0))  # 2 = discharge alive
  d$t_fg_h <- ifelse(d$st == 1, d$time_landmark, pmin(600, d$hours_to_disch - 72))
  d$t_fg <- pmax(d$t_fg_h, 1) / 24                  # days, strictly positive
  cc <- d[complete.cases(d[, c('t_fg', 'st', 'r_w', 'n_meas', 'log_baseline',
                               cov_adm, cov_lm)]), ]
  n_postdisch <- sum(cc$st == 1 & cc$hosp_mort == 0, na.rm = TRUE)
  X <- model.matrix(as.formula(paste('~ r_w + n_meas + log_baseline +',
                                     paste(c(cov_adm, cov_lm), collapse = ' + '))),
                    data = cc)[, -1, drop = FALSE]
  fg <- crr(cc$t_fg, cc$st, X)
  i <- which(colnames(X) == 'r_w')
  se <- sqrt(fg$var[i, i])
  cat(sprintf('S3 Fine-Gray: n = %d, death = %d (of which %d post-discharge), discharge-alive = %d | sHR(r_w) = %.4f (%.4f-%.4f), P = %.4f\n',
              nrow(cc), sum(cc$st == 1), n_postdisch, sum(cc$st == 2),
              exp(fg$coef[i]), exp(fg$coef[i] - 1.96 * se), exp(fg$coef[i] + 1.96 * se),
              2 * pnorm(-abs(fg$coef[i] / se))))
  write.csv(data.frame(analysis = 'S3 Fine-Gray at landmark (discharge alive = competing)',
                       n = nrow(cc), events_death = sum(cc$st == 1),
                       events_discharge = sum(cc$st == 2), term = 'r_w',
                       sHR = exp(fg$coef[i]), lo = exp(fg$coef[i] - 1.96 * se),
                       hi = exp(fg$coef[i] + 1.96 * se),
                       p = 2 * pnorm(-abs(fg$coef[i] / se)),
                       note = sprintf('28-day all-cause death timed at recorded death; %d/%d deaths occurred after hospital discharge',
                                      n_postdisch, sum(cc$st == 1))),
            file.path(OUT, 'sens_finegray_landmark.csv'), row.names = FALSE)
}

n_ckd_missing <- sum(is.na(d$ckd))
d6 <- d[is.na(d$ckd) | d$ckd != 1, ]          # ckd coding missing -> presumed non-CKD (conservative)
d6 <- d6[d6$rrt != 1, ]
n_excluded <- nrow(d) - nrow(d6)
addrow('S6 excluding CKD or RRT',
       sprintf('excluded %d with documented CKD/RRT; %d with missing CKD coding retained as presumed non-CKD (conservative)',
               n_excluded, n_ckd_missing),
       coxph(f_main, data = d6, x = TRUE, y = TRUE))

# ---- S7: in-hospital death as the outcome -----------------------------------
# Use actual death time (time_landmark) for in-hospital deaths and restrict to
# deaths within the 28-day window (mortality_28 == 1); hours_to_disch would
# misclassify patients whose recorded discharge time was later than their death.
d$t_hosp <- pmax(d$time_landmark / 24, 1 / 24)
d$ev_hosp <- as.integer(d$hosp_mort == 1 & d$mortality_28 == 1)
f_hosp <- as.formula(paste('Surv(t_hosp, ev_hosp) ~ r_w + n_meas + log_baseline +',
                           paste(c(cov_adm, cov_lm), collapse = ' + ')))
m_s7 <- coxph(f_hosp, data = d, x = TRUE, y = TRUE)
addrow('S7 in-hospital death (not 28-day all-cause)',
       sprintf('%d in-hospital deaths within 28 days in the full landmark cohort; %d in the fitted complete-case set (see events column)',
               sum(d$ev_hosp, na.rm = TRUE), m_s7$nevent),
       m_s7)

# ---- S9: echo subgroup (M3b) ------------------------------------------------
dl <- d[!is.na(d$min_lvef_72h), ]
f_lv <- as.formula(paste('Surv(time_days, death28) ~ r_w + n_meas + log_baseline + min_lvef_72h +',
                         paste(c(cov_adm, cov_lm), collapse = ' + ')))
addrow('S9 echo subgroup (M3b)', 'nadir LVEF added as a covariate', coxph(f_lv, data = dl, x = TRUE, y = TRUE))

# ---- S4: complete case vs multiple imputation (m = 20, Rubin pooling) -------
if (requireNamespace('mice', quietly = TRUE)) {
  library(mice)
  vars <- c('time_days', 'death28', 'r_w', 'log_baseline', 'n_meas', cov_adm, cov_lm)
  md <- d[, vars]
  cat('\nS4 missingness in the main-model variables:\n'); print(colSums(is.na(md)))
  # Methodological note: r_w is a derived variable (function of first_tn, last_tn,
  # dt_hours). mice imputes it independently, which does not preserve the functional
  # relationship. However, r_w itself is virtually complete (see count above);
  # missingness is driven by sofa72/vaso72/lactate_max. Independent imputation of r_w
  # is therefore expected to have negligible impact.
  imp <- mice(md, m = 20, printFlag = FALSE, seed = 20261005)
  co <- se <- numeric(20)
  for (k in 1:20) {
    mk <- coxph(f_main, data = complete(imp, k))
    co[k] <- coef(mk)[['r_w']]
    se[k] <- sqrt(vcov(mk)[['r_w', 'r_w']])
  }
  qbar <- mean(co); ubar <- mean(se^2); b <- var(co)
  tot <- ubar + (1 + 1 / 20) * b
  dfm <- (19) * (1 + ubar / ((1 + 1 / 20) * b))^2
  lo_i <- qbar - qt(0.975, dfm) * sqrt(tot); hi_i <- qbar + qt(0.975, dfm) * sqrt(tot)
  p_i <- 2 * pt(-abs(qbar / sqrt(tot)), dfm)
  cat(sprintf('S4 multiple imputation (m=20): HR(r_w) = %.4f (%.4f-%.4f), P = %.4f\n',
              exp(qbar), exp(lo_i), exp(hi_i), p_i))
  rows[[length(rows) + 1]] <- data.frame(
    analysis = 'S4 multiple imputation (m = 20, Rubin)',
    note = sprintf('complete-case n = %d; imputed n = %d', m_ref$n, nrow(md)),
    n = nrow(md), events = sum(md$death28), term = 'r_w',
    HR = exp(qbar), lo = exp(lo_i), hi = exp(hi_i), p = p_i)
  write.csv(data.frame(analysis = 'S4 MI (m=20)', m = 20, qbar = qbar,
                       between_var = b, total_var = tot, df = dfm,
                       HR = exp(qbar), lo = exp(lo_i), hi = exp(hi_i), p = p_i),
            file.path(OUT, 'sens_multiple_imputation.csv'), row.names = FALSE)
}

# ---- S5: time-dependent (IPCW) AUC, reference vs +rate ----------------------
if (requireNamespace('riskRegression', quietly = TRUE)) {
  library(riskRegression)
  need <- c('time_days', 'death28', 'r_w', 'log_baseline', 'sofa_24h', 'age', 'female')
  dd <- d[complete.cases(d[, need]), ]
  m0 <- coxph(Surv(time_days, death28) ~ sofa_24h + log_baseline + age + female, data = dd, x = TRUE, y = TRUE)
  m1 <- coxph(Surv(time_days, death28) ~ sofa_24h + log_baseline + age + female + r_w, data = dd, x = TRUE, y = TRUE)
  # Follow-up ends at day 25 from the landmark, so t = 28
  # is beyond the risk set and would be silently dropped by riskRegression.
  # Using t = 24 aligns with the package-wide TAU = 24.
  sc <- Score(list(reference = m0, 'reference+rate' = m1),
              formula = Surv(time_days, death28) ~ 1, data = dd,
              times = c(7, 14, 21, 24), metrics = c('AUC'), conf.int = TRUE)
  au <- as.data.frame(sc$AUC$score)
  print(au)
  stopifnot(nrow(au) == 8)   # 4 timepoints x 2 models
  write.csv(au, file.path(OUT, 'sens_timedependent_auc.csv'), row.names = FALSE)
}

# ---- S10: LOD/sqrt(2) substitution sensitivity ------------------------------
LOD <- 0.01
d$first_lod <- ifelse(d$first_tn <= LOD, LOD / sqrt(2), d$first_tn)
d$last_lod  <- ifelse(d$last_tn  <= LOD, LOD / sqrt(2), d$last_tn)
d$r_lod <- (log(d$last_lod) - log(d$first_lod)) / d$dt_hours * 72
q3 <- quantile(d$r_lod, c(0.01, 0.99), na.rm = TRUE)
d$r_w_lod <- pmin(pmax(d$r_lod, q3[1]), q3[2])
d$log_baseline_lod <- log(d$first_lod)
f_lod <- as.formula(paste('Surv(time_days, death28) ~ r_w_lod + n_meas + log_baseline_lod +',
                          paste(c(cov_adm, cov_lm), collapse = ' + ')))
m_lod <- coxph(f_lod, data = d, x = TRUE, y = TRUE)
g <- grab(m_lod, 'r_w_lod')
cat(sprintf('S10 LOD/sqrt(2) substitution (%d values at/below %g affected in first or last): ',
            sum(d$first_tn <= LOD | d$last_tn <= LOD), LOD))
cat(sprintf('HR = %.4f (%.4f-%.4f), P = %.4f\n', g$HR, g$lo, g$hi, g$p))
rows[[length(rows) + 1]] <- data.frame(
  analysis = 'S10 LOD/sqrt(2) substitution',
  note = sprintf('%d patients with first or last value <= %g ng/mL',
                 sum(d$first_tn <= LOD | d$last_tn <= LOD), LOD),
  n = m_lod$n, events = m_lod$nevent, g)

# ---- S11: Schoenfeld PH diagnostics for every rung of the ladder ------------
ladder <- list(
  M1  = Surv(time_days, death28) ~ r_w,
  M2  = as.formula(paste('Surv(time_days, death28) ~ r_w +', paste(cov_adm, collapse = ' + '))),
  M2n = as.formula(paste('Surv(time_days, death28) ~ r_w + n_meas +', paste(cov_adm, collapse = ' + '))),
  M3a = f_main)
ph <- do.call(rbind, lapply(names(ladder), function(nm) {
  m <- coxph(ladder[[nm]], data = d, x = TRUE, y = TRUE)
  z <- cox.zph(m)
  data.frame(model = nm, term = rownames(z$table), chisq = z$table[, 'chisq'],
             df = z$table[, 'df'], p = z$table[, 'p'])
}))
write.csv(ph, file.path(OUT, 'sens_ph_by_model.csv'), row.names = FALSE)
cat('\nS11 global Schoenfeld P by model:\n')
print(ph[ph$term == 'GLOBAL', ])

# ---- write the consolidated sensitivity table -------------------------------
tab <- do.call(rbind, rows)
tab$HR_per_SD <- tab$HR^sd(d$r_w, na.rm = TRUE)
tab$lo_per_SD <- tab$lo^sd(d$r_w, na.rm = TRUE)
tab$hi_per_SD <- tab$hi^sd(d$r_w, na.rm = TRUE)
tab[, c('HR', 'lo', 'hi', 'p', 'HR_per_SD', 'lo_per_SD', 'hi_per_SD')] <-
  lapply(tab[, c('HR', 'lo', 'hi', 'p', 'HR_per_SD', 'lo_per_SD', 'hi_per_SD')],
         function(x) signif(as.numeric(x), 5))
print(tab, row.names = FALSE)
write.csv(tab, file.path(OUT, 'sens_main_table.csv'), row.names = FALSE)
cat('\nDONE. outputs in', OUT, '\n')

