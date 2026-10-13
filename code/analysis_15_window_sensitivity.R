# Sensitivity analyses from raw MIMIC-IV: 72-h-window organ support,
# admission year, and last measurement >=48 h. All numbers code-derived.
library(survival)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT <- file.path(BASE, 'output')
df <- read.csv(file.path(BASE, 'data', 'v4_cohort_2plus_true_time.csv'), stringsAsFactors = FALSE)
wv <- read.csv(file.path(OUT, 'mimic_window_vars.csv'), stringsAsFactors = FALSE)
df <- merge(df, wv, by = 'stay_id', all.x = TRUE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group,
                                levels = c('Stable-low', 'Falling', 'Stable-high', 'Rising'))
df$gender <- factor(df$gender)
S <- function(d) Surv(d$time_28d_true, d$mortality_28)
fmt <- function(hr, lo, hi, p) sprintf('%.2f (%.2f-%.2f), P=%s', hr, lo, hi,
                                       ifelse(p < 0.001, '<0.001', sprintf('%.3f', p)))
hr_rising <- function(fit) {
  ci <- summary(fit)$conf.int['trop_dynamic_groupRising', ]
  p <- summary(fit)$coefficients['trop_dynamic_groupRising', 5]
  fmt(ci[1], ci[3], ci[4], p)
}
cc2 <- complete.cases(df[, c('myocardial_infarct', 'congestive_heart_failure', 'lactate_max', 'sofa_24h')])
d2 <- df[cc2, ]
d3 <- d2[!is.na(d2$min_lvef_72h), ]

cat('--- flag comparison (whole-stay vs 72h window), full cohort ---\n')
cat('vent:', sum(df$invasive_vent), '->', sum(df$vent72),
    '| norepi:', sum(df$norepinephrine), '->', sum(df$norepi72),
    '| rrt:', sum(df$rrt), '->', sum(df$rrt72), '\n')

# ---------- B4: 72h-window organ support ----------
covs_common <- ~ trop_dynamic_group + age + gender + sofa_24h + creatinine_max +
  myocardial_infarct + congestive_heart_failure + lactate_max
f2w <- coxph(update(covs_common, S(d2) ~ . + vent72 + norepi72 + rrt72), data = d2)
f3w <- coxph(update(covs_common, S(d3) ~ . + vent72 + norepi72 + rrt72 +
                      log(baseline_troponin + 0.001) + min_lvef_72h), data = d3)
cat('--- B4 window Model 2 (n=1,435):', hr_rising(f2w), '\n')
cat('--- B4 window Model 3 (n=902):  ', hr_rising(f3w), '\n')

# ---------- F4: admission year ----------
f3y <- coxph(update(covs_common, S(d3) ~ . + invasive_vent + norepinephrine + rrt +
                      log(baseline_troponin + 0.001) + min_lvef_72h + anchor_year), data = d3)
cat('--- F4 Model 3 + anchor_year (n=902):', hr_rising(f3y), '\n')
cy <- summary(f3y)$conf.int['anchor_year', ]
py <- summary(f3y)$coefficients['anchor_year', 5]
cat('    anchor_year HR per year:', fmt(cy[1], cy[3], cy[4], py), '\n')
# era-stratified (within anchor_year_group, M3 complete-case)
d3$era <- factor(d3$anchor_year_group)
for (e in levels(d3$era)) {
  de <- d3[d3$era == e, ]
  if (nrow(de) < 80 || sum(de$mortality_28) < 20) { cat('    era', e, ': too small (n=', nrow(de), ')\n'); next }
  fe <- tryCatch(coxph(update(covs_common, S(de) ~ . + invasive_vent + norepinephrine + rrt +
                                log(baseline_troponin + 0.001) + min_lvef_72h), data = de),
                 error = function(x) NULL)
  if (!is.null(fe)) cat('    era', e, sprintf('(n=%d, deaths=%d):', nrow(de), sum(de$mortality_28)),
                        hr_rising(fe), '\n')
}

# ---------- F2: last-measurement timing (LOCF concern) ----------
cat('--- last_trop_hours quantiles by group (full cohort) ---\n')
print(tapply(df$last_trop_hours, df$trop_dynamic_group,
             function(x) round(quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE), 1)))
d3l <- d3[!is.na(d3$last_trop_hours) & d3$last_trop_hours >= 48, ]
cat(sprintf('strict >=48h restriction: n=%d (of 902) -- underpowered, reported descriptively only\n',
            nrow(d3l)))
# Preferred sensitivity: exclude Falling patients whose last measurement was <48 h
# (classification-uncertain), keep all other groups intact
d3x <- d3[!(d3$trop_dynamic_group == 'Falling' & !is.na(d3$last_trop_hours) &
              d3$last_trop_hours < 48), ]
d3x$trop_dynamic_group <- droplevels(d3x$trop_dynamic_group)
f3x <- coxph(update(covs_common, S(d3x) ~ . + invasive_vent + norepinephrine + rrt +
                      log(baseline_troponin + 0.001) + min_lvef_72h), data = d3x)
n_excl <- nrow(d3) - nrow(d3x)
cat(sprintf('--- F2 exclude early-last Falling (excluded n=%d; model n=%d, deaths=%d): %s\n',
            n_excl, nrow(d3x), sum(d3x$mortality_28), hr_rising(f3x)))
# Falling estimate within the reclassified sensitivity (secondary contrast of interest)
ci_f <- summary(f3x)$conf.int['trop_dynamic_groupFalling', ]
p_f <- summary(f3x)$coefficients['trop_dynamic_groupFalling', 5]
cat('    Falling vs Stable-low in this sensitivity:', fmt(ci_f[1], ci_f[3], ci_f[4], p_f), '\n')

# ---------- export summary ----------
out <- data.frame(
  analysis = c('72h-window organ support, Model 2 covariates (n=1,435)',
               '72h-window organ support, Model 3 covariates (n=902)',
               'Model 3 + admission year (n=902)',
               'Model 3 excluding Falling patients with last measurement <48 h'),
  hr_txt = c(hr_rising(f2w), hr_rising(f3w), hr_rising(f3y), hr_rising(f3x)))
write.csv(out, file.path(OUT, 'tableS12_window_sens.csv'), row.names = FALSE)
cat('--- exported tableS12_window_sens.csv ---\n')
