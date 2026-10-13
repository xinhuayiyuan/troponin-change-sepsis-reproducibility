# Supplementary analyses, part 2: multiple imputation of missing covariates for the
# landmark model (m = 20, chained equations, Rubin pooling).
library(survival)
library(mice)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT <- file.path(BASE, 'output')
df <- read.csv(file.path(BASE, 'data', 'v4_cohort_2plus_true_time.csv'), stringsAsFactors = FALSE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group,
                                levels = c('Stable-low', 'Falling', 'Stable-high', 'Rising'))
df$gender <- factor(df$gender)
S <- function(d) Surv(d$time_28d_true, d$mortality_28)
m2rhs <- ~ trop_dynamic_group + age + gender + sofa_24h + invasive_vent + norepinephrine +
  creatinine_max + myocardial_infarct + congestive_heart_failure + lactate_max
m3rhs <- update(m2rhs, . ~ . + log(baseline_troponin + 0.001) + min_lvef_72h)
cc2 <- complete.cases(df[, c('myocardial_infarct', 'congestive_heart_failure', 'lactate_max', 'sofa_24h')])
d2 <- df[cc2, ]
d3 <- d2[!is.na(d2$min_lvef_72h), ]

fmt <- function(hr, lo, hi, p) sprintf('%.2f (%.2f-%.2f), P=%s', hr, lo, hi,
                                       ifelse(p < 0.001, '<0.001', sprintf('%.3f', p)))

# ---------- 1.3 Multiple imputation sensitivity (Table S8) ----------
impvars <- c('trop_dynamic_group', 'age', 'gender', 'sofa_24h', 'invasive_vent', 'norepinephrine',
             'creatinine_max', 'myocardial_infarct', 'congestive_heart_failure', 'lactate_max',
             'baseline_troponin', 'min_lvef_72h', 'time_28d_true', 'mortality_28')
dm <- df[, impvars]
dm$log_base <- log(dm$baseline_troponin + 0.001)
dm$baseline_troponin <- NULL
meth <- make.method(dm)
meth[c('myocardial_infarct', 'congestive_heart_failure')] <- 'logreg'
set.seed(20261001)
imp <- mice(dm, m = 20, maxit = 10, method = meth, printFlag = FALSE)
fitmi <- with(imp, coxph(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender +
                           sofa_24h + invasive_vent + norepinephrine + creatinine_max +
                           myocardial_infarct + congestive_heart_failure + lactate_max +
                           log_base + min_lvef_72h))
pooled <- summary(pool(fitmi), conf.int = TRUE)
row_r <- pooled[pooled$term == 'trop_dynamic_groupRising', ]
mi_txt <- fmt(exp(row_r$estimate), exp(row_r$`2.5 %`), exp(row_r$`97.5 %`), row_r$p.value)
cat('--- Table S8: MI pooled Rising vs Stable-low ---\n', mi_txt, '\n')

# export Table S8 as CSV (complete-case comparator recomputed from the same data)
fit_cc <- coxph(update(m3rhs, Surv(time_28d_true, mortality_28) ~ .), data = d3)
cc_ci <- summary(fit_cc)$conf.int['trop_dynamic_groupRising', ]
cc_p <- summary(fit_cc)$coefficients['trop_dynamic_groupRising', 5]
s8 <- data.frame(
  approach = c('Complete-case (Model 3)', 'Multiple imputation (m=20)'),
  n = c(nrow(d3), nrow(df)),
  hr = c(round(cc_ci[1], 3), round(exp(row_r$estimate), 3)),
  ci_low = c(round(cc_ci[3], 3), round(exp(row_r$`2.5 %`), 3)),
  ci_high = c(round(cc_ci[4], 3), round(exp(row_r$`97.5 %`), 3)),
  p = c(signif(cc_p, 3), signif(row_r$p.value, 3)))
write.csv(s8, file.path(OUT, 'tableS8_mi.csv'), row.names = FALSE)
cat('--- exported tableS8_mi.csv ---\n')
