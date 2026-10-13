library(survival)

df <- read.csv('v4_cohort_2plus_clean.csv', stringsAsFactors=FALSE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group, levels=c('Stable-low','Falling','Stable-high','Rising'))
df$time_28d <- pmin(ifelse(df$mortality_28==1, df$los_icu*24, 672), 672)

cat('=== TEST 1: PROPORTIONAL HAZARDS ASSUMPTION ===\n\n')
m3 <- coxph(Surv(time_28d, mortality_28) ~ trop_dynamic_group + age + gender +
              sofa_24h + invasive_vent + norepinephrine + creatinine_max +
              myocardial_infarct + congestive_heart_failure + lactate_max +
              log(baseline_troponin + 0.001) + min_lvef_72h, data=df)
ph <- cox.zph(m3)
print(ph)
cat('\nGlobal PH test p-value:', ph$table['GLOBAL','p'], '\n')

cat('\n=== TEST 2: VIF (collinearity diagnostics) ===\n\n')
# Approximate VIF using linear model
vars <- c('age','sofa_24h','creatinine_max','lactate_max','baseline_troponin','min_lvef_72h')
df_cc <- df[complete.cases(df[,c(vars,'trop_dynamic_group')]),]
df_cc$log_base <- log(df_cc$baseline_troponin + 0.001)
df_cc$delta <- df_cc$delta_troponin
fit_lm <- lm(log_base ~ delta, data=df_cc)
r2 <- summary(fit_lm)$r.squared
cat('Collinearity between delta and log(baseline): R2 =', round(r2,3),
    ' -> VIF approx =', round(1/(1-r2),2), '\n')
# Correlation between baseline and group
cat('Mean baseline troponin by group:\n')
print(tapply(df_cc$baseline_troponin, df_cc$trop_dynamic_group, median))

cat('\n=== TEST 3: 72-HOUR LANDMARK ANALYSIS (immortal time bias correction) ===\n\n')
# Patients alive at 72h, follow-up starts at 72h
df_lm <- df[df$time_28d >= 72 | df$mortality_28 == 0, ]
df_lm$time_lm <- df_lm$time_28d - 72
df_lm$time_lm <- pmax(df_lm$time_lm, 0.5)
cat('Landmark cohort (alive at 72h): n =', nrow(df_lm), '\n')
cat('Deaths after 72h:', sum(df_lm$mortality_28), '\n')
cat('Early deaths excluded (<72h):', nrow(df) - nrow(df_lm), '\n\n')

m1_lm <- coxph(Surv(time_lm, mortality_28) ~ trop_dynamic_group, data=df_lm)
cat('Landmark Model 1 (crude):\n')
print(round(summary(m1_lm)$coefficients, 4))

m3_lm <- coxph(Surv(time_lm, mortality_28) ~ trop_dynamic_group + age + gender +
                 sofa_24h + invasive_vent + norepinephrine + creatinine_max +
                 myocardial_infarct + congestive_heart_failure + lactate_max +
                 log(baseline_troponin + 0.001) + min_lvef_72h, data=df_lm)
cat('\nLandmark Model 3 (full):\n')
s <- summary(m3_lm)
for(i in 1:3) {
  cat(sprintf('  %s: HR=%.2f (%.2f-%.2f), P=%.4f\n', rownames(s$coefficients)[i],
              s$coefficients[i,2], s$conf.int[i,3], s$conf.int[i,4], s$coefficients[i,5]))
}

cat('\n=== TEST 4: EQUIVALENCE TEST (TOST) for Rising vs Stable-low ===\n\n')
# Is the 95% CI entirely within equivalence margin HR (1/1.25, 1.25)?
hr <- s$coefficients['trop_dynamic_groupRising', 2]
lo <- s$conf.int['trop_dynamic_groupRising', 3]
hi <- s$conf.int['trop_dynamic_groupRising', 4]
cat(sprintf('Rising vs Stable-low (landmark M3): HR=%.2f (%.2f-%.2f)\n', hr, lo, hi))
margin <- 1.25
cat('Equivalence margin: HR within (', round(1/margin,2), ',', margin, ')\n')
if(lo > 1/margin & hi < margin) {
  cat('-> EQUIVALENCE ESTABLISHED (CI entirely within margin)\n')
} else {
  cat('-> EQUIVALENCE NOT established: CI (', round(lo,2), '-', round(hi,2),
      ') crosses margin -> cannot rule out HR up to', round(hi,2), '\n')
}

# Post-hoc power for main analysis
cat('\n=== TEST 5: POST-HOC POWER (approximate) ===\n')
df_cc2 <- df[complete.cases(df[,c('trop_dynamic_group','age','gender','sofa_24h','invasive_vent',
                                   'norepinephrine','creatinine_max','myocardial_infarct',
                                   'congestive_heart_failure','lactate_max','baseline_troponin','min_lvef_72h')]),]
n_rising <- sum(df_cc2$trop_dynamic_group=='Rising')
n_ref <- sum(df_cc2$trop_dynamic_group=='Stable-low')
d_events <- sum(df_cc2$mortality_28[df_cc2$trop_dynamic_group %in% c('Rising','Stable-low')])
cat('Rising n =', n_rising, ', Stable-low n =', n_ref, ', deaths (both) =', d_events, '\n')
# Rough power calc: with ~140 events in reference and ~160 in rising, SE for log-HR ~ 0.11
# Minimum detectable HR at 80% power, alpha 0.05 (2-sided): exp(2.8 * SE)
se_approx <- sqrt(1/(d_events * n_rising/(n_rising+n_ref) * n_ref/(n_rising+n_ref)))
cat('Approx SE(log HR):', round(se_approx,3), '\n')
cat('Minimum detectable HR (80% power):', round(exp(2.8*se_approx),2), '\n')
cat('-> We can rule out HR > ~', round(exp(2.8*se_approx),2),
    'but observed CI upper bound was 1.73, wider than this threshold\n')
