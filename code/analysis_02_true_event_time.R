library(survival)

df <- read.csv('v4_cohort_2plus_true_time.csv', stringsAsFactors=FALSE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group, levels=c('Stable-low','Falling','Stable-high','Rising'))

cat('=== RE-ANALYSIS WITH TRUE EVENT TIMES ===\n\n')

# Model 1 (crude) with true time
m1 <- coxph(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group, data=df)
s1 <- summary(m1)
cat('Model 1 (true event time):\n')
for(i in 1:3) {
  cat(sprintf('  %s: HR=%.2f (%.2f-%.2f), P=%.4f\n', rownames(s1$coefficients)[i],
              s1$coefficients[i,2], s1$conf.int[i,3], s1$conf.int[i,4], s1$coefficients[i,5]))
}

# Model 3 with true time
m3 <- coxph(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender +
              sofa_24h + invasive_vent + norepinephrine + creatinine_max +
              myocardial_infarct + congestive_heart_failure + lactate_max +
              log(baseline_troponin + 0.001) + min_lvef_72h, data=df)
s3 <- summary(m3)
cat('\nModel 3 (true event time):\n')
for(i in 1:3) {
  cat(sprintf('  %s: HR=%.2f (%.2f-%.2f), P=%.4f\n', rownames(s3$coefficients)[i],
              s3$coefficients[i,2], s3$conf.int[i,3], s3$conf.int[i,4], s3$coefficients[i,5]))
}

# PH test on corrected model
cat('\nPH test (true time model):\n')
ph <- cox.zph(m3)
print(ph$table[,3])

# 90-day with true time
m90 <- coxph(Surv(time_90d_true, mortality_90) ~ trop_dynamic_group + age + gender +
               sofa_24h + invasive_vent + norepinephrine + creatinine_max +
               myocardial_infarct + congestive_heart_failure + lactate_max +
               log(baseline_troponin + 0.001) + min_lvef_72h, data=df)
s90 <- summary(m90)
cat('\n90-day (true event time):\n')
for(i in 1:3) {
  cat(sprintf('  %s: HR=%.2f (%.2f-%.2f), P=%.4f\n', rownames(s90$coefficients)[i],
              s90$coefficients[i,2], s90$conf.int[i,3], s90$conf.int[i,4], s90$coefficients[i,5]))
}

# Landmark 72h with true time
df_lm <- df[df$time_28d_true >= 72 | df$mortality_28 == 0, ]
df_lm$time_lm <- pmax(df_lm$time_28d_true - 72, 0.5)
m3_lm <- coxph(Surv(time_lm, mortality_28) ~ trop_dynamic_group + age + gender +
                 sofa_24h + invasive_vent + norepinephrine + creatinine_max +
                 myocardial_infarct + congestive_heart_failure + lactate_max +
                 log(baseline_troponin + 0.001) + min_lvef_72h, data=df_lm)
s_lm <- summary(m3_lm)
cat('\nLandmark 72h Model 3 (true time, n=', nrow(df_lm), '):\n')
for(i in 1:3) {
  cat(sprintf('  %s: HR=%.2f (%.2f-%.2f), P=%.4f\n', rownames(s_lm$coefficients)[i],
              s_lm$coefficients[i,2], s_lm$conf.int[i,3], s_lm$conf.int[i,4], s_lm$coefficients[i,5]))
}
