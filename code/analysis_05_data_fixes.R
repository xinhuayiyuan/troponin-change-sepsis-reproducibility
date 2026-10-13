library(survival)

df <- read.csv('v4_cohort_2plus_true_time.csv', stringsAsFactors=FALSE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group, levels=c('Stable-low','Falling','Stable-high','Rising'))
covs <- 'age + gender + sofa_24h + invasive_vent + norepinephrine + creatinine_max + myocardial_infarct + congestive_heart_failure + lactate_max'
f3 <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ trop_dynamic_group +', covs, '+ log(baseline_troponin + 0.001) + min_lvef_72h'))

# ===== D7: 90-day landmark =====
cat('=== D7: 90-DAY LANDMARK (alive at 72h, outcome=90d) ===\n')
df_lm <- df[df$mortality_28==0 | df$time_28d_true >= 72, ]
df_lm$time_lm90 <- pmax(df_lm$time_90d_true - 72, 0.5)
flm90 <- as.formula(paste('Surv(time_lm90, mortality_90) ~ trop_dynamic_group +', covs, '+ log(baseline_troponin + 0.001) + min_lvef_72h'))
m90lm <- coxph(flm90, data=df_lm)
s <- summary(m90lm)
cat('n =', s$n, ', events =', s$nevent, '\n')
for(i in 1:3) cat(sprintf('  %s: HR=%.3f (%.3f-%.3f), P=%.4f\n', rownames(s$coefficients)[i],
            s$coefficients[i,2], s$conf.int[i,3], s$conf.int[i,4], s$coefficients[i,5]))

# ===== D8: delta nonlinearity (log1p + dichotomised) =====
cat('\n=== D8: DELTA SPECIFICATION CHECKS ===\n')
df$delta_log <- log(df$delta_troponin + 0.01)
df$delta_pos <- as.integer(df$delta_troponin > 0)
fc1 <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ delta_log + log(baseline_troponin + 0.001) +', covs, '+ min_lvef_72h'))
fc2 <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ delta_pos + log(baseline_troponin + 0.001) +', covs, '+ min_lvef_72h'))
mc1 <- summary(coxph(fc1, data=df))
mc2 <- summary(coxph(fc2, data=df))
cat('log(delta+0.01): HR=', round(mc1$coefficients['delta_log',2],3), ', P=', round(mc1$coefficients['delta_log',5],3), '\n')
cat('delta>0 (any rise): HR=', round(mc2$coefficients['delta_pos',2],3), ', P=', round(mc2$coefficients['delta_pos',5],3), '\n')

# ===== D11: included vs excluded (indeterminate) comparison =====
cat('\n=== D11: need indeterminate group from raw trop table (SQL) ===\n')

# ===== D2: Model 3 full covariate table =====
cat('\n=== D2: MODEL 3 FULL COVARIATE RESULTS ===\n')
m3 <- coxph(f3, data=df)
s3 <- summary(m3)
out <- data.frame(term=rownames(s3$coefficients), HR=round(s3$coefficients[,2],3),
                  CI_low=round(s3$conf.int[,3],3), CI_high=round(s3$conf.int[,4],3),
                  P=signif(s3$coefficients[,5],3))
print(out)
