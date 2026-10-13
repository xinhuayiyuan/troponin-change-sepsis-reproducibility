library(survival)
library(cmprsk)

df <- read.csv('v4_cohort_2plus_true_time.csv', stringsAsFactors=FALSE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group, levels=c('Stable-low','Falling','Stable-high','Rising'))

pr <- function(s, label, idx=1:3) {
  cat(label, '\n')
  for(i in idx) cat(sprintf('  %s: HR=%.2f (%.2f-%.2f), P=%.4f\n', rownames(s$coefficients)[i], s$coefficients[i,2], s$conf.int[i,3], s$conf.int[i,4], s$coefficients[i,5]))
}

covs <- 'age + gender + sofa_24h + invasive_vent + norepinephrine + creatinine_max + myocardial_infarct + congestive_heart_failure + lactate_max'
f1 <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ trop_dynamic_group'))
f2 <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ trop_dynamic_group +', covs))
f3 <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ trop_dynamic_group +', covs, '+ log(baseline_troponin + 0.001) + min_lvef_72h'))

pr(summary(coxph(f1, data=df)), 'M1 crude:')
pr(summary(coxph(f2, data=df)), 'M2 main:')
m3 <- coxph(f3, data=df); s3 <- summary(m3)
pr(s3, 'M3 full:')
cat('  log(baseline): HR=', round(s3$coefficients['log(baseline_troponin + 0.001)',2],2),
    ' P=', round(s3$coefficients['log(baseline_troponin + 0.001)',5],4), '\n')

df$trop_ref_sh <- relevel(df$trop_dynamic_group, ref='Stable-high')
f3sh <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ trop_ref_sh +', covs, '+ log(baseline_troponin + 0.001) + min_lvef_72h'))
s3sh <- summary(coxph(f3sh, data=df))
cat(sprintf('Rising vs Stable-high: HR=%.2f (%.2f-%.2f), P=%.4f\n', s3sh$coefficients[1,2], s3sh$conf.int[1,3], s3sh$conf.int[1,4], s3sh$coefficients[1,5]))

fc <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ delta_troponin + log(baseline_troponin + 0.001) + delta_troponin:log(baseline_troponin + 0.001) +', covs, '+ min_lvef_72h'))
sc <- summary(coxph(fc, data=df))
cat('\nContinuous model:\n')
for(v in c('delta_troponin','log(baseline_troponin + 0.001)','delta_troponin:log(baseline_troponin + 0.001)'))
  cat(sprintf('  %s: HR=%.2f (%.2f-%.2f), P=%.4f\n', v, sc$coefficients[v,2], sc$conf.int[v,3], sc$conf.int[v,4], sc$coefficients[v,5]))

f90 <- as.formula(paste('Surv(time_90d_true, mortality_90) ~ trop_dynamic_group +', covs, '+ log(baseline_troponin + 0.001) + min_lvef_72h'))
pr(summary(coxph(f90, data=df)), '\n90-day M3:')

# Landmark 72h
df_lm <- df[df$time_28d_true >= 72 | df$mortality_28 == 0, ]
df_lm$time_lm <- pmax(df_lm$time_28d_true - 72, 0.5)
flm <- as.formula(paste('Surv(time_lm, mortality_28) ~ trop_dynamic_group +', covs, '+ log(baseline_troponin + 0.001) + min_lvef_72h'))
slm <- summary(coxph(flm, data=df_lm))
cat('\nLandmark 72h M3 (n=', nrow(df_lm), ', deaths=', sum(df_lm$mortality_28), '):\n', sep='')
pr(slm, '')

# PH-stratified
df$sofa_tert <- cut(df$sofa_24h, breaks=c(-Inf,4,8,Inf))
df$lact_tert <- cut(df$lactate_max, breaks=c(-Inf,2,4,Inf))
fstrat <- as.formula('Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender + creatinine_max + myocardial_infarct + congestive_heart_failure + log(baseline_troponin + 0.001) + min_lvef_72h + strata(sofa_tert) + strata(invasive_vent) + strata(lact_tert)')
mstrat <- coxph(fstrat, data=df)
sst <- summary(mstrat)
cat('\nPH-stratified M3:\n')
pr(sst, '')
cat('  GLOBAL PH P =', round(cox.zph(mstrat)$table['GLOBAL','p'], 3), '\n')

# >=4 measurements
df4 <- df[df$trop_valid_count>=4, ]
s4 <- summary(coxph(f3, data=df4))
cat('\n>=4 meas M3 (n=', nrow(df4), '):\n', sep='')
cat(sprintf('  Rising: HR=%.2f (%.2f-%.2f), P=%.4f\n', s4$coefficients[3,2], s4$conf.int[3,3], s4$conf.int[3,4], s4$coefficients[3,5]))

# Excluding prior MI
dfnomi <- df[!is.na(df$myocardial_infarct) & df$myocardial_infarct==0, ]
fnomi <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender + sofa_24h + invasive_vent + norepinephrine + creatinine_max + congestive_heart_failure + lactate_max + log(baseline_troponin + 0.001) + min_lvef_72h'))
snomi <- summary(coxph(fnomi, data=dfnomi))
cat('\nExcl prior MI M3 (n=', nrow(dfnomi), '):\n', sep='')
cat(sprintf('  Rising: HR=%.2f (%.2f-%.2f), P=%.4f\n', snomi$coefficients[3,2], snomi$conf.int[3,3], snomi$conf.int[3,4], snomi$coefficients[3,5]))

# Fine-Gray
df$event_fg <- ifelse(df$mortality_28==1, 1, ifelse(df$hosp_mort==0 & df$los_icu<28, 2, 0))
fg1 <- crr(df$time_28d_true, df$event_fg, cov1=model.matrix(~trop_dynamic_group, data=df)[,-1])
sfg1 <- summary(fg1)
cat('\nFine-Gray crude:\n')
for(i in 1:3) cat(sprintf('  %s: SHR=%.2f (%.2f-%.2f), P=%.4f\n', rownames(sfg1$coef)[i], sfg1$coef[i,2], sfg1$conf.int[i,3], sfg1$conf.int[i,4], sfg1$coef[i,5]))

vars_needed <- c('age','gender','sofa_24h','invasive_vent','norepinephrine',
                 'creatinine_max','myocardial_infarct','congestive_heart_failure',
                 'lactate_max','baseline_troponin','min_lvef_72h')
dfc <- df[complete.cases(df[,vars_needed]), ]
fg3 <- crr(dfc$time_28d_true, dfc$event_fg,
           cov1=model.matrix(~trop_dynamic_group + age + gender + sofa_24h +
                               invasive_vent + norepinephrine + creatinine_max +
                               myocardial_infarct + congestive_heart_failure + lactate_max +
                               log(baseline_troponin + 0.001) + min_lvef_72h, data=dfc)[,-1])
sfg3 <- summary(fg3)
cat('\nFine-Gray adjusted (n=', nrow(dfc), '):\n', sep='')
for(i in 1:3) cat(sprintf('  %s: SHR=%.2f (%.2f-%.2f), P=%.4f\n', rownames(sfg3$coef)[i], sfg3$coef[i,2], sfg3$conf.int[i,3], sfg3$conf.int[i,4], sfg3$coef[i,5]))

cat('\n=== DONE ===\n')

