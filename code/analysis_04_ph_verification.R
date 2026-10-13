library(survival)

df <- read.csv('v4_cohort_2plus_true_time.csv', stringsAsFactors=FALSE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group, levels=c('Stable-low','Falling','Stable-high','Rising'))

# THE model used in the manuscript (true event time, Model 3)
m3 <- coxph(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender +
              sofa_24h + invasive_vent + norepinephrine + creatinine_max +
              myocardial_infarct + congestive_heart_failure + lactate_max +
              log(baseline_troponin + 0.001) + min_lvef_72h, data=df)

cat('=== DEFINITIVE Schoenfeld test on MANUSCRIPT MODEL (true event time, Model 3) ===\n')
ph <- cox.zph(m3)
res <- data.frame(term=rownames(ph$table), chisq=round(ph$table[,'chisq'],3),
                  df=ph$table[,'df'], p=signif(ph$table[,'p'],4))
print(res)
cat('\nViolated (P<0.05):', paste(res$term[res$p<0.05], collapse=', '), '\n')

# Stratify on the ACTUALLY violated variables
df$sofa_tert <- cut(df$sofa_24h, breaks=c(-Inf,4,8,Inf))
df$lact_tert <- cut(df$lactate_max, breaks=c(-Inf,2,4,Inf))
mstrat <- coxph(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender +
                  creatinine_max + myocardial_infarct + congestive_heart_failure +
                  log(baseline_troponin + 0.001) + min_lvef_72h +
                  strata(sofa_tert) + strata(norepinephrine) + strata(lact_tert), data=df)
ph2 <- cox.zph(mstrat)
cat('\n=== Stratified on sofa_tert + norepinephrine + lact_tert (the TRUE violators) ===\n')
cat('GLOBAL P after stratification:', signif(ph2$table['GLOBAL','p'],4), '\n')
s <- summary(mstrat)
for(i in 1:3) cat(sprintf('  %s: HR=%.3f (%.3f-%.3f), P=%.4f\n', rownames(s$coefficients)[i],
            s$coefficients[i,2], s$conf.int[i,3], s$conf.int[i,4], s$coefficients[i,5]))

# Compare with the previous (wrong-variable) stratification
mstrat_old <- coxph(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender +
                      creatinine_max + myocardial_infarct + congestive_heart_failure +
                      log(baseline_troponin + 0.001) + min_lvef_72h +
                      strata(sofa_tert) + strata(invasive_vent) + strata(lact_tert), data=df)
cat('\nPrevious stratification (sofa+invasive_vent+lact): GLOBAL P =',
    signif(cox.zph(mstrat_old)$table['GLOBAL','p'],4), '\n')
s_old <- summary(mstrat_old)
for(i in 1:3) cat(sprintf('  %s: HR=%.3f, P=%.4f\n', rownames(s_old$coefficients)[i],
            s_old$coefficients[i,2], s_old$coefficients[i,5]))
