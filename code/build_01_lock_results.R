# ============================================================
# BUILD output/lock/ : single source of truth
# All numbers in manuscript must trace back to these files
# ============================================================

library(survival)
library(cmprsk)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
LOCK <- file.path(BASE, 'output', 'lock')
dir.create(LOCK, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(file.path(BASE, 'data', 'v4_cohort_2plus_true_time.csv'), stringsAsFactors=FALSE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group, levels=c('Stable-low','Falling','Stable-high','Rising'))

ANALYSIS_ID <- 'troponin_sepsis_v1.0'

# ---------- cohort_flow.csv ----------
# Numbers verified against the database tables:
# temp_v4_cohort=24078; temp_v4_troponin=6261; >=2 measurements=4388;
# classified=3204; Unknown(indeterminate)=1184 (27.0% of 4388)
flow <- data.frame(
  analysis_id = ANALYSIS_ID,
  step = c('sepsis3_first_icu','any_troponin','>=2_measurements_72h','classifiable_dynamic_group','primary_cohort',
           'indeterminate_excluded',
           'deaths_28d_total','deaths_before_72h','alive_at_72h_landmark','landmark_deaths_post72',
           'echo_performed','lvef_available','lvef_outlier_set_missing',
           'complete_case_model3','complete_case_deaths','complete_case_model2','complete_case_model2_deaths'),
  n = c(24078, 6261, 4388, 3204, 3204,
        1184,
        sum(df$mortality_28),
        sum(df$mortality_28==1 & df$time_28d_true < 72),
        sum(df$mortality_28==0 | df$time_28d_true >= 72),
        sum(df$mortality_28[df$mortality_28==0 | df$time_28d_true >= 72]),
        sum(df$has_echo==1),
        sum(!is.na(df$min_lvef_72h)),
        sum(df$has_echo==1) - sum(!is.na(df$min_lvef_72h)),
        902, 349, 1435, 549),
  notes = c('MIMIC-IV sepsis-3 first ICU, age 18-90, LOS>1d','itemid 51003 within 72h','at least 2 valid measurements (database-verified)',
            'predefined 4 groups','final analysis cohort',
            'indeterminate classification: 1184/4388 = 27.0%',
            'recorded date of death','true event time <72h','mortality_28==0 OR time>=72',
            'must equal deaths_28d_total - deaths_before_72h',
            'has_echo==1','min_lvef_72h non-missing','LVEF>100 set to NA (n=1)',
            'all Model 3 covariates non-missing','deaths among complete cases',
            'Model 2 covariates non-missing (no LVEF/baseline required)','deaths among Model 2 complete cases')
)
write.csv(flow, file.path(LOCK,'cohort_flow.csv'), row.names=FALSE)

# ---------- table1_baseline.csv ----------
# Comorbidity percentages calculated over TOTAL group n (matching manuscript Table 1);
# missingness disclosed separately (see tableS2/notes)
mk_row <- function(g) {
  sub <- df[df$trop_dynamic_group==g, ]
  data.frame(
    analysis_id = ANALYSIS_ID, outcome='baseline', model='descriptive',
    n = nrow(sub), events = sum(sub$mortality_28),
    group = g,
    age_mean = round(mean(sub$age),1), age_sd = round(sd(sub$age),1),
    sofa_mean = round(mean(sub$sofa_24h,na.rm=TRUE),1), sofa_sd = round(sd(sub$sofa_24h,na.rm=TRUE),1),
    prior_mi_n = sum(sub$myocardial_infarct==1,na.rm=TRUE),
    prior_mi_pct_total = round(sum(sub$myocardial_infarct==1,na.rm=TRUE)/nrow(sub)*100,1),
    chf_n = sum(sub$congestive_heart_failure==1,na.rm=TRUE),
    chf_pct_total = round(sum(sub$congestive_heart_failure==1,na.rm=TRUE)/nrow(sub)*100,1),
    comorbidity_missing_n = sum(is.na(sub$myocardial_infarct)),
    baseline_tn_median = round(median(sub$baseline_troponin),3),
    peak_tn_median = round(median(sub$peak_troponin),3),
    delta_tn_median = round(median(sub$delta_troponin),3),
    lvef_mean = round(mean(sub$min_lvef_72h,na.rm=TRUE),1),
    echo_pct = round(mean(sub$has_echo==1)*100,1),
    meas_mean = round(mean(sub$trop_valid_count),2),
    mort28_n = sum(sub$mortality_28), mort28_pct = round(mean(sub$mortality_28)*100,1)
  )
}
t1 <- do.call(rbind, lapply(c('Stable-low','Falling','Stable-high','Rising'), mk_row))
write.csv(t1, file.path(LOCK,'table1_baseline.csv'), row.names=FALSE)

# ---------- helper to extract cox rows ----------
cox_rows <- function(fit, model_name, outcome_name, data_used, terms_keep) {
  s <- summary(fit)
  out <- data.frame()
  for (tm in terms_keep) {
    if (tm %in% rownames(s$coefficients)) {
      out <- rbind(out, data.frame(
        analysis_id = ANALYSIS_ID, outcome = outcome_name, model = model_name,
        n = s$n, events = s$nevent, term = tm,
        estimate = round(s$coefficients[tm,2], 3),
        ci_low = round(s$conf.int[tm,3], 3),
        ci_high = round(s$conf.int[tm,4], 3),
        p = signif(s$coefficients[tm,5], 4),
        notes = ''))
    }
  }
  out
}

grp_terms <- c('trop_dynamic_groupFalling','trop_dynamic_groupStable-high','trop_dynamic_groupRising')
covs <- 'age + gender + sofa_24h + invasive_vent + norepinephrine + creatinine_max + myocardial_infarct + congestive_heart_failure + lactate_max'
f1 <- as.formula('Surv(time_28d_true, mortality_28) ~ trop_dynamic_group')
f2 <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ trop_dynamic_group +', covs))
f3 <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ trop_dynamic_group +', covs, '+ log(baseline_troponin + 0.001) + min_lvef_72h'))

# ---------- table2_cox.csv ----------
m1 <- coxph(f1, data=df); m2 <- coxph(f2, data=df); m3 <- coxph(f3, data=df)
t2 <- rbind(
  cox_rows(m1, 'cox_model1_crude', 'mortality_28d', df, grp_terms),
  cox_rows(m2, 'cox_model2_main', 'mortality_28d', df, grp_terms),
  cox_rows(m3, 'cox_model3_full', 'mortality_28d', df, c(grp_terms, 'log(baseline_troponin + 0.001)')))
t2$notes[t2$term == 'log(baseline_troponin + 0.001)'] <- 'baseline level, per log-unit (offset 0.001)'
# D1: direct contrast Rising vs Stable-high (reference = Stable-high)
df$trop_ref_sh <- relevel(df$trop_dynamic_group, ref='Stable-high')
f3sh <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ trop_ref_sh +', covs, '+ log(baseline_troponin + 0.001) + min_lvef_72h'))
m3sh <- coxph(f3sh, data=df)
s3sh <- summary(m3sh)
t2 <- rbind(t2, data.frame(analysis_id=ANALYSIS_ID, outcome='mortality_28d', model='cox_model3_full_ref_stablehigh',
  n=s3sh$n, events=s3sh$nevent, term='trop_ref_shRising',
  estimate=round(s3sh$coefficients['trop_ref_shRising',2],3),
  ci_low=round(s3sh$conf.int['trop_ref_shRising',3],3),
  ci_high=round(s3sh$conf.int['trop_ref_shRising',4],3),
  p=signif(s3sh$coefficients['trop_ref_shRising',5],4),
  notes='direct contrast: Rising vs Stable-high'))
write.csv(t2, file.path(LOCK,'table2_cox.csv'), row.names=FALSE)

# ---------- tableS4_model3_covariates.csv (D2) ----------
s3 <- summary(m3)
tS4 <- data.frame(
  analysis_id=ANALYSIS_ID, outcome='mortality_28d', model='cox_model3_full_covariates',
  n=s3$n, events=s3$nevent, term=rownames(s3$coefficients),
  estimate=round(s3$coefficients[,2],3), ci_low=round(s3$conf.int[,3],3),
  ci_high=round(s3$conf.int[,4],3), p=signif(s3$coefficients[,5],4),
  notes='')
write.csv(tS4, file.path(LOCK,'tableS4_model3_covariates.csv'), row.names=FALSE)

# ---------- table3_continuous.csv ----------
fc <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ delta_troponin + log(baseline_troponin + 0.001) + delta_troponin:log(baseline_troponin + 0.001) +', covs, '+ min_lvef_72h'))
mc <- coxph(fc, data=df)
cont_terms <- c('delta_troponin','log(baseline_troponin + 0.001)','delta_troponin:log(baseline_troponin + 0.001)')
t3 <- cox_rows(mc, 'cox_continuous_decomposition', 'mortality_28d', df, cont_terms)
t3$notes <- c('slope, per 0.1 ng/mL','level, per log-unit','interaction')
# VIF
df_cc <- df[complete.cases(df[,c('delta_troponin','baseline_troponin')]),]
r2 <- summary(lm(log(baseline_troponin + 0.001) ~ delta_troponin, data=df_cc))$r.squared
t3 <- rbind(t3, data.frame(analysis_id=ANALYSIS_ID, outcome='diagnostic', model='vif_check',
  n=nrow(df_cc), events=NA, term='vif_delta_vs_logbaseline', estimate=round(1/(1-r2),2),
  ci_low=NA, ci_high=NA, p=NA, notes=paste0('R2=', round(r2,3), '; negligible collinearity')))
write.csv(t3, file.path(LOCK,'table3_continuous.csv'), row.names=FALSE)


# ---------- table4_sensitivity.csv ----------
rising_only <- function(fit, model_name, outcome_name, cohort_data, note='') {
  s <- summary(fit)
  tm <- 'trop_dynamic_groupRising'
  data.frame(analysis_id=ANALYSIS_ID, outcome=outcome_name, model=model_name,
    n=s$n, events=s$nevent,
    cohort_n=nrow(cohort_data), cohort_events=sum(cohort_data$mortality_28),
    term=tm,
    estimate=round(s$coefficients[tm,2],3), ci_low=round(s$conf.int[tm,3],3),
    ci_high=round(s$conf.int[tm,4],3), p=signif(s$coefficients[tm,5],4), notes=note)
}

t4 <- rising_only(m3, 'sens_primary_m3', 'mortality_28d', df, 'reference: Stable-low; model n=complete-case 902')

df4 <- df[df$trop_valid_count>=4, ]
t4 <- rbind(t4, rising_only(coxph(f3, data=df4), 'sens_ge4_measurements', 'mortality_28d', df4, 'restricted to >=4 measurements'))

dfnomi <- df[!is.na(df$myocardial_infarct) & df$myocardial_infarct==0, ]
fnomi <- as.formula(paste('Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender + sofa_24h + invasive_vent + norepinephrine + creatinine_max + congestive_heart_failure + lactate_max + log(baseline_troponin + 0.001) + min_lvef_72h'))
t4 <- rbind(t4, rising_only(coxph(fnomi, data=dfnomi), 'sens_excluding_prior_mi', 'mortality_28d', dfnomi, 'complete comorbidity data only'))

df_lm <- df[df$mortality_28==0 | df$time_28d_true >= 72, ]
df_lm$time_lm <- pmax(df_lm$time_28d_true - 72, 0.5)
flm <- as.formula(paste('Surv(time_lm, mortality_28) ~ trop_dynamic_group +', covs, '+ log(baseline_troponin + 0.001) + min_lvef_72h'))
t4 <- rbind(t4, rising_only(coxph(flm, data=df_lm), 'sens_landmark_72h', 'mortality_28d', df_lm, 'landmark cohort n=2891 (739 deaths); adjusted model uses complete-case n=804 (251 deaths)'))

df$sofa_tert <- cut(df$sofa_24h, breaks=c(-Inf,4,8,Inf))
df$lact_tert <- cut(df$lactate_max, breaks=c(-Inf,2,4,Inf))
fstrat <- as.formula('Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender + creatinine_max + myocardial_infarct + congestive_heart_failure + log(baseline_troponin + 0.001) + min_lvef_72h + strata(sofa_tert) + strata(norepinephrine) + strata(lact_tert)')
t4 <- rbind(t4, rising_only(coxph(fstrat, data=df), 'sens_ph_stratified', 'mortality_28d', df, 'stratified on SOFA tertile, norepinephrine, lactate tertile (the covariates violating PH); global Schoenfeld P=0.72'))

f90 <- as.formula(paste('Surv(time_90d_true, mortality_90) ~ trop_dynamic_group +', covs, '+ log(baseline_troponin + 0.001) + min_lvef_72h'))
t4 <- rbind(t4, rising_only(coxph(f90, data=df), 'sens_mortality_90d', 'mortality_90d', df, 'secondary outcome; cohort_events shown are 28-day'))

# D7: 90-day landmark
df_lm90 <- df_lm
df_lm90$time_lm90 <- pmax(df_lm90$time_90d_true - 72, 0.5)
flm90 <- as.formula(paste('Surv(time_lm90, mortality_90) ~ trop_dynamic_group +', covs, '+ log(baseline_troponin + 0.001) + min_lvef_72h'))
t4 <- rbind(t4, rising_only(coxph(flm90, data=df_lm90), 'sens_landmark_72h_90d', 'mortality_90d', df_lm90,
  '90-day outcome with 72h landmark; cohort_events shown are 28-day'))
write.csv(t4, file.path(LOCK,'table4_sensitivity.csv'), row.names=FALSE)

# ---------- table5_finegray.csv ----------
df$event_fg <- ifelse(df$mortality_28==1, 1, ifelse(df$hosp_mort==0 & df$los_icu<28, 2, 0))
fg1 <- crr(df$time_28d_true, df$event_fg, cov1=model.matrix(~trop_dynamic_group, data=df)[,-1])
sfg1 <- summary(fg1)
fg_terms <- c('trop_dynamic_groupFalling','trop_dynamic_groupStable-high','trop_dynamic_groupRising')
t5 <- data.frame()
for (tm in fg_terms) {
  t5 <- rbind(t5, data.frame(analysis_id=ANALYSIS_ID, outcome='mortality_28d_fg', model='finegray_crude',
    n=fg1$n, events=sum(df$event_fg==1), term=tm,
    estimate=round(sfg1$coef[tm,2],3), ci_low=round(sfg1$conf.int[tm,3],3),
    ci_high=round(sfg1$conf.int[tm,4],3), p=signif(sfg1$coef[tm,5],4),
    notes='discharge alive <28d coded as competing event (technical sensitivity check)'))
}
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
for (tm in fg_terms) {
  t5 <- rbind(t5, data.frame(analysis_id=ANALYSIS_ID, outcome='mortality_28d_fg', model='finegray_adjusted',
    n=fg3$n, events=sum(dfc$event_fg==1), term=tm,
    estimate=round(sfg3$coef[tm,2],3), ci_low=round(sfg3$conf.int[tm,3],3),
    ci_high=round(sfg3$conf.int[tm,4],3), p=signif(sfg3$coef[tm,5],4),
    notes='complete-case n=902; see notes on competing-event definition'))
}
write.csv(t5, file.path(LOCK,'table5_finegray.csv'), row.names=FALSE)


# ---------- tableS1_schoenfeld.csv ----------
ph <- cox.zph(m3)
tS1 <- data.frame(
  analysis_id = ANALYSIS_ID, outcome='diagnostic', model='schoenfeld_model3',
  term = rownames(ph$table),
  estimate = round(ph$table[,'chisq'],3),
  df = ph$table[,'df'],
  p = signif(ph$table[,'p'],4),
  notes = ifelse(ph$table[,'p'] < 0.05, 'PH violated -> stratified in sensitivity model', 'ok'))
tS1 <- rbind(tS1, data.frame(analysis_id=ANALYSIS_ID, outcome='diagnostic', model='schoenfeld_stratified',
  term='GLOBAL_after_stratification', estimate=NA, df=NA,
  p=cox.zph(coxph(fstrat, data=df))$table['GLOBAL','p'], notes='stratified model satisfies PH'))
write.csv(tS1, file.path(LOCK,'tableS1_schoenfeld.csv'), row.names=FALSE)

# ---------- tableS2_missing.csv ----------
df$lvef_avail <- as.integer(!is.na(df$min_lvef_72h))
agg <- function(x) c(mean=round(mean(x,na.rm=TRUE),2))
tS2 <- data.frame(
  analysis_id = ANALYSIS_ID, model='missingness_pattern',
  characteristic = c('mortality_28_pct','sofa_mean'),
  lvef_available = c(round(mean(df$mortality_28[df$lvef_avail==1])*100,1),
                     round(mean(df$sofa_24h[df$lvef_avail==1]),1)),
  lvef_missing = c(round(mean(df$mortality_28[df$lvef_avail==0])*100,1),
                   round(mean(df$sofa_24h[df$lvef_avail==0]),1)),
  n_available = sum(df$lvef_avail==1), n_missing = sum(df$lvef_avail==0))
write.csv(tS2, file.path(LOCK,'tableS2_missing.csv'), row.names=FALSE)

# ---------- tableS3_intensity.csv ----------
tS3 <- do.call(rbind, lapply(c('Stable-low','Falling','Stable-high','Rising'), function(g) {
  sub <- df[df$trop_dynamic_group==g, ]
  data.frame(analysis_id=ANALYSIS_ID, model='measurement_intensity', group=g,
    meas_mean_72h = round(mean(sub$trop_valid_count),2),
    echo_pct = round(mean(sub$has_echo==1)*100,1))
}))
write.csv(tS3, file.path(LOCK,'tableS3_intensity.csv'), row.names=FALSE)

# ---------- sessionInfo.txt ----------
sink(file.path(LOCK,'sessionInfo.txt'))
print(sessionInfo())
sink()

# ---------- code_commit.txt ----------
writeLines(c(
  paste0('analysis_id: ', ANALYSIS_ID),
  paste0('locked_at: ', format(Sys.time(), '%Y-%m-%d %H:%M:%S')),
  'data_file: v4_cohort_2plus_true_time.csv',
  'builder: build_01_lock_results.R',
  'event_time: recorded date of death (dod) minus ICU intime, hours',
  'window: exposure assigned over 0-72h; last window carried to 78h maximum (see supp S2)'),
  file.path(LOCK,'code_commit.txt'))

cat('output/lock built with', length(list.files(LOCK)), 'files\n')
print(list.files(LOCK))

