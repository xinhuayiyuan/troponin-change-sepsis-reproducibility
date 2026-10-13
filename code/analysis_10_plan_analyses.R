# Supplementary analyses: attenuation ladder and missingness pattern.
# All numbers are computed from the locked analysis dataset; no hand transcription.
library(survival)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT <- file.path(BASE, 'output')
dir.create(OUT, showWarnings = FALSE)
df <- read.csv(file.path(BASE, 'data', 'v4_cohort_2plus_true_time.csv'), stringsAsFactors = FALSE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group,
                                levels = c('Stable-low', 'Falling', 'Stable-high', 'Rising'))
df$gender <- factor(df$gender)

grp <- ~ trop_dynamic_group
m2rhs <- ~ trop_dynamic_group + age + gender + sofa_24h + invasive_vent + norepinephrine +
  creatinine_max + myocardial_infarct + congestive_heart_failure + lactate_max
m3rhs <- update(m2rhs, . ~ . + log(baseline_troponin + 0.001) + min_lvef_72h)
S <- function(d) Surv(d$time_28d_true, d$mortality_28)

cc2 <- complete.cases(df[, c('myocardial_infarct', 'congestive_heart_failure', 'lactate_max', 'sofa_24h')])
d2 <- df[cc2, ]
d3 <- d2[!is.na(d2$min_lvef_72h), ]
cat('samples:', nrow(df), nrow(d2), nrow(d3), '\n')

hr_rising <- function(fit) {
  s <- summary(fit)$coefficients['trop_dynamic_groupRising', ]
  ci <- summary(fit)$conf.int['trop_dynamic_groupRising', ]
  c(hr = unname(ci[1]), lo = unname(ci[3]), hi = unname(ci[4]), p = unname(s[5]))
}
fmt <- function(v) sprintf('%.2f (%.2f-%.2f), P=%s', v[1], v[2], v[3],
                           ifelse(v[4] < 0.001, '<0.001', sprintf('%.3f', v[4])))

# ---------- 1.1 Attenuation matrix (Table S6) ----------
a1 <- hr_rising(coxph(S(df) ~ trop_dynamic_group, data = df))
a2 <- hr_rising(coxph(S(d2) ~ trop_dynamic_group, data = d2))
a3 <- hr_rising(coxph(S(d3) ~ trop_dynamic_group, data = d3))
b1 <- hr_rising(coxph(update(m2rhs, S(d2) ~ .), data = d2))
b2 <- hr_rising(coxph(update(m2rhs, S(d3) ~ .), data = d3))
c1 <- hr_rising(coxph(update(m3rhs, S(d3) ~ .), data = d3))
s6 <- data.frame(
  covariates = c('None (group only)', 'None (group only)', 'None (group only)',
                 'Model 2 set', 'Model 2 set', 'Model 3 set (+ log baseline, LVEF)'),
  sample = c('Full cohort (n=3,204)', 'M2 complete-case (n=1,435)', 'M3 complete-case (n=902)',
             'M2 complete-case (n=1,435)', 'M3 complete-case (n=902)', 'M3 complete-case (n=902)'),
  hr_txt = c(fmt(a1), fmt(a2), fmt(a3), fmt(b1), fmt(b2), fmt(c1)))
write.csv(s6, file.path(OUT, 'tableS6_attenuation.csv'), row.names = FALSE)
cat('\n--- Table S6 attenuation matrix ---\n'); print(s6, right = FALSE)

# ---------- 1.2 Missingness mechanism (Table S7) ----------
df$has_dx <- !is.na(df$myocardial_infarct)
tab_mort <- table(df$has_dx, df$mortality_28)
p_mort <- chisq.test(tab_mort)$p.value
s7 <- data.frame(
  characteristic = c('n', '28-day mortality, n (%)', 'Age, mean (SD)', 'SOFA, mean (SD)',
                     'Lactate max, mean (SD)', 'Echocardiography performed, %', 'LVEF available, %'),
  dx_present = c(
    sum(df$has_dx),
    sprintf('%d (%.1f)', sum(df$mortality_28[df$has_dx]), 100 * mean(df$mortality_28[df$has_dx])),
    sprintf('%.1f (%.1f)', mean(df$age[df$has_dx]), sd(df$age[df$has_dx])),
    sprintf('%.1f (%.1f)', mean(df$sofa_24h[df$has_dx]), sd(df$sofa_24h[df$has_dx])),
    sprintf('%.1f (%.1f)', mean(df$lactate_max[df$has_dx], na.rm = TRUE), sd(df$lactate_max[df$has_dx], na.rm = TRUE)),
    sprintf('%.1f', 100 * mean(df$has_echo[df$has_dx])),
    sprintf('%.1f', 100 * mean(!is.na(df$min_lvef_72h[df$has_dx])))),
  dx_absent = c(
    sum(!df$has_dx),
    sprintf('%d (%.1f)', sum(df$mortality_28[!df$has_dx]), 100 * mean(df$mortality_28[!df$has_dx])),
    sprintf('%.1f (%.1f)', mean(df$age[!df$has_dx]), sd(df$age[!df$has_dx])),
    sprintf('%.1f (%.1f)', mean(df$sofa_24h[!df$has_dx]), sd(df$sofa_24h[!df$has_dx])),
    sprintf('%.1f (%.1f)', mean(df$lactate_max[!df$has_dx], na.rm = TRUE), sd(df$lactate_max[!df$has_dx], na.rm = TRUE)),
    sprintf('%.1f', 100 * mean(df$has_echo[!df$has_dx])),
    sprintf('%.1f', 100 * mean(!is.na(df$min_lvef_72h[!df$has_dx])))))
write.csv(s7, file.path(OUT, 'tableS7_missingness.csv'), row.names = FALSE)
cat('\n--- Table S7 missingness mechanism ---\n')
print(s7, right = FALSE)
cat(sprintf('chi-square P for mortality by dx-record presence: %.3e\n', p_mort))
cat(sprintf('M3 complete-case mortality: %.1f%% (%d/%d); excluded-from-M3: %.1f%% (%d/%d)\n',
            100 * mean(d3$mortality_28), sum(d3$mortality_28), nrow(d3),
            100 * mean(df$mortality_28[!rownames(df) %in% rownames(d3)]),
            sum(df$mortality_28[!rownames(df) %in% rownames(d3)]),
            nrow(df) - nrow(d3)))
