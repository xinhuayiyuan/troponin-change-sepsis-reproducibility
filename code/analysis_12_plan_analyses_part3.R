# Supplementary analyses, part 3: last-baseline trend test, specification checks,
# group sizes and events, attrition, number-at-risk table, standardised differences.
library(survival)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
OUT <- file.path(BASE, 'output')
df <- read.csv(file.path(BASE, 'data', 'v4_cohort_2plus_true_time.csv'), stringsAsFactors = FALSE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group,
                                levels = c('Stable-low', 'Falling', 'Stable-high', 'Rising'))
df$gender <- factor(df$gender)
S <- function(d) Surv(d$time_28d_true, d$mortality_28)
cc2 <- complete.cases(df[, c('myocardial_infarct', 'congestive_heart_failure', 'lactate_max', 'sofa_24h')])
d2 <- df[cc2, ]
d3 <- d2[!is.na(d2$min_lvef_72h), ]
fmt <- function(hr, lo, hi, p) sprintf('%.2f (%.2f-%.2f), P=%s', hr, lo, hi,
                                       ifelse(p < 0.001, '<0.001', sprintf('%.3f', p)))
m3cov <- ~ age + gender + sofa_24h + invasive_vent + norepinephrine + creatinine_max +
  myocardial_infarct + congestive_heart_failure + lactate_max + min_lvef_72h

# ---------- 1.5 Declining-change continuous test (Table S9) ----------
d3$change_last <- (d3$last_troponin - d3$baseline_troponin) / 0.1
fit_l <- coxph(update(m3cov, S(d3) ~ change_last + log(baseline_troponin + 0.001) +
                        change_last:log(baseline_troponin + 0.001) + .), data = d3)
sc <- summary(fit_l)$conf.int
sp <- summary(fit_l)$coefficients
t9 <- data.frame(
  term = c('Change (last - baseline), per 0.1 ng/mL', 'Log(baseline troponin + 0.001)',
           'Change x baseline interaction'),
  txt = c(fmt(sc['change_last', 1], sc['change_last', 3], sc['change_last', 4], sp['change_last', 5]),
          fmt(sc['log(baseline_troponin + 0.001)', 1], sc['log(baseline_troponin + 0.001)', 3],
              sc['log(baseline_troponin + 0.001)', 4], sp['log(baseline_troponin + 0.001)', 5]),
          fmt(sc['change_last:log(baseline_troponin + 0.001)', 1],
              sc['change_last:log(baseline_troponin + 0.001)', 3],
              sc['change_last:log(baseline_troponin + 0.001)', 4],
              sp['change_last:log(baseline_troponin + 0.001)', 5])))
write.csv(t9, file.path(OUT, 'tableS9_change_last.csv'), row.names = FALSE)
cat('--- Table S9 ---\n'); print(t9, right = FALSE)

# ---------- specification checks with full CI ----------
d3$logdelta <- log(d3$delta_troponin + 0.01)
f1 <- coxph(update(m3cov, S(d3) ~ logdelta + log(baseline_troponin + 0.001) + .), data = d3)
s1c <- summary(f1)$conf.int['logdelta', ]; s1p <- summary(f1)$coefficients['logdelta', 5]
d3$any_incr <- as.integer(d3$delta_troponin > 0)
f2 <- coxph(update(m3cov, S(d3) ~ any_incr + log(baseline_troponin + 0.001) + .), data = d3)
s2c <- summary(f2)$conf.int['any_incr', ]; s2p <- summary(f2)$coefficients['any_incr', 5]
cat('--- C1 log(delta+0.01):', fmt(s1c[1], s1c[3], s1c[4], s1p), '\n')
cat('--- C2 delta>0:', fmt(s2c[1], s2c[3], s2c[4], s2p),
    sprintf('(n=%d, events=%d; n with increase=%d)\n', nrow(d3), sum(d3$mortality_28), sum(d3$any_incr)))

# ---------- B10 / Table S13 group n and deaths per model ----------
g13 <- function(d, name) {
  as.data.frame(t(sapply(levels(d$trop_dynamic_group), function(g) {
    sub <- d[d$trop_dynamic_group == g, ]
    c(sample = name, group = g, n = nrow(sub), deaths = sum(sub$mortality_28))
  })))
}
s13 <- rbind(g13(df, 'Model 1 (n=3,204)'), g13(d2, 'Model 2 (n=1,435)'), g13(d3, 'Model 3 (n=902)'))
write.csv(s13, file.path(OUT, 'tableS13_group_n.csv'), row.names = FALSE)
cat('--- Table S13 ---\n'); print(s13, right = FALSE)

# ---------- C6 attrition ----------
miss_dx <- sum(is.na(df$myocardial_infarct))
rem1 <- df[!is.na(df$myocardial_infarct), ]
miss_lac <- sum(is.na(rem1$lactate_max) | is.na(rem1$sofa_24h))
rem2 <- rem1[!is.na(rem1$lactate_max) & !is.na(rem1$sofa_24h), ]
miss_lvef <- sum(is.na(rem2$min_lvef_72h))
cat(sprintf('--- C6: 3,204 -> -%d (no dx records) -> %d -> -%d (lactate/SOFA missing) -> %d -> -%d (LVEF missing) -> %d\n',
            miss_dx, nrow(rem1), miss_lac, nrow(rem2), miss_lvef, nrow(rem2) - miss_lvef))

# ---------- C5 at-risk table export ----------
fit_km <- survfit(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group, data = df)
nar <- summary(fit_km, times = c(0, 168, 336, 504, 672), extend = TRUE)
s15 <- as.data.frame(t(matrix(nar$n.risk, nrow = 5,
                              dimnames = list(paste0('t=', c(0, 168, 336, 504, 672)),
                                              levels(df$trop_dynamic_group)))))
s15 <- cbind(group = rownames(s15), s15)
write.csv(s15, file.path(OUT, 'tableS15_atrisk.csv'), row.names = FALSE)
cat('--- Table S15 exported ---\n')

# ---------- C9 SMD vs Stable-low ----------
numvars <- c('age', 'sofa_24h', 'los_icu', 'creatinine_max', 'lactate_max', 'platelets_min', 'min_lvef_72h')
binvars <- c('male', 'invasive_vent', 'norepinephrine', 'rrt', 'myocardial_infarct',
             'congestive_heart_failure', 'ischemic_heart_disease', 'hypertension', 'diabetes', 'ckd')
df$male <- as.integer(df$gender == 'M')
smd_num <- function(v, g) {
  m <- tapply(v, g, mean, na.rm = TRUE); s <- tapply(v, g, sd, na.rm = TRUE)
  sp <- sqrt((s['Stable-low']^2 + s^2) / 2)
  round(abs(m - m['Stable-low']) / sp, 2)
}
smd_bin <- function(v, g) {
  p <- tapply(v, g, mean, na.rm = TRUE)
  pp <- (p + p['Stable-low']) / 2
  round(abs(p - p['Stable-low']) / sqrt(pp * (1 - pp)), 2)
}
smd_tab <- data.frame(variable = c(numvars, binvars))
for (v in numvars) smd_tab[smd_tab$variable == v, levels(df$trop_dynamic_group)] <- as.list(smd_num(df[[v]], df$trop_dynamic_group))
for (v in binvars) smd_tab[smd_tab$variable == v, levels(df$trop_dynamic_group)] <- as.list(smd_bin(df[[v]], df$trop_dynamic_group))
write.csv(smd_tab, file.path(OUT, 'table1_smd.csv'), row.names = FALSE)
cat('--- SMD table exported ---\n')
