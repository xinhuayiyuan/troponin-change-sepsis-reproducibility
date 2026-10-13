# Supplementary figures: S1 baseline cTnT dose-response (RCS), S2 90-day KM,
# S3 Schoenfeld residual plots for PH-violating covariates.
library(survival)
library(rms)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
FIG <- file.path(BASE, 'figures')
df <- read.csv(file.path(BASE, 'data', 'v4_cohort_2plus_true_time.csv'), stringsAsFactors = FALSE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group,
                                levels = c('Stable-low', 'Falling', 'Stable-high', 'Rising'))
df$gender <- factor(df$gender)
cc2 <- complete.cases(df[, c('myocardial_infarct', 'congestive_heart_failure', 'lactate_max', 'sofa_24h')])
d3 <- df[cc2 & !is.na(df$min_lvef_72h), ]
d3$log_base <- log(d3$baseline_troponin + 0.001)
colors <- c('#1a9850', '#2166ac', '#d08c00', '#d73027')
grps <- levels(df$trop_dynamic_group)

# ---------- Figure S1: baseline cTnT dose-response (RCS) ----------
dd <- datadist(d3); options(datadist = 'dd')
f_rcs <- cph(Surv(time_28d_true, mortality_28) ~ rcs(log_base, 4) + age + gender + sofa_24h +
               invasive_vent + norepinephrine + creatinine_max + myocardial_infarct +
               congestive_heart_failure + lactate_max + min_lvef_72h, data = d3, x = TRUE, y = TRUE)
png(file.path(FIG, 'figS1_rcs_baseline.png'), width = 2800, height = 1900, res = 300)
par(mar = c(5, 5, 3, 1), cex.axis = 1.15, cex.lab = 1.25, cex.main = 1.25)
p <- as.data.frame(Predict(f_rcs, log_base, fun = exp, ref.zero = TRUE))
xlb <- 'Log(baseline troponin + 0.001), ng/mL'
yl <- range(c(p$lower, p$upper), finite = TRUE)
plot(p$log_base, p$yhat, type = 'l', lwd = 2.5, col = '#2166ac', ylim = yl,
     xlab = xlb, ylab = 'Adjusted hazard ratio (95% CI)',
     main = 'Dose-response of baseline troponin (Model 3, n=902)')
lines(p$log_base, p$lower, lty = 2, col = '#2166ac', lwd = 1.5)
lines(p$log_base, p$upper, lty = 2, col = '#2166ac', lwd = 1.5)
abline(h = 1, lty = 3, col = '#333333')
rug(d3$log_base, side = 1, col = '#666666')
dev.off()
cat('saved figS1_rcs_baseline.png\n')

# ---------- Figure S2: 90-day KM ----------
fit90 <- survfit(Surv(time_90d_true, mortality_90) ~ trop_dynamic_group, data = df)
png(file.path(FIG, 'figS2_km90.png'), width = 2800, height = 1900, res = 300)
par(mar = c(4.6, 5, 3, 1), cex.axis = 1.15, cex.lab = 1.25, cex.main = 1.3)
plot(fit90, col = colors, lwd = 2.5, xlab = '', ylab = '90-day survival probability',
     main = 'Unadjusted Kaplan-Meier curves, 90-day mortality', xlim = c(0, 2160), xaxs = 'i')
title(xlab = 'Hours from ICU admission', line = 2.6)
legend('bottomleft', legend = grps, col = colors, lwd = 2.5, bty = 'n', cex = 1.05, text.col = colors)
dev.off()
cat('saved figS2_km90.png\n')

# ---------- Figure S3: Schoenfeld residuals ----------
m3 <- coxph(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender + sofa_24h +
              invasive_vent + norepinephrine + creatinine_max + myocardial_infarct +
              congestive_heart_failure + lactate_max + log(baseline_troponin + 0.001) +
              min_lvef_72h, data = df)
zph <- cox.zph(m3)
png(file.path(FIG, 'figS3_schoenfeld.png'), width = 3000, height = 1250, res = 300)
par(mfrow = c(1, 3), mar = c(4.6, 4.6, 3, 1), cex.axis = 1.1, cex.lab = 1.15, cex.main = 1.2)
for (v in c('sofa_24h', 'norepinephrine', 'lactate_max')) {
  nice <- c(sofa_24h = 'SOFA score', norepinephrine = 'Norepinephrine', lactate_max = 'Lactate')[[v]]
  pv <- zph$table[v, 'p']
  plab <- if (pv < 0.001) 'P<0.001' else sprintf('P=%.3f', pv)
  # 1x3 layout avoids an empty panel; the horizontal axis is labelled explicitly as
  # the 1-KM scale used by plot.cox.zph (tick values are times in hours).
  plot(zph, var = v, xlab = 'Time (hours; axis on the 1-KM scale)',
       ylab = paste('Scaled Schoenfeld residuals:', nice),
       main = paste0(nice, ' (', plab, ')'))
  abline(h = 0, lty = 2, col = '#666666')
}
dev.off()
cat('saved figS3_schoenfeld.png\n')
