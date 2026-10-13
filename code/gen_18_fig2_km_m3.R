# Figure 2: (A) unadjusted Kaplan-Meier, (B) Model 3-conditional survival,
# with a shared number-at-risk strip. Also emits Supplementary Figure SF5 (Model 2-conditional).
library(survival)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
FIG  <- file.path(BASE, 'figures')

df <- read.csv(file.path(BASE, 'data', 'v4_cohort_2plus_true_time.csv'), stringsAsFactors = FALSE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group,
                                levels = c('Stable-low', 'Falling', 'Stable-high', 'Rising'))

fit_km <- survfit(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group, data = df)

m3 <- coxph(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender +
              sofa_24h + invasive_vent + norepinephrine + creatinine_max +
              myocardial_infarct + congestive_heart_failure + lactate_max +
              log(baseline_troponin + 0.001) + min_lvef_72h, data = df)

cc2 <- complete.cases(df[, c('myocardial_infarct', 'congestive_heart_failure', 'lactate_max', 'sofa_24h')])
d2 <- df[cc2, ]
m2 <- coxph(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender +
              sofa_24h + invasive_vent + norepinephrine + creatinine_max +
              myocardial_infarct + congestive_heart_failure + lactate_max, data = d2)

nd_base <- data.frame(
  trop_dynamic_group = factor(c('Stable-low', 'Falling', 'Stable-high', 'Rising'),
                              levels = c('Stable-low', 'Falling', 'Stable-high', 'Rising')),
  age = mean(df$age, na.rm = TRUE), gender = 'M',
  sofa_24h = mean(df$sofa_24h, na.rm = TRUE), invasive_vent = 1, norepinephrine = 1,
  creatinine_max = mean(df$creatinine_max, na.rm = TRUE),
  myocardial_infarct = 0, congestive_heart_failure = 0,
  lactate_max = mean(df$lactate_max, na.rm = TRUE))
nd_full <- cbind(nd_base,
  baseline_troponin = exp(mean(log(df$baseline_troponin + 0.001), na.rm = TRUE)) - 0.001,
  min_lvef_72h = mean(df$min_lvef_72h, na.rm = TRUE))
fit_m2 <- survfit(m2, newdata = nd_base)
fit_m3 <- survfit(m3, newdata = nd_full)

colors <- c('#1a9850', '#2166ac', '#d08c00', '#d73027')
grps   <- c('Stable-low', 'Falling', 'Stable-high', 'Rising')
times_eval <- c(0, 168, 336, 504, 672)
nar <- summary(fit_km, times = times_eval, extend = TRUE)
nar_mat <- t(matrix(nar$n.risk, nrow = length(times_eval),
                    dimnames = list(paste0('t=', times_eval), grps)))

# ---------------- Figure 2: two panels + shared number-at-risk strip ----------------
png(file.path(FIG, 'fig2_km_m3.png'), width = 4400, height = 2450, res = 300)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(4.6, 1.25))

par(mar = c(4.8, 7, 3.2, 1.2), cex.axis = 1.5, cex.lab = 1.6, cex.main = 1.5)
plot(fit_km, col = colors, lwd = 3, xlab = '', ylab = 'Survival probability',
     main = 'a  Unadjusted (n=3,204)', xlim = c(0, 672), xaxs = 'i')
title(xlab = 'Hours from ICU admission', line = 2.6)
legend('bottomleft', legend = grps, col = colors, lwd = 3, bty = 'n', cex = 1.3, text.col = colors)
text(300, 1.0, expression(paste('log-rank ', italic(P), '<2×10'^'-7')), cex = 1.3, adj = 0, col = '#111111')

plot(fit_m3, col = colors, lwd = 3, xlab = '', ylab = 'Adjusted survival probability',
     main = 'b  Model 3-conditional (n=902)', xlim = c(0, 672), xaxs = 'i')
title(xlab = 'Hours from ICU admission', line = 2.6)
legend('bottomleft', legend = grps, col = colors, lwd = 3, bty = 'n', cex = 1.3, text.col = colors)

par(mar = c(0.8, 7, 0.2, 1.2))
plot.new()
plot.window(xlim = c(0, 672), ylim = c(0, 6.4))
text(-8, 5.6, 'Number at risk (unadjusted cohort)', adj = 1, font = 2, cex = 1.25, xpd = NA, col = '#111111')
text(8, 5.6, times_eval[1], adj = 0, cex = 1.25, font = 2, col = '#111111')
text(times_eval[-1], rep(5.6, length(times_eval) - 1), times_eval[-1], adj = 0.5,
     cex = 1.25, font = 2, col = '#111111')
for (g in seq_along(grps)) {
  y <- 5.6 - g * 1.15
  text(-8, y, grps[g], adj = 1, cex = 1.25, col = colors[g], xpd = NA, font = 2)
  text(8, y, format(nar_mat[g, 1], big.mark = ','), adj = 0, cex = 1.25, col = colors[g])
  text(times_eval[-1], rep(y, length(times_eval) - 1),
       format(nar_mat[g, -1], big.mark = ','), adj = 0.5, cex = 1.25, col = colors[g])
}
dev.off()

# ---------------- Supplementary Figure SF5: Model 2-conditional ----------------
png(file.path(FIG, 'figS5_m2.png'), width = 2400, height = 2000, res = 300)
par(mar = c(4.8, 7, 3.2, 1.2), cex.axis = 1.5, cex.lab = 1.6, cex.main = 1.5)
plot(fit_m2, col = colors, lwd = 3, xlab = '', ylab = 'Adjusted survival probability',
     main = 'Model 2-conditional (n=1,435)', xlim = c(0, 672), xaxs = 'i')
title(xlab = 'Hours from ICU admission', line = 2.6)
legend('bottomleft', legend = grps, col = colors, lwd = 3, bty = 'n', cex = 1.3, text.col = colors)
dev.off()

cat('saved fig2_km_m3.png and figS5_m2.png\n')
print(nar_mat)
