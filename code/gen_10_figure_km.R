# Regenerate KM + covariate-conditional survival figure with clean number-at-risk
# strip (no overlaps), larger fonts, darker colours. Output: fig2_survival.png
library(survival)

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
df <- read.csv(file.path(BASE, 'data', 'v4_cohort_2plus_true_time.csv'), stringsAsFactors = FALSE)
df$trop_dynamic_group <- factor(df$trop_dynamic_group,
                                levels = c('Stable-low', 'Falling', 'Stable-high', 'Rising'))

fit_km <- survfit(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group, data = df)

m3 <- coxph(Surv(time_28d_true, mortality_28) ~ trop_dynamic_group + age + gender +
              sofa_24h + invasive_vent + norepinephrine + creatinine_max +
              myocardial_infarct + congestive_heart_failure + lactate_max +
              log(baseline_troponin + 0.001) + min_lvef_72h, data = df)

# Model 2 fit on its own complete-case sample (n=1,435) for the attenuation triptych
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
newdata <- cbind(nd_base,
  baseline_troponin = exp(mean(log(df$baseline_troponin + 0.001), na.rm = TRUE)) - 0.001,
  min_lvef_72h = mean(df$min_lvef_72h, na.rm = TRUE))
fit_adj2 <- survfit(m2, newdata = nd_base)
fit_adj <- survfit(m3, newdata = newdata)

colors <- c('#1a9850', '#2166ac', '#d08c00', '#d73027')  # darker green/blue/orange/red
grps <- c('Stable-low', 'Falling', 'Stable-high', 'Rising')

times_eval <- c(0, 168, 336, 504, 672)
nar <- summary(fit_km, times = times_eval, extend = TRUE)
nar_mat <- t(matrix(nar$n.risk, nrow = length(times_eval),
                    dimnames = list(paste0('t=', times_eval), grps)))

png(file.path(BASE, 'figures', 'fig2_survival.png'), width = 4200, height = 1900, res = 300)
layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE), heights = c(4.2, 1.35))

par(mar = c(4.6, 9, 3, 1), cex.axis = 1.55, cex.lab = 1.65, cex.main = 1.6)

# Panel A
plot(fit_km, col = colors, lwd = 3, xlab = '', ylab = 'Survival probability',
     main = 'A  Unadjusted (n=3,204)', xlim = c(0, 672), xaxs = 'i')
title(xlab = 'Hours from ICU admission', line = 2.6)
legend('bottomleft', legend = grps, col = colors, lwd = 3, bty = 'n', cex = 1.45,
       text.col = colors)
text(380, 0.99, expression(paste('log-rank ', italic(P), '<2×10'^'-7')), cex = 1.4, adj = 0, col = '#111111')

# Panel B: Model 2 conditional
plot(fit_adj2, col = colors, lwd = 3, xlab = '', ylab = 'Adjusted survival probability',
     main = 'B  Model 2-conditional (n=1,435)', xlim = c(0, 672), xaxs = 'i')
title(xlab = 'Hours from ICU admission', line = 2.6)
legend('bottomleft', legend = grps, col = colors, lwd = 3, bty = 'n', cex = 1.45,
       text.col = colors)

# Panel C: Model 3 conditional
plot(fit_adj, col = colors, lwd = 3, xlab = '', ylab = 'Adjusted survival probability',
     main = 'C  Model 3-conditional (n=902)', xlim = c(0, 672), xaxs = 'i')
title(xlab = 'Hours from ICU admission', line = 2.6)
legend('bottomleft', legend = grps, col = colors, lwd = 3, bty = 'n', cex = 1.45,
       text.col = colors)

# Number-at-risk strip: same left margin (8 lines) as panels so x=0 aligns;
# labels right-aligned in the margin ending at x=-8 (never reach the numbers)
par(mar = c(0.8, 9, 0.2, 1))
plot.new()
plot.window(xlim = c(0, 672), ylim = c(0, 6.4))
text(-8, 5.6, 'Number at risk', adj = 1, font = 2, cex = 1.3, xpd = NA, col = '#111111')
text(8, 5.6, times_eval[1], adj = 0, cex = 1.3, font = 2, col = '#111111')
text(times_eval[-1], rep(5.6, length(times_eval) - 1), times_eval[-1], adj = 0.5,
     cex = 1.3, font = 2, col = '#111111')
for (g in seq_along(grps)) {
  y <- 5.6 - g * 1.15
  text(-8, y, grps[g], adj = 1, cex = 1.3, col = colors[g], xpd = NA, font = 2)
  # first column left-adjusted at x=8; remaining columns centred on time points
  text(8, y, format(nar_mat[g, 1], big.mark = ','), adj = 0, cex = 1.3, col = colors[g])
  text(times_eval[-1], rep(y, length(times_eval) - 1),
       format(nar_mat[g, -1], big.mark = ','), adj = 0.5, cex = 1.3, col = colors[g])
}
dev.off()

cat('saved fig2_survival.png\n')
print(nar_mat)
