# =============================================================================
# eICU-CRD site-gradient external validation
#   Exposure : 0-72 h log-fold change RATE of cTnI, r (same estimand as MIMIC track A)
#   Landmark : 72 h (still in hospital, hence alive)
#   Outcome  : in-hospital death within 28 days after the landmark
#   Design   : sites stratified into tertiles of troponin MEASUREMENT INTENSITY;
#              r x intensity interaction with cluster-robust SE by hospitalid
#
# Inputs (produced by the operator running 07_Scripts/eicu_extract_01_cohort.sql):
#   data/eicu_extract_cohort_stays.csv   (required)
#   data/eicu_extract_site_unit_lod.csv  (required - unit/LOD audit)
#   data/eicu_extract_site_hospital.csv  (optional)
#   data/eicu_extract_apache.csv         (optional - severity)
# Output: output/lock/eicu/*.csv
# No value is imputed or invented; every number below is produced by these models.
# =============================================================================
suppressPackageStartupMessages({
  library(survival)
})

BASE <- Sys.getenv("PROJECT_ROOT", unset = ".")
DAT  <- file.path(BASE, 'data')
OUT  <- file.path(BASE, 'output', 'lock', 'eicu')
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

f_stays <- file.path(DAT, 'eicu_extract_cohort_stays.csv')
f_lod   <- file.path(DAT, 'eicu_extract_site_unit_lod.csv')
f_hosp  <- file.path(DAT, 'eicu_extract_site_hospital.csv')
f_apach <- file.path(DAT, 'eicu_extract_apache.csv')
f_win   <- file.path(DAT, 'eicu_extract_window72.csv')

if (!file.exists(f_stays)) {
  stop('Missing ', f_stays, '\nRun 07_Scripts/eicu_extract_01_cohort.sql PART 1 in Navicat/pgAdmin ',
       'and export to data/ (CSV, UTF-8, header, comma).')
}

d <- read.csv(f_stays, stringsAsFactors = FALSE)
cat('rows read:', nrow(d), '\n')

# ---- minimal, documented cleaning ------------------------------------------
num <- function(x) suppressWarnings(as.numeric(x))
for (v in c('age_num', 'unitdischargeoffset', 'hospitaldischargeoffset', 'n_meas',
            'first_offset_min', 'last_offset_min', 'first_val', 'last_val',
            'n_rows_72h', 'n_null', 'n_nonpos', 'n_unit_implausible',
            'n_below_detect', 'max_val', 'in_hosp_death')) {
  if (v %in% names(d)) d[[v]] <- num(d[[v]])
}
d$female <- as.integer(tolower(d$gender) %in% c('female', 'f'))
d$stays_any <- 1L                       # every exported row has >=1 cTnI (n_meas >= 1)

cat('n_meas distribution:\n'); print(table(d$n_meas, useNA = 'ifany'))
cat('hospitals:', length(unique(d$hospitalid)), '\n')

# ---- flow (real counts) ----------------------------------------------------
land <- d[which(d$n_meas >= 2 & !is.na(d$first_val) & !is.na(d$last_val) &
                  d$first_val > 0 & d$last_val > 0 &
                  !is.na(d$first_offset_min) & !is.na(d$last_offset_min) &
                  d$last_offset_min > d$first_offset_min), ]
flow <- data.frame(
  step = c('eligible stays with >=1 cTnI in 0-72h (exported)',
           'with >=2 cTnI in 0-72h (landmark cohort)',
           'with a positive measurement span (dt > 0)'),
  n = c(nrow(d), sum(d$n_meas >= 2), nrow(land)))
print(flow)

# ---- exposure --------------------------------------------------------------
land$log_baseline <- log(land$first_val)
land$dt_hours <- (land$last_offset_min - land$first_offset_min) / 60
land$r <- (log(land$last_val) - log(land$first_val)) / land$dt_hours * 72
q <- quantile(land$r, c(0.01, 0.99), na.rm = TRUE)
land$r_w <- pmin(pmax(land$r, q[1]), q[2])
SD_R <- sd(land$r_w, na.rm = TRUE)
cat('rate winsor limits:', round(q, 4), '| sd(r_w) =', round(SD_R, 4), '\n')

# ---- outcome (landmark-relative, in-hospital death) ------------------------
land$t_days <- (land$hospitaldischargeoffset - 4320) / 1440
land <- land[which(!is.na(land$t_days) & land$t_days > 0), ]
land$event <- as.integer(land$in_hosp_death == 1 & land$t_days <= 28)
land$t <- pmin(land$t_days, 28)
n_landmark_valid <- nrow(land)
n_death_landmark <- sum(land$event)
n_late_censored <- sum(land$in_hosp_death == 1 & land$t_days > 28)
cat('landmark cohort n =', nrow(land), ' events =', sum(land$event),
    ' late (>28d) censored =', n_late_censored, '\n')

# ---- site-level measurement intensity (over ALL exported stays) ------------
site <- aggregate(cbind(n_stays = d$stays_any, n_ge2 = as.integer(d$n_meas >= 2)) ~ hospitalid,
                  data = d, FUN = sum)
site$prop_ge2 <- site$n_ge2 / site$n_stays
mm <- aggregate(n_meas ~ hospitalid, data = d, FUN = function(x) mean(x, na.rm = TRUE))
site <- merge(site, mm, by = 'hospitalid', all.x = TRUE)
names(site)[names(site) == 'n_meas'] <- 'mean_meas'

MIN_SITE <- 20                                   # pre-specified: drop low-volume sites
site <- site[order(-site$n_stays), ]
site_all <- site
site <- site[which(site$n_stays >= MIN_SITE), ]
cat('sites total =', nrow(site_all), '| retained (>=', MIN_SITE, 'stays) =', nrow(site), '\n')

br <- quantile(site$prop_ge2, c(1/3, 2/3))
site$tertile <- cut(site$prop_ge2, c(-Inf, br, Inf), labels = c('low', 'mid', 'high'))
cat('tertile cut points (prop with >=2 cTnI):', round(br, 4), '\n')
print(table(site$tertile))

land$intensity <- site$prop_ge2[match(land$hospitalid, site$hospitalid)]
land$tertile <- site$tertile[match(land$hospitalid, site$hospitalid)]

# low-volume sites are dropped from EVERY model (pre-specified, not per-analysis)
n_before <- nrow(land)
land <- land[which(!is.na(land$intensity)), ]
n_analysis <- nrow(land)
n_events_analysis <- sum(land$event)
cat('dropped rows in low-volume sites:', n_before - n_analysis,
    '| analysis cohort n =', n_analysis, ' events =', n_events_analysis, '\n')

# ---- optional APACHE severity ----------------------------------------------
have_sev <- FALSE
if (file.exists(f_apach)) {
  ap <- read.csv(f_apach, stringsAsFactors = FALSE)
  ap$pred_hosp_mortality <- num(ap$pred_hosp_mortality)
  if (max(ap$pred_hosp_mortality, na.rm = TRUE) > 1) ap$pred_hosp_mortality <- ap$pred_hosp_mortality / 100
  land$pred_mort <- ap$pred_hosp_mortality[match(land$patientunitstayid, ap$patientunitstayid)]
  have_sev <- sum(!is.na(land$pred_mort)) > 0.5 * nrow(land)
  cat('APACHE predicted mortality joined:', sum(!is.na(land$pred_mort)), '/', nrow(land), '\n')
}

# ---- optional 72 h clinical-state proxy (vasopressor / lactate) -------------
# eICU infusiondrug has only infusionoffset (no start/stop), so vaso72 is a
# snapshot proxy: a vasopressor entry charted inside the 48-72 h window.
have_win <- FALSE
ok_lac <- FALSE
if (file.exists(f_win)) {
  w <- read.csv(f_win, stringsAsFactors = FALSE)
  w$vaso72 <- num(w$vaso72); w$lac72 <- num(w$lac72)
  land$vaso72 <- w$vaso72[match(land$patientunitstayid, w$patientunitstayid)]
  land$lac72  <- w$lac72[match(land$patientunitstayid, w$patientunitstayid)]
  land$vaso72[is.na(land$vaso72)] <- 0        # no infusion record == no vasopressor
  ok_lac <- mean(!is.na(land$lac72)) >= 0.5    # only use lactate if coverage >= 50%
  have_win <- TRUE
  cat('72h window proxy joined: vaso72 = 1 in', sum(land$vaso72 == 1),
      '| lac72 available', sum(!is.na(land$lac72)), '/', nrow(land),
      '| lactate used in model:', ok_lac, '\n')
}

# ---- models (cluster-robust SE by hospitalid) -------------------------------
grab <- function(m, term) {
  cf <- coef(m); se <- sqrt(diag(vcov(m)))
  i <- which(names(cf) == term)
  if (!length(i)) return(data.frame(term = term, HR = NA, lo = NA, hi = NA, p = NA))
  data.frame(term = term, HR = exp(cf[i]),
             lo = exp(cf[i] - 1.96 * se[i]), hi = exp(cf[i] + 1.96 * se[i]),
             p = summary(m)$coefficients[i, ncol(summary(m)$coefficients)])
}

adj_cov <- c('age_num', 'female', 'log_baseline', 'n_meas')
if (have_sev) adj_cov <- c(adj_cov, 'pred_mort')
# 72 h clinical state, mirroring the stage at which the MIMIC association vanished
state_cov <- character(0)
if (have_win) {
  state_cov <- 'vaso72'
  if (ok_lac) state_cov <- c(state_cov, 'lac72')
}
full_cov <- c(adj_cov, state_cov)
f_cru  <- Surv(t, event) ~ r_w
f_adj  <- as.formula(paste('Surv(t, event) ~ r_w +', paste(adj_cov, collapse = ' + ')))
f_full <- as.formula(paste('Surv(t, event) ~ r_w +', paste(full_cov, collapse = ' + ')))

mk <- function(dat, f, label, cluster = TRUE) {
  if (nrow(dat) < 50 || sum(dat$event) < 20) return(NULL)
  m <- if (cluster) coxph(f, data = dat, cluster = hospitalid, x = TRUE) else coxph(f, data = dat, x = TRUE)
  g <- grab(m, 'r_w')
  # per-1-SD effect on the same scale used for the pre-specified MCID (HR/SD >= 1.20)
  out <- data.frame(model = label, n = m$n, events = m$nevent, g)
  out$HR_per_SD  <- out$HR^SD_R
  out$lo_per_SD  <- out$lo^SD_R
  out$hi_per_SD  <- out$hi^SD_R
  out
}

res <- list()
res[['crude']] <- mk(land, f_cru,  'eICU crude (cluster by hospital)')
res[['adj']]   <- mk(land, f_adj,  'eICU + admission-level covariates (cluster)')
res[['full']]  <- mk(land, f_full, 'eICU + admission + 72h clinical state (cluster)')
for (g in levels(site$tertile)) {
  sub <- land[which(land$tertile == g), ]
  res[[paste0('tert_', g)]] <- mk(sub, f_full, paste0('eICU tertile ', g, ' (full adj, cluster)'))
}
main <- do.call(rbind, res[!vapply(res, is.null, TRUE)])
print(main)

# ---- interaction: rate x site measurement intensity -------------------------
# NOTE: anova.coxph() refuses models with robust (cluster) variances, so the
# joint test on the interaction block is a Wald test using the robust vcov.
int_rows <- list()
land$intensity_c <- land$intensity - mean(land$intensity, na.rm = TRUE)
f_int <- as.formula(paste('Surv(t, event) ~ r_w * intensity_c +', paste(full_cov, collapse = ' + ')))
mi <- coxph(f_int, data = land, cluster = hospitalid, x = TRUE)
mkint <- function(g, spec, n, ev, chi2, df) {
  cbind(data.frame(specification = spec, n = n, events = ev), g,
        data.frame(chi2 = chi2, df = df))
}
se_ci <- function(g) (log(g$hi) - log(g$lo)) / (2 * 1.96)
g_cont <- grab(mi, 'r_w:intensity_c')
zi <- log(g_cont$HR) / se_ci(g_cont)
int_rows[['cont']] <- rbind(
  mkint(grab(mi, 'r_w'),         'r_w x intensity (continuous): main effect of rate',       mi$n, mi$nevent, NA_real_, NA_real_),
  mkint(grab(mi, 'intensity_c'), 'r_w x intensity (continuous): main effect of intensity',  mi$n, mi$nevent, NA_real_, NA_real_),
  mkint(g_cont,                  'r_w x intensity (continuous): INTERACTION',               mi$n, mi$nevent, zi^2, 1))

f_intf <- as.formula(paste('Surv(t, event) ~ r_w * tertile +', paste(full_cov, collapse = ' + ')))
mif <- coxph(f_intf, data = land, cluster = hospitalid, x = TRUE)
idx <- grep('^r_w:tertile', names(coef(mif)))
if (length(idx)) {
  b <- coef(mif)[idx]; V <- vcov(mif)[idx, idx, drop = FALSE]
  chi2 <- as.numeric(t(b) %*% solve(V) %*% b)
  p_fac <- pchisq(chi2, length(idx), lower.tail = FALSE)
} else {
  chi2 <- NA_real_; p_fac <- NA_real_
}
int_rows[['fac']] <- data.frame(
  specification = 'r_w x intensity tertile (robust Wald on interaction block)',
  n = mif$n, events = mif$nevent, term = 'r_w:tertile (joint)',
  HR = NA_real_, lo = NA_real_, hi = NA_real_, p = p_fac, chi2 = chi2, df = length(idx))
interaction_out <- do.call(rbind, int_rows)
print(interaction_out)

# ---- unit / LOD audit -------------------------------------------------------
if (file.exists(f_lod)) {
  lod <- read.csv(f_lod, stringsAsFactors = FALSE)
  lod$implausible_units <- grepl('mg/ml|mg/dl', tolower(lod$interface_units))
  lod$prop_nonpos <- lod$n_nonpos / lod$ctni_rows
  if (file.exists(f_hosp)) {
    hs <- read.csv(f_hosp, stringsAsFactors = FALSE)
    lod <- merge(lod, hs, by = 'hospitalid', all.x = TRUE)
  }
  write.csv(lod, file.path(OUT, 'eicu_units_lod_by_site.csv'), row.names = FALSE)
  cat('sites with any implausible (mg/mL) interface unit:',
      sum(lod$implausible_units), '/', nrow(lod), '\n')
}

# ---- write outputs ----------------------------------------------------------
# complete flow (the preliminary `flow` above only covered the first three steps)
flow_full <- data.frame(
  step = c('eligible stays with >=1 cTnI in 0-72h (exported)',
           'with >=2 cTnI in 0-72h and a positive span (landmark cohort)',
           '- in-hospital death within 28 days after the landmark',
           '- censored more than 28 days after the landmark',
           'sites with any cTnI in the 0-72h window',
           paste0('sites retained (pre-specified minimum ', MIN_SITE, ' stays)'),
           'analysis cohort after dropping low-volume sites',
           '- events in the analysis cohort',
           '- stays dropped with low-volume sites'),
  n = c(nrow(d), n_landmark_valid, n_death_landmark, n_late_censored,
        nrow(site_all), nrow(site), n_analysis, n_events_analysis,
        n_before - n_analysis))
print(flow_full)
write.csv(flow_full, file.path(OUT, 'eicu_flow.csv'), row.names = FALSE)
write.csv(site, file.path(OUT, 'eicu_site_intensity.csv'), row.names = FALSE)
write.csv(main, file.path(OUT, 'eicu_site_gradient_models.csv'), row.names = FALSE)
write.csv(interaction_out, file.path(OUT, 'eicu_intensity_interaction.csv'), row.names = FALSE)
cat('\nDONE. outputs in', OUT, '\n')

