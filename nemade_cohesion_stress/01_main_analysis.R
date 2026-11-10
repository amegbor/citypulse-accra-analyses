## 01_main_analysis.R ----------------------------------------------------------------------------
## Social capital, neighbourhood problems and perceived stress: CityPulse Accra survey (2024-25).
## Scales and reliability, reconstructed sampling clusters, leave-self-out neighbourhood means, multilevel
## models A-G with model-based (REML, Wald) inference, model-based sensitivity analyses and descriptives.
## CR2 cluster-robust inference, used for the reported results, is computed in 03_robust_inference.R.
## Outputs are aggregate tables (outputs/tables) and figures (outputs/figures); the only participant-level
## output is the participant ID -> cluster lookup in outputs/derived. No coordinates are written.
##
## Sampling clusters: no enumeration-area, household or interviewer identifier is available, so clusters
## are reconstructed from home coordinates (complete linkage within district; cut height giving a median
## cluster size close to 10; singletons merged into the nearest same-district cluster within 1 km).
## Between-neighbourhood components use leave-self-out means within 500 m; within-cluster means are a
## sensitivity check.

source("config.R")
set.seed(20241001)
pk <- c("lme4", "performance", "psych", "lavaan", "ordinal", "sandwich", "lmtest", "spdep", "sf", "ggplot2")
miss <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install first: ", paste(miss, collapse = ", "))
library(lme4); library(ggplot2)        # others called with :: (ordinal masks lme4's ranef/VarCorr)

tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
der_dir <- file.path(out_dir, "derived")
for (d in c(tab_dir, fig_dir, der_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

covars_ind <- c("age", "female", "educ_f", "wealth_quintile", "Marital", "Emp", "YrsInHome")
crs_ll <- 4326     # WGS84 longitude/latitude (survey coordinates)
crs_m  <- 32630    # WGS84 / UTM zone 30N (metres) for all distances
near_m <- 500      # radius (m) for leave-self-out neighbour means

p <- readRDS(data_file)

## 1. Items and variable mapping ------------------------------------------------------------------
## Problem items: "Concerned with ..." (0-4). Trust items: agree (1-5). Cooperation items: likely (1-5).
pn_items <- c(PN_Noise = "Concerned with Noise", PN_Air = "Concerned with Air Quality",
              PN_Unemp = "Concerned with Unemployment", PN_Litter = "Concerned with Littering and Sanitation",
              PN_Robbery = "Concerned with Theft and Robbery", PN_Drugs = "Concerned with Drugs and Alcohol",
              PN_Pov = "Concerned with Poverty", PN_SexViol = "Concerned with Sexual Assault",
              PN_Vandalism = "Concerned with Vandalism", PN_Housing = "Concerned with Housing and Rent",
              PN_Disaster = "Concerned with Disasters (e.g. floods)", PN_Homeless = "Concerned with Homelessness",
              PN_Prost = "Concerned with Prostitution", PN_Safety = "Concerned about Safety")
sc_items <- c(PN_TS1  = "Most people who live in this neighbourhood can be trusted.",
              PN_TS2  = "In this neighbourhood, one has to be alert or someone is likely to take advantage of you.",
              PN_TS3  = "Most people in this neighbourhood are willing to help if you need it.",
              PN_TS4  = "In this neighbourhood, people generally do not trust each other in matters of lending and borrowing money.",
              PN_CAC1 = "If there is a problem, how likely is it that people in your neighbourhood will get together to help: problem in the community, e.g. water supply",
              PN_CAC2 = "If there is a problem, how likely is it that people in your neighbourhood will get together to help: unfortunate events, e.g. serious illness or death of a relative")
coh_items <- c("PN_TS1", "PN_TS3", "PN_CAC1", "PN_CAC2")
dis_items <- c("PN_TS2", "PN_TS4")
pss_items <- c("PSS_1", "PSS_2", "PSS_3", "revPSS_4", "revPSS_5", "PSS_6", "revPSS_7", "revPSS_8", "PSS_9", "PSS_10")
subscales <- list(crime = c("PN_Robbery", "PN_Drugs", "PN_SexViol", "PN_Vandalism", "PN_Prost", "PN_Safety"),
                  envir = c("PN_Noise", "PN_Air", "PN_Litter", "PN_Disaster"),
                  socio = c("PN_Unemp", "PN_Pov", "PN_Housing", "PN_Homeless"))

rest_coh <- function(v) rowMeans(p[setdiff(coh_items, v)])   # cohesion items other than v
mapping <- do.call(rbind, lapply(names(c(pn_items, sc_items)), function(v) {
  x <- as.numeric(p[[v]])
  data.frame(variable = v, wording = c(pn_items, sc_items)[[v]],
             role = if (v %in% names(pn_items)) "problems" else if (v %in% coh_items) "cohesion" else "distrust (negatively worded)",
             min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE), n_missing = sum(is.na(x)),
             r_with_cohesion_rest = if (v %in% names(sc_items)) cor(x, rest_coh(v), use = "complete.obs") else NA)
}))
write.csv(mapping, file.path(tab_dir, "variable_mapping.csv"), row.names = FALSE)

## 2. Analytic sample ------------------------------------------------------------------------------
## Ethnicity: groups with < 50 respondents (Grusi, Gurma, Mande) folded into Other.
p$ethn5 <- factor(ifelse(as.character(p$Ethn) %in% c("Grusi", "Gurma", "Mande"), "Other", as.character(p$Ethn)))
p$ethn5 <- relevel(p$ethn5, ref = grep("^Akan", levels(p$ethn5), value = TRUE))
p$DISTRICT <- factor(p$DISTRICT)
steps <- list("PSS-10 score" = "PSS_score", "Neighbourhood problem items" = names(pn_items),
              "Trust/cooperation items" = names(sc_items), "Covariates" = covars_ind,
              "Ethnicity" = "ethn5", "District" = "DISTRICT", "Home coordinates" = c("LONG", "LAT"))
keep <- rep(TRUE, nrow(p)); excl <- NULL
for (s in names(steps)) {
  drop <- keep & !complete.cases(p[steps[[s]]])
  excl <- rbind(excl, data.frame(section = "exclusions", variable = s, n = sum(drop)))
  keep <- keep & !drop
}
A <- p[keep, ]
z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
A$problems <- rowMeans(A[names(pn_items)])
A$cohesion <- rowMeans(A[coh_items])
A$distrust <- rowMeans(A[dis_items])                # higher = more distrust (items not reversed)
for (s in names(subscales)) A[[s]] <- rowMeans(A[subscales[[s]]])
A$pss_z <- z(A$PSS_score); A$prob_z <- z(A$problems); A$coh_z <- z(A$cohesion); A$dis_z <- z(A$distrust)
A$crime_z <- z(A$crime); A$envir_z <- z(A$envir); A$socio_z <- z(A$socio)
A$poor   <- as.numeric(A$wealth_quintile %in% c("Poorest", "Poorer"))
A$hiedu  <- as.numeric(as.character(A$educ_f) %in% c("4", "5"))   # completed secondary or post-secondary
A$pss_hi <- as.numeric(A$PSS_score > 13)

## 3. Reconstructed sampling clusters (full sample) ------------------------------------------------
ok <- complete.cases(p[c("LONG", "LAT", "DISTRICT")])
xy <- matrix(NA_real_, nrow(p), 2)
xy[ok, ] <- sf::st_coordinates(sf::st_transform(sf::st_as_sf(p[ok, c("LONG", "LAT")], coords = c("LONG", "LAT"), crs = crs_ll), crs_m))
by_d  <- split(which(ok), droplevels(p$DISTRICT[ok]))
trees <- lapply(by_d, function(i) if (length(i) > 1) hclust(dist(xy[i, ]), method = "complete"))
cut_at <- function(h) {
  cl <- rep(NA_integer_, nrow(p)); nxt <- 0L
  for (d in names(by_d)) {
    i <- by_d[[d]]; k <- if (length(i) > 1) cutree(trees[[d]], h = h) else 1L
    cl[i] <- k + nxt; nxt <- nxt + max(k)
  }
  cl
}
grid <- t(sapply(seq(100, 3000, by = 10), function(h) { sz <- table(cut_at(h)); c(h = h, k = length(sz), med = median(sz)) }))
## cut height: smallest whose median cluster size is closest to 10 (the design's 300 enumeration areas are
## only a rough guide to the number of clusters)
h_star <- grid[which.min(abs(grid[, "med"] - 10)), "h"]
cl0 <- cut_at(h_star)
single <- as.integer(names(which(table(cl0) == 1)))
## singleton -> cluster of the nearest member of a non-singleton cluster in the same district, if within max_m
merge_singletons <- function(max_m) {
  cl <- cl0; d_m <- c()
  for (s in single) {
    i <- which(cl0 == s)
    cand <- which(ok & p$DISTRICT == p$DISTRICT[i] & !(cl0 %in% single))
    if (!length(cand)) next
    d2 <- colSums((t(xy[cand, , drop = FALSE]) - xy[i, ])^2)
    d_m <- c(d_m, sqrt(min(d2)))
    if (sqrt(min(d2)) <= max_m) cl[i] <- cl0[cand[which.min(d2)]]
  }
  list(cl = as.integer(factor(cl)), d_m = d_m)
}
mg  <- merge_singletons(1000)                      # main definition
mgA <- merge_singletons(Inf)                       # alternative definition (all singletons merged): sensitivity analysis
cl  <- mg$cl
diam <- function(g) sapply(split(which(ok), g[ok]), function(i) if (length(i) > 1) max(dist(xy[i, ])) else 0)
fmt3 <- function(x) paste(min(x), median(x), max(x), sep = " / ")
fmt2 <- function(x) paste(round(median(x)), round(max(x)), sep = " / ")
sz0 <- table(cl0); sz1 <- table(cl); szA <- table(mgA$cl); dm0 <- diam(cl0); dm1 <- diam(cl); dmA <- diam(mgA$cl)
design <- data.frame(section = "clusters", variable = c(
  "cut height (m, complete linkage = max diameter before merging)",
  "BEFORE merging: clusters", "BEFORE merging: singletons", "BEFORE merging: size min / median / max",
  "BEFORE merging: diameter median / max, clusters with >= 2 members (m)",
  "singleton distance to nearest same-district cluster member: median / max (m)",
  "singletons within 1 km (merged)", "singletons beyond 1 km (kept as own cluster)",
  "MAIN definition (merge within 1 km): clusters", "MAIN: singletons remaining", "MAIN: size min / median / max",
  "MAIN: diameter median / max, clusters with >= 2 members (m)",
  "alternative definition (all singletons merged): clusters", "alternative definition: size min / median / max",
  "alternative definition: diameter median / max (m)"),
  value = c(h_star, length(sz0), sum(sz0 == 1), fmt3(sz0), fmt2(dm0[sz0 > 1]), fmt2(mg$d_m),
            sum(mg$d_m <= 1000), sum(mg$d_m > 1000), length(sz1), sum(sz1 == 1), fmt3(sz1), fmt2(dm1[sz1 > 1]),
            length(szA), fmt3(szA), fmt2(dmA)))
saveRDS(data.frame(Respondent_ID = p$Respondent_ID[ok], cluster = cl[ok]), file.path(der_dir, "clusters.rds"))
iA <- match(A$Respondent_ID, p$Respondent_ID)
A$cluster <- factor(cl[iA])
A$cluster_all <- factor(mgA$cl[iA])

## 4. Leave-self-out collective measures -----------------------------------------------------------
W <- as.matrix(dist(xy[iA, ])) <= near_m; diag(W) <- FALSE
A$nb500 <- rowSums(W)
lso_500 <- function(x) ifelse(A$nb500 >= 3, as.vector(W %*% x) / A$nb500, NA)
n_cl <- ave(rep(1, nrow(A)), A$cluster, FUN = sum) - 1
lso_cl <- function(x) ifelse(n_cl >= 3, (ave(x, A$cluster, FUN = sum) - x) / n_cl, NA)
A$prob_b <- lso_500(A$prob_z); A$prob_w <- A$prob_z - A$prob_b
A$coh_b  <- lso_500(A$coh_z);  A$coh_w  <- A$coh_z - A$coh_b
A$prob_bc <- lso_cl(A$prob_z); A$prob_wc <- A$prob_z - A$prob_bc
A$coh_bc  <- lso_cl(A$coh_z);  A$coh_wc  <- A$coh_z - A$coh_bc
rm(W)
q <- quantile(A$nb500, c(.25, .5, .75)); qc <- quantile(n_cl + 1, c(.25, .5, .75))
design <- rbind(design, data.frame(section = "neighbours", variable = c(
  "other respondents within 500 m: median (IQR)", "respondents with < 3 others within 500 m (between component missing)",
  "respondents per reconstructed cluster (analytic sample): median (IQR)", "respondents with < 3 others in cluster (cluster sensitivity missing)"),
  value = c(sprintf("%g (%g-%g)", q[2], q[1], q[3]), sum(A$nb500 < 3), sprintf("%g (%g-%g)", qc[2], qc[1], qc[3]), sum(n_cl < 3))))

## 5. Reliability and factor structure -------------------------------------------------------------
rel <- NULL
add_rel <- function(scale, items, data) {
  X <- data[complete.cases(data[items]), items]
  a <- suppressWarnings(suppressMessages(psych::alpha(X, check.keys = FALSE, warnings = FALSE)))
  om <- if (length(items) >= 3) tryCatch(suppressWarnings(suppressMessages(
    psych::omega(X, nfactors = 1, plot = FALSE)))$omega.tot, error = function(e) NA) else NA
  rel <<- rbind(rel,
    data.frame(scale = scale, item = "(scale)", statistic = c("n", "n_items", "alpha", "std_alpha", "omega_total", "mean_inter_item_r"),
               value = c(nrow(X), length(items), a$total$raw_alpha, a$total$std.alpha, om, a$total$average_r)),
    data.frame(scale = scale, item = rep(items, 2), statistic = rep(c("item_total_r_corrected", "alpha_if_dropped"), each = length(items)),
               value = c(a$item.stats[items, "r.drop"], a$alpha.drop[items, "raw_alpha"])))
}
add_rel("problems (14 items)", names(pn_items), A)
add_rel("cohesion (TS1, TS3, CAC1, CAC2)", coh_items, A)
add_rel("distrust (TS2, TS4)", dis_items, A)
add_rel("all 6 trust/cooperation items, TS2/TS4 not reversed", names(sc_items), A)
for (s in names(subscales)) add_rel(paste0("problems subscale: ", s), subscales[[s]], A)
add_rel("PSS-10", pss_items, A)
efa <- psych::fa(A[names(pn_items)], nfactors = 1, fm = "ml")
ev  <- eigen(cor(A[names(pn_items)]))$values
rel <- rbind(rel,
  data.frame(scale = "problems EFA (1 factor, ML)", item = names(pn_items), statistic = "loading", value = as.numeric(efa$loadings[, 1])),
  data.frame(scale = "problems EFA (1 factor, ML)", item = "(scale)", statistic = c("variance_explained", "eigen_1", "eigen_2", "eigen_3", "TLI", "RMSEA"),
             value = c(efa$Vaccounted["Proportion Var", 1], ev[1:3], efa$TLI, efa$RMSEA[1])))
cfa <- lavaan::cfa(paste("prob =~", paste(names(pn_items), collapse = " + ")), data = A, estimator = "MLR")
fm  <- lavaan::fitMeasures(cfa, c("chisq.scaled", "df", "cfi.robust", "tli.robust", "rmsea.robust", "srmr"))
std <- lavaan::standardizedSolution(cfa); std <- std[std$op == "=~", ]
rel <- rbind(rel,
  data.frame(scale = "problems CFA (1 factor, MLR)", item = "(scale)", statistic = names(fm), value = as.numeric(fm)),
  data.frame(scale = "problems CFA (1 factor, MLR)", item = std$rhs, statistic = "std_loading", value = std$est.std))
write.csv(rel, file.path(tab_dir, "reliability.csv"), row.names = FALSE)

## 6. Main models --------------------------------------------------------------------------------
cov_rhs <- "age + female + educ_f + wealth_quintile + Marital + Emp + YrsInHome + ethn5 + DISTRICT"
fit_lmm <- function(rhs, y = "pss_z", data = A, re = "cluster")
  lmer(as.formula(paste(y, "~", rhs, "+", cov_rhs, "+ (1 |", re, ")")), data = data, REML = TRUE)
## poor / hiedu main effects are absorbed by wealth_quintile / educ_f, so only their interactions enter
rhs <- c(A = "prob_z", B = "coh_z", C = "prob_z + coh_z", D = "prob_z * coh_z",
         E = "prob_b + prob_w + coh_b + coh_w",
         F = "prob_z * coh_z + prob_z:poor + coh_z:poor + prob_z:coh_z:poor",
         G = "prob_z * coh_z + prob_z:hiedu + coh_z:hiedu + prob_z:coh_z:hiedu")
tidy_lmm <- function(m, model) {
  b <- fixef(m); se <- sqrt(diag(as.matrix(vcov(m)))); r2 <- performance::r2_nakagawa(m)
  vc <- as.data.frame(lme4::VarCorr(m))
  data.frame(model = model, term = names(b), estimate = b, se = se, ci_low = b - qnorm(.975) * se,
             ci_high = b + qnorm(.975) * se, p = 2 * pnorm(-abs(b / se)), n = nobs(m), n_clusters = ngrps(m)[[1]],
             r2_marginal = unname(r2$R2_marginal), r2_conditional = unname(r2$R2_conditional),
             icc = vc$vcov[1] / sum(vc$vcov), row.names = NULL)
}
fits <- lapply(rhs, fit_lmm)
m0  <- lmer(pss_z ~ 1 + (1 | cluster), data = A)
m0d <- lmer(pss_z ~ DISTRICT + (1 | cluster), data = A)
icc0 <- sapply(list(m0, m0d), function(m) { vc <- as.data.frame(lme4::VarCorr(m)); vc$vcov[1] / sum(vc$vcov) })
mods <- do.call(rbind, c(
  list(data.frame(model = c("Null", "Null + district FE"), term = "(Intercept)", estimate = NA, se = NA, ci_low = NA,
                  ci_high = NA, p = NA, n = nobs(m0), n_clusters = ngrps(m0)[["cluster"]], r2_marginal = NA,
                  r2_conditional = NA, icc = icc0)),
  Map(tidy_lmm, fits, names(fits))))

## variance components (REML) for the null models and models A-G
vcomp <- do.call(rbind, Map(function(m, model) {
  vc <- as.data.frame(lme4::VarCorr(m))
  data.frame(model = model, cluster_variance = vc$vcov[vc$grp == "cluster"], residual_variance = vc$vcov[vc$grp == "Residual"],
             icc = vc$vcov[vc$grp == "cluster"] / sum(vc$vcov), n_clusters = ngrps(m)[["cluster"]], n = nobs(m))
}, c(list(m0, m0d), fits), c("Null", "Null + district FE", names(fits))))
write.csv(vcomp, file.path(tab_dir, "variance_components.csv"), row.names = FALSE)

## simple slopes of problems at cohesion -1, 0, +1 SD (model D)
mD <- fits$D; bD <- fixef(mD); VD <- as.matrix(vcov(mD))
slopes <- do.call(rbind, lapply(c(-1, 0, 1), function(cz) {
  g <- setNames(numeric(length(bD)), names(bD)); g["prob_z"] <- 1; g["prob_z:coh_z"] <- cz
  est <- sum(g * bD); se <- sqrt(drop(t(g) %*% VD %*% g))
  data.frame(model = "D simple slopes", term = sprintf("prob_z at coh_z = %+d", cz), estimate = est, se = se,
             ci_low = est - qnorm(.975) * se, ci_high = est + qnorm(.975) * se, p = 2 * pnorm(-abs(est / se)),
             n = nobs(mD), n_clusters = ngrps(mD)[["cluster"]], r2_marginal = NA, r2_conditional = NA, icc = NA)
}))
mods <- rbind(mods, slopes)
write.csv(mods, file.path(tab_dir, "models.csv"), row.names = FALSE)

## 7. Predicted PSS across problems at cohesion -1, 0, +1 SD (model D, model-based CIs) -------------
## Average of fixed-effect predictions over the analytic sample (linear, so = mean covariate row).
X <- getME(mD, "X"); xb <- colMeans(X)
pred <- do.call(rbind, lapply(c(-1, 0, 1), function(cz) do.call(rbind, lapply(seq(-2, 2, by = 0.1), function(pz) {
  x <- xb; x["prob_z"] <- pz; x["coh_z"] <- cz; x["prob_z:coh_z"] <- pz * cz
  data.frame(prob_z = pz, coh_z = cz, fit = sum(x * bD), se = sqrt(drop(t(x) %*% VD %*% x)))
}))))
mu <- mean(A$PSS_score); sdp <- sd(A$PSS_score)
pred <- transform(pred, pss = mu + sdp * fit, lo = mu + sdp * (fit - qnorm(.975) * se), hi = mu + sdp * (fit + qnorm(.975) * se),
                  cohesion = factor(coh_z, levels = c(-1, 0, 1), labels = c("Low cohesion (-1 SD)", "Mean cohesion", "High cohesion (+1 SD)")))
write.csv(pred, file.path(tab_dir, "predicted_pss_modelbased.csv"), row.names = FALSE)
seq_blue <- c("#86b6ef", "#2a78d6", "#104281")      # sequential blue steps (ordered levels)
ends <- pred[pred$prob_z == 2, ]
g <- ggplot(pred, aes(prob_z, pss, colour = cohesion, fill = cohesion, linetype = cohesion)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_text(data = ends, aes(label = cohesion), hjust = 0, nudge_x = 0.06, size = 3.2, colour = "#52514e", show.legend = FALSE) +
  scale_colour_manual(values = seq_blue) + scale_fill_manual(values = seq_blue) +
  scale_linetype_manual(values = c("dashed", "solid", "dotdash")) +
  scale_x_continuous(breaks = -2:2, expand = expansion(mult = c(0.02, 0.38))) +
  labs(x = "Perceived neighbourhood problems (SD from mean)", y = "Predicted PSS-10 score (95% CI)",
       colour = NULL, fill = NULL, linetype = NULL,
       caption = "Model D: linear mixed model, reconstructed-cluster random intercept, district fixed effects, model-based CIs, covariates averaged over the sample.") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
        axis.text = element_text(colour = "#52514e"), plot.caption = element_text(colour = "#52514e", size = 8, hjust = 0),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))
ggsave(file.path(fig_dir, "predicted_pss_modelbased.png"), g, width = 7, height = 4.8, dpi = 300)

## 8. Sensitivity analyses (model-based inference) --------------------------------------------------
key <- "prob|coh|dis_z|crime|envir|socio"
pick <- function(d, analysis) { d <- d[grepl(key, d$term), ]; cbind(analysis = analysis, d) }
sens_lmm <- function(m, analysis, model) {
  d <- tidy_lmm(m, model); pick(transform(d[, c("model", "term", "estimate", "se", "ci_low", "ci_high", "p", "n", "n_clusters")],
                                         scale = "SD of outcome"), analysis)
}
sens <- list()

## ordinal PSS: cumulative logit with SEs clustered by reconstructed cluster (a cumulative link mixed
## model with 34 thresholds does not converge in reasonable time)
A$pss_ord <- factor(A$PSS_score, ordered = TRUE)
for (k in c("C", "D")) {
  mo <- ordinal::clm(as.formula(paste("pss_ord ~", rhs[[k]], "+", cov_rhs)), data = A, link = "logit")
  b  <- coef(mo); se <- sqrt(diag(sandwich::vcovCL(mo, cluster = A$cluster, type = "HC0")))
  b  <- b[!grepl("\\|", names(b))]; se <- se[names(b)]
  sens[[paste0("a", k)]] <- pick(data.frame(model = k, term = names(b), estimate = b, se = se,
    ci_low = b - qnorm(.975) * se, ci_high = b + qnorm(.975) * se, p = 2 * pnorm(-abs(b / se)),
    n = nrow(A), n_clusters = nlevels(A$cluster), scale = "log cumulative odds", row.names = NULL),
    "Ordinal PSS-10 (clm, cluster-robust SEs)")
}
## binary PSS > 13. Four districts have 100% PSS > 13, so their fixed effects are quasi-separated
## (not interpretable); the bobyqa optimizer stalls on them, nloptwrap converges.
for (k in c("C", "D")) {
  mb <- glmer(as.formula(paste("pss_hi ~", rhs[[k]], "+", cov_rhs, "+ (1 | cluster)")), data = A, family = binomial,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
  b <- fixef(mb); se <- sqrt(diag(as.matrix(vcov(mb))))
  sens[[paste0("b", k)]] <- pick(data.frame(model = k, term = names(b), estimate = b, se = se, ci_low = b - qnorm(.975) * se,
    ci_high = b + qnorm(.975) * se, p = 2 * pnorm(-abs(b / se)), n = nobs(mb), n_clusters = ngrps(mb)[["cluster"]],
    scale = "log odds", row.names = NULL), "Binary PSS-10 > 13 (glmer)")
}
## problem subscales
sens$c1 <- sens_lmm(fit_lmm("crime_z + envir_z + coh_z"), "Problem subscales", "C: crime + environmental")
sens$c2 <- sens_lmm(fit_lmm("crime_z + envir_z + socio_z + coh_z"), "Problem subscales", "C: crime + environmental + socioeconomic")
for (s in names(subscales)) {
  v <- paste0(s, "_z")
  sens[[paste0("c_", s, "C")]] <- sens_lmm(fit_lmm(paste(v, "+ coh_z")), "Problem subscales", paste0("C: ", s, " only"))
  sens[[paste0("c_", s, "D")]] <- sens_lmm(fit_lmm(paste(v, "* coh_z")), "Problem subscales", paste0("D: ", s, " only"))
}
## distrust added
sens$dC <- sens_lmm(fit_lmm(paste(rhs[["C"]], "+ dis_z")), "Distrust added", "C + distrust")
sens$dD <- sens_lmm(fit_lmm(paste(rhs[["D"]], "+ dis_z")), "Distrust added", "D + distrust")
## clustering: OLS + district FE with cluster-robust SEs (no cluster random effect)
ols <- lm(as.formula(paste("pss_z ~", rhs[["D"]], "+", cov_rhs)), data = A)
for (cv in c("DISTRICT", "cluster")) {
  ct <- lmtest::coeftest(ols, vcov = sandwich::vcovCL(ols, cluster = A[[cv]], type = "HC1"))
  sens[[paste0("g", cv)]] <- pick(data.frame(model = paste0("D, OLS, SEs clustered by ", cv), term = rownames(ct), estimate = ct[, 1],
    se = ct[, 2], ci_low = ct[, 1] - qnorm(.975) * ct[, 2], ci_high = ct[, 1] + qnorm(.975) * ct[, 2], p = ct[, 4],
    n = nobs(ols), n_clusters = length(unique(A[[cv]])), scale = "SD of outcome", row.names = NULL), "OLS with cluster-robust SEs")
}
## model E with leave-self-out means within the reconstructed cluster
sens$h <- sens_lmm(fit_lmm("prob_bc + prob_wc + coh_bc + coh_wc"), "Model E with within-cluster leave-self-out means", "E (cluster)")

## Cluster-definition comparison (models C, D, E): main clusters vs alternative (all singletons merged)
## vs OLS with district-clustered SEs. Flag: shift > 0.5 SE of the main estimate, or 95% CI crosses 0 differently.
cmp <- do.call(rbind, lapply(c("C", "D", "E"), function(k) {
  fin <- tidy_lmm(fits[[k]], k); old <- tidy_lmm(fit_lmm(rhs[[k]], re = "cluster_all"), k)
  o   <- lm(as.formula(paste("pss_z ~", rhs[[k]], "+", cov_rhs)), data = A)
  ct  <- lmtest::coeftest(o, vcov = sandwich::vcovCL(o, cluster = A$DISTRICT, type = "HC1"))
  terms <- fin$term[grepl(key, fin$term)]
  do.call(rbind, lapply(terms, function(tm) {
    f <- fin[fin$term == tm, ]; a <- old[old$term == tm, ]
    ci <- function(e, s) c(e - qnorm(.975) * s, e + qnorm(.975) * s)
    oc <- ci(ct[tm, 1], ct[tm, 2]); ac <- c(a$ci_low, a$ci_high); fc <- c(f$ci_low, f$ci_high)
    excl0 <- function(x) x[1] > 0 | x[2] < 0
    data.frame(model = k, term = tm, n = f$n,
      main_est = f$estimate, main_ci = sprintf("%.3f to %.3f", fc[1], fc[2]), main_n_clusters = f$n_clusters,
      allmerged_est = a$estimate, allmerged_ci = sprintf("%.3f to %.3f", ac[1], ac[2]), allmerged_n_clusters = a$n_clusters,
      ols_district_cl_est = ct[tm, 1], ols_district_cl_ci = sprintf("%.3f to %.3f", oc[1], oc[2]),
      flag = paste(c(if (abs(a$estimate - f$estimate) > 0.5 * f$se || excl0(ac) != excl0(fc)) "all-merged differs",
                     if (abs(ct[tm, 1] - f$estimate) > 0.5 * f$se || excl0(oc) != excl0(fc)) "OLS/district-SE differs"),
                   collapse = "; "))
  }))
}))
write.csv(cmp, file.path(tab_dir, "cluster_comparison.csv"), row.names = FALSE)

## Moran's I of model D residuals at cluster level
cent <- apply(xy[iA, ], 2, function(v) tapply(v, A$cluster, mean))
blup <- ranef(mD)$cluster[levels(A$cluster), 1]
mres <- tapply(A$pss_z - as.vector(X %*% bD), A$cluster, mean)[levels(A$cluster)]
lw_k <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(cent, k = 6)), style = "W")
lw_d <- spdep::nb2listw(spdep::dnearneigh(cent, 0, 2000), style = "W", zero.policy = TRUE)
moran <- do.call(rbind, lapply(list(list("cluster random intercept (BLUP)", blup), list("cluster mean marginal residual", mres)), function(r)
  do.call(rbind, lapply(list(list("6 nearest clusters", lw_k), list("clusters within 2 km", lw_d)), function(w) {
    mt <- spdep::moran.test(r[[2]], w[[2]], zero.policy = TRUE)
    mc <- spdep::moran.mc(r[[2]], w[[2]], nsim = 999, zero.policy = TRUE)
    data.frame(analysis = "Moran's I, model D", model = "D", term = paste0(r[[1]], "; ", w[[1]]),
               estimate = unname(mt$estimate[1]), se = sqrt(unname(mt$estimate[3])), ci_low = NA, ci_high = NA,
               p = mc$p.value, n = nrow(A), n_clusters = length(r[[2]]),
               scale = sprintf("Moran's I (E[I] = %.4f; p from 999 permutations)", mt$estimate[2]))
  }))))
sens$f <- moran
sensitivity <- do.call(rbind, lapply(sens, function(d) d[, c("analysis", "model", "term", "estimate", "se", "ci_low", "ci_high", "p", "n", "n_clusters", "scale")]))
write.csv(sensitivity, file.path(tab_dir, "sensitivity_modelbased.csv"), row.names = FALSE)

## 9. Descriptives (analytic sample) --------------------------------------------------------------
cont <- c(PSS_score = "PSS-10 score", problems = "Neighbourhood problems (mean, 0-4)", cohesion = "Cohesion (mean, 1-5)",
          distrust = "Distrust (mean, 1-5)", crime = "Problems: crime/disorder", envir = "Problems: environmental",
          socio = "Problems: socioeconomic", age = "Age, years", YrsInHome = "Years in home",
          nb500 = "Other respondents within 500 m")
num_row <- function(v) { x <- A[[v]]; x <- x[!is.na(x)]
  data.frame(section = "continuous", variable = cont[[v]], level = NA, n = length(x), pct = NA, mean = mean(x), sd = sd(x),
             median = median(x), q1 = quantile(x, .25), q3 = quantile(x, .75), min = min(x), max = max(x), row.names = NULL) }
educ_labs <- c("0" = "No formal education", "1" = "Some primary school", "2" = "Completed primary school",
               "3" = "Some secondary school", "4" = "Completed secondary school", "5" = "Post-secondary")
cat_rows <- function(x, name) { t <- table(x, useNA = "ifany")
  data.frame(section = "categorical", variable = name, level = names(t), n = as.integer(t), pct = 100 * as.numeric(t) / length(x),
             mean = NA, sd = NA, median = NA, q1 = NA, q3 = NA, min = NA, max = NA, row.names = NULL) }
cors <- cor(A[c("PSS_score", "problems", "cohesion", "distrust", "crime", "envir", "socio")])
pr <- combn(colnames(cors), 2)
cor_rows <- data.frame(section = "correlation (Pearson)", variable = pr[1, ], level = pr[2, ], n = nrow(A), pct = NA,
                       mean = cors[t(pr)], sd = NA, median = NA, q1 = NA, q3 = NA, min = NA, max = NA)
blank <- function(d) data.frame(section = d$section, variable = d$variable, level = if (is.null(d$value)) NA else d$value,
                                n = if (is.null(d$n)) NA else d$n, pct = NA, mean = NA, sd = NA, median = NA, q1 = NA, q3 = NA, min = NA, max = NA)
desc <- rbind(
  data.frame(section = "sample", variable = c("Survey respondents", "Analytic sample"), level = NA, n = c(nrow(p), nrow(A)),
             pct = NA, mean = NA, sd = NA, median = NA, q1 = NA, q3 = NA, min = NA, max = NA),
  blank(excl), blank(design),
  do.call(rbind, lapply(names(cont), num_row)),
  cat_rows(ifelse(A$female == 1, "Female", "Male/other"), "Sex"),
  cat_rows(educ_labs[as.character(A$educ_f)], "Highest education"),
  cat_rows(A$wealth_quintile, "Wealth quintile"), cat_rows(A$Marital, "Marital status"), cat_rows(A$Emp, "Employment"),
  cat_rows(A$ethn5, "Ethnicity (Grusi, Gurma, Mande in Other)"), cat_rows(A$DISTRICT, "District"),
  cat_rows(ifelse(A$pss_hi == 1, "PSS > 13", "PSS <= 13"), "Binary PSS-10 (> 13)"),
  cat_rows(ifelse(A$poor == 1, "Poor (Q1-2)", "Non-poor (Q3-5)"), "Wealth group (model F)"),
  cat_rows(ifelse(A$hiedu == 1, "Completed secondary or higher", "Lower"), "Education group (model G)"),
  cor_rows)
write.csv(desc, file.path(tab_dir, "descriptives.csv"), row.names = FALSE)

## 10. Short summary ------------------------------------------------------------------------------
f2 <- function(d) sprintf("%-26s %6.3f (%6.3f, %6.3f) p=%.3f", d$term, d$estimate, d$ci_low, d$ci_high, d$p)
cat("\n==== Sample:", nrow(A), "of", nrow(p), "respondents;", nlevels(A$cluster), "reconstructed clusters\n")
print(design[, c("variable", "value")], row.names = FALSE, right = FALSE)
cat(sprintf("\nICC (null): %.3f; with district FE: %.3f\n", icc0[1], icc0[2]))
r <- rel[rel$statistic %in% c("alpha", "omega_total") & rel$item == "(scale)", ]
cat("\n==== Reliability\n"); print(reshape(r[, c("scale", "statistic", "value")], idvar = "scale", timevar = "statistic", direction = "wide"), row.names = FALSE)
cat(sprintf("CFA problems: CFI %.3f, TLI %.3f, RMSEA %.3f, SRMR %.3f\n", fm["cfi.robust"], fm["tli.robust"], fm["rmsea.robust"], fm["srmr"]))
cat("\nTrust-item correlation with other cohesion items (negative = negatively worded):\n")
print(mapping[mapping$variable %in% names(sc_items), c("variable", "r_with_cohesion_rest")], row.names = FALSE)
cat("\n==== Main models (outcome: PSS-10, SD units; model-based Wald 95% CI)\n")
for (k in names(fits)) {
  d <- mods[mods$model == k & grepl(key, mods$term), ]
  cat(sprintf("Model %s (n = %d; R2m %.3f, R2c %.3f)\n", k, d$n[1], d$r2_marginal[1], d$r2_conditional[1])); cat(paste(" ", f2(d)), sep = "\n")
}
cat("Simple slopes of problems (model D):\n"); cat(paste(" ", f2(slopes)), sep = "\n")
cat("\n==== Sensitivity (key terms)\n")
for (a in unique(sensitivity$analysis)) {
  d <- sensitivity[sensitivity$analysis == a, ]; cat(a, "\n")
  for (mm in unique(d$model)) { dd <- d[d$model == mm, ]; cat("  ", mm, "(n =", dd$n[1], ")\n"); cat(paste("    ", f2(dd)), sep = "\n") }
}
cat("\n==== Cluster-definition comparison (PSS-10, SD units)\n")
print(cmp[, c("model", "term", "main_est", "main_ci", "allmerged_est", "allmerged_ci", "ols_district_cl_est", "ols_district_cl_ci", "flag")],
      row.names = FALSE, digits = 3)
if (any(moran$estimate > 0.1 & moran$p < 0.01))
  cat("\n*** Moran's I of model D cluster residuals is clearly non-zero: check before fitting a spatial model. ***\n")
