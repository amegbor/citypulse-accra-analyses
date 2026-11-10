## 03_robust_inference.R -------------------------------------------------------------------------
## CR2 cluster-robust inference for every reported term. Same data, models and point estimates as
## 01_main_analysis.R; only the variance estimator changes.
##   Main inference    : CR2 (Bell-McCaffrey) standard errors clustered by reconstructed cluster,
##                       with Satterthwaite degrees of freedom (clubSandwich).
##   Sensitivity bound : CR2 standard errors clustered by district (22 clusters).
## Why: residual diagnostics (02_diagnostics.R) show non-constant residual variance, so model-based (REML,
## Wald) standard errors understate uncertainty.
## Data preparation and fits are taken verbatim from 01_main_analysis.R (everything up to the main fits),
## with its output folders redirected to a temporary directory. Run 01_main_analysis.R first.
## Outputs: outputs/tables/{robust_estimates, modelC_robust, predicted_pss, moran}.csv (aggregates only).

source("config.R")
res_dir <- file.path(out_dir, "tables"); dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
for (pk in c("clubSandwich", "sandwich", "ordinal")) if (!requireNamespace(pk, quietly = TRUE)) stop("Install first: ", pk)

## 1. Data preparation and model fits, verbatim from 01_main_analysis.R --------------------------------
src <- readLines("01_main_analysis.R")
stop_at <- grep("^fits <- lapply\\(rhs, fit_lmm\\)", src); stopifnot(length(stop_at) == 1)
prep <- function() {
  s <- src[1:stop_at]
  tmp <- gsub("\\\\", "/", tempfile("scs_")); dir.create(tmp)
  dl <- grep("^(tab_dir|fig_dir|der_dir) +<- ", s); stopifnot(length(dl) == 3)
  s[dl] <- sprintf('%s <- "%s"', sub(" .*", "", s[dl]), tmp)
  e <- new.env()
  suppressWarnings(suppressMessages(eval(parse(text = s), envir = e)))
  unlink(tmp, recursive = TRUE)
  e
}
E0 <- prep()
A <- E0$A; fits <- E0$fits; rhs <- E0$rhs; cov_rhs <- E0$cov_rhs; fit_lmm <- E0$fit_lmm
library(lme4)

pub <- read.csv(file.path(out_dir, "tables", "models.csv"))
for (k in names(fits)) stopifnot(isTRUE(all.equal(unname(fixef(fits[[k]])[pub$term[pub$model == k]]),
                                                  pub$estimate[pub$model == k], tolerance = 1e-6)))
cat("Fixed effects of models A-G reproduce models.csv\n")

## 2. Robust inference helpers ----------------------------------------------------------------------
key <- "prob|coh|crime|envir|socio|dis_z"
fe  <- function(m) if (inherits(m, "merMod")) fixef(m) else coef(m)
cl_lab <- c(cluster = "CR2, reconstructed cluster", cluster_all = "CR2, alternative cluster (all singletons merged)",
            DISTRICT = "CR2, district")
## CR2 + Satterthwaite; falls back to CR1S (df = G - 1) if CR2 cannot be computed (reported in se_type)
robust_v <- function(m, cl) {
  V <- tryCatch(clubSandwich::vcovCR(m, cluster = cl, type = "CR2"), error = function(e) NULL)
  if (!is.null(V) && all(is.finite(diag(as.matrix(V))))) return(list(V = V, type = "CR2"))
  list(V = clubSandwich::vcovCR(m, cluster = cl, type = "CR1S"), type = "CR1S")
}
cl_for <- function(m, cv, data) { cl <- data[[cv]][match(rownames(model.frame(m)), rownames(data))]; stopifnot(!anyNA(cl)); cl }
rob <- function(m, analysis, model, terms = NULL, data = A, clv = c("cluster", "DISTRICT"), mb = TRUE, scale = "SD of outcome") {
  b <- fe(m); if (is.null(terms)) terms <- names(b)[grepl(key, names(b))]
  n <- nobs(m); out <- list()
  if (mb) {
    se <- sqrt(diag(as.matrix(vcov(m))))[terms]
    out$mb <- data.frame(se_type = "model-based (Wald)", estimate = b[terms], se = se, df = NA, ci_low = b[terms] - qnorm(.975) * se,
                         ci_high = b[terms] + qnorm(.975) * se, p = 2 * pnorm(-abs(b[terms] / se)),
                         n_clusters = if (inherits(m, "merMod")) ngrps(m)[[1]] else NA)
  }
  for (cv in clv) {
    cl <- cl_for(m, cv, data); rv <- robust_v(m, cl)
    ct <- clubSandwich::coef_test(m, vcov = rv$V, test = if (rv$type == "CR2") "Satterthwaite" else "naive-tp", coefs = terms)
    df <- if (rv$type == "CR2") ct$df_Satt else rep(length(unique(cl)) - 1, nrow(ct))
    pv <- if (rv$type == "CR2") ct$p_Satt else 2 * pt(-abs(ct$tstat), df)
    lab <- if (rv$type == "CR2") cl_lab[[cv]] else sub("CR2", "CR1S (CR2 not computable)", cl_lab[[cv]])
    out[[cv]] <- data.frame(se_type = lab, estimate = ct$beta, se = ct$SE, df = df, ci_low = ct$beta - qt(.975, df) * ct$SE,
                            ci_high = ct$beta + qt(.975, df) * ct$SE, p = pv, n_clusters = length(unique(cl)))
  }
  r <- do.call(rbind, out)
  cbind(analysis = analysis, model = model, term = rep(terms, length(out)), r[, c("se_type", "estimate", "se", "df", "ci_low", "ci_high", "p")],
        n = n, n_clusters = r$n_clusters, scale = scale, row.names = NULL)
}
## linear combinations (simple slopes, group-specific slopes), model-based and CR2
lincom <- function(m, L, analysis, model, data = A, clv = c("cluster", "DISTRICT")) {
  b <- fe(m); L <- L[, names(b), drop = FALSE]
  Vm <- as.matrix(vcov(m)); est <- drop(L %*% b); se <- sqrt(diag(L %*% Vm %*% t(L)))
  out <- list(data.frame(term = rownames(L), se_type = "model-based (Wald)", estimate = est, se = se, df = NA,
                         ci_low = est - qnorm(.975) * se, ci_high = est + qnorm(.975) * se, p = 2 * pnorm(-abs(est / se)),
                         n_clusters = ngrps(m)[[1]]))
  for (cv in clv) {
    cl <- cl_for(m, cv, data); rv <- robust_v(m, cl); stopifnot(rv$type == "CR2")
    lc <- clubSandwich::linear_contrast(m, vcov = rv$V, contrasts = L, test = "Satterthwaite")
    out[[cv]] <- data.frame(term = lc$Coef, se_type = cl_lab[[cv]], estimate = lc$Est, se = lc$SE, df = lc$df,
                            ci_low = lc$CI_L, ci_high = lc$CI_U, p = 2 * pt(-abs(lc$Est / lc$SE), lc$df), n_clusters = length(unique(cl)))
  }
  r <- do.call(rbind, out)
  cbind(analysis = analysis, model = model, r[, c("term", "se_type", "estimate", "se", "df", "ci_low", "ci_high", "p")],
        n = nobs(m), n_clusters = r$n_clusters, scale = "SD of outcome", row.names = NULL)
}
Lrow <- function(m, ...) { b <- fe(m); x <- setNames(numeric(length(b)), names(b)); v <- list(...); x[names(v)] <- unlist(v); x }

res <- list()
## 3. Main models A-G ------------------------------------------------------------------------------
for (k in names(fits)) res[[paste0("main_", k)]] <- rob(fits[[k]], "Main models", k,
  terms = names(fixef(fits[[k]]))[grepl(paste0(key, "|poor|hiedu"), names(fixef(fits[[k]])))])
## Model C: every fixed effect (S1 Table), CR2 by cluster and district
full_C <- rob(fits$C, "Model C, all terms", "C", terms = names(fixef(fits$C)))

mD <- fits$D
res$slopes <- lincom(mD, rbind(`prob_z at coh_z = -1` = Lrow(mD, prob_z = 1, `prob_z:coh_z` = -1),
                               `prob_z at coh_z = 0`  = Lrow(mD, prob_z = 1),
                               `prob_z at coh_z = +1` = Lrow(mD, prob_z = 1, `prob_z:coh_z` = 1)), "Simple slopes", "D")
mF <- fits$F; mG <- fits$G
res$F_groups <- lincom(mF, rbind(
  `prob_z, poorer households (Q1-2)` = Lrow(mF, prob_z = 1, `prob_z:poor` = 1),
  `coh_z, poorer households (Q1-2)`  = Lrow(mF, coh_z = 1, `coh_z:poor` = 1),
  `prob_z:coh_z, poorer households (Q1-2)` = Lrow(mF, `prob_z:coh_z` = 1, `prob_z:coh_z:poor` = 1)), "Group-specific slopes", "F")
res$G_groups <- lincom(mG, rbind(
  `prob_z, completed secondary or higher` = Lrow(mG, prob_z = 1, `prob_z:hiedu` = 1),
  `coh_z, completed secondary or higher`  = Lrow(mG, coh_z = 1, `coh_z:hiedu` = 1),
  `prob_z:coh_z, completed secondary or higher` = Lrow(mG, `prob_z:coh_z` = 1, `prob_z:coh_z:hiedu` = 1)), "Group-specific slopes", "G")
## cohesion (and problems) slope in each wealth quintile: model D + exposure x quintile interactions
mQ <- fit_lmm("prob_z * coh_z + prob_z:wealth_quintile + coh_z:wealth_quintile")
qs <- levels(A$wealth_quintile)
Lq <- do.call(rbind, lapply(qs, function(q) {
  a <- list(coh_z = 1); if (q != qs[1]) a[[paste0("coh_z:wealth_quintile", q)]] <- 1
  b <- list(prob_z = 1); if (q != qs[1]) b[[paste0("prob_z:wealth_quintile", q)]] <- 1
  rbind(do.call(Lrow, c(list(mQ), a)), do.call(Lrow, c(list(mQ), b)))
}))
rownames(Lq) <- c(rbind(paste0("coh_z, wealth quintile: ", qs), paste0("prob_z, wealth quintile: ", qs)))
res$quintile <- lincom(mQ, Lq, "Slopes by wealth quintile", "D + exposure x wealth quintile")

## predicted PSS-10 (Fig 1): sample-averaged fixed-effect predictions, CIs from CR2 (cluster)
X <- getME(mD, "X"); xb <- colMeans(X); grid <- expand.grid(prob_z = seq(-2, 2, by = 0.1), coh_z = c(-1, 0, 1))
Lp <- t(apply(grid, 1, function(g) { x <- xb; x["prob_z"] <- g[1]; x["coh_z"] <- g[2]; x["prob_z:coh_z"] <- g[1] * g[2]; x }))
rownames(Lp) <- sprintf("p%+.1f_c%+d", grid$prob_z, grid$coh_z)
pr <- lincom(mD, Lp, "Predicted", "D", clv = "cluster")
mu <- mean(A$PSS_score); sdp <- sd(A$PSS_score)
pred <- cbind(grid[rep(seq_len(nrow(grid)), 2), ], pr[, c("se_type", "estimate", "se", "df", "ci_low", "ci_high")])
pred <- transform(pred, pss = mu + sdp * estimate, lo = mu + sdp * ci_low, hi = mu + sdp * ci_high)
write.csv(pred, file.path(res_dir, "predicted_pss.csv"), row.names = FALSE)

## 4. Sensitivity analyses with CR2 ---------------------------------------------------------------
## ordinal PSS-10: clm is not supported by clubSandwich; cluster-robust sandwich (HC0) by cluster and district
A$pss_ord <- factor(A$PSS_score, ordered = TRUE)
for (k in c("C", "D")) {
  mo <- ordinal::clm(as.formula(paste("pss_ord ~", rhs[[k]], "+", cov_rhs)), data = A, link = "logit")
  b <- coef(mo); tm <- names(b)[!grepl("\\|", names(b)) & grepl(key, names(b))]
  res[[paste0("a", k)]] <- do.call(rbind, lapply(c("cluster", "DISTRICT"), function(cv) {
    se <- sqrt(diag(sandwich::vcovCL(mo, cluster = A[[cv]], type = "HC0")))[tm]
    data.frame(analysis = "Ordinal PSS-10 (clm)", model = k, term = tm,
               se_type = paste0("sandwich HC0, ", if (cv == "cluster") "reconstructed cluster" else "district", " (clm: CR2 not available)"),
               estimate = b[tm], se = se, df = NA, ci_low = b[tm] - qnorm(.975) * se, ci_high = b[tm] + qnorm(.975) * se,
               p = 2 * pnorm(-abs(b[tm] / se)), n = nrow(A), n_clusters = length(unique(A[[cv]])), scale = "log cumulative odds", row.names = NULL)
  }))
}
## PSS-10 > 13: glmer (model-based) and logistic GLM with district FE + CR2 (glmer not supported by clubSandwich)
for (k in c("C", "D")) {
  mb <- glmer(as.formula(paste("pss_hi ~", rhs[[k]], "+", cov_rhs, "+ (1 | cluster)")), data = A, family = binomial,
              control = glmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
  res[[paste0("b", k, "_glmer")]] <- rob(mb, "PSS-10 > 13", paste(k, "(glmer)"), clv = character(0), scale = "log odds")
  mg <- glm(as.formula(paste("pss_hi ~", rhs[[k]], "+", cov_rhs)), data = A, family = binomial)
  res[[paste0("b", k, "_glm")]] <- rob(mg, "PSS-10 > 13", paste(k, "(logistic GLM, district FE)"), mb = FALSE, scale = "log odds")
}
## problem subscales
res$c1 <- rob(fit_lmm("crime_z + envir_z + coh_z"), "Problem subscales", "C: crime + environmental")
res$c2 <- rob(fit_lmm("crime_z + envir_z + socio_z + coh_z"), "Problem subscales", "C: crime + environmental + socioeconomic")
for (s in c("crime", "envir", "socio")) {
  v <- paste0(s, "_z")
  res[[paste0("c_", s, "C")]] <- rob(fit_lmm(paste(v, "+ coh_z")), "Problem subscales", paste0("C: ", s, " only"))
  res[[paste0("c_", s, "D")]] <- rob(fit_lmm(paste(v, "* coh_z")), "Problem subscales", paste0("D: ", s, " only"))
}
## distrust
res$dC <- rob(fit_lmm(paste(rhs[["C"]], "+ dis_z")), "Distrust added", "C + distrust")
res$dD <- rob(fit_lmm(paste(rhs[["D"]], "+ dis_z")), "Distrust added", "D + distrust")
## cluster definition: alternative clusters (RE and CR2 by the alternative cluster); OLS with CR2
for (k in c("C", "D", "E")) {
  res[[paste0("g_alt", k)]] <- rob(fit_lmm(rhs[[k]], re = "cluster_all"), "Alternative cluster definition", k,
                                   clv = c("cluster_all", "DISTRICT"))
  o <- lm(as.formula(paste("pss_z ~", rhs[[k]], "+", cov_rhs)), data = A)
  res[[paste0("g_ols", k)]] <- rob(o, "OLS, district FE, no random effect", k, mb = FALSE)
}
## model E with within-cluster leave-self-out means
res$h <- rob(fit_lmm("prob_bc + prob_wc + coh_bc + coh_wc"), "Model E, within-cluster leave-self-out means", "E (cluster)")

## 5. Write -----------------------------------------------------------------------------------------
allr <- do.call(rbind, res)
allr <- allr[, c("analysis", "model", "term", "se_type", "estimate", "se", "df", "ci_low", "ci_high", "p", "n", "n_clusters", "scale")]
write.csv(allr, file.path(res_dir, "robust_estimates.csv"), row.names = FALSE)
write.csv(full_C, file.path(res_dir, "modelC_robust.csv"), row.names = FALSE)
moran <- read.csv(file.path(out_dir, "tables", "sensitivity_modelbased.csv"))
write.csv(moran[moran$analysis == "Moran's I, model D", ], file.path(res_dir, "moran.csv"), row.names = FALSE)

options(width = 220)
print(within(allr[allr$analysis %in% c("Main models", "Simple slopes", "Group-specific slopes"), c("model", "term", "se_type", "estimate", "ci_low", "ci_high", "p")],
             { estimate <- round(estimate, 3); ci_low <- round(ci_low, 3); ci_high <- round(ci_high, 3); p <- signif(p, 2) }), row.names = FALSE)
cat("clubSandwich", as.character(packageVersion("clubSandwich")), ";", R.version.string, "\n")
