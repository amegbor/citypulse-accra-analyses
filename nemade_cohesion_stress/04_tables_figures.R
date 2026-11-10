## 04_tables_figures.R ---------------------------------------------------------------------------
## Manuscript tables and figures. Reported inference uses CR2 cluster-robust standard errors clustered by
## reconstructed cluster (Satterthwaite degrees of freedom), with CR2 clustered by district as a bound.
## Reads only aggregate outputs of 01_main_analysis.R and 03_robust_inference.R in outputs/tables.
## Writes outputs/tables/{Table1, Table2, Table3, S1_Table, S2_Table}.docx and
## outputs/figures/{Fig1, Fig2}.tif (PLOS ONE format: 300 dpi, LZW, RGB) with PNG previews.

source("config.R")
pk <- c("ggplot2", "ragg", "flextable", "officer")
miss <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install first: ", paste(miss, collapse = ", "))
library(ggplot2); library(flextable)

tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures"); dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

R     <- read.csv(file.path(tab_dir, "robust_estimates.csv"), check.names = FALSE)
fC    <- read.csv(file.path(tab_dir, "modelC_robust.csv"))
pred  <- read.csv(file.path(tab_dir, "predicted_pss.csv"))
moran <- read.csv(file.path(tab_dir, "moran.csv"))
mods  <- read.csv(file.path(tab_dir, "models.csv"))
vcmp  <- read.csv(file.path(tab_dir, "variance_components.csv"))
desc  <- read.csv(file.path(tab_dir, "descriptives.csv"))

ty <- function(x) ifelse(grepl("^model-based", x), "mb", ifelse(grepl("district", x), "dist", "cl"))
R$ty <- ty(R$se_type); fC$ty <- ty(fC$se_type)
get1 <- function(analysis, model, term, type = "cl") {
  r <- R[R$analysis == analysis & R$model == model & R$term == term & R$ty == type, ]
  if (nrow(r) != 1) stop("lookup failed: ", paste(analysis, model, term, type, sep = " | "))
  r
}

## formatting (2 decimals; p to 3 decimals or <0.001)
f2  <- function(x) ifelse(is.na(x), "", sub("^-", "−", sub("^-(0\\.00)$", "\\1", formatC(x, format = "f", digits = 2))))
fp  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
fci <- function(lo, hi) ifelse(is.na(lo), "", paste0(f2(lo), ", ", f2(hi)))
eci <- function(r) paste0(f2(r$estimate), " (", fci(r$ci_low, r$ci_high), ")")
b0  <- officer::fp_border(color = "black", width = 0.75)
style <- function(ft, first_left = TRUE) {
  ft <- font(ft, fontname = "Times New Roman", part = "all"); ft <- fontsize(ft, size = 10, part = "all")
  ft <- border_remove(ft); ft <- hline_top(ft, border = b0, part = "header"); ft <- hline_bottom(ft, border = b0, part = "header")
  ft <- hline_bottom(ft, border = b0, part = "body")
  ft <- align(ft, j = -1, align = "center", part = "all"); if (first_left) ft <- align(ft, j = 1, align = "left", part = "all")
  ft <- padding(ft, padding.top = 1, padding.bottom = 1, part = "all"); autofit(ft)
}
footer <- function(ft, lines) { ft <- add_footer_lines(ft, lines); ft <- font(ft, fontname = "Times New Roman", part = "footer"); fontsize(ft, size = 9, part = "footer") }
write_doc <- function(ft, title, file, landscape = FALSE) {
  doc <- officer::read_docx()
  doc <- officer::body_add_fpar(doc, officer::fpar(officer::ftext(title, officer::fp_text(font.family = "Times New Roman", font.size = 10, bold = TRUE))))
  doc <- flextable::body_add_flextable(doc, ft, align = "left")
  if (landscape) doc <- officer::body_set_default_section(doc, officer::prop_section(page_size = officer::page_size(orient = "landscape")))
  print(doc, target = file)
}
inf_note <- paste("95% CIs and p values use CR2 (Bell-McCaffrey) cluster-robust standard errors clustered by reconstructed",
                  "cluster, with Satterthwaite degrees of freedom; CR2 bounds clustered by district are in S2 Table.")

## Table 1: participant characteristics ----------------------------------------------------------------
d_cont <- function(lab, digits) { r <- desc[desc$section == "continuous" & desc$variable == lab, ]
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), r$mean, r$sd) }
d_cat  <- function(var, level) { r <- desc[desc$section == "categorical" & desc$variable == var & desc$level == level, ]
  stopifnot(nrow(r) == 1); sprintf("%s (%.1f)", format(r$n, big.mark = ","), r$pct) }
d_med  <- function(lab) { r <- desc[desc$section == "continuous" & desc$variable == lab, ]
  sprintf("%g (%g–%g)", r$median, r$q1, r$q3) }
akan  <- grep("^Akan", unique(desc$level[desc$variable == "Ethnicity (Grusi, Gurma, Mande in Other)"]), value = TRUE)
ind <- "   "
t1 <- rbind(
  c("Sociodemographic characteristics", ""),
  c("Age, years", d_cont("Age, years", 1)),
  c("Female", d_cat("Sex", "Female")),
  c("Highest level of education", ""),
  cbind(paste0(ind, c("No formal education", "Some primary school", "Completed primary school", "Some secondary school",
                      "Completed secondary school", "Post-secondary")),
        sapply(c("No formal education", "Some primary school", "Completed primary school", "Some secondary school",
                 "Completed secondary school", "Post-secondary"), function(l) d_cat("Highest education", l))),
  c("Household wealth quintile", ""),
  cbind(paste0(ind, c("Poorest", "Poorer", "Middle", "Richer", "Richest")),
        sapply(c("Poorest", "Poorer", "Middle", "Richer", "Richest"), function(l) d_cat("Wealth quintile", l))),
  c("Marital status", ""),
  cbind(paste0(ind, c("Married or living with a partner", "Single", "Divorced or separated", "Widowed")),
        sapply(c("Married or in a domestic partnership", "Single", "Divorced/Separated", "Widowed"), function(l) d_cat("Marital status", l))),
  c("Employment status", ""),
  cbind(paste0(ind, c("Employed or self-employed", "Unemployed", "Student", "Retired")),
        sapply(c("Employed/Self-employed", "Unemployed", "Student", "Retired"), function(l) d_cat("Employment", l))),
  c("Ethnicity", ""),
  cbind(paste0(ind, c("Akan", "Ga-Adangbe", "Ewe", "Mole-Dagbani", "Other")),
        sapply(c(akan, "Ga-Adangbe", "Ewe", "Mole-Dagbani", "Other"), function(l) d_cat("Ethnicity (Grusi, Gurma, Mande in Other)", l))),
  c("Years at current dwelling", d_cont("Years in home", 1)),
  c("Perceived stress", ""),
  c("PSS-10 score (0–40)", d_cont("PSS-10 score", 1)),
  c("PSS-10 score above 13", d_cat("Binary PSS-10 (> 13)", "PSS > 13")),
  c("Neighborhood measures", ""),
  c("Neighborhood problems (0–4)", d_cont("Neighbourhood problems (mean, 0-4)", 2)),
  c(paste0(ind, "Crime and disorder"), d_cont("Problems: crime/disorder", 2)),
  c(paste0(ind, "Environmental problems"), d_cont("Problems: environmental", 2)),
  c(paste0(ind, "Socioeconomic problems"), d_cont("Problems: socioeconomic", 2)),
  c("Social cohesion (1–5)", d_cont("Cohesion (mean, 1-5)", 2)),
  c("Distrust (1–5)", d_cont("Distrust (mean, 1-5)", 2)),
  c("Other participants within 500 m, median (IQR)", d_med("Other respondents within 500 m")))
t1 <- data.frame(Characteristic = t1[, 1], `n (%) or mean (SD)` = t1[, 2], check.names = FALSE, row.names = NULL)
n_all <- desc$n[desc$section == "sample" & desc$variable == "Analytic sample"]
h1 <- which(t1[[2]] == "" & !startsWith(t1[[1]], ind) & t1[[1]] %in% c("Sociodemographic characteristics", "Perceived stress", "Neighborhood measures"))
ft <- style(flextable(t1)); ft <- bold(ft, i = h1, part = "body")
ft <- footer(ft, "PSS-10, 10-item Perceived Stress Scale; SD, standard deviation; IQR, interquartile range.")
write_doc(ft, sprintf("Table 1. Characteristics of participants (n = %s).", format(n_all, big.mark = ",")), file.path(tab_dir, "Table1.docx"))

## Table 2: models A-D -----------------------------------------------------------------------------
cellv <- function(k, term) { r <- R[R$analysis == "Main models" & R$model == k & R$term == term & R$ty == "cl", ]; if (nrow(r)) eci(r) else "—" }
vv <- function(k, col) f2(vcmp[vcmp$model == k, col]); rr <- function(k, col) f2(mods[mods$model == k, col][1])
t2 <- data.frame(check.names = FALSE,
  ` ` = c("Neighborhood problems (per SD)", "Social cohesion (per SD)", "Problems × cohesion",
          "Cluster variance", "Residual variance", "ICC (cluster)", "Marginal R²", "Conditional R²"),
  `Model A` = c(cellv("A", "prob_z"), "—", "—", vv("A", "cluster_variance"), vv("A", "residual_variance"), vv("A", "icc"), rr("A", "r2_marginal"), rr("A", "r2_conditional")),
  `Model B` = c("—", cellv("B", "coh_z"), "—", vv("B", "cluster_variance"), vv("B", "residual_variance"), vv("B", "icc"), rr("B", "r2_marginal"), rr("B", "r2_conditional")),
  `Model C` = c(cellv("C", "prob_z"), cellv("C", "coh_z"), "—", vv("C", "cluster_variance"), vv("C", "residual_variance"), vv("C", "icc"), rr("C", "r2_marginal"), rr("C", "r2_conditional")),
  `Model D` = c(cellv("D", "prob_z"), cellv("D", "coh_z"), cellv("D", "prob_z:coh_z"), vv("D", "cluster_variance"), vv("D", "residual_variance"), vv("D", "icc"), rr("D", "r2_marginal"), rr("D", "r2_conditional")))
ft <- style(flextable(t2)); ft <- hline(ft, i = 3, border = officer::fp_border(color = "black", width = 0.5), part = "body")
ft <- footer(ft, c(paste("Coefficients (95% CI) are differences in standardized PSS-10 scores per SD difference in each exposure, from linear mixed models",
  "with a random intercept for reconstructed clusters and fixed effects for district, adjusted for age, sex, education, household wealth quintile,",
  "marital status, employment status, years at the current dwelling, and ethnicity (REML)."), inf_note,
  "One SD of the PSS-10 corresponds to 4.2 points. ICC, intraclass correlation coefficient; SD, standard deviation."))
write_doc(ft, "Table 2. Associations of perceived neighborhood problems and social cohesion with perceived stress (n = 3,150).",
          file.path(tab_dir, "Table2.docx"))

## Table 3: models E, F, G ----------------------------------------------------------------------------
row3 <- function(label, analysis, model, term) { r <- get1(analysis, model, term); data.frame(Term = label, `Estimate (95% CI)` = eci(r), p = fp(r$p), check.names = FALSE) }
hdr3 <- function(label) data.frame(Term = label, `Estimate (95% CI)` = "", p = "", check.names = FALSE)
nE <- format(get1("Main models", "E", "prob_b")$n, big.mark = ",")
t3 <- rbind(
  hdr3(paste0("Model E: own perceptions and neighbors’ reports (n = ", nE, ")")),
  row3("Neighborhood problems, neighbors’ reports (between)", "Main models", "E", "prob_b"),
  row3("Neighborhood problems, own perception (within)", "Main models", "E", "prob_w"),
  row3("Social cohesion, neighbors’ reports (between)", "Main models", "E", "coh_b"),
  row3("Social cohesion, own perception (within)", "Main models", "E", "coh_w"),
  hdr3("Model F: household wealth (poorest two quintiles vs rest; n = 3,150)"),
  row3("Neighborhood problems (wealthier households)", "Main models", "F", "prob_z"),
  row3("Social cohesion (wealthier households)", "Main models", "F", "coh_z"),
  row3("Problems × cohesion", "Main models", "F", "prob_z:coh_z"),
  row3("Problems × poorer household", "Main models", "F", "prob_z:poor"),
  row3("Cohesion × poorer household", "Main models", "F", "coh_z:poor"),
  row3("Problems × cohesion × poorer household", "Main models", "F", "prob_z:coh_z:poor"),
  row3("   Neighborhood problems, poorer households (derived)", "Group-specific slopes", "F", "prob_z, poorer households (Q1-2)"),
  row3("   Social cohesion, poorer households (derived)", "Group-specific slopes", "F", "coh_z, poorer households (Q1-2)"),
  hdr3("Model G: education (completed secondary or higher vs lower; n = 3,150)"),
  row3("Neighborhood problems (lower education)", "Main models", "G", "prob_z"),
  row3("Social cohesion (lower education)", "Main models", "G", "coh_z"),
  row3("Problems × cohesion", "Main models", "G", "prob_z:coh_z"),
  row3("Problems × higher education", "Main models", "G", "prob_z:hiedu"),
  row3("Cohesion × higher education", "Main models", "G", "coh_z:hiedu"),
  row3("Problems × cohesion × higher education", "Main models", "G", "prob_z:coh_z:hiedu"),
  row3("   Neighborhood problems, higher education (derived)", "Group-specific slopes", "G", "prob_z, completed secondary or higher"),
  row3("   Social cohesion, higher education (derived)", "Group-specific slopes", "G", "coh_z, completed secondary or higher"))
h3 <- which(t3$`Estimate (95% CI)` == "")
ft <- style(flextable(t3)); ft <- bold(ft, i = h3, part = "body"); ft <- merge_h_range(ft, i = h3, j1 = 1, j2 = 3); ft <- align(ft, i = h3, j = 1, align = "left")
ft <- footer(ft, c(paste("Coefficients are differences in standardized PSS-10 scores per SD difference in each exposure, from linear mixed models adjusted as in Table 2.",
  "In Model E, neighbors’ reports are the mean scores of other participants within 500 m of each home (excluding the participant), and own perception",
  "is the participant’s score minus that mean; participants with fewer than three neighbors within 500 m were excluded. In Models F and G, main",
  "effects refer to the reference group shown in parentheses; derived rows are the slopes in the other group (main effect plus interaction)."), inf_note))
ft <- width(ft, j = 1, width = 4.3)
write_doc(ft, "Table 3. Individual and collective components of neighborhood problems and cohesion, and differences by socioeconomic position.",
          file.path(tab_dir, "Table3.docx"))

## S1 Table: full model C -----------------------------------------------------------------------------
title_case <- function(x) gsub("(^|[ -])([a-z])", "\\1\\U\\2", tolower(x), perl = TRUE)
lab_term <- function(t) {
  edu <- c("1" = "Some primary school", "2" = "Completed primary school", "3" = "Some secondary school", "4" = "Completed secondary school", "5" = "Post-secondary")
  out <- t
  out[t == "(Intercept)"] <- "Intercept"; out[t == "prob_z"] <- "Perceived neighborhood problems (per SD)"; out[t == "coh_z"] <- "Social cohesion (per SD)"
  out[t == "age"] <- "Age (per 10 years)"; out[t == "female"] <- "Female"; out[t == "YrsInHome"] <- "Years in current home (per 10 years)"
  i <- grepl("^educ_f", t); out[i] <- paste0("Education: ", edu[sub("educ_f", "", t[i])])
  i <- grepl("^wealth_quintile", t); out[i] <- paste0("Wealth quintile: ", sub("wealth_quintile", "", t[i]))
  i <- grepl("^Marital", t); out[i] <- paste0("Marital status: ", sub("Marital", "", t[i]))
  i <- grepl("^Emp", t); out[i] <- paste0("Employment: ", sub("Emp", "", t[i]))
  i <- grepl("^ethn5", t); out[i] <- paste0("Ethnicity: ", sub("ethn5", "", t[i]))
  i <- grepl("^DISTRICT", t); out[i] <- paste0("District: ", title_case(sub("DISTRICT", "", t[i])))
  out
}
C1 <- fC[fC$ty == "cl", ]
per10 <- C1$term %in% c("age", "YrsInHome")
C1[per10, c("estimate", "ci_low", "ci_high")] <- 10 * C1[per10, c("estimate", "ci_low", "ci_high")]
s1 <- data.frame(Term = lab_term(C1$term), Estimate = f2(C1$estimate), `95% CI` = fci(C1$ci_low, C1$ci_high), p = fp(C1$p), check.names = FALSE)
v <- vcmp[vcmp$model == "C", ]; r2 <- mods[mods$model == "C", ][1, ]
s1 <- rbind(s1, data.frame(Term = c("Cluster variance (random intercept)", "Residual variance", "ICC (cluster variance / (cluster + residual variance))",
                                    "R² marginal", "R² conditional", "Respondents (n)", "Reconstructed clusters (n)"),
                           Estimate = c(f2(v$cluster_variance), f2(v$residual_variance), f2(v$icc), f2(r2$r2_marginal), f2(r2$r2_conditional),
                                        format(v$n, big.mark = ","), v$n_clusters), `95% CI` = "", p = "", check.names = FALSE))
ft <- style(flextable(s1)); ft <- hline(ft, i = nrow(C1), border = officer::fp_border(color = "black", width = 0.5), part = "body")
ft <- footer(ft, c(
  "Outcome: PSS-10 score standardized in the analytic sample; estimates are differences in SD. Linear mixed model with a reconstructed-cluster random intercept and district fixed effects (REML).",
  inf_note,
  "Reference categories: sex, male (includes 7 respondents reporting another gender); education, no formal education; wealth quintile, poorest; marital status, married or in a domestic partnership; employment, employed or self-employed; ethnicity, Akan (Grusi, Gurma and Mande grouped with Other); district, Ablekuma Central Municipal.",
  "Neighborhood problems and social cohesion are standardized (mean 0, SD 1). CI, confidence interval; ICC, intraclass correlation; PSS-10, 10-item Perceived Stress Scale; SD, standard deviation."))
ft <- width(ft, j = 1, width = 4.2)
write_doc(ft, "S1 Table. Full results of Model C: perceived neighborhood problems, social cohesion, and perceived stress (PSS-10).",
          file.path(tab_dir, "S1_Table.docx"))

## S2 Table: sensitivity analyses (sections follow the order of the Methods) ---------------------------
term_lab <- c(prob_z = "Neighborhood problems", coh_z = "Social cohesion", `prob_z:coh_z` = "Problems × cohesion",
  prob_b = "Problems, neighbors within 500 m (between)", prob_w = "Problems, own perception (within)",
  coh_b = "Cohesion, neighbors within 500 m (between)", coh_w = "Cohesion, own perception (within)",
  prob_bc = "Problems, same-cluster neighbors (between)", prob_wc = "Problems, own perception (within)",
  coh_bc = "Cohesion, same-cluster neighbors (between)", coh_wc = "Cohesion, own perception (within)",
  crime_z = "Problems: crime and disorder", envir_z = "Problems: environmental", socio_z = "Problems: socioeconomic",
  `crime_z:coh_z` = "Crime and disorder × cohesion", `envir_z:coh_z` = "Environmental × cohesion", `socio_z:coh_z` = "Socioeconomic × cohesion",
  `prob_z:poor` = "Problems × poorer household", `coh_z:poor` = "Cohesion × poorer household", `prob_z:coh_z:poor` = "Problems × cohesion × poorer",
  `prob_z:hiedu` = "Problems × higher education", `coh_z:hiedu` = "Cohesion × higher education", `prob_z:coh_z:hiedu` = "Problems × cohesion × higher education",
  `prob_z at coh_z = -1` = "Problems at low cohesion (−1 SD)", `prob_z at coh_z = 0` = "Problems at mean cohesion", `prob_z at coh_z = +1` = "Problems at high cohesion (+1 SD)")
tl <- function(t) ifelse(t %in% names(term_lab), term_lab[t], sub("^coh_z, ", "Social cohesion, ", sub("^prob_z, ", "Neighborhood problems, ", sub("^prob_z:coh_z, ", "Problems × cohesion, ", t))))
mlab <- function(m) {
  m <- sub("^([A-G])$", "Model \\1", m)
  m <- sub("^C: ", "Model C, ", sub("^D: ", "Model D, ", m)); m <- sub(" only$", " subscale only", m)
  m <- sub("\\bcrime\\b", "crime and disorder", sub("\\benvir\\b", "environmental", sub("\\bsocio\\b", "socioeconomic", m)))
  m <- sub("^([CD]) \\+ distrust", "Model \\1 + distrust", m); m <- sub("^([CD]) \\((.*)\\)$", "Model \\1, \\2", m)
  m <- sub("^E \\(cluster\\)", "Model E, within-cluster leave-self-out means", m)
  m
}
rows_of <- function(d, block, keep_type = "cl", note_type = FALSE) {
  d <- d[d$ty %in% keep_type & !(d$term == "dis_z"), ]
  data.frame(block = block, Model = paste0(mlab(d$model), if (note_type) paste0(" [", d$se_type, "]") else ""),
             Term = tl(d$term), n = format(d$n, big.mark = ","), Estimate = f2(d$estimate), `95% CI` = fci(d$ci_low, d$ci_high),
             p = fp(d$p), check.names = FALSE)
}
sub_a <- function(a) R[R$analysis == a, ]
bin <- sub_a("PSS-10 > 13")
blocks <- list(
  a = list("(a) Ordinal PSS-10: cumulative logit; sandwich SEs (CR2 not available for clm), clustered by reconstructed cluster and by district", "log cumulative odds per SD",
           rows_of(sub_a("Ordinal PSS-10 (clm)"), "a", c("cl", "dist"), TRUE)),
  b = list("(b) Binary PSS-10 > 13: logistic mixed model (model-based SEs) and logistic GLM with district fixed effects (CR2 by reconstructed cluster)", "log odds per SD",
           rbind(rows_of(bin[grepl("glmer", bin$model), ], "b", "mb"), rows_of(bin[grepl("GLM", bin$model), ], "b", "cl"))),
  c = list("(c) Problem subscales", "SD of PSS-10 per SD", rows_of(sub_a("Problem subscales"), "c")),
  d = list("(d) Distrust added (distrust estimate not shown)", "SD of PSS-10 per SD", rows_of(sub_a("Distrust added"), "d")),
  e = list("(e) Cluster definition: alternative clusters (all singletons merged; random intercept and CR2 by that cluster) and OLS with district fixed effects (CR2 by reconstructed cluster and by district)", "SD of PSS-10 per SD",
           rbind(transform(rows_of(sub_a("Alternative cluster definition"), "e"), Model = paste0(Model, ", alternative clusters")),
                 transform(rows_of(sub_a("OLS, district FE, no random effect"), "e", c("cl", "dist"), TRUE), Model = paste0(Model, ", OLS")))),
  f = list("(f) Model E with leave-self-out means within the reconstructed cluster", "SD of PSS-10 per SD", rows_of(sub_a("Model E, within-cluster leave-self-out means"), "f")),
  g = list("(g) Moran’s I of model D cluster-level residuals", "Moran’s I",
           data.frame(block = "g", Model = "Model D", Term = sub("cluster random intercept \\(BLUP\\)", "Cluster random intercept",
                      sub("cluster mean marginal residual", "Cluster mean marginal residual", moran$term)),
                      n = format(moran$n, big.mark = ","), Estimate = formatC(moran$estimate, format = "f", digits = 3), `95% CI` = "",
                      p = fp(moran$p), check.names = FALSE)),
  h = list("(h) CR2 standard errors clustered by district (22 clusters): sensitivity bound for the main models", "SD of PSS-10 per SD",
           rbind(rows_of(sub_a("Main models"), "h", "dist"),
                 transform(rows_of(sub_a("Simple slopes"), "h", "dist"), Model = "Model D, simple slopes"),
                 transform(rows_of(sub_a("Group-specific slopes"), "h", "dist"), Model = paste0(Model, ", derived group-specific slopes")))),
  i = list("(i) Slopes by household wealth quintile (model D plus exposure × wealth quintile interactions)", "SD of PSS-10 per SD",
           transform(rows_of(sub_a("Slopes by wealth quintile"), "i"), Model = "Model D + exposure × wealth quintile")))
out <- do.call(rbind, lapply(names(blocks), function(b) rbind(
  data.frame(block = b, Model = paste0(blocks[[b]][[1]], " [scale: ", blocks[[b]][[2]], "]"), Term = "", n = "", Estimate = "", `95% CI` = "", p = "", check.names = FALSE),
  blocks[[b]][[3]])))
hdr <- which(out$Term == "")
out$Model[-hdr] <- ifelse(duplicated(paste(out$block, out$Model))[-hdr], "", out$Model[-hdr])
ft <- flextable(out[, -1]); ft <- style(ft)
ft <- merge_h_range(ft, i = hdr, j1 = 1, j2 = 6); ft <- bold(ft, i = hdr, part = "body")
ft <- align(ft, j = 1:2, align = "left", part = "all")
ft <- footer(ft, c(
  paste("Unless stated otherwise, models are linear mixed models with a reconstructed-cluster random intercept, district fixed effects and the covariates of Model C,",
        "with CR2 (Bell-McCaffrey) cluster-robust SEs clustered by reconstructed cluster and Satterthwaite degrees of freedom; exposures are standardized (per SD)."),
  "(b) clubSandwich does not support logistic mixed models, so the CR2 version is a logistic GLM with the same fixed effects and no random intercept; four districts with 100% PSS-10 > 13 have quasi-separated fixed effects.",
  "(g) Moran’s I with permutation p values (999 permutations); expected value under no autocorrelation about −0.005.",
  "CI, confidence interval; GLM, generalized linear model; OLS, ordinary least squares; PSS-10, 10-item Perceived Stress Scale; SD, standard deviation; SE, standard error."))
ft <- width(ft, j = 1, width = 2.7); ft <- width(ft, j = 2, width = 2.5)
write_doc(ft, "S2 Table. Sensitivity analyses.", file.path(tab_dir, "S2_Table.docx"), landscape = TRUE)

## Figures --------------------------------------------------------------------------------------------
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e3e2de"; pt <- ggplot2::.pt
theme_plos <- theme_classic(base_size = 10, base_family = "Arial") +
  theme(plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3), axis.line = element_line(colour = ink2, linewidth = 0.4),
        axis.ticks = element_line(colour = ink2, linewidth = 0.4), axis.text = element_text(colour = ink2, size = 9),
        axis.title = element_text(colour = ink, size = 10), strip.background = element_blank(),
        strip.text = element_text(colour = ink, size = 10, face = "bold", hjust = 0), legend.position = "none")
save_fig <- function(g, name, w, h) {
  ragg::agg_tiff(file.path(fig_dir, paste0(name, ".tif")), width = w, height = h, units = "in", res = 300, compression = "lzw", background = "white")
  print(g); invisible(dev.off())
  ragg::agg_png(file.path(fig_dir, paste0(name, "_preview.png")), width = w, height = h, units = "in", res = 150, background = "white")
  print(g); invisible(dev.off())
}
## Fig 1: predicted PSS-10 (model D), CIs from CR2 by cluster
P <- pred[grepl("reconstructed cluster", pred$se_type), ]
labs1 <- c("Low cohesion (−1 SD)", "Mean cohesion", "High cohesion (+1 SD)")
P$cohesion <- factor(P$coh_z, levels = c(-1, 0, 1), labels = labs1)
blues <- setNames(c("#86b6ef", "#2a78d6", "#104281"), labs1); ltys <- setNames(c("dashed", "solid", "dotdash"), labs1)
ends <- P[abs(P$prob_z - 2) < 1e-9, ]
g1 <- ggplot(P, aes(prob_z, pss, group = cohesion)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = cohesion), alpha = 0.10) +
  geom_line(aes(colour = cohesion, linetype = cohesion), linewidth = 0.9) +
  geom_text(data = ends, aes(label = cohesion), hjust = 0, nudge_x = 0.06, size = 9 / pt, family = "Arial", colour = ink) +
  scale_colour_manual(values = blues) + scale_fill_manual(values = blues) + scale_linetype_manual(values = ltys) +
  scale_x_continuous(breaks = -2:2, labels = c("−2", "−1", "0", "1", "2")) +
  coord_cartesian(xlim = c(-2, 2), expand = FALSE, clip = "off") +
  labs(x = "Perceived neighborhood problems (SD from mean)", y = "Predicted PSS-10 score") +
  theme_plos + theme(plot.margin = margin(8, 125, 6, 6))
save_fig(g1, "Fig1", 7.5, 5)
## Fig 2: model E coefficients, CR2 by cluster (filled: 500 m; hollow: same reconstructed cluster)
row_lab <- c(w = "Own perception (within)", b = "Neighbors’ reports (between)")
mk <- function(d, src) data.frame(panel = ifelse(grepl("^prob", d$term), "Neighborhood problems", "Social cohesion"),
                                  row = row_lab[ifelse(grepl("_w", d$term), "w", "b")], source = src, estimate = d$estimate, lo = d$ci_low, hi = d$ci_high)
src_lab <- c("Neighbors within 500 m (main)", "Neighbors in the same reconstructed cluster (sensitivity)")
F2 <- rbind(mk(R[R$analysis == "Main models" & R$model == "E" & R$ty == "cl", ], src_lab[1]),
            mk(R[R$analysis == "Model E, within-cluster leave-self-out means" & R$ty == "cl", ], src_lab[2]))
F2$row <- factor(F2$row, levels = rev(row_lab)); F2$source <- factor(F2$source, levels = src_lab)
dodge <- position_dodge(width = 0.45)
g2 <- ggplot(F2, aes(estimate, row, shape = source, linetype = source)) +
  geom_vline(xintercept = 0, colour = ink2, linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.15, orientation = "y", colour = "#2a78d6", linewidth = 0.6, position = dodge) +
  geom_point(size = 2.8, colour = "#2a78d6", fill = "white", stroke = 0.9, position = dodge) +
  scale_shape_manual(values = c(16, 21), name = NULL) + scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_x_continuous(labels = function(x) sub("^-", "−", format(x))) +
  facet_wrap(~ panel, nrow = 1) +
  labs(x = "Difference in perceived stress (SD) per SD of exposure", y = NULL) +
  theme_plos + theme(legend.position = "top", legend.justification = "left", legend.text = element_text(size = 9, colour = ink),
                     panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.3),
                     axis.text.y = element_text(colour = ink, size = 9), panel.spacing = unit(14, "pt"), legend.key.width = unit(22, "pt"))
save_fig(g2, "Fig2", 7.5, 3.5)
cat("Tables and figures written to", tab_dir, "and", fig_dir, "\n")
