## 02_diagnostics.R ------------------------------------------------------------------------------
## Residual diagnostics for models C-E: residuals vs fitted, Q-Q plots of residuals and of cluster random
## intercepts, and summary statistics (skewness, kurtosis, residual variance by fifth of fitted values,
## Breusch-Pagan-type test). Data preparation and model fits are taken verbatim from 01_main_analysis.R
## (everything up to the main model fits), with its output folders redirected to a temporary directory.
## Fixed effects are checked against outputs/tables/models.csv, so run 01_main_analysis.R first.
## Outputs: outputs/figures/diagnostics/ (PNG plots and diagnostics_summary.csv, aggregates only).

src  <- readLines("01_main_analysis.R")
stop_at <- grep("^fits <- lapply\\(rhs, fit_lmm\\)", src); stopifnot(length(stop_at) == 1)
tmp <- gsub("\\\\", "/", tempfile("scs_")); dir.create(tmp)
dir_lines <- grep("^(tab_dir|fig_dir|der_dir) +<- ", src[1:stop_at]); stopifnot(length(dir_lines) == 3)
src[dir_lines] <- sprintf('%s <- "%s"', sub(" .*", "", src[dir_lines]), tmp)
eval(parse(text = src[1:stop_at]))

diag_dir <- file.path(out_dir, "figures", "diagnostics")
dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

## fixed effects must match the main analysis
pub <- read.csv(file.path(out_dir, "tables", "models.csv"))
for (k in c("C", "D", "E")) {
  b <- fixef(fits[[k]]); pb <- pub[pub$model == k, ]
  stopifnot(isTRUE(all.equal(unname(b[pb$term]), pb$estimate, tolerance = 1e-6)))
}
cat("Fixed effects of models C-E reproduce models.csv\n")

ink2 <- "#52514e"; pt_col <- "#2a78d6"; ln_col <- "#0b0b0b"
theme_d <- theme_classic(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA), axis.text = element_text(colour = ink2),
        plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8.5, colour = ink2))
save_png <- function(g, f) ggsave(file.path(diag_dir, f), g, width = 5, height = 4.2, dpi = 200, bg = "white")

summ <- NULL
for (k in c("C", "D", "E")) {
  m  <- fits[[k]]
  d  <- data.frame(fitted = fitted(m), resid = resid(m))            # conditional (within-cluster) residuals
  re <- ranef(m)$cluster[, 1]
  lab <- paste("Model", k)
  sub_n <- sprintf("n = %s respondents, %d clusters", format(nobs(m), big.mark = ","), ngrps(m)[["cluster"]])

  g1 <- ggplot(d, aes(fitted, resid)) +
    geom_hline(yintercept = 0, colour = ink2, linewidth = 0.4) +
    geom_point(colour = pt_col, alpha = 0.25, size = 0.9, stroke = 0) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = ln_col, linewidth = 0.7) +
    labs(title = paste0(lab, ": residuals vs fitted"), subtitle = paste0(sub_n, "; conditional residuals; black line = loess smooth"),
         x = "Fitted value (PSS-10, SD units)", y = "Residual") + theme_d
  g2 <- ggplot(d, aes(sample = resid)) + stat_qq(colour = pt_col, alpha = 0.35, size = 0.9, stroke = 0) +
    stat_qq_line(colour = ln_col, linewidth = 0.6) +
    labs(title = paste0(lab, ": Q-Q plot of residuals"), subtitle = sub_n, x = "Theoretical normal quantile", y = "Sample quantile") + theme_d
  g3 <- ggplot(data.frame(re = re), aes(sample = re)) + stat_qq(colour = pt_col, size = 1.4) +
    stat_qq_line(colour = ln_col, linewidth = 0.6) +
    labs(title = paste0(lab, ": Q-Q plot of cluster intercepts"), subtitle = sub_n, x = "Theoretical normal quantile", y = "Cluster intercept (BLUP)") + theme_d
  save_png(g1, sprintf("model%s_resid_vs_fitted.png", k))
  save_png(g2, sprintf("model%s_qq_residuals.png", k))
  save_png(g3, sprintf("model%s_qq_cluster_intercepts.png", k))

  ## aggregate numbers behind the plots
  sk <- function(x) mean((x - mean(x))^3) / sd(x)^3
  ku <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3
  qf <- cut(d$fitted, quantile(d$fitted, 0:5 / 5), include.lowest = TRUE)
  v5 <- tapply(d$resid, qf, var)
  bp <- summary(lm(I(resid^2) ~ fitted, data = d))                 # Breusch-Pagan-type (studentized: n * R2)
  summ <- rbind(summ, data.frame(model = k, n = nobs(m), n_clusters = ngrps(m)[["cluster"]],
    resid_skewness = sk(d$resid), resid_excess_kurtosis = ku(d$resid),
    resid_pct_beyond_3sd = 100 * mean(abs(d$resid / sd(d$resid)) > 3),
    var_resid_by_fitted_quintile = paste(sprintf("%.2f", v5), collapse = " / "),
    var_ratio_max_min_quintile = max(v5) / min(v5),
    bp_stat_nR2 = nobs(m) * bp$r.squared, bp_p = pchisq(nobs(m) * bp$r.squared, 1, lower.tail = FALSE),
    spearman_absresid_fitted = cor(abs(d$resid), d$fitted, method = "spearman"),
    re_skewness = sk(re), re_excess_kurtosis = ku(re), re_shapiro_p = shapiro.test(re)$p.value))
}
write.csv(summ, file.path(diag_dir, "diagnostics_summary.csv"), row.names = FALSE)
print(summ, digits = 3, row.names = FALSE)
unlink(tmp, recursive = TRUE)
