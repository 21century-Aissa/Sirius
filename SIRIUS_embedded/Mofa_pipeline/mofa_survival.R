#!/usr/bin/env Rscript
# =============================================================================
# MOFA SURVIVAL — Modèles prédictifs de résultats cliniques
# =============================================================================
# Cox Proportional Hazards Models + Kaplan-Meier plots
# Basé sur les facteurs latents inférés par MOFA2
#
# Sections :
#   SECTION 0  — Lecture des arguments
#   SECTION 1  — Chargement du modèle & construction de l'objet Surv
#   SECTION 2  — Ajustement du modèle de Cox (coxph)
#   SECTION 3  — Plot des Hazard Ratios (forest plot)
#   SECTION 4  — Plot Kaplan-Meier (un ou plusieurs facteurs)
#   SECTION 5  — Tableau de résumé des p-values (tous facteurs)
# =============================================================================

suppressMessages({
  library(MOFA2)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
})

# ── Source des utilitaires partagés (lib/) ───────────────────────────────────
.script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
source(file.path(.script_dir, "lib", "mofa_utils.R"))
source(file.path(.script_dir, "lib", "mofa_io.R"))

# =============================================================================
# SECTION 0 — Lecture des arguments
# =============================================================================

args   <- commandArgs(trailingOnly = TRUE)
params <- parse_args(args)

model_path   <- params[["model_path"]]
work_dir     <- params[["work_dir"]]
analysis     <- params[["analysis"]]    # "cox_hr" | "kaplan_meier" | "cox_summary"
out_name     <- params[["out_name"]]

# Paramètres survival (colonnes des metadata)
col_time     <- params[["col_time"]]    # ex: "TTT"
col_event    <- params[["col_event"]]   # ex: "treatedAfter"

# Paramètres spécifiques selon le type
factor_index <- if (!is.null(params[["factor_index"]])) as.integer(params[["factor_index"]]) else 1
group        <- params[["group"]]       # pour KM : colonne metadata à utiliser comme group (optionnel)
km_label_low <- if (!is.null(params[["km_label_low"]])  && params[["km_label_low"]]  != "") params[["km_label_low"]]  else paste("Low LF",  factor_index)
km_label_high<- if (!is.null(params[["km_label_high"]]) && params[["km_label_high"]] != "") params[["km_label_high"]] else paste("High LF", factor_index)
xlab_surv    <- if (!is.null(params[["xlab_surv"]])     && params[["xlab_surv"]]     != "") params[["xlab_surv"]]     else "Time"
title_surv   <- if (!is.null(params[["title_surv"]])    && params[["title_surv"]]    != "") params[["title_surv"]]    else paste("Factor", factor_index)
conf_int     <- if (!is.null(params[["conf_int"]]))     as.logical(params[["conf_int"]])    else TRUE
pval_show    <- if (!is.null(params[["pval_show"]]))    as.logical(params[["pval_show"]])   else TRUE
palette_km   <- if (!is.null(params[["palette_km"]])    && params[["palette_km"]]    != "") params[["palette_km"]]    else "hue"
preview_only <- isTRUE(params[["preview_only"]] == "TRUE")  # PNG only, no PDF

setwd(work_dir)

mofa_header("MOFA SURVIVAL", analysis, list(
  "Temps"      = col_time,
  "Événement" = col_event
))

# =============================================================================
# SECTION 1 — Chargement du modèle & objet Survie
# =============================================================================

cat("[ÉTAPE 1] Chargement du modèle MOFA...\n")
mofa_model <- load_model_auto(model_path)

meta <- samples_metadata(mofa_model)

# Vérification des colonnes
for (col in c(col_time, col_event)) {
  if (!col %in% colnames(meta)) {
    cat(sprintf(
      "  [ERREUR] Colonne '%s' absente des métadonnées.\n  Colonnes disponibles : %s\n",
      col, paste(colnames(meta), collapse = ", ")
    ))
    quit(status = 1)
  }
}

cat(sprintf("[ÉTAPE 1b] Construction de l'objet Surv (%s, %s)...\n", col_time, col_event))
SurvObject <- Surv(meta[[col_time]], as.numeric(meta[[col_event]]))

# Facteurs latents — groupe 1 (première group/view)
Z <- get_factors(mofa_model)[[1]]
cat(sprintf("  ✓ Matrice de facteurs : %d samples × %d facteurs\n\n", nrow(Z), ncol(Z)))

outfile <- file.path(work_dir, paste0(out_name, ".pdf"))

# =============================================================================
# SECTION 2 & 3 — Cox Model + Forest Plot des Hazard Ratios
# =============================================================================

if (analysis %in% c("cox_hr", "cox_summary")) {

  cat("[ÉTAPE 2] Ajustement du modèle Cox (tous facteurs)...\n")

  fit_cox <- coxph(SurvObject ~ Z)
  s       <- summary(fit_cox)
  coef_mx <- s[["coefficients"]]

  cat("  ✓ Modèle Cox ajusté.\n")
  cat(sprintf("  Concordance : %.3f\n", s$concordance["C"]))
  cat(sprintf("  Likelihood ratio test p-value : %.4f\n\n", s$logtest["pvalue"]))

  # ── Tableau de résumé exporté en CSV
  df_coef <- data.frame(
    factor    = rownames(coef_mx),
    coef      = coef_mx[, "coef"],
    exp_coef  = coef_mx[, "exp(coef)"],
    se_coef   = coef_mx[, "se(coef)"],
    z         = coef_mx[, "z"],
    p_value   = coef_mx[, "Pr(>|z|)"],
    lower_95  = s[["conf.int"]][, "lower .95"],
    upper_95  = s[["conf.int"]][, "upper .95"],
    signif    = ifelse(coef_mx[, "Pr(>|z|)"] < 0.001, "***",
               ifelse(coef_mx[, "Pr(>|z|)"] < 0.01,  "**",
               ifelse(coef_mx[, "Pr(>|z|)"] < 0.05,  "*",
               ifelse(coef_mx[, "Pr(>|z|)"] < 0.1,   ".", ""))))
  )

  csv_path <- file.path(work_dir, paste0(out_name, "_cox_summary.csv"))
  write.csv(df_coef, csv_path, row.names = FALSE)
  cat(sprintf("  ✓ Tableau Cox exporté : %s\n\n", csv_path))

  if (analysis == "cox_hr") {
    cat("[ÉTAPE 3] Forest plot des Hazard Ratios...\n")

    df_plot <- data.frame(
      factor = factor(rownames(coef_mx), levels = rev(rownames(coef_mx))),
      p      = coef_mx[, "Pr(>|z|)"],
      coef   = coef_mx[, "exp(coef)"],
      lower  = s[["conf.int"]][, "lower .95"],
      higher = s[["conf.int"]][, "upper .95"]
    )

    # Couleur selon significativité
    df_plot$sig_color <- ifelse(df_plot$p < 0.05, "#FF6B6B", "#619CFF")

    p_forest <- ggplot(df_plot,
        aes(x = factor, y = coef, ymin = lower, ymax = higher, color = sig_color)) +
      geom_pointrange(size = 0.8, linewidth = 0.9) +
      geom_hline(aes(yintercept = 1), linetype = "dotted", color = "grey50", linewidth = 0.8) +
      scale_color_identity(
        guide  = "legend",
        labels = c("#FF6B6B" = "p < 0.05", "#619CFF" = "p ≥ 0.05"),
        name   = "Significance"
      ) +
      coord_flip() +
      scale_x_discrete() +
      labs(
        title    = "Cox Model — Hazard Ratios per Factor",
        subtitle = paste("Outcome:", col_time, "·", col_event),
        y        = "Hazard Ratio (95% CI)",
        x        = ""
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title    = element_text(face = "bold"),
        legend.position = "bottom"
      )

    save_gg(p_forest, outfile, width = 9, height = max(5, nrow(df_plot) * 0.5 + 2))
  }

  if (analysis == "cox_summary") {
    cat("[ÉTAPE 3] Tableau p-values tous facteurs (bubble plot)...\n")

    df_plot <- df_coef
    df_plot$factor <- factor(df_plot$factor, levels = rev(df_plot$factor))
    df_plot$neg_log_p <- -log10(df_plot$p_value)

    p_bubble <- ggplot(df_plot,
        aes(x = factor, y = neg_log_p, size = abs(coef), color = coef)) +
      geom_point(alpha = 0.85) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.7) +
      scale_color_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                            midpoint = 0, name = "Coefficient") +
      scale_size_continuous(range = c(3, 12), name = "|coef|") +
      coord_flip() +
      labs(
        title    = "Cox Model — Factor significance overview",
        subtitle = paste("Outcome:", col_time, "·", col_event),
        y        = "-log10(p-value)  [dashed = 0.05]",
        x        = ""
      ) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(face = "bold"))

    save_gg(p_bubble, outfile, width = 9, height = max(5, nrow(df_coef) * 0.5 + 2))
  }
}

# =============================================================================
# SECTION 4 — Kaplan-Meier
# =============================================================================

if (analysis == "kaplan_meier") {

  cat(sprintf("[ÉTAPE 4] Kaplan-Meier pour le Facteur %d...\n", factor_index))

  # Vérification que le facteur existe
  if (factor_index > ncol(Z)) {
    cat(sprintf("  [ERREUR] Facteur %d inexistant (max = %d).\n", factor_index, ncol(Z)))
    quit(status = 1)
  }

  factor_col <- paste0("Factor", factor_index)
  df_km <- data.frame(
    time  = SurvObject[, 1],
    event = SurvObject[, 2],
    Z_col = Z[, factor_index]
  )
  names(df_km)[3] <- factor_col

  # ── Coupure optimale par maximally selected rank statistics
  if (!requireNamespace("maxstat", quietly = TRUE)) {
    cat("  [INFO] Package maxstat non disponible. Utilisation de la médiane comme coupure.\n")
    cutpoint_val <- median(df_km[[factor_col]], na.rm = TRUE)
  } else {
    suppressMessages(library(maxstat))
    cut_res      <- surv_cutpoint(df_km, time = "time", event = "event", variables = factor_col)
    cutpoint_val <- cut_res$cutpoint$cutpoint
    cat(sprintf("  → Coupure optimale (maxstat) : %.4f\n", cutpoint_val))
  }

  df_km$FactorGroup <- ifelse(df_km[[factor_col]] > cutpoint_val, "HIGH", "LOW")
  df_km$FactorGroup <- factor(df_km$FactorGroup, levels = c("LOW", "HIGH"))

  n_low  <- sum(df_km$FactorGroup == "LOW",  na.rm = TRUE)
  n_high <- sum(df_km$FactorGroup == "HIGH", na.rm = TRUE)
  cat(sprintf("  → Groupe LOW  : %d samples\n", n_low))
  cat(sprintf("  → Groupe HIGH : %d samples\n", n_high))

  fit_km <- survfit(Surv(time, event) ~ FactorGroup, data = df_km)

  p_km <- ggsurvplot(
    fit_km,
    data          = df_km,
    conf.int      = conf_int,
    pval          = pval_show,
    pval.method   = TRUE,
    fun           = function(y) y * 100,
    palette       = if (palette_km == "hue") c("#619CFF", "#FF6B6B") else palette_km,
    legend        = "top",
    legend.title  = paste("Factor", factor_index),
    legend.labs   = c(
      paste0(km_label_low,  " (n=", n_low,  ")"),
      paste0(km_label_high, " (n=", n_high, ")")
    ),
    xlab          = xlab_surv,
    ylab          = "Survival probability (%)",
    title         = title_surv,
    subtitle      = paste("Cutpoint =", round(cutpoint_val, 4)),
    risk.table    = TRUE,
    risk.table.height = 0.25,
    ggtheme       = theme_bw(base_size = 13)
  )

  if (!preview_only) {
    pdf(outfile, width = 10, height = 8)
    print(p_km)
    dev.off()
    cat(sprintf("  ✓ Kaplan-Meier PDF : %s\n", outfile))
  }
  png_path <- sub("\\.pdf$", "_preview.png", outfile)
  grDevices::png(png_path, width = 1000, height = 800, res = 100)
  print(p_km)
  dev.off()
  cat(sprintf("  ✓ Preview PNG : %s\n\n", png_path))
}

# =============================================================================
# SECTION 5 — Multi-factor KM (tous facteurs, une page par facteur)
# =============================================================================

if (analysis == "kaplan_meier_all") {

  cat("[ÉTAPE 5] Kaplan-Meier — tous les facteurs...\n")

  # Pour preview: uniquement le premier facteur en PNG
  if (preview_only) {
    png_path <- sub("\\.pdf$", "_preview.png", outfile)
    grDevices::png(png_path, width = 1000, height = 800, res = 100)
  } else {
    pdf(outfile, width = 10, height = 8)
  }

  for (fi in seq_len(ncol(Z))) {
    factor_col <- paste0("Factor", fi)
    df_km <- data.frame(
      time  = SurvObject[, 1],
      event = SurvObject[, 2],
      Z_col = Z[, fi]
    )
    names(df_km)[3] <- factor_col

    if (!requireNamespace("maxstat", quietly = TRUE)) {
      cutpoint_val <- median(df_km[[factor_col]], na.rm = TRUE)
    } else {
      suppressMessages(library(maxstat))
      tryCatch({
        cut_res      <- surv_cutpoint(df_km, time = "time", event = "event", variables = factor_col)
        cutpoint_val <- cut_res$cutpoint$cutpoint
      }, error = function(e) {
        cutpoint_val <<- median(df_km[[factor_col]], na.rm = TRUE)
      })
    }

    df_km$FactorGroup <- factor(
      ifelse(df_km[[factor_col]] > cutpoint_val, "HIGH", "LOW"),
      levels = c("LOW", "HIGH")
    )

    fit_km <- tryCatch(
      survfit(Surv(time, event) ~ FactorGroup, data = df_km),
      error = function(e) NULL
    )

    if (!is.null(fit_km)) {
      p_km <- ggsurvplot(
        fit_km, data = df_km,
        conf.int    = conf_int,
        pval        = pval_show,
        fun         = function(y) y * 100,
        palette     = c("#619CFF", "#FF6B6B"),
        legend      = "top",
        legend.labs = c(paste("Low LF", fi), paste("High LF", fi)),
        xlab        = xlab_surv,
        ylab        = "Survival probability (%)",
        title       = paste("Factor", fi),
        subtitle    = paste("Cutpoint =", round(cutpoint_val, 4)),
        risk.table  = TRUE, risk.table.height = 0.25,
        ggtheme     = theme_bw(base_size = 12)
      )
      print(p_km)
      cat(sprintf("  ✓ Facteur %d plotté\n", fi))
    } else {
      cat(sprintf("  [AVERTISSEMENT] Facteur %d ignoré (erreur survfit).\n", fi))
    }
  }
  dev.off()
  cat(sprintf("\n  ✓ Tous les KM enregistrés : %s\n\n", outfile))
}

# =============================================================================
# FIN
# =============================================================================

cat("========================================================\n")
cat(sprintf("  SURVIVAL ANALYSIS TERMINÉE : %s.pdf\n", out_name))
cat("========================================================\n")

# =============================================================================
# SECTION 6 — Survival sur clusters (kaplan_meier_clusters / cox_clusters /
#             cox_diagnostics)
# =============================================================================

kaplan_meier_clusters <- function(surv_data, output_dir) {
  cat("[kaplan_meier_clusters] Lecture clusters.rds...\n")
  clusters <- load_clusters(output_dir)

  df_cl <- data.frame(
    sample  = names(clusters),
    cluster = factor(clusters),
    stringsAsFactors = FALSE
  )
  surv_data$sample <- rownames(surv_data)
  merged <- merge(surv_data, df_cl, by = "sample")
  cat(sprintf("  ✓ %d samples appariés\n", nrow(merged)))

  fit <- survfit(Surv(time, status) ~ cluster, data = merged)
  p_km <- ggsurvplot(
    fit, data = merged,
    pval       = TRUE,
    conf.int   = TRUE,
    risk.table = TRUE,
    palette    = "nejm"
  )

  outfile <- file.path(output_dir, "km_clusters.pdf")
  save_gg(p_km, outfile, width = 10, height = 8)
}

cox_clusters <- function(surv_data, covariates = NULL, output_dir) {
  cat("[cox_clusters] Lecture clusters.rds...\n")
  clusters <- load_clusters(output_dir)

  df_cl <- data.frame(
    sample  = names(clusters),
    cluster = factor(clusters),
    stringsAsFactors = FALSE
  )
  surv_data$sample <- rownames(surv_data)
  merged <- merge(surv_data, df_cl, by = "sample")
  cat(sprintf("  ✓ %d samples appariés\n", nrow(merged)))

  rhs <- "cluster"
  if (!is.null(covariates) && length(covariates) > 0) {
    cov_present <- intersect(covariates, colnames(merged))
    if (length(cov_present) > 0) {
      rhs <- paste(c("cluster", cov_present), collapse = " + ")
    }
  }
  fml <- as.formula(paste0("Surv(time, status) ~ ", rhs))
  cat(sprintf("  ✓ Formule : %s\n", deparse(fml)))

  fit <- coxph(fml, data = merged)

  p_forest <- ggforest(fit, data = merged)
  outfile <- file.path(output_dir, "cox_clusters.pdf")
  save_gg(p_forest, outfile, width = 10, height = 8)

  sum_path <- file.path(output_dir, "cox_clusters_summary.txt")
  sink(sum_path)
  print(summary(fit))
  sink()
  cat(sprintf("  ✓ Résumé enregistré : %s\n", sum_path))

  invisible(fit)
}

cox_diagnostics <- function(cox_model, output_dir) {
  cat("[cox_diagnostics] cox.zph...\n")
  zph <- cox.zph(cox_model)

  p_zph <- ggcoxzph(zph)
  outfile <- file.path(output_dir, "cox_zph.pdf")
  save_gg(p_zph, outfile, width = 10, height = 8)

  sum_path <- file.path(output_dir, "cox_zph_summary.txt")
  sink(sum_path)
  print(zph)
  sink()
  cat(sprintf("  ✓ Résumé enregistré : %s\n", sum_path))
}

# --- Dispatch additionnel : analyses basées sur clusters ---
if (analysis %in% c("kaplan_meier_clusters", "cox_clusters", "cox_diagnostics")) {
  surv_data <- data.frame(
    time   = SurvObject[, 1],
    status = as.numeric(SurvObject[, 2]),
    row.names = rownames(Z),
    stringsAsFactors = FALSE
  )
  # Joindre les colonnes meta (covariables potentielles)
  meta_join <- meta
  if (!is.null(rownames(meta_join))) {
    common_cols <- setdiff(colnames(meta_join), c("time", "status", "sample"))
    for (cc in common_cols) {
      surv_data[[cc]] <- meta_join[[cc]][match(rownames(surv_data), rownames(meta_join))]
    }
  }

  cov_str <- params[["covariates"]]
  covariates <- NULL
  if (!is.null(cov_str) && nzchar(cov_str)) {
    covariates <- trimws(strsplit(cov_str, ",", fixed = TRUE)[[1]])
    covariates <- covariates[nzchar(covariates)]
  }

  if (analysis == "kaplan_meier_clusters") {
    kaplan_meier_clusters(surv_data, work_dir)
  } else if (analysis == "cox_clusters") {
    cox_clusters(surv_data, covariates, work_dir)
  } else if (analysis == "cox_diagnostics") {
    fit <- cox_clusters(surv_data, covariates, work_dir)
    cox_diagnostics(fit, work_dir)
  }

  cat("========================================================\n")
  cat(sprintf("  CLUSTER SURVIVAL ANALYSIS TERMINÉE : %s\n", analysis))
  cat("========================================================\n")
}
