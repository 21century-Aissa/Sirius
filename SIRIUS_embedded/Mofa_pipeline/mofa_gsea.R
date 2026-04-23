#!/usr/bin/env Rscript
# =============================================================================
# MOFA GSEA — Gene Set Enrichment Analysis
# =============================================================================
# Analyse d'enrichissement de gene sets sur les poids des facteurs MOFA.
# Basé sur run_enrichment() de MOFA2 (implémentation PCGSE modifiée).
#
# SECTIONS :
#   SECTION 0  — Parsing des arguments
#   SECTION 1  — Chargement du modèle + gene set matrix
#   SECTION 2  — Préparation & harmonisation des noms de gènes
#   SECTION 3  — run_enrichment() (un ou plusieurs tests)
#   SECTION 4  — plot_enrichment_heatmap  (vue globale)
#   SECTION 5  — plot_enrichment          (top pathways par facteur)
#   SECTION 6  — plot_enrichment_detailed (gènes contributeurs)
#   SECTION 7  — Comparaison deux tests   (histogrammes + scatter)
#   SECTION 8  — Export CSV des résultats (p-values, p-adj, statistiques)
# =============================================================================

suppressMessages({
  library(MOFA2)
  library(ggplot2)
  library(data.table)
  library(purrr)
})

# ── cowplot optionnel (comparaison de tests)
has_cowplot <- requireNamespace("cowplot", quietly = TRUE)
if (has_cowplot) suppressMessages(library(cowplot))

# =============================================================================
# UTILITAIRES
# =============================================================================

parse_args <- function(args) {
  params <- list()
  for (arg in args) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    if (length(parts) >= 2) {
      key   <- parts[1]
      value <- paste(parts[-1], collapse = "=")
      params[[key]] <- value
    }
  }
  return(params)
}

save_gg <- function(p, filepath, width = 10, height = 7) {
  if (!preview_only) {
    pdf(filepath, width = width, height = height)
    if (inherits(p, "gg") || inherits(p, "ggplot")) print(p) else print(p)
    dev.off()
    cat(sprintf("  ✓ PDF enregistré : %s\n", filepath))
  }
  png_path <- sub("\\.pdf$", "_preview.png", filepath)
  grDevices::png(png_path, width = round(width * 100), height = round(height * 100), res = 100)
  print(p)
  dev.off()
  cat(sprintf("  ✓ Preview PNG : %s\n", png_path))
}

# =============================================================================
# SECTION 0 — Lecture des arguments
# =============================================================================

args   <- commandArgs(trailingOnly = TRUE)
params <- parse_args(args)

model_path     <- params[["model_path"]]
work_dir       <- params[["work_dir"]]
analysis       <- params[["analysis"]]      # voir liste ci-dessous
out_name       <- params[["out_name"]]

# Paramètres gene sets
geneset_source <- params[["geneset_source"]]  # "reactome" | "msigdb_c2_human" |
                                               # "msigdb_c5_human" | "msigdb_c2_mouse" |
                                               # "msigdb_c5_mouse" | "custom"
custom_gs_path <- params[["custom_gs_path"]]   # si geneset_source == "custom"

# Paramètres run_enrichment
view_gsea      <- params[["view_gsea"]]        # ex: "mRNA"
factors_gsea   <- params[["factors_gsea"]]     # ex: "1:3" ou "1,2,3"
sign_gsea      <- params[["sign_gsea"]]        # "positive" | "negative" | "all"
stat_test      <- params[["stat_test"]]        # "parametric" | "cor.adj.parametric" | "permutation"
stat_test2     <- params[["stat_test2"]]       # second test (pour comparaison, section 7)
gene_case      <- params[["gene_case"]]        # "upper" | "lower" | "as_is"

# Paramètres pour les nouveaux plots
result_type    <- params[["result_type"]]      # "positive" | "negative" (pour charger un résultat sauvegardé)
load_result    <- params[["load_result"]]      # "TRUE" | "FALSE" (pour charger un résultat depuis mémoire)
factor_plot    <- if (!is.null(params[["factor_plot"]]))  as.integer(params[["factor_plot"]])  else 5
max_pathways   <- if (!is.null(params[["max_pathways"]])) as.integer(params[["max_pathways"]]) else 15
preview_only   <- isTRUE(params[["preview_only"]] == "TRUE")  # PNG only, no PDF

# Paramètres de plot
max_genes      <- if (!is.null(params[["max_genes"]]))    as.integer(params[["max_genes"]])    else 8
fdr_threshold  <- if (!is.null(params[["fdr_threshold"]])) as.numeric(params[["fdr_threshold"]]) else 0.05

setwd(work_dir)

cat("========================================================\n")
cat(sprintf("  MOFA GSEA — %s\n", analysis))
cat(sprintf("  Modèle      : %s\n", model_path))
cat(sprintf("  Vue         : %s\n", view_gsea))
cat(sprintf("  Facteurs    : %s\n", factors_gsea))
cat(sprintf("  Signe       : %s\n", sign_gsea))
cat(sprintf("  Test stat.  : %s\n", stat_test))
cat(sprintf("  Gene sets   : %s\n", geneset_source))
cat("========================================================\n\n")

# =============================================================================
# SECTION 1 — Chargement du modèle
# =============================================================================

cat("[ÉTAPE 1] Chargement du modèle MOFA...\n")
# Support both RDS (metadata included) and HDF5
if (tolower(tools::file_ext(model_path)) == "rds") {
  mofa_model <- readRDS(model_path)
  cat("  ✓ Modèle chargé depuis RDS.\n\n")
} else {
  mofa_model <- load_model(model_path)
  cat("  ✓ Modèle chargé depuis HDF5.\n\n")
}

# =============================================================================
# SECTION 2 — Chargement du gene set matrix
# =============================================================================

cat("[ÉTAPE 2] Chargement de la matrice de gene sets...\n")

if (geneset_source == "custom") {

  # ── Matrice personnalisée fournie par l'utilisateur
  if (is.null(custom_gs_path) || !file.exists(custom_gs_path)) {
    cat("  [ERREUR] Fichier gene set personnalisé introuvable :", custom_gs_path, "\n")
    cat("  Format attendu : matrice binaire CSV/RDS (lignes = pathways, colonnes = gènes)\n")
    quit(status = 1)
  }
  ext <- tolower(tools::file_ext(custom_gs_path))
  gs_matrix <- if (ext == "rds") {
    readRDS(custom_gs_path)
  } else {
    m <- read.csv(custom_gs_path, row.names = 1, check.names = FALSE)
    as.matrix(m)
  }
  cat(sprintf("  ✓ Gene set personnalisé : %d pathways × %d gènes\n\n",
              nrow(gs_matrix), ncol(gs_matrix)))

} else {

  # ── Gene sets pré-définis (MOFAdata)
  if (!requireNamespace("MOFAdata", quietly = TRUE)) {
    cat("  [ERREUR] Package MOFAdata requis pour les gene sets pré-définis.\n")
    cat("  Installez-le avec : devtools::install_github('bioFAM/MOFAdata')\n")
    quit(status = 1)
  }
  suppressMessages(library(MOFAdata))

  gs_matrix <- switch(geneset_source,
    "reactome"        = { utils::data("reactomeGS");            get("reactomeGS") },
    "msigdb_c2_human" = { data("MSigDB_v6.0_C2_human"); get("MSigDB_v6.0_C2_human") },
    "msigdb_c5_human" = { data("MSigDB_v6.0_C5_human"); get("MSigDB_v6.0_C5_human") },
    "msigdb_c2_mouse" = { data("MSigDB_v6.0_C2_mouse"); get("MSigDB_v6.0_C2_mouse") },
    "msigdb_c5_mouse" = { data("MSigDB_v6.0_C5_mouse"); get("MSigDB_v6.0_C5_mouse") },
    {
      cat(sprintf("  [ERREUR] Source de gene set inconnue : %s\n", geneset_source))
      quit(status = 1)
    }
  )
  cat(sprintf("  ✓ Gene set '%s' chargé : %d pathways × %d gènes\n\n",
              geneset_source, nrow(gs_matrix), ncol(gs_matrix)))
}

# =============================================================================
# SECTION 3 — Harmonisation des noms de gènes
# =============================================================================

cat("[ÉTAPE 3] Harmonisation des noms de gènes (vue :", view_gsea, ")...\n")

current_features <- features_names(mofa_model)[[view_gsea]]
if (is.null(current_features)) {
  cat(sprintf("  [ERREUR] Vue '%s' non trouvée dans le modèle.\n", view_gsea))
  cat("  Vues disponibles :", paste(names(features_names(mofa_model)), collapse=", "), "\n")
  quit(status = 1)
}

if (gene_case == "upper") {
  features_names(mofa_model)[[view_gsea]] <- toupper(current_features)
  cat("  → Noms de gènes convertis en MAJUSCULES\n")
} else if (gene_case == "lower") {
  features_names(mofa_model)[[view_gsea]] <- tolower(current_features)
  cat("  → Noms de gènes convertis en minuscules\n")
} else {
  cat("  → Noms de gènes conservés tels quels\n")
}

# Résumé overlap
mofa_genes <- features_names(mofa_model)[[view_gsea]]
gs_genes   <- colnames(gs_matrix)
n_overlap  <- length(intersect(mofa_genes, gs_genes))
cat(sprintf("  → Overlap gènes MOFA ∩ gene sets : %d / %d\n\n", n_overlap, length(mofa_genes)))

if (n_overlap == 0) {
  cat("  [ERREUR] Aucun gène en commun entre le modèle et la matrice de gene sets.\n")
  cat("  Vérifiez le paramètre 'Conversion des noms de gènes'.\n")
  quit(status = 1)
}

# =============================================================================
# SECTION 4 — run_enrichment (test principal)
# =============================================================================

cat(sprintf("[ÉTAPE 4] run_enrichment (test = %s, signe = %s)...\n", stat_test, sign_gsea))

# Parsing des facteurs : "1:3" → 1:3, "1,2,3" → c(1,2,3)
factors_parsed <- tryCatch({
  if (grepl(":", factors_gsea)) {
    parts <- as.integer(strsplit(factors_gsea, ":")[[1]])
    seq(parts[1], parts[2])
  } else {
    as.integer(trimws(strsplit(factors_gsea, ",")[[1]]))
  }
}, error = function(e) {
  cat("  [ERREUR] Format de facteurs invalide. Utilisez '1:3' ou '1,2,3'.\n")
  quit(status = 1)
})

cat(sprintf("  → Facteurs analysés : %s\n", paste(factors_parsed, collapse = ", ")))

# ── Cache: on sauvegarde/recharge res.positive ou res.negative selon sign_gsea
cache_file   <- file.path(work_dir, sprintf("res_%s.rds", sign_gsea))
fp_file      <- file.path(work_dir, sprintf("res_%s.fp",  sign_gsea))
current_fp   <- paste(view_gsea, paste(factors_parsed, collapse = ","),
                      sign_gsea, stat_test, geneset_source, sep = "|")

use_cache <- FALSE
if (file.exists(cache_file) && file.exists(fp_file)) {
  cached_fp <- tryCatch(readLines(fp_file, warn = FALSE)[1], error = function(e) "")
  if (identical(cached_fp, current_fp)) {
    use_cache <- TRUE
  }
}

if (use_cache) {
  cat(sprintf("  ✓ Cache trouvé, rechargement depuis %s\n", basename(cache_file)))
  enrichment_result <- readRDS(cache_file)
} else {
  enrichment_result <- tryCatch({
    run_enrichment(
      mofa_model,
      view             = view_gsea,
      factors          = factors_parsed,
      feature.sets     = gs_matrix,
      sign             = sign_gsea,
      statistical.test = stat_test
    )
  }, error = function(e) {
    cat(sprintf("  [ERREUR run_enrichment] %s\n", e$message))
    quit(status = 1)
  })
  cat("  ✓ Enrichissement calculé.\n")
  # Sauvegarde sur disque
  tryCatch({
    saveRDS(enrichment_result, cache_file)
    writeLines(current_fp, fp_file)
    cat(sprintf("  ✓ Sauvegardé : %s\n", basename(cache_file)))
  }, error = function(e) cat(sprintf("  [AVERT.] Cache non sauvé : %s\n", e$message)))
}
cat(sprintf("  → %d gene sets testés\n", nrow(enrichment_result$pval)))

# ── Pathways significatifs par facteur
for (fi in seq_along(factors_parsed)) {
  sig <- enrichment_result$sigPathways[[fi]]
  if (!is.null(sig) && length(sig) > 0) {
    cat(sprintf("  → Facteur %d : %d pathways significatifs (FDR < %.2f)\n",
                factors_parsed[fi], length(sig), fdr_threshold))
  } else {
    cat(sprintf("  → Facteur %d : aucun pathway significatif\n", factors_parsed[fi]))
  }
}
cat("\n")

# =============================================================================
# SECTION 5 — Export CSV des résultats
# =============================================================================

cat("[ÉTAPE 5] Export des résultats en CSV...\n")

# p-values nominales
pval_df <- as.data.frame(enrichment_result$pval)
pval_df$pathway <- rownames(pval_df)
pval_df <- pval_df[, c("pathway", setdiff(names(pval_df), "pathway"))]
write.csv(pval_df, file.path(work_dir, paste0(out_name, "_pvalues.csv")), row.names = FALSE)

# p-values ajustées (FDR)
padj_df <- as.data.frame(enrichment_result$pval.adj)
padj_df$pathway <- rownames(padj_df)
padj_df <- padj_df[, c("pathway", setdiff(names(padj_df), "pathway"))]
write.csv(padj_df, file.path(work_dir, paste0(out_name, "_pvalues_adj.csv")), row.names = FALSE)

# Statistiques de gene sets
setstats_df <- as.data.frame(enrichment_result$set.statistics)
setstats_df$pathway <- rownames(setstats_df)
setstats_df <- setstats_df[, c("pathway", setdiff(names(setstats_df), "pathway"))]
write.csv(setstats_df, file.path(work_dir, paste0(out_name, "_set_statistics.csv")), row.names = FALSE)

cat(sprintf("  ✓ CSV exportés : %s_[pvalues|pvalues_adj|set_statistics].csv\n\n", out_name))

# ═══════════════════════════════════════════════════════════════════════════
# BLOCS DE PLOTS — déclenchés selon analysis=
# ═══════════════════════════════════════════════════════════════════════════
outfile <- file.path(work_dir, paste0(out_name, ".pdf"))

# Variable pour stocker les résultats en mémoire
enrichment_positive <- NULL
enrichment_negative <- NULL

# Sauvegarder les résultats en mémoire
save_enrichment_results <- function() {
  saveRDS(enrichment_positive, file.path(work_dir, "res_positive.rds"))
  saveRDS(enrichment_negative, file.path(work_dir, "res_negative.rds"))
  cat("  ✓ Résultats sauvegardés en mémoire (res_positive.rds, res_negative.rds)\n")
}

# Charger les résultats depuis mémoire
load_enrichment_results <- function() {
  if (file.exists(file.path(work_dir, "res_positive.rds"))) {
    enrichment_positive <<- readRDS(file.path(work_dir, "res_positive.rds"))
    cat("  ✓ res_positive chargé depuis mémoire\n")
  }
  if (file.exists(file.path(work_dir, "res_negative.rds"))) {
    enrichment_negative <<- readRDS(file.path(work_dir, "res_negative.rds"))
    cat("  ✓ res_negative chargé depuis mémoire\n")
  }
}

# =============================================================================
# NEW ANALYSIS — GSEA Positive/Negative with reactomeGS
# analysis == "gsea_positive" or "gsea_negative"
# =============================================================================

if (analysis == "gsea_positive") {
  cat("[ANALYSIS] GSEA on positive weights...\n")
  
  enrichment_positive <- tryCatch({
    run_enrichment(
      mofa_model,
      view             = view_gsea,
      factors          = factors_parsed,
      feature.sets     = gs_matrix,
      sign             = "positive",
      statistical.test = stat_test
    )
  }, error = function(e) {
    cat(sprintf("  [ERREUR run_enrichment] %s\n", e$message))
    quit(status = 1)
  })
  
  cat("  ✓ Enrichissement positive calculé.\n")
  save_enrichment_results()
  cat(sprintf("  ✓ Résultat sauvegardé : %s\n\n", file.path(work_dir, "res_positive.rds")))
}

if (analysis == "gsea_negative") {
  cat("[ANALYSIS] GSEA on negative weights...\n")
  
  enrichment_negative <- tryCatch({
    run_enrichment(
      mofa_model,
      view             = view_gsea,
      factors          = factors_parsed,
      feature.sets     = gs_matrix,
      sign             = "negative",
      statistical.test = stat_test
    )
  }, error = function(e) {
    cat(sprintf("  [ERREUR run_enrichment] %s\n", e$message))
    quit(status = 1)
  })
  
  cat("  ✓ Enrichissement negative calculé.\n")
  save_enrichment_results()
  cat(sprintf("  ✓ Résultat sauvegardé : %s\n\n", file.path(work_dir, "res_negative.rds")))
}

# =============================================================================
# PLOT A — Heatmap globale (pathways × facteurs)
# analysis == "gsea_heatmap"
# =============================================================================

if (analysis == "gsea_heatmap") {
  cat("[PLOT A] Heatmap globale des enrichissements...\n")
  
  # Charger le résultat approprié
  if (isTRUE(load_result == "TRUE") && !is.null(result_type)) {
    load_enrichment_results()
    enrichment_result <- if (result_type == "positive") enrichment_positive else enrichment_negative
  }

  draw_heatmap <- function() {
    tryCatch(
      plot_enrichment_heatmap(enrichment_result),
      error = function(e) cat(sprintf("  [AVERTISSEMENT] %s\n", e$message))
    )
  }

  if (!preview_only) {
    pdf(outfile, width = 12, height = 9)
    draw_heatmap()
    dev.off()
    cat(sprintf("  ✓ Heatmap GSEA PDF : %s\n", outfile))
  }
  png_path <- sub("\\.pdf$", "_preview.png", outfile)
  grDevices::png(png_path, width = 1200, height = 900, res = 100)
  draw_heatmap()
  dev.off()
  cat(sprintf("  ✓ Preview PNG : %s\n\n", png_path))
}

# =============================================================================
# PLOT B — Top pathways pour un facteur donné
# analysis == "gsea_top_pathways"
# =============================================================================

if (analysis == "gsea_top_pathways") {
  # Charger le résultat approprié
  if (isTRUE(load_result == "TRUE") && !is.null(result_type)) {
    load_enrichment_results()
    enrichment_result <- if (result_type == "positive") enrichment_positive else enrichment_negative
  }

  if (!factor_plot %in% factors_parsed) {
    cat(sprintf("  [INFO] Facteur %d non présent dans l'analyse, ajustement au facteur %d.\n", factor_plot, factors_parsed[1]))
    factor_plot <- factors_parsed[1]
  }
  cat(sprintf("[PLOT B] Top %d pathways — Facteur %d...\n", max_pathways, factor_plot))

  p <- tryCatch(
    plot_enrichment(enrichment_result, factor = factor_plot, max.pathways = max_pathways),
    error = function(e) { cat(sprintf("  [ERREUR] %s\n", e$message)); NULL }
  )

  if (!is.null(p)) {
    save_gg(p, outfile, width = 12, height = max(6, max_pathways * 0.45 + 2))
  }
}

# =============================================================================
# PLOT C — Détail des gènes contributeurs
# analysis == "gsea_detailed"
# =============================================================================

if (analysis == "gsea_detailed") {
  # Charger le résultat approprié
  if (isTRUE(load_result == "TRUE") && !is.null(result_type)) {
    load_enrichment_results()
    enrichment_result <- if (result_type == "positive") enrichment_positive else enrichment_negative
  }

  if (!factor_plot %in% factors_parsed) {
    cat(sprintf("  [INFO] Facteur %d non présent dans l'analyse, ajustement au facteur %d.\n", factor_plot, factors_parsed[1]))
    factor_plot <- factors_parsed[1]
  }
  cat(sprintf("[PLOT C] Détail gènes — Facteur %d (%d pathways, %d gènes)...\n",
              factor_plot, max_pathways, max_genes))

  p <- tryCatch(
    plot_enrichment_detailed(
      enrichment_result,
      factor       = factor_plot,
      max.genes    = max_genes,
      max.pathways = max_pathways
    ),
    error = function(e) { cat(sprintf("  [ERREUR] %s\n", e$message)); NULL }
  )

  if (!is.null(p)) {
    save_gg(p, outfile, width = 13, height = max(7, max_pathways * 0.6 + 2))
  }
}

# =============================================================================
# PLOT D — Comparaison deux tests statistiques
# analysis == "gsea_compare_tests"
# =============================================================================

if (analysis == "gsea_compare_tests") {

  if (is.null(stat_test2) || stat_test2 == "") {
    cat("  [ERREUR] Un second test (stat_test2) est requis pour la comparaison.\n")
    quit(status = 1)
  }

  cat(sprintf("[PLOT D] Comparaison tests : %s vs %s...\n", stat_test, stat_test2))

  enrichment2 <- tryCatch({
    run_enrichment(
      mofa_model,
      view             = view_gsea,
      factors          = factors_parsed,
      feature.sets     = gs_matrix,
      sign             = sign_gsea,
      statistical.test = stat_test2
    )
  }, error = function(e) {
    cat(sprintf("  [ERREUR run_enrichment test2] %s\n", e$message))
    quit(status = 1)
  })

  cat(sprintf("  ✓ Second enrichissement calculé (%s).\n\n", stat_test2))

  # ── Histogrammes des p-values
  cat("  → Histogrammes p-values...\n")
  n_factors_plot <- min(length(factors_parsed), 3)
  fi_sel <- factors_parsed[seq_len(n_factors_plot)]
  fi_names <- paste0("Factor", fi_sel)

  dt1 <- as.data.table(enrichment_result$pval[, fi_names, drop = FALSE])
  dt1[, test := stat_test][, pathway := seq_len(.N)]

  dt2 <- as.data.table(enrichment2$pval[, fi_names, drop = FALSE])
  dt2[, test := stat_test2][, pathway := seq_len(.N)]

  dt_hist <- rbind(dt1, dt2) %>%
    melt(id.vars = c("test", "pathway"), variable.name = "factor")

  p_hist <- ggplot(dt_hist, aes(x = value, fill = test)) +
    facet_wrap(~factor, scales = "free_y", nrow = 1) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c("#619CFF", "#FF6B6B"), name = "Test") +
    labs(
      title    = "Distribution des p-values par test statistique",
      subtitle = paste("Signe :", sign_gsea, "| Vue :", view_gsea),
      x        = "p-value",
      y        = "Fréquence"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", legend.title = element_blank())

  # ── Scatter corrélation entre les deux tests
  cat("  → Scatter corrélation p-values...\n")
  dt_s1 <- as.data.table(enrichment_result$pval[, fi_names, drop = FALSE])
  dt_s1[, pathway := seq_len(.N)]
  dt_s2 <- as.data.table(enrichment2$pval[, fi_names, drop = FALSE])
  dt_s2[, pathway := seq_len(.N)]

  dt_scatter <- merge(
    melt(dt_s1, id.vars = "pathway", variable.name = "factor", value.name = stat_test),
    melt(dt_s2, id.vars = "pathway", variable.name = "factor", value.name = stat_test2),
    by = c("pathway", "factor")
  )

  p_scatter <- ggplot(dt_scatter,
      aes_string(x = stat_test, y = stat_test2)) +
    geom_point(size = 0.6, alpha = 0.5, color = "#555555") +
    geom_abline(slope = 1, intercept = 0, color = "orange", linewidth = 0.8) +
    facet_wrap(~factor, scales = "free", nrow = 1) +
    labs(
      title    = paste("Corrélation des p-values :", stat_test, "vs", stat_test2),
      x        = paste("p-value —", stat_test),
      y        = paste("p-value —", stat_test2)
    ) +
    theme_bw(base_size = 12)

  # ── Assemblage
  pdf(outfile, width = 14, height = 10)
  if (has_cowplot) {
    print(cowplot::plot_grid(p_hist, p_scatter, ncol = 1, labels = c("A", "B")))
  } else {
    print(p_hist)
    print(p_scatter)
  }
  dev.off()
  cat(sprintf("  ✓ Comparaison enregistrée : %s\n\n", outfile))
}

# =============================================================================
# PLOT E — Tous les plots en un seul PDF multi-pages
# analysis == "gsea_full_report"
# =============================================================================

if (analysis == "gsea_full_report") {
  cat("[PLOT E] Rapport complet GSEA (PDF multi-pages)...\n")

  pdf(outfile, width = 12, height = 9)

  # Page 1 : heatmap globale
  cat("  → Heatmap globale...\n")
  tryCatch(
    plot_enrichment_heatmap(enrichment_result),
    error = function(e) cat(sprintf("  [AVERTISSEMENT heatmap] %s\n", e$message))
  )

  # Pages suivantes : top pathways + détail pour chaque facteur
  for (fi in factors_parsed) {
    cat(sprintf("  → Facteur %d : top pathways...\n", fi))
    p_top <- tryCatch(
      plot_enrichment(enrichment_result, factor = fi, max.pathways = max_pathways),
      error = function(e) NULL
    )
    if (!is.null(p_top)) print(p_top)

    cat(sprintf("  → Facteur %d : détail gènes...\n", fi))
    p_det <- tryCatch(
      plot_enrichment_detailed(
        enrichment_result,
        factor       = fi,
        max.genes    = max_genes,
        max.pathways = min(5, max_pathways)
      ),
      error = function(e) NULL
    )
    if (!is.null(p_det)) print(p_det)
  }

  dev.off()
  cat(sprintf("  ✓ Rapport complet enregistré : %s\n\n", outfile))
}

# =============================================================================
# FIN
# =============================================================================

cat("========================================================\n")
cat(sprintf("  GSEA TERMINÉ : %s\n", out_name))
cat("  Fichiers générés :\n")
cat(sprintf("  • %s.pdf\n", out_name))
cat(sprintf("  • %s_pvalues.csv\n", out_name))
cat(sprintf("  • %s_pvalues_adj.csv\n", out_name))
cat(sprintf("  • %s_set_statistics.csv\n", out_name))
cat("========================================================\n")
