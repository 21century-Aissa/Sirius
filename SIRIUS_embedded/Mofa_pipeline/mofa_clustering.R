#!/usr/bin/env Rscript
# =============================================================================
# MOFA CLUSTERING — Clustering non supervisé sur les facteurs latents MOFA
# =============================================================================
# Méthodes : kmeans / pam / hclust_ward
# Sortie  : clusters.rds, clusters.csv, choose_k.pdf, cluster_scatter.pdf
# =============================================================================

suppressPackageStartupMessages(library(MOFA2))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(NbClust))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

# ── Source des utilitaires partagés (lib/) ───────────────────────────────────
.script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
source(file.path(.script_dir, "lib", "mofa_utils.R"))
source(file.path(.script_dir, "lib", "mofa_io.R"))

set.seed(42)

# =============================================================================
# Lecture des arguments
# =============================================================================
args   <- commandArgs(trailingOnly = TRUE)
params <- parse_args(args)

model_path     <- params[["model"]]
output_dir     <- params[["output_dir"]]
analysis       <- params[["analysis"]]
method         <- if (!is.null(params[["method"]])) params[["method"]] else "kmeans"
k              <- if (!is.null(params[["k"]])) as.integer(params[["k"]]) else 3L
var_threshold  <- if (!is.null(params[["var_threshold"]])) as.numeric(params[["var_threshold"]]) else 0
factor_x       <- if (!is.null(params[["factor_x"]])) as.integer(params[["factor_x"]]) else 1L
factor_y       <- if (!is.null(params[["factor_y"]])) as.integer(params[["factor_y"]]) else 2L

if (is.null(model_path) || is.null(output_dir) || is.null(analysis)) {
  stop("Arguments requis : --model, --output_dir, --analysis")
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

mofa_header("MOFA CLUSTERING", analysis, list(
  "Méthode"  = method,
  "k"        = k,
  "Filtrage" = var_threshold
))

# =============================================================================
# Chargement du modèle et extraction des facteurs filtrés
# =============================================================================
cat("[ÉTAPE 1] Chargement du modèle et extraction des facteurs...\n")
mofa_model <- load_model_auto(model_path)

factors_all <- get_factors(mofa_model)[[1]]
r2_per_factor <- get_variance_explained(mofa_model)$r2_per_factor[[1]]

# r2_per_factor : matrice (views x factors) — somme/maximum sur les vues
if (is.matrix(r2_per_factor)) {
  factor_var <- colSums(r2_per_factor, na.rm = TRUE)
} else {
  factor_var <- as.numeric(r2_per_factor)
}

keep <- factor_var >= var_threshold
cat(sprintf("  ✓ %d/%d facteurs retenus (seuil = %.4f)\n",
            sum(keep), length(keep), var_threshold))

if (!any(keep)) {
  stop(sprintf(
    "Aucun facteur ne dépasse le seuil de variance expliquée (%.4f).",
    var_threshold))
}

factors <- factors_all[, keep, drop = FALSE]
cat(sprintf("  ✓ Matrice finale : %d samples × %d facteurs\n\n",
            nrow(factors), ncol(factors)))

# =============================================================================
# Function 1 — choose_k
# =============================================================================
choose_k <- function(factors, output_dir) {
  cat("[choose_k] Calcul WSS / Silhouette / Gap...\n")
  set.seed(42)
  p_wss <- fviz_nbclust(factors, kmeans, method = "wss") +
    ggtitle("Within Sum of Squares")
  p_sil <- fviz_nbclust(factors, kmeans, method = "silhouette") +
    ggtitle("Silhouette")
  p_gap <- fviz_nbclust(factors, kmeans, method = "gap_stat",
                        nboot = 50, verbose = FALSE) +
    ggtitle("Gap statistic")

  outfile <- file.path(output_dir, "choose_k.pdf")
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    arr <- gridExtra::arrangeGrob(p_wss, p_sil, p_gap, ncol = 3)
    save_gg(arr, outfile, width = 14, height = 5)
  } else {
    save_gg(p_wss, outfile, width = 14, height = 5)
  }
}

# =============================================================================
# Function 2 — run_clustering
# =============================================================================
run_clustering <- function(factors, method, k, factor_x, factor_y, output_dir) {
  cat(sprintf("[run_clustering] %s, k=%d\n", method, k))
  set.seed(42)

  clusters <- if (method == "kmeans") {
    kmeans(factors, centers = k, nstart = 25)$cluster
  } else if (method == "pam") {
    cluster::pam(factors, k = k)$clustering
  } else if (method == "hclust_ward") {
    cutree(hclust(dist(factors), method = "ward.D2"), k = k)
  } else {
    stop(sprintf("Méthode inconnue : %s", method))
  }

  names(clusters) <- rownames(factors)

  save_clusters(clusters, output_dir)

  fx_col <- paste0("Factor", factor_x)
  fy_col <- paste0("Factor", factor_y)
  available <- colnames(factors)
  if (!(fx_col %in% available)) fx_col <- available[min(factor_x, length(available))]
  if (!(fy_col %in% available)) fy_col <- available[min(factor_y, length(available))]

  df <- data.frame(
    x       = factors[, fx_col],
    y       = factors[, fy_col],
    cluster = factor(clusters)
  )

  p <- ggplot(df, aes(x = x, y = y, color = cluster)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(
      title = sprintf("MOFA clusters (%s, k=%d)", method, k),
      x     = fx_col,
      y     = fy_col,
      color = "Cluster"
    ) +
    theme_minimal()

  outfile <- file.path(output_dir, "cluster_scatter.pdf")
  save_gg(p, outfile, width = 7, height = 6)
}

# =============================================================================
# Dispatch
# =============================================================================
if (analysis == "choose_k") {
  choose_k(factors, output_dir)
} else if (analysis == "run_clustering") {
  run_clustering(factors, method, k, factor_x, factor_y, output_dir)
} else {
  stop(sprintf("Analyse inconnue : %s (attendu : choose_k | run_clustering)",
               analysis))
}

cat("========================================================\n")
cat(sprintf("  CLUSTERING TERMINÉ : %s\n", analysis))
cat("========================================================\n")
