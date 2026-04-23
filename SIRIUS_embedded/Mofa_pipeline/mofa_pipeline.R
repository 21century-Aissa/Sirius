#!/usr/bin/env Rscript
# =============================================================================
# MOFA PIPELINE — Multi-Omics Factor Analysis
# =============================================================================
# Script auto-généré par l'interface MOFA Desktop App
# Chaque section est clairement délimitée et nommée
# =============================================================================

suppressMessages({
  library(MOFA2)
  library(ggplot2)
  library(data.table)
})

# =============================================================================
# SECTION 0 — Lecture des arguments passés par l'interface Python
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

# Les arguments sont passés sous forme de paires clé=valeur
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

params <- parse_args(args)

work_dir      <- params[["work_dir"]]
num_factors   <- as.integer(params[["num_factors"]])
max_r2        <- if (!is.null(params[["max_r2"]]) && params[["max_r2"]] != "") as.numeric(params[["max_r2"]]) else NULL
has_metadata  <- as.logical(params[["has_metadata"]])
metadata_file <- params[["metadata_file"]]
python_path   <- params[["python_path"]]

# Fichiers omiques : passés comme omic_NAME=path
omic_files <- list()
for (key in names(params)) {
  if (startsWith(key, "omic_")) {
    omic_name <- sub("^omic_", "", key)
    omic_files[[omic_name]] <- params[[key]]
  }
}

cat("========================================================\n")
cat("  MOFA PIPELINE — Démarrage\n")
cat("========================================================\n")
cat(sprintf("  Dossier de travail : %s\n", work_dir))
cat(sprintf("  Nombre de facteurs : %d\n", num_factors))
cat(sprintf("  Max R² : %s\n", if (!is.null(max_r2)) as.character(max_r2) else "NULL (default)"))
cat(sprintf("  Nombre de vues omiques : %d\n", length(omic_files)))
cat("========================================================\n\n")

setwd(work_dir)

# =============================================================================
# SECTION 1 — Chargement et validation des données omiques
# =============================================================================

cat("[ÉTAPE 1] Chargement des matrices omiques...\n")

load_matrix <- function(filepath) {
  ext <- tolower(tools::file_ext(filepath))
  mat <- if (ext %in% c("csv")) {
    read.csv(filepath, row.names = 1, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    read.delim(filepath, row.names = 1, check.names = FALSE)
  } else if (ext %in% c("rds")) {
    readRDS(filepath)
  } else {
    stop(paste("Format non supporté :", ext))
  }
  return(as.matrix(mat))
}

omics.list <- list()
for (omic_name in names(omic_files)) {
  filepath <- omic_files[[omic_name]]
  cat(sprintf("  → Chargement : %s depuis %s\n", omic_name, filepath))
  tryCatch({
    mat <- load_matrix(filepath)
    # Vérification : features en lignes, samples en colonnes
    cat(sprintf("     Dimensions : %d features × %d samples\n", nrow(mat), ncol(mat)))
    omics.list[[omic_name]] <- mat
  }, error = function(e) {
    cat(sprintf("  [ERREUR] Impossible de charger %s : %s\n", omic_name, e$message))
    quit(status = 1)
  })
}

cat(sprintf("\n  ✓ %d vues omiques chargées avec succès.\n\n", length(omics.list)))

# =============================================================================
# SECTION 2 — Création de l'objet MOFA
# =============================================================================

cat("[ÉTAPE 2] Création de l'objet MOFA...\n")

MOFAobject <- create_mofa(omics.list)
cat("  ✓ Objet MOFA créé.\n\n")

# =============================================================================
# SECTION 3 — Vue d'ensemble des données (plot_data_overview)
# =============================================================================

cat("[ÉTAPE 3] Génération du plot d'aperçu des données...\n")

data_overview_path <- file.path(work_dir, "data_overview.pdf")
pdf(data_overview_path, width = 10, height = 6)
  plot_data_overview(MOFAobject)
dev.off()

cat(sprintf("  ✓ Aperçu enregistré : %s\n\n", data_overview_path))

# =============================================================================
# SECTION 4 — Définition des options des données
# =============================================================================

cat("[ÉTAPE 4] Définition des options des données...\n")

data_opts <- get_default_data_options(MOFAobject)
cat("  ✓ Options des données définies (valeurs par défaut).\n\n")

# =============================================================================
# SECTION 5 — Définition des options du modèle
# =============================================================================

cat("[ÉTAPE 5] Définition des options du modèle...\n")

model_opts                <- get_default_model_options(MOFAobject)
model_opts$num_factors    <- num_factors

cat(sprintf("  → Nombre de facteurs : %d\n", model_opts$num_factors))
cat("  ✓ Options du modèle définies.\n\n")

# =============================================================================
# SECTION 6 — Définition des options d'entraînement
# =============================================================================

cat("[ÉTAPE 6] Définition des options d'entraînement...\n")

train_opts                    <- get_default_training_options(MOFAobject)
train_opts$convergence_mode   <- "slow"
train_opts$seed               <- 42

cat("  → Mode de convergence : slow\n")
cat("  → Seed : 42\n")
cat("  ✓ Options d'entraînement définies.\n\n")

# =============================================================================
# SECTION 7 — Préparation et entraînement du modèle MOFA
# =============================================================================

cat("[ÉTAPE 7] Préparation du modèle MOFA...\n")

MOFAobject <- prepare_mofa(
  MOFAobject,
  data_options     = data_opts,
  model_options    = model_opts,
  training_options = train_opts
)

cat("  ✓ Modèle préparé.\n\n")
cat("[ÉTAPE 7b] Entraînement du modèle MOFA (peut prendre plusieurs minutes)...\n")

outfile <- normalizePath(file.path(work_dir, "model.hdf5"), winslash = "/", mustWork = FALSE)

cat(sprintf("  → Fichier de sortie : %s\n", outfile))

if (!is.null(python_path) && nchar(python_path) > 0) {
  reticulate::use_python(python_path, required = TRUE)
}
if (!reticulate::py_module_available("mofapy2")) {
  cat("  → mofapy2 absent, installation via pip...\n")
  reticulate::py_install("mofapy2", pip = TRUE)
}

MOFAobject <- run_mofa(MOFAobject, outfile, use_basilisk = FALSE)

cat(sprintf("  ✓ Modèle entraîné et enregistré : %s\n\n", outfile))

# =============================================================================
# SECTION 8 — Rechargement du modèle entraîné
# =============================================================================

cat("[ÉTAPE 8] Rechargement du modèle entraîné...\n")

mofa_model <- load_model(outfile)

cat("  ✓ Modèle rechargé.\n\n")

# =============================================================================
# SECTION 9 — Ajout des métadonnées (si fournies)
# =============================================================================

cat("[ÉTAPE 9] Gestion des métadonnées...\n")

if (has_metadata && !is.null(metadata_file) && file.exists(metadata_file)) {
  cat(sprintf("  → Chargement des métadonnées : %s\n", metadata_file))
  ext      <- tolower(tools::file_ext(metadata_file))
  metadata <- if (ext == "rds") {
    readRDS(metadata_file)
  } else if (ext == "csv") {
    read.csv(metadata_file, check.names = FALSE)
  } else {
    read.delim(metadata_file, check.names = FALSE)
  }
  samples_metadata(mofa_model) <- metadata
  cat("  ✓ Metadata attached to model.\n")
  # Save model with metadata as RDS so plots can reload it with metadata intact
  rds_path <- file.path(work_dir, "model_with_metadata.rds")
  saveRDS(mofa_model, rds_path)
  cat(sprintf("  ✓ Model with metadata saved: %s\n\n", rds_path))
} else {
  cat("  → No metadata provided, step skipped.\n\n")
}

# =============================================================================
# SECTION 10 — Corrélation entre facteurs (plot_factor_cor)
# =============================================================================

cat("[ÉTAPE 10] Génération du plot de corrélation entre facteurs...\n")

factor_cor_path <- file.path(work_dir, "factor_correlation_plot.pdf")
pdf(factor_cor_path, width = 8, height = 8)
  plot_factor_cor(mofa_model)
dev.off()

cat(sprintf("  ✓ Corrélation des facteurs enregistrée : %s\n\n", factor_cor_path))

# =============================================================================
# SECTION 11 — Variance expliquée par facteur (plot_variance_explained)
# =============================================================================

cat("[ÉTAPE 11] Génération du plot de décomposition de la variance...\n")

variance_by_factor_path <- file.path(work_dir, "variance_decomposition_by_factor.pdf")
pdf(variance_by_factor_path, width = 10, height = 6)
  if (!is.null(max_r2)) plot_variance_explained(mofa_model, max_r2 = max_r2) else plot_variance_explained(mofa_model)
dev.off()

cat(sprintf("  ✓ Variance par facteur enregistrée : %s\n\n", variance_by_factor_path))

# =============================================================================
# SECTION 12 — Variance totale expliquée par vue (plot_total variance)
# =============================================================================

cat("[ÉTAPE 12] Génération du plot de variance totale expliquée par vue...\n")

total_variance_path <- file.path(work_dir, "Total_variance_explained_per_view.pdf")
pdf(total_variance_path, width = 8, height = 6)
  print(plot_variance_explained(mofa_model, plot_total = TRUE)[[2]])
dev.off()

cat(sprintf("  ✓ Variance totale par vue enregistrée : %s\n\n", total_variance_path))

# =============================================================================
# SECTION FINALE — Résumé des sorties
# =============================================================================

cat("========================================================\n")
cat("  PIPELINE MOFA TERMINÉ AVEC SUCCÈS\n")
cat("========================================================\n")
cat("  Fichiers générés :\n")
cat(sprintf("  • data_overview.pdf\n"))
cat(sprintf("  • model.hdf5\n"))
cat(sprintf("  • factor_correlation_plot.pdf\n"))
cat(sprintf("  • variance_decomposition_by_factor.pdf\n"))
cat(sprintf("  • Total_variance_explained_per_view.pdf\n"))
cat("========================================================\n")
