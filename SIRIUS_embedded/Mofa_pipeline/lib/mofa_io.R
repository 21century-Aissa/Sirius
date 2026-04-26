# =============================================================================
# lib/mofa_io.R
# -----------------------------------------------------------------------------
# Entrées/sorties partagées : sauvegarde et chargement des assignations de
# clusters produites par mofa_clustering.R.
# =============================================================================

# Sauvegarde un vecteur de clusters (nommé par sample) en RDS + CSV.
save_clusters <- function(clusters, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(clusters, file.path(output_dir, "clusters.rds"))
  write.csv(
    data.frame(sample = names(clusters), cluster = clusters),
    file.path(output_dir, "clusters.csv"),
    row.names = FALSE
  )
  cat(sprintf("  ✓ clusters.rds / clusters.csv enregistrés dans %s\n", output_dir))
}

# Recharge le vecteur de clusters depuis output_dir.
# Échoue avec un message clair si clusters.rds n'existe pas.
load_clusters <- function(output_dir) {
  rds <- file.path(output_dir, "clusters.rds")
  if (!file.exists(rds)) {
    stop(sprintf(
      "clusters.rds introuvable dans %s. Lancez d'abord run_clustering.",
      output_dir
    ), call. = FALSE)
  }
  readRDS(rds)
}
