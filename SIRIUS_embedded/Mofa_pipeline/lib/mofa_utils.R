# =============================================================================
# lib/mofa_utils.R
# -----------------------------------------------------------------------------
# Utilitaires partagés par les scripts du pipeline MOFA :
#   - parse_args()      : parsing canonique des arguments CLI
#   - load_model_auto() : chargement modèle MOFA (RDS ou HDF5)
#   - save_gg()         : sauvegarde PDF + preview PNG d'un ggplot
#   - mofa_header()     : impression d'un en-tête de section
#   - check_required()  : validation des arguments obligatoires
# =============================================================================

# Canonical CLI parser.
# Format principal : "--key value".
#   - Si le token suivant commence par "--", la clé courante est traitée
#     comme un flag booléen (TRUE).
#   - Pour rétrocompatibilité, les paires "key=value" (sans --) sont aussi
#     reconnues, ce qui évite de casser les scripts hérités.
parse_args <- function(args) {
  params <- list()
  i <- 1
  n <- length(args)
  while (i <= n) {
    a <- args[[i]]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (i + 1 <= n && !startsWith(args[[i + 1]], "--")) {
        params[[key]] <- args[[i + 1]]
        i <- i + 2
      } else {
        params[[key]] <- TRUE
        i <- i + 1
      }
    } else if (grepl("=", a, fixed = TRUE)) {
      parts <- strsplit(a, "=", fixed = TRUE)[[1]]
      key   <- parts[1]
      value <- paste(parts[-1], collapse = "=")
      params[[key]] <- value
      i <- i + 1
    } else {
      i <- i + 1
    }
  }
  params
}

# Charge un modèle MOFA depuis un fichier .rds ou .hdf5
load_model_auto <- function(model_path) {
  ext <- tolower(tools::file_ext(model_path))
  if (ext == "rds") {
    m <- readRDS(model_path)
    cat("  ✓ Modèle chargé depuis RDS\n")
  } else {
    m <- MOFA2::load_model(model_path)
    cat("  ✓ Modèle chargé depuis HDF5\n")
  }
  m
}

# Enregistre un ggplot en PDF et écrit toujours un preview PNG associé.
save_gg <- function(p, filepath, width = 10, height = 7, preview_only = FALSE) {
  if (!isTRUE(preview_only)) {
    pdf(filepath, width = width, height = height)
    print(p)
    dev.off()
    cat(sprintf("  ✓ PDF enregistré : %s\n", filepath))
  }
  png_path <- sub("\\.pdf$", "_preview.png", filepath)
  grDevices::png(png_path,
                 width  = round(width  * 100),
                 height = round(height * 100),
                 res    = 100)
  print(p)
  dev.off()
  cat(sprintf("  ✓ Preview PNG : %s\n", png_path))
}

# En-tête imprimé en début de chaque script.
# extra_lines : named list label → value (chaque ligne devient "  label : value")
mofa_header <- function(script_name, analysis, extra_lines = list()) {
  cat("========================================================\n")
  cat(sprintf("  %s — %s\n", script_name, analysis))
  if (length(extra_lines) > 0) {
    labels <- names(extra_lines)
    pad    <- max(nchar(labels))
    for (k in labels) {
      val <- extra_lines[[k]]
      cat(sprintf("  %-*s : %s\n", pad, k, format(val)))
    }
  }
  cat("========================================================\n\n")
}

# Vérifie la présence d'une liste d'arguments requis dans `params`.
check_required <- function(params, keys) {
  for (k in keys) {
    if (is.null(params[[k]])) {
      stop(sprintf("Argument requis manquant : --%s", k), call. = FALSE)
    }
  }
  invisible(TRUE)
}
