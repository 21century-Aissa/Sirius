#!/usr/bin/env Rscript
# =============================================================================
# MOFA PLOTS — Visualisation avancée du modèle entraîné
# =============================================================================
# Ce script est appelé par l'interface Python avec des arguments ciblés.
# Chaque type de plot est une fonction indépendante déclenchée par plot_type=...
# =============================================================================

suppressMessages({
  library(MOFA2)
  library(ggplot2)
  library(pheatmap)
})

# =============================================================================
# UTILITAIRES — Parsing des arguments & helpers
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

# Convertit "1,2,3" → c(1,2,3) entiers
parse_int_vector <- function(s) {
  as.integer(trimws(strsplit(s, ",")[[1]]))
}

# Convertit "Feature1,Feature2" → c("Feature1","Feature2")
parse_str_vector <- function(s) {
  trimws(strsplit(s, ",")[[1]])
}

# Saves PNG preview (always) + PDF (only when preview_only is FALSE)
save_plot <- function(p, filepath, width = 10, height = 7) {
  if (!preview_only) {
    pdf(filepath, width = width, height = height)
    print(p)
    dev.off()
    cat(sprintf("  ✓ PDF saved: %s\n", filepath))
  }
  png_path <- sub("\\.pdf$", "_preview.png", filepath)
  grDevices::png(png_path, width = round(width * 100), height = round(height * 100), res = 100)
  print(p)
  dev.off()
  cat(sprintf("  ✓ Preview saved: %s\n", png_path))
}

# =============================================================================
# LECTURE DES ARGUMENTS
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
params  <- parse_args(args)

model_path   <- params[["model_path"]]
work_dir     <- params[["work_dir"]]
plot_type    <- params[["plot_type"]]
out_name     <- params[["out_name"]]   # nom du fichier de sortie sans extension
preview_only <- isTRUE(params[["preview_only"]] == "TRUE")  # PNG only, no PDF

setwd(work_dir)

cat("========================================================\n")
cat(sprintf("  MOFA PLOTS — %s\n", plot_type))
cat(sprintf("  Modèle : %s\n", model_path))
cat(sprintf("  Sortie : %s/%s.pdf\n", work_dir, out_name))
cat("========================================================\n\n")

# =============================================================================
# CHARGEMENT DU MODÈLE
# =============================================================================

cat("[INIT] Loading MOFA model...\n")
# Prefer RDS (includes metadata) over raw HDF5 if available
rds_path <- file.path(work_dir, "model_with_metadata.rds")
if (file.exists(rds_path)) {
  mofa_model <- readRDS(rds_path)
  cat("  Model loaded from RDS (metadata included).\n\n")
} else {
  mofa_model <- load_model(model_path)
  cat("  Model loaded from HDF5 (no metadata).\n\n")
}

outfile <- file.path(work_dir, paste0(out_name, ".pdf"))

# =============================================================================
# PLOT 1 — Factor Values  (plot_factor)
# Paramètres : by factors
# =============================================================================

if (plot_type == "factor_values") {
  cat("[PLOT] Factor Values...\n")

  factors    <- parse_int_vector(params[["factors"]])
  color_by   <- paste0("Factor", factors[1])

  p <- plot_factor(
    mofa_model,
    factors    = factors,
    color_by   = color_by
  )
  save_plot(p, outfile)
}

# =============================================================================
# PLOT 2 — Top Weights  (plot_top_weights)
# Paramètres : view, factor, nfeatures, scale
# =============================================================================

if (plot_type == "top_weights") {
  cat("[PLOT] Top Weights...\n")

  view_raw  <- params[["view"]]
  factor    <- as.integer(params[["factor"]])
  nfeatures <- as.integer(params[["nfeatures"]])
  scale     <- as.logical(params[["scale"]])

  # view accepte :
  #   - un nom precis            : "mRNA"
  #   - plusieurs noms (CSV)     : "mRNA,Methylation"
  #   - toutes les vues          : "ALL" (ou champ vide)
  all_views <- views_names(mofa_model)
  if (is.null(view_raw) || !nzchar(trimws(view_raw)) ||
      toupper(trimws(view_raw)) == "ALL") {
    view <- all_views
    cat(sprintf("  View = ALL (%s)\n", paste(view, collapse = ", ")))
  } else {
    view    <- parse_str_vector(view_raw)
    unknown <- setdiff(view, all_views)
    if (length(unknown) > 0) {
      cat(sprintf("  [ERREUR] Vue(s) inconnue(s) : %s. Disponibles : %s\n",
                  paste(unknown, collapse = ", "),
                  paste(all_views, collapse = ", ")))
      quit(status = 1)
    }
    cat(sprintf("  View = %s\n", paste(view, collapse = ", ")))
  }

  p <- plot_top_weights(
    mofa_model,
    view      = view,
    factor    = factor,
    nfeatures = nfeatures,
    scale     = scale
  )
  save_plot(p, outfile)
}

# =============================================================================
# PLOT 3 — Factor avec violin  (plot_factor + violin)
# Paramètres : factors, color_by, add_violin=TRUE, dodge=TRUE
# =============================================================================

if (plot_type == "factor_violin") {
  cat("[PLOT] Factor Violin...\n")

  factors  <- parse_int_vector(params[["factors"]])
  color_by <- params[["color_by"]]

  p <- plot_factor(
    mofa_model,
    factors    = factors,
    color_by   = color_by,
    add_violin = TRUE,
    dodge      = TRUE
  )
  save_plot(p, outfile)
}

# =============================================================================
# PLOT 4 — Feature Weights  (plot_weights)
# Paramètres : view, factor, nfeatures (scale=TRUE)
# =============================================================================

if (plot_type == "feature_weights") {
  cat("[PLOT] Feature Weights...\n")

  view      <- params[["view"]]
  factor    <- as.integer(params[["factor"]])
  nfeatures <- as.integer(params[["nfeatures"]])

  p <- plot_weights(
    mofa_model,
    view      = view,
    factor    = factor,
    nfeatures = nfeatures,
    scale     = TRUE
  )
  save_plot(p, outfile)
}

# =============================================================================
# PLOT 5 — Data Scatter  (plot_data_scatter)
# Paramètres : view, factor, features (virgule-séparés), sign, color_by, y_label
# =============================================================================

if (plot_type == "data_scatter") {
  cat("[PLOT] Data Scatter...\n")

  view     <- params[["view"]]
  factor   <- as.integer(params[["factor"]])
  features <- parse_str_vector(params[["features"]])
  sign     <- params[["sign"]]        # "positive" ou "negative"
  color_by <- params[["color_by"]]
  y_label  <- if (!is.null(params[["y_label"]]) && params[["y_label"]] != "")
                params[["y_label"]] else view

  p <- plot_data_scatter(
    mofa_model,
    view     = view,
    factor   = factor,
    features = features,
    sign     = sign,
    color_by = color_by
  ) + labs(y = y_label)

  save_plot(p, outfile, width = 12, height = 8)
}

# =============================================================================
# PLOT 6 — Data Heatmap  (plot_data_heatmap)
# Paramètres : view, factor, features (int), cluster_rows, cluster_cols,
#              show_rownames, show_colnames, scale
# =============================================================================

if (plot_type == "data_heatmap") {
  cat("[PLOT] Data Heatmap...\n")

  view          <- params[["view"]]
  factor        <- as.integer(params[["factor"]])
  features      <- as.integer(params[["features"]])
  cluster_rows  <- as.logical(params[["cluster_rows"]])
  cluster_cols  <- as.logical(params[["cluster_cols"]])
  show_rownames <- as.logical(params[["show_rownames"]])
  show_colnames <- as.logical(params[["show_colnames"]])
  scale_mode    <- params[["scale_mode"]]  # "row", "column", "none"

  draw_heatmap <- function() {
    plot_data_heatmap(
      mofa_model,
      view          = view,
      factor        = factor,
      features      = features,
      cluster_rows  = cluster_rows,
      cluster_cols  = cluster_cols,
      show_rownames = show_rownames,
      show_colnames = show_colnames,
      scale         = scale_mode
    )
  }

  if (!preview_only) {
    pdf(outfile, width = 12, height = 9)
    draw_heatmap()
    dev.off()
    cat(sprintf("  \u2713 PDF saved: %s\n", outfile))
  }
  png_path <- sub("[.]pdf$", "_preview.png", outfile)
  grDevices::png(png_path, width = 1200, height = 900, res = 100)
  draw_heatmap()
  dev.off()
  cat(sprintf("  \u2713 Preview saved: %s\n", png_path))
}

# =============================================================================
# PLOT 7 — Factors Scatter (2 facteurs en XY)  (plot_factors)
# Paramètres : factor_x, factor_y, color_by, shape_by
# =============================================================================

if (plot_type == "factors_scatter") {
  cat("[PLOT] Factors Scatter (2D)...\n")

  factor_x <- as.integer(params[["factor_x"]])
  factor_y <- as.integer(params[["factor_y"]])
  color_by <- params[["color_by"]]
  shape_by <- if (!is.null(params[["shape_by"]]) && params[["shape_by"]] != "")
                params[["shape_by"]] else NULL

  p <- plot_factors(
    mofa_model,
    factors  = c(factor_x, factor_y),
    color_by = color_by,
    shape_by = shape_by
  )
  save_plot(p, outfile)
}

# =============================================================================
# PLOT 8 — Variance Explained par facteur  (plot_variance_explained)
# Paramètres : max_r2
# =============================================================================

if (plot_type == "variance_by_factor") {
  cat("[PLOT] Variance Explained by Factor...\n")

  min_r2  <- if (!is.null(params[["min_r2"]]) && params[["min_r2"]] != "") as.numeric(params[["min_r2"]]) else 0
  max_r2  <- if (!is.null(params[["max_r2"]]) && params[["max_r2"]] != "") as.numeric(params[["max_r2"]]) else NULL
  x_axis  <- if (!is.null(params[["x"]]) && params[["x"]] != "") params[["x"]] else "view"
  y_axis  <- if (!is.null(params[["y"]]) && params[["y"]] != "") params[["y"]] else "factor"
  palette <- if (!is.null(params[["palette"]]) && params[["palette"]] != "") params[["palette"]] else NULL

  p <- if (!is.null(max_r2)) {
    plot_variance_explained(mofa_model, x = x_axis, y = y_axis, min_r2 = min_r2, max_r2 = max_r2)
  } else {
    plot_variance_explained(mofa_model, x = x_axis, y = y_axis, min_r2 = min_r2)
  }
  if (!is.null(palette)) p <- p + ggplot2::scale_fill_distiller(palette = palette, direction = 1, na.value = "grey90")
  save_plot(p, outfile, width = 10, height = 6)
}

# =============================================================================
# PLOT 9 — Variance totale par vue  (plot_variance_explained total)
# =============================================================================

if (plot_type == "variance_total") {
  cat("[PLOT] Total Variance Explained per View...\n")

  factors_raw  <- if (!is.null(params[["factors"]]) && params[["factors"]] != "") params[["factors"]] else "all"
  factors_sel  <- if (factors_raw == "all") "all" else parse_int_vector(factors_raw)
  palette      <- if (!is.null(params[["palette"]]) && params[["palette"]] != "") params[["palette"]] else NULL
  label_thresh <- if (!is.null(params[["label_threshold"]]) && params[["label_threshold"]] != "") as.numeric(params[["label_threshold"]]) else 1

  # Extract variance explained per factor per view (average across groups)
  r2_list <- get_variance_explained(mofa_model)$r2_per_factor
  if (length(r2_list) > 1) {
    r2_mat <- Reduce("+", r2_list) / length(r2_list)
  } else {
    r2_mat <- r2_list[[1]]
  }

  # Build long data frame (r2_mat: rows=factors, cols=views)
  df <- do.call(rbind, lapply(rownames(r2_mat), function(f) {
    data.frame(View     = colnames(r2_mat),
               Factor   = f,
               Variance = as.numeric(r2_mat[f, ]),
               stringsAsFactors = FALSE)
  }))

  # Filter factors if requested
  if (!identical(factors_sel, "all")) {
    valid_factors    <- paste0("Factor", factors_sel)
    existing_factors <- intersect(valid_factors, unique(df$Factor))
    if (length(existing_factors) == 0) {
      warning("No valid factors selected. Using all factors.")
    } else {
      df <- df[df$Factor %in% existing_factors, ]
    }
  }

  # Order levels: Factor1 at the base (bottom) of the stack
  factor_nums   <- sort(as.integer(gsub("Factor", "", unique(df$Factor))))
  factor_levels <- rev(paste0("Factor", factor_nums))
  df$Factor <- factor(df$Factor, levels = factor_levels)

  # Compute label positions (with proper numeric ordering)
  df_labels <- df %>%
    dplyr::group_by(View) %>%
    dplyr::arrange(as.numeric(gsub("Factor", "", Factor)), .by_group = TRUE) %>%
    dplyr::mutate(
      cumsum_before  = dplyr::lag(cumsum(Variance), default = 0),
      label_position = cumsum_before + Variance / 2,
      label_text     = ifelse(Variance > label_thresh,
                              paste0(round(Variance, 1), "%"), "")
    ) %>%
    dplyr::ungroup()

  # Default: saturated orange gradient (medium → very dark) — no washed-out yellows
  n_f <- length(levels(df$Factor))
  if (n_f == 1) {
    def_colors <- c("Factor1" = "#CC5500")
  } else {
    def_colors <- colorRampPalette(c("#7A2500", "#F5A623"))(n_f)
    names(def_colors) <- levels(df$Factor)
  }

  # Dynamic text color: dark text on light fill, white on dark fill
  active_colors <- if (!is.null(palette)) {
    tryCatch(
      RColorBrewer::brewer.pal(max(3, n_f), palette)[seq_len(n_f)],
      error = function(e) def_colors
    )
  } else def_colors

  text_colors <- sapply(active_colors, function(hex) {
    r <- strtoi(substr(hex, 2, 3), 16L)
    g <- strtoi(substr(hex, 4, 5), 16L)
    b <- strtoi(substr(hex, 6, 7), 16L)
    lum <- (0.299 * r + 0.587 * g + 0.114 * b) / 255
    if (lum > 0.45) "#222222" else "white"
  })
  df_labels$text_color <- text_colors[as.character(df_labels$Factor)]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Variance, y = View, fill = Factor)) +
    ggplot2::geom_col(position = "stack", width = 0.7) +
    ggplot2::geom_text(
      data = df_labels,
      ggplot2::aes(x = label_position, y = View, label = label_text,
                   color = I(text_color)),
      inherit.aes = FALSE,
      size = 3.8, fontface = "bold"
    ) +
    ggplot2::labs(x = "% of Variance", y = "", fill = "Factor") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "right"
    )

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_fill_brewer(palette = palette, direction = -1)
  } else {
    p <- p + ggplot2::scale_fill_manual(values = def_colors)
  }

  save_plot(p, outfile, width = 8, height = 6)
}

# =============================================================================
# PLOT 10 — Corrélation entre facteurs  (plot_factor_cor)
# =============================================================================

if (plot_type == "factor_correlation") {
  cat("[PLOT] Factor Correlation...\n")

  if (!preview_only) {
    pdf(outfile, width = 8, height = 8)
    plot_factor_cor(mofa_model)
    dev.off()
    cat(sprintf("  PDF saved: %s\n", outfile))
  }
  png_path <- sub("\\.pdf$", "_preview.png", outfile)
  grDevices::png(png_path, width = 800, height = 800, res = 100)
  plot_factor_cor(mofa_model)
  dev.off()
  cat(sprintf("  Preview saved: %s\n", png_path))
}

# =============================================================================
# PLOT 11 — Enrichissement / annotation des facteurs  (plot_enrichment_heatmap)
# Paramètres : view, factor, omic_features_col (optionnel)
# Note : nécessite que les features aient des annotations GSEAbase
# =============================================================================

if (plot_type == "enrichment_heatmap") {
  cat("[PLOT] Enrichment Heatmap...\n")

  if (requireNamespace("GSEABase", quietly = TRUE)) {
    view   <- params[["view"]]
    factor <- as.integer(params[["factor"]])

    if (!preview_only) {
      pdf(outfile, width = 12, height = 9)
      plot_enrichment_heatmap(mofa_model, view = view, factor = factor)
      dev.off()
      cat(sprintf("  \u2713 PDF saved: %s\n", outfile))
    }
    png_path <- sub("[.]pdf$", "_preview.png", outfile)
    grDevices::png(png_path, width = 1200, height = 900, res = 100)
    plot_enrichment_heatmap(mofa_model, view = view, factor = factor)
    dev.off()
    cat(sprintf("  \u2713 Preview saved: %s\n", png_path))
  } else {
    cat("  [AVERTISSEMENT] Package GSEABase non installé. Plot ignoré.\n")
  }
}

# =============================================================================
# PLOT 12 — Loadings bubbles / weights overview  (plot_weights_overview)
# Paramètres : view
# =============================================================================

if (plot_type == "weights_overview") {
  cat("[PLOT] Weights Overview...\n")

  view   <- params[["view"]]
  factor <- as.integer(params[["factor"]])

  p <- plot_weights(
    mofa_model,
    view      = view,
    factor    = factor,
    nfeatures = 0,      # 0 = tous
    scale     = FALSE,
    abs       = FALSE
  )
  save_plot(p, outfile, width = 12, height = 7)
}

# =============================================================================
# PLOT 13 — UMAP / tSNE sur les facteurs  (run_umap + plot_dimred)
# Paramètres : method ("UMAP" ou "TSNE"), color_by, n_neighbors, min_dist
# =============================================================================

if (plot_type == "dimred") {
  cat("[PLOT] Dimensionality Reduction (UMAP/tSNE)...\n")

  method     <- toupper(params[["method"]])   # "UMAP" ou "TSNE"
  color_by   <- params[["color_by"]]
  n_neighbors <- as.integer(params[["n_neighbors"]])
  min_dist    <- as.numeric(params[["min_dist"]])

  if (method == "UMAP") {
    if (!requireNamespace("uwot", quietly = TRUE)) {
      cat("  [ERREUR] Package uwot requis pour UMAP. Installez-le avec install.packages('uwot').\n")
      quit(status = 1)
    }
    mofa_model <- run_umap(mofa_model, n_neighbors = n_neighbors, min_dist = min_dist)
  } else if (method == "TSNE") {
    if (!requireNamespace("Rtsne", quietly = TRUE)) {
      cat("  [ERREUR] Package Rtsne requis pour tSNE. Installez-le avec install.packages('Rtsne').\n")
      quit(status = 1)
    }
    mofa_model <- run_tsne(mofa_model)
  }

  p <- plot_dimred(mofa_model, method = method, color_by = color_by)
  save_plot(p, outfile, width = 9, height = 7)
}

# =============================================================================
# PLOT 14 — Overview des données manquantes par vue/sample
# =============================================================================

if (plot_type == "data_overview") {
  cat("[PLOT] Data Overview...\n")

  if (!preview_only) {
    pdf(outfile, width = 10, height = 6)
    plot_data_overview(mofa_model)
    dev.off()
    cat(sprintf("  PDF saved: %s\n", outfile))
  }
  png_path <- sub("\\.pdf$", "_preview.png", outfile)
  grDevices::png(png_path, width = 1000, height = 600, res = 100)
  plot_data_overview(mofa_model)
  dev.off()
  cat(sprintf("  Preview saved: %s\n", png_path))
}

# =============================================================================
# PLOT 15 — Association Analysis  (correlate_factors_with_covariates)
# Parameters: covariates (comma-separated metadata column names), plot
# =============================================================================

if (plot_type == "association_analysis") {
  cat("[PLOT] Association Analysis...\n")

  covariates_raw <- params[["covariates"]]
  if (is.null(covariates_raw) || covariates_raw == "") {
    cat("  [ERROR] No covariates provided. Specify comma-separated metadata column names.\n")
    quit(status = 1)
  }

  # Check metadata is attached and show available columns
  meta <- tryCatch(samples_metadata(mofa_model), error = function(e) NULL)
  if (is.null(meta) || ncol(meta) == 0) {
    cat("  [ERROR] No metadata attached to this model. Run the pipeline with a metadata file first.\n")
    quit(status = 1)
  }
  cat(sprintf("  Available metadata columns: %s\n", paste(colnames(meta), collapse = ", ")))

  # Validate and filter covariates against actual column names
  covariates_req <- trimws(strsplit(covariates_raw, ",")[[1]])
  cat(sprintf("  Covariates requested      : %s\n", paste(covariates_req, collapse = ", ")))
  valid   <- covariates_req[covariates_req %in% colnames(meta)]
  invalid <- covariates_req[!covariates_req %in% colnames(meta)]
  if (length(invalid) > 0)
    cat(sprintf("  [WARNING] Not found and ignored: %s\n", paste(invalid, collapse = ", ")))
  if (length(valid) == 0) {
    cat("  [ERROR] None of the specified covariates match metadata columns. Aborting.\n")
    quit(status = 1)
  }
  cat(sprintf("  Using covariates: %s\n", paste(valid, collapse = ", ")))

  plot_stat  <- if (!is.null(params[["plot"]]) && params[["plot"]] != "") params[["plot"]] else "log_pval"

  p <- correlate_factors_with_covariates(
    mofa_model,
    covariates = valid,
    plot       = plot_stat
  )

  if (!preview_only) {
    ggplot2::ggsave(outfile, plot = p, width = 9, height = 6)
    cat(sprintf("  \u2713 PDF saved: %s\n", outfile))
  }
  png_path <- sub("\\.pdf$", "_preview.png", outfile)
  ggplot2::ggsave(png_path, plot = p, width = 9, height = 6, dpi = 100)
  cat(sprintf("  \u2713 Preview saved: %s\n", png_path))
}

# =============================================================================
# PLOT 16 — Variance Decomposition by Factor with Total  (plot_variance_explained + plot_total)
# Parameters: max_r2 (optional)
# =============================================================================

if (plot_type == "variance_detailed") {
  cat("[PLOT] Variance Decomposition by Factor (with total)...\n")

  min_r2  <- if (!is.null(params[["min_r2"]]) && params[["min_r2"]] != "") as.numeric(params[["min_r2"]]) else 0
  max_r2  <- if (!is.null(params[["max_r2"]]) && params[["max_r2"]] != "") as.numeric(params[["max_r2"]]) else NULL
  x_axis  <- if (!is.null(params[["x"]]) && params[["x"]] != "") params[["x"]] else "view"
  y_axis  <- if (!is.null(params[["y"]]) && params[["y"]] != "") params[["y"]] else "factor"
  palette <- if (!is.null(params[["palette"]]) && params[["palette"]] != "") params[["palette"]] else NULL
  cat(sprintf("  max_r2: %s\n", if (!is.null(max_r2)) as.character(max_r2) else "NULL (default)"))

  plots <- if (!is.null(max_r2)) {
    plot_variance_explained(mofa_model, x = x_axis, y = y_axis, min_r2 = min_r2, max_r2 = max_r2, plot_total = TRUE)
  } else {
    plot_variance_explained(mofa_model, x = x_axis, y = y_axis, min_r2 = min_r2, plot_total = TRUE)
  }
  if (!is.null(palette)) {
    plots[[1]] <- plots[[1]] + ggplot2::scale_fill_distiller(palette = palette, direction = 1, na.value = "grey90")
    plots[[2]] <- plots[[2]] + ggplot2::scale_fill_distiller(palette = palette, direction = 1, na.value = "grey90")
  }

  if (!preview_only) {
    pdf(outfile, width = 10, height = 6)
    print(plots[[1]])
    print(plots[[2]])
    dev.off()
    cat(sprintf("  PDF saved: %s\n", outfile))
  }
  png_path <- sub("\\.pdf$", "_preview.png", outfile)
  ggplot2::ggsave(png_path, plot = plots[[1]], width = 10, height = 6, dpi = 100)
  cat(sprintf("  Preview saved (page 1): %s\n", png_path))
}

# =============================================================================
# FIN
# =============================================================================

cat("\n========================================================\n")
cat(sprintf("  PLOT TERMINÉ : %s.pdf\n", out_name))
cat("========================================================\n")
