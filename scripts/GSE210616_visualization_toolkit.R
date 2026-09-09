# ============================================================
# GSE210616 Visium Pipeline — VISUALIZATION TOOLKIT
#
# Companion script to your Stage 1 (QC/clustering) and Stage 2
# (RCTD deconvolution) scripts. Nothing here re-runs QC or RCTD;
# every function just reads objects you've already saved and
# plots them.
#
# HOW TO USE
#   1. source() this file once per R session.
#   2. Load a processed sample:
#        result <- readRDS(".../1_data_preprocessing/GSM6433587_093A_processed.rds")
#        obj <- result$object
#   3. Call whichever plot function you need, e.g.:
#        plot_qc_violin(obj, sample_id = "GSM6433587_093A")
#        plot_pca(obj)
#        plot_umap_clusters(obj)
#
# Each SECTION below is independent — copy just the function(s)
# you need if you don't want the whole file.
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Output directory for saved figures (created if missing)
FIG_DIR <- "C:/Users/Sravani_123/Downloads/GSE210616/figures"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# SECTION 1 — QC VIOLIN / SCATTER PLOTS (per sample)
# ============================================================
# Reproduces the standard Seurat spatial-QC panel: nCount,
# nFeature, percent.mt, shown as violins, plus a scatter of
# nCount vs percent.mt to eyeball outliers.

plot_qc_violin <- function(seurat_obj, sample_id = NULL, save = TRUE) {
  sample_id <- sample_id %||% "sample"

  p1 <- VlnPlot(seurat_obj, features = "nCount_Spatial", pt.size = 0.1) +
    NoLegend() + ggtitle(paste0(sample_id, " — nCount_Spatial"))
  p2 <- VlnPlot(seurat_obj, features = "nFeature_Spatial", pt.size = 0.1) +
    NoLegend() + ggtitle(paste0(sample_id, " — nFeature_Spatial"))
  p3 <- VlnPlot(seurat_obj, features = "percent.mt", pt.size = 0.1) +
    NoLegend() + ggtitle(paste0(sample_id, " — percent.mt"))

  p4 <- FeatureScatter(seurat_obj, feature1 = "nCount_Spatial",
                        feature2 = "percent.mt") +
    ggtitle(paste0(sample_id, " — nCount vs %mt"))

  combined <- (p1 | p2 | p3) / p4

  if (save) {
    ggsave(file.path(FIG_DIR, paste0(sample_id, "_QC_violin.png")),
           combined, width = 12, height = 8, dpi = 150)
  }
  combined
}

# Pooled version across ALL samples, reusing the per-spot table
# already written by qc_diagnostics() in Stage 1 — no re-reading
# of raw h5 files needed.
plot_qc_pooled <- function(qc_per_spot_csv =
                              "C:/Users/Sravani_123/Downloads/GSE210616/0_qc_diagnostics/qc_per_spot_raw.csv",
                            save = TRUE) {
  per_spot <- read.csv(qc_per_spot_csv)

  p <- ggplot(per_spot, aes(x = reorder(sample_id, nCount, median),
                             y = nCount)) +
    geom_violin(scale = "width", fill = "steelblue", alpha = 0.6) +
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    labs(title = "nCount_Spatial across all samples (pre-filter)",
         x = NULL, y = "nCount_Spatial (log10)")

  if (save) {
    ggsave(file.path(FIG_DIR, "pooled_QC_nCount.png"),
           p, width = 14, height = 5, dpi = 150)
  }
  p
}


# ============================================================
# SECTION 2 — PCA PLOTS
# ============================================================
# Elbow plot (how many PCs carry signal), PCA scatter colored
# by cluster, and a loadings heatmap for the first few PCs.

plot_pca <- function(seurat_obj, sample_id = NULL, ndims = 30, save = TRUE) {
  sample_id <- sample_id %||% "sample"

  p_elbow <- ElbowPlot(seurat_obj, ndims = ndims) +
    ggtitle(paste0(sample_id, " — Elbow plot"))

  p_scatter <- DimPlot(seurat_obj, reduction = "pca",
                        group.by = "seurat_clusters", label = TRUE) +
    ggtitle(paste0(sample_id, " — PCA (PC1 vs PC2)"))

  combined <- p_elbow | p_scatter

  if (save) {
    ggsave(file.path(FIG_DIR, paste0(sample_id, "_PCA.png")),
           combined, width = 10, height = 5, dpi = 150)
  }
  combined
}

# Gene loading heatmap for the top PCs — useful for sanity-
# checking what biology each PC captures.
plot_pca_loadings <- function(seurat_obj, dims = 1:6, save = TRUE,
                               sample_id = NULL) {
  sample_id <- sample_id %||% "sample"
  png(file.path(FIG_DIR, paste0(sample_id, "_PCA_loadings.png")),
      width = 1000, height = 800, res = 120)
  DimHeatmap(seurat_obj, dims = dims, cells = 500, balanced = TRUE)
  dev.off()
  message("Saved PCA loadings heatmap for ", sample_id)
}


# ============================================================
# SECTION 3 — UMAP PLOTS
# ============================================================

# Standard cluster UMAP (matches Fig. 1B/1C style: one dot per
# spot, colored by cluster/LCP label).
plot_umap_clusters <- function(seurat_obj, sample_id = NULL, save = TRUE) {
  sample_id <- sample_id %||% "sample"

  p <- DimPlot(seurat_obj, reduction = "umap",
               group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
    ggtitle(paste0(sample_id, " — UMAP by cluster")) +
    coord_fixed()

  if (save) {
    ggsave(file.path(FIG_DIR, paste0(sample_id, "_UMAP_clusters.png")),
           p, width = 6, height = 5, dpi = 150)
  }
  p
}

# UMAP colored by a continuous QC metric — good for spotting
# clusters that are just "low quality spots" rather than biology.
plot_umap_qc <- function(seurat_obj, sample_id = NULL, save = TRUE) {
  sample_id <- sample_id %||% "sample"

  p1 <- FeaturePlot(seurat_obj, features = "nCount_Spatial", reduction = "umap") +
    ggtitle("nCount_Spatial")
  p2 <- FeaturePlot(seurat_obj, features = "percent.mt", reduction = "umap") +
    ggtitle("percent.mt")

  combined <- p1 | p2
  if (save) {
    ggsave(file.path(FIG_DIR, paste0(sample_id, "_UMAP_QC.png")),
           combined, width = 10, height = 5, dpi = 150)
  }
  combined
}

# UMAP colored by expression of specific marker genes (e.g. from
# your cluster markers table, or canonical markers you choose).
plot_umap_markers <- function(seurat_obj, genes, sample_id = NULL, save = TRUE) {
  sample_id <- sample_id %||% "sample"
  DefaultAssay(seurat_obj) <- "SCT"
  genes <- genes[genes %in% rownames(seurat_obj)]
  if (length(genes) == 0) {
    message("None of the requested genes are in this object.")
    return(invisible(NULL))
  }
  p <- FeaturePlot(seurat_obj, features = genes, reduction = "umap", ncol = 3)
  if (save) {
    ggsave(file.path(FIG_DIR, paste0(sample_id, "_UMAP_markers.png")),
           p, width = 4 * min(3, length(genes)), height = 4 * ceiling(length(genes) / 3),
           dpi = 150)
  }
  p
}


# ============================================================
# SECTION 4 — SPATIAL MAPS (cluster + gene expression, in tissue)
# ============================================================
# Works whether or not the Seurat object has an attached
# Images() slot: if it does, we use Seurat's SpatialDimPlot /
# SpatialFeaturePlot (nicer, shows the H&E image); if not
# (e.g. the manually-built GSM6433585/586 objects), we fall
# back to a manual ggplot using GetTissueCoordinates().

.has_image <- function(seurat_obj) length(Images(seurat_obj)) > 0

plot_spatial_clusters <- function(seurat_obj, sample_id = NULL, save = TRUE) {
  sample_id <- sample_id %||% "sample"

  if (.has_image(seurat_obj)) {
    p <- SpatialDimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE) +
      ggtitle(paste0(sample_id, " — spatial clusters"))
  } else {
    coords <- GetTissueCoordinates(seurat_obj)
    coords$cluster <- seurat_obj$seurat_clusters[rownames(coords)]
    p <- ggplot(coords, aes(x = x, y = -y, color = cluster)) +
      geom_point(size = 0.8) +
      coord_fixed() +
      theme_void() +
      ggtitle(paste0(sample_id, " — spatial clusters (no H&E image attached)"))
  }

  if (save) {
    ggsave(file.path(FIG_DIR, paste0(sample_id, "_spatial_clusters.png")),
           p, width = 7, height = 6, dpi = 150)
  }
  p
}

plot_spatial_gene <- function(seurat_obj, gene, sample_id = NULL, save = TRUE) {
  sample_id <- sample_id %||% "sample"
  DefaultAssay(seurat_obj) <- "SCT"
  if (!(gene %in% rownames(seurat_obj))) {
    message(gene, " not found in this object.")
    return(invisible(NULL))
  }

  if (.has_image(seurat_obj)) {
    p <- SpatialFeaturePlot(seurat_obj, features = gene) +
      ggtitle(paste0(sample_id, " — ", gene))
  } else {
    coords <- GetTissueCoordinates(seurat_obj)
    expr <- FetchData(seurat_obj, vars = gene)
    coords$expr <- expr[rownames(coords), 1]
    p <- ggplot(coords, aes(x = x, y = -y, color = expr)) +
      geom_point(size = 0.8) +
      scale_color_viridis_c() +
      coord_fixed() +
      theme_void() +
      ggtitle(paste0(sample_id, " — ", gene, " (no H&E image attached)"))
  }

  if (save) {
    ggsave(file.path(FIG_DIR, paste0(sample_id, "_spatial_", gene, ".png")),
           p, width = 7, height = 6, dpi = 150)
  }
  p
}


# ============================================================
# SECTION 5 — RCTD DECONVOLUTION PLOTS
# ============================================================
# Assumes `rctd_obj` is the RCTD S4 object AFTER run.RCTD()
# (i.e. rctd_obj@results is populated), matching your Stage 2c
# / console-history workflow.

# 5a. Spot-class breakdown (reject / singlet / doublet_certain /
#     doublet_uncertain) — quick sanity check of deconvolution quality.
plot_rctd_spotclass <- function(rctd_obj, sample_id = NULL, save = TRUE) {
  sample_id <- sample_id %||% "sample"
  df <- rctd_obj@results$results_df

  p <- ggplot(df, aes(x = spot_class, fill = spot_class)) +
    geom_bar() +
    theme_minimal() +
    labs(title = paste0(sample_id, " — RCTD spot classification"),
         x = NULL, y = "Number of spots") +
    theme(legend.position = "none")

  if (save) {
    ggsave(file.path(FIG_DIR, paste0(sample_id, "_RCTD_spotclass.png")),
           p, width = 5, height = 4, dpi = 150)
  }
  p
}

# 5b. Cell-type composition, stacked bar per spot-class group
#     (mirrors Figure 1D/1E-style "mean composition" summaries).
plot_rctd_composition_bar <- function(rctd_obj, sample_id = NULL, save = TRUE) {
  sample_id <- sample_id %||% "sample"

  weights_norm <- normalize_weights(rctd_obj@results$weights)
  comp <- as.data.frame(as.matrix(weights_norm))
  comp$spot_class <- rctd_obj@results$results_df[rownames(comp), "spot_class"]

  comp_long <- comp %>%
    tidyr::pivot_longer(-spot_class, names_to = "cell_type", values_to = "fraction") %>%
    group_by(spot_class, cell_type) %>%
    summarise(mean_fraction = mean(fraction), .groups = "drop")

  p <- ggplot(comp_long, aes(x = spot_class, y = mean_fraction, fill = cell_type)) +
    geom_col(position = "stack") +
    theme_minimal() +
    labs(title = paste0(sample_id, " — mean cell-type composition by spot class"),
         x = NULL, y = "Mean fraction") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (save) {
    ggsave(file.path(FIG_DIR, paste0(sample_id, "_RCTD_composition_bar.png")),
           p, width = 8, height = 5, dpi = 150)
  }
  p
}

# 5c. Spatial map of the DOMINANT cell type per spot (the
#     "first_type" RCTD already assigns) — cheap way to see the
#     tissue architecture before doing full LCP/niche clustering.
plot_rctd_spatial_dominant <- function(rctd_obj, sample_id = NULL, save = TRUE) {
  sample_id <- sample_id %||% "sample"

  coords <- rctd_obj@spatialRNA@coords
  coords$dominant_type <- rctd_obj@results$results_df[rownames(coords), "first_type"]

  p <- ggplot(coords, aes(x = x, y = -y, color = dominant_type)) +
    geom_point(size = 0.8) +
    coord_fixed() +
    theme_void() +
    ggtitle(paste0(sample_id, " — dominant cell type per spot (RCTD first_type)"))

  if (save) {
    ggsave(file.path(FIG_DIR, paste0(sample_id, "_RCTD_spatial_dominant.png")),
           p, width = 7, height = 6, dpi = 150)
  }
  p
}

# 5d. Spatial "pie chart" map showing full cell-type mixture per
#     spot (closest visual match to the paper's spot-composition
#     pie plots in Figure 2C). Requires the `scatterpie` package
#     (install.packages("scatterpie") once); subsampled by
#     default since plotting a pie per spot is slow for >1000 spots.
plot_rctd_spatial_piechart <- function(rctd_obj, sample_id = NULL,
                                        n_spots_subsample = 300,
                                        pie_radius = 40, save = TRUE) {
  if (!requireNamespace("scatterpie", quietly = TRUE)) {
    stop("Run install.packages('scatterpie') first.")
  }
  sample_id <- sample_id %||% "sample"

  weights_norm <- normalize_weights(rctd_obj@results$weights)
  comp <- as.data.frame(as.matrix(weights_norm))
  coords <- rctd_obj@spatialRNA@coords
  comp <- cbind(coords[rownames(comp), ], comp)

  if (nrow(comp) > n_spots_subsample) {
    set.seed(1)
    comp <- comp[sample(nrow(comp), n_spots_subsample), ]
  }

  cell_types <- setdiff(colnames(comp), c("x", "y"))

  p <- ggplot() +
    scatterpie::geom_scatterpie(data = comp, aes(x = x, y = -y, r = pie_radius),
                                 cols = cell_types, color = NA) +
    coord_fixed() +
    theme_void() +
    ggtitle(paste0(sample_id, " — spatial cell-type composition (subsampled spots)"))

  if (save) {
    ggsave(file.path(FIG_DIR, paste0(sample_id, "_RCTD_spatial_piechart.png")),
           p, width = 8, height = 7, dpi = 150)
  }
  p
}


# ============================================================
# SECTION 6 — BATCH EXPORT ACROSS ALL SAMPLES
# ============================================================
# Loops through every processed .rds in Stage 1 output, saving
# QC + PCA + UMAP + spatial-cluster figures per sample into one
# multi-page PDF. RCTD plots are batched separately since they
# live in a different folder/object type.

batch_export_stage1_plots <- function(
    stage1_dir = "C:/Users/Sravani_123/Downloads/GSE210616/1_data_preprocessing",
    out_pdf = file.path(FIG_DIR, "all_samples_stage1_plots.pdf")) {

  files <- list.files(stage1_dir, pattern = "_processed\\.rds$", full.names = TRUE)
  message("Found ", length(files), " processed samples to plot.")

  pdf(out_pdf, width = 12, height = 8)
  for (f in files) {
    sample_id <- sub("_processed\\.rds$", "", basename(f))
    message("Plotting: ", sample_id)

    result <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(result)) next
    obj <- if (is.list(result) && "object" %in% names(result)) result$object else result

    print(plot_qc_violin(obj, sample_id, save = FALSE))
    print(plot_pca(obj, sample_id, save = FALSE))
    print(plot_umap_clusters(obj, sample_id, save = FALSE))
    print(plot_spatial_clusters(obj, sample_id, save = FALSE))
  }
  dev.off()
  message("\nSaved multi-page PDF to: ", out_pdf)
}

batch_export_rctd_plots <- function(
    rctd_dir = "C:/Users/Sravani_123/Downloads/GSE210616/2_rctd",
    out_pdf = file.path(FIG_DIR, "all_samples_RCTD_plots.pdf")) {

  files <- list.files(rctd_dir, pattern = "_rctd.*\\.rds$", full.names = TRUE)
  message("Found ", length(files), " RCTD objects to plot.")

  pdf(out_pdf, width = 10, height = 7)
  for (f in files) {
    sample_id <- sub("_rctd.*\\.rds$", "", basename(f))
    message("Plotting RCTD: ", sample_id)

    rctd_obj <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(rctd_obj) || !("results" %in% slotNames(rctd_obj))) next

    print(plot_rctd_spotclass(rctd_obj, sample_id, save = FALSE))
    print(plot_rctd_composition_bar(rctd_obj, sample_id, save = FALSE))
    print(plot_rctd_spatial_dominant(rctd_obj, sample_id, save = FALSE))
  }
  dev.off()
  message("\nSaved multi-page PDF to: ", out_pdf)
}


# ============================================================
# QUICK USAGE EXAMPLES (not run automatically)
# ============================================================
# result <- readRDS("C:/Users/Sravani_123/Downloads/GSE210616/1_data_preprocessing/GSM6433587_093A_processed.rds")
# obj <- result$object
#
# plot_qc_violin(obj, "GSM6433587_093A")
# plot_pca(obj, "GSM6433587_093A")
# plot_umap_clusters(obj, "GSM6433587_093A")
# plot_spatial_clusters(obj, "GSM6433587_093A")
# plot_spatial_gene(obj, "EPCAM", "GSM6433587_093A")
#
# rctd_obj <- readRDS("C:/Users/Sravani_123/Downloads/GSE210616/2_rctd/GSM6433586_092B_rctd.rds")
# plot_rctd_spotclass(rctd_obj, "GSM6433586_092B")
# plot_rctd_spatial_dominant(rctd_obj, "GSM6433586_092B")
#
# # Everything, for every sample, in one go:
# batch_export_stage1_plots()
# batch_export_rctd_plots()
