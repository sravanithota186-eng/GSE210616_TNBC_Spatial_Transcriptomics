# ============================================================
# GSE210616 Visium Pipeline — Stage 1 (v2)
# Single-sample processing, wrapped for batch use across all
# samples that have complete spatial files.
#
# CHANGES FROM v1:
#   - Added qc_diagnostics(): profiles nCount/nFeature/pct.mt for
#     EVERY sample BEFORE any filtering, so thresholds are chosen
#     from the actual pooled distribution instead of one sample.
#   - Replaced the fixed thresholds (1500/3000/5%) — which were
#     calibrated on the single deepest sample, GSM6433586, and
#     wiped out 85-99% of spots in shallower samples — with gentler
#     defaults that mirror the paper's own STAR Methods numbers:
#     the paper's sample-level filter is spot number>500,
#     median UMI>500, median mito%<20%. We use the same UMI/mito
#     numbers as spot-level floors (500 / 20%) plus a standard
#     Seurat spatial nFeature floor (200), rather than inventing
#     new absolute cutoffs.
#
# This reproduces (in simplified, marker-based form) the first
# stage of Li et al. 2026's pipeline:
#   raw counts -> QC diagnostics -> QC -> normalize -> PCA -> cluster -> annotate
#
# NOT yet implemented (later stages, see notes at bottom):
#   - RCTD deconvolution against a reference scRNA-seq dataset
#   - Malignant/Stromal/Benign compartment calling via inferCNV
#   - LCP / niche definition
#   - Ligand-receptor, NCDEG, survival analyses
# ============================================================

library(Seurat)
library(SeuratObject)
library(Matrix)
library(R.utils)
library(png)
library(jsonlite)
library(magick)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 0. CONFIG — your machine's paths
# ------------------------------------------------------------
base_dir <- "C:/Users/Sravani_123/Downloads/GSE210616"   # note: forward slashes, not backslashes

setwd(base_dir)   # makes all relative paths below resolve inside this folder

extracted_dir <- file.path(base_dir, "extracted")            # where the raw GSM_*.gz files live
staging_dir   <- file.path(base_dir, "spatial_staging")      # where we'll build standard Visium folders
output_dir    <- file.path(base_dir, "1_data_preprocessing") # where processed .rds files go
qc_dir        <- file.path(base_dir, "0_qc_diagnostics")     # NEW: pre-filtering QC plots/tables

dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir,      recursive = TRUE, showWarnings = FALSE)

# --- Gentler QC thresholds, aligned with the paper's own numbers ---
# Paper (STAR Methods): sample-level filter is
#   spot number > 500, median UMI > 500, median mito% < 20%
# We apply the UMI/mito values as spot-level floors (not the
# stricter values calibrated from one deep sample), plus the
# standard Seurat spatial-vignette nFeature floor.
QC_MIN_FEATURES <- 200
QC_MIN_COUNTS   <- 500
QC_MAX_PCT_MT   <- 20

# ------------------------------------------------------------
# 1. Build spatial_info table: one row per sample with complete files
#    (unchanged from v1)
# ------------------------------------------------------------
build_spatial_info <- function(extracted_dir) {

  h5_files <- list.files(
    extracted_dir,
    pattern = "filtered_feature_bc_matrix\\.h5$",
    recursive = TRUE, full.names = TRUE
  )
  h5_files <- h5_files[grepl("^GSM", basename(h5_files))]

  gsm_ids <- sub("_filtered_feature_bc_matrix\\.h5$", "", basename(h5_files))

  sf_files  <- list.files(extracted_dir, pattern = "scalefactors_json\\.json(\\.gz)?$",
                           recursive = TRUE, full.names = TRUE)
  pos_files <- list.files(extracted_dir, pattern = "tissue_positions.*\\.csv(\\.gz)?$",
                           recursive = TRUE, full.names = TRUE)
  img_files <- list.files(extracted_dir, pattern = "tissue_hires_image.*\\.png(\\.gz)?$",
                           recursive = TRUE, full.names = TRUE)

  info <- data.frame(sample = gsm_ids, h5 = h5_files,
                      scalefactors = NA_character_,
                      hires_image  = NA_character_,
                      positions    = NA_character_,
                      stringsAsFactors = FALSE)

  for (i in seq_len(nrow(info))) {
    s <- info$sample[i]
    sf  <- sf_files[grepl(paste0("^", s, "_scalefactors_json\\.json(\\.gz)?$"), basename(sf_files))]
    pos <- pos_files[grepl(paste0("^", s, "_tissue_positions.*\\.csv(\\.gz)?$"), basename(pos_files))]
    img <- img_files[grepl(paste0("^", s, "_tissue_hires_image.*\\.png(\\.gz)?$"), basename(img_files))]
    if (length(sf) > 0)  info$scalefactors[i] <- sf[1]
    if (length(img) > 0) info$hires_image[i]  <- img[1]
    if (length(pos) > 0) info$positions[i]    <- pos[1]
  }

  complete <- info[complete.cases(info), ]
  incomplete <- setdiff(info$sample, complete$sample)
  if (length(incomplete) > 0) {
    message("Skipping samples missing spatial files: ", paste(incomplete, collapse = ", "))
  }
  complete
}

# ------------------------------------------------------------
# 2. Stage one sample into a standard Visium folder structure
#    (unchanged from v1)
# ------------------------------------------------------------
stage_sample <- function(row, staging_dir) {
  sample_dir  <- file.path(staging_dir, row$sample)
  spatial_dir <- file.path(sample_dir, "spatial")
  dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

  h5_dest <- file.path(sample_dir, "filtered_feature_bc_matrix.h5")
  if (!file.exists(h5_dest)) file.copy(row$h5, h5_dest, overwrite = TRUE)

  decompress_if_needed <- function(src, dest) {
    if (grepl("\\.gz$", src)) {
      R.utils::gunzip(src, destname = dest, overwrite = TRUE, remove = FALSE)
    } else {
      file.copy(src, dest, overwrite = TRUE)
    }
  }

  decompress_if_needed(row$scalefactors, file.path(spatial_dir, "scalefactors_json.json"))
  decompress_if_needed(row$hires_image,  file.path(spatial_dir, "tissue_hires_image.png"))
  decompress_if_needed(row$positions,    file.path(spatial_dir, "tissue_positions_list.csv"))

  sample_dir
}

# ------------------------------------------------------------
# 3. NEW — QC diagnostics: profile every sample BEFORE filtering
# ------------------------------------------------------------
# Loads each sample's raw (unfiltered) counts, computes nCount/
# nFeature/percent.mt per spot, and returns both the pooled
# per-spot table and a per-sample summary. Use this to actually
# LOOK at the distributions before picking any cutoff — don't
# guess from one sample again.
qc_diagnostics <- function(spatial_info, staging_dir, qc_dir) {

  per_spot_list <- list()

  for (i in seq_len(nrow(spatial_info))) {
    row <- spatial_info[i, ]
    message("QC scan: ", row$sample)

    sample_dir <- stage_sample(row, staging_dir)

    obj <- tryCatch({
      Load10X_Spatial(
        data.dir   = sample_dir,
        filename   = "filtered_feature_bc_matrix.h5",
        assay      = "Spatial",
        slice      = row$sample,
        image.name = "tissue_hires_image.png"
      )
    }, error = function(e) {
      message("  Load10X_Spatial failed for ", row$sample, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(obj)) next

    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

    per_spot_list[[row$sample]] <- data.frame(
      sample_id  = row$sample,
      nCount     = obj$nCount_Spatial,
      nFeature   = obj$nFeature_Spatial,
      percent.mt = obj$percent.mt
    )
  }

  per_spot <- do.call(rbind, per_spot_list)
  rownames(per_spot) <- NULL

  summary_tbl <- per_spot %>%
    group_by(sample_id) %>%
    summarise(
      n_spots_raw      = n(),
      median_nCount    = median(nCount),
      q25_nCount       = quantile(nCount, 0.25),
      q75_nCount       = quantile(nCount, 0.75),
      median_nFeature  = median(nFeature),
      q25_nFeature     = quantile(nFeature, 0.25),
      q75_nFeature     = quantile(nFeature, 0.75),
      median_pct_mt    = median(percent.mt),
      pct_spots_would_fail_v1_thresholds = mean(
        nFeature < 1500 | nCount < 3000 | percent.mt > 5
      ) * 100,
      pct_spots_would_fail_new_thresholds = mean(
        nFeature < QC_MIN_FEATURES | nCount < QC_MIN_COUNTS | percent.mt > QC_MAX_PCT_MT
      ) * 100,
      .groups = "drop"
    ) %>%
    arrange(median_nCount)

  write.csv(summary_tbl, file.path(qc_dir, "qc_summary_by_sample.csv"), row.names = FALSE)
  write.csv(per_spot,   file.path(qc_dir, "qc_per_spot_raw.csv"),      row.names = FALSE)

  # Order samples by depth so the plots are readable
  per_spot$sample_id <- factor(per_spot$sample_id, levels = summary_tbl$sample_id)

  p_count <- ggplot(per_spot, aes(x = sample_id, y = nCount)) +
    geom_violin(scale = "width", fill = "steelblue", alpha = 0.6) +
    scale_y_log10() +
    geom_hline(yintercept = QC_MIN_COUNTS, linetype = "dashed", color = "red") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    labs(title = "nCount_Spatial per sample (pre-filter, log10 scale)",
         y = "nCount_Spatial", x = NULL)

  p_feature <- ggplot(per_spot, aes(x = sample_id, y = nFeature)) +
    geom_violin(scale = "width", fill = "darkorange", alpha = 0.6) +
    scale_y_log10() +
    geom_hline(yintercept = QC_MIN_FEATURES, linetype = "dashed", color = "red") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    labs(title = "nFeature_Spatial per sample (pre-filter, log10 scale)",
         y = "nFeature_Spatial", x = NULL)

  p_mt <- ggplot(per_spot, aes(x = sample_id, y = percent.mt)) +
    geom_violin(scale = "width", fill = "forestgreen", alpha = 0.6) +
    geom_hline(yintercept = QC_MAX_PCT_MT, linetype = "dashed", color = "red") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
    labs(title = "percent.mt per sample (pre-filter)", y = "% mitochondrial", x = NULL)

  ggsave(file.path(qc_dir, "qc_nCount_by_sample.png"),   p_count,   width = 14, height = 5, dpi = 150)
  ggsave(file.path(qc_dir, "qc_nFeature_by_sample.png"), p_feature, width = 14, height = 5, dpi = 150)
  ggsave(file.path(qc_dir, "qc_pct_mt_by_sample.png"),   p_mt,      width = 14, height = 5, dpi = 150)

  message("\nQC diagnostics written to: ", qc_dir)
  message("  - qc_summary_by_sample.csv  (per-sample medians/quartiles + %fail under old vs new thresholds)")
  message("  - qc_per_spot_raw.csv       (every spot, every sample, unfiltered)")
  message("  - 3 PNGs                    (nCount / nFeature / pct.mt violins across samples)")
  message("\nInspect qc_summary_by_sample.csv now, BEFORE running the batch pipeline.")
  message("If any sample's median_nCount or median_pct_mt already looks pathological")
  message("(e.g. median_nCount << 500, or median_pct_mt > 20), that sample is a candidate")
  message("for exclusion at the SAMPLE level (as the paper does), not for being force-fit")
  message("through spot-level filtering.")

  list(per_spot = per_spot, summary = summary_tbl)
}

# ------------------------------------------------------------
# 4. Full per-sample pipeline: load -> QC -> normalize -> cluster -> annotate
# ------------------------------------------------------------
run_visium_sample <- function(sample_id, sample_dir,
                               min_features = QC_MIN_FEATURES,
                               min_counts   = QC_MIN_COUNTS,
                               max_pct_mt   = QC_MAX_PCT_MT) {

  message("=== Processing ", sample_id, " ===")

  obj <- tryCatch({
    Load10X_Spatial(
      data.dir  = sample_dir,
      filename  = "filtered_feature_bc_matrix.h5",
      assay     = "Spatial",
      slice     = sample_id,
      image.name = "tissue_hires_image.png"   # avoids the lowres-image lookup failure
    )
  }, error = function(e) {
    message("  Load10X_Spatial failed for ", sample_id, ": ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(obj)) return(NULL)

  obj$sample_id <- sample_id

  # --- QC ---
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  n_before <- ncol(obj)
  obj <- subset(obj, subset = nFeature_Spatial >= min_features &
                              nCount_Spatial   >= min_counts &
                              percent.mt       <= max_pct_mt)
  pct_kept <- round(100 * ncol(obj) / n_before, 1)
  message("  QC: ", n_before, " -> ", ncol(obj), " spots retained (", pct_kept, "%)")
  if (pct_kept < 30) {
    message("  WARNING: retained <30% of spots for ", sample_id,
            " even under gentler thresholds — check qc_summary_by_sample.csv,",
            " this sample may need sample-level exclusion rather than spot-level QC.")
  }
  if (ncol(obj) < 50) {
    message("  Too few spots after QC, skipping ", sample_id)
    return(NULL)
  }

  # --- gene filter: keep genes detected in >=3 spots ---
  counts <- GetAssayData(obj, assay = "Spatial", layer = "counts")
  keep_genes <- Matrix::rowSums(counts > 0) >= 3
  obj <- obj[keep_genes, ]

  # --- normalize + SCT + PCA + cluster ---
  obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)
  obj <- RunPCA(obj, assay = "SCT", npcs = 30, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)

  # --- markers per cluster (provisional; RCTD deconvolution is the
  #     rigorous replacement for this in the paper's actual pipeline) ---
  DefaultAssay(obj) <- "SCT"
  markers <- tryCatch(
    FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE),
    error = function(e) { message("  FindAllMarkers failed: ", conditionMessage(e)); NULL }
  )

  list(object = obj, markers = markers)
}

# ------------------------------------------------------------
# 5. Batch driver
# ------------------------------------------------------------
run_pipeline_batch <- function(extracted_dir, staging_dir, output_dir, qc_dir, overwrite = FALSE) {

  spatial_info <- build_spatial_info(extracted_dir)
  message("Found ", nrow(spatial_info), " samples with complete spatial files.")

  # --- NEW: run diagnostics on every sample BEFORE any filtering ---
  qc <- qc_diagnostics(spatial_info, staging_dir, qc_dir)

  results_summary <- list()

  for (i in seq_len(nrow(spatial_info))) {
    row <- spatial_info[i, ]
    out_rds <- file.path(output_dir, paste0(row$sample, "_processed.rds"))

    if (file.exists(out_rds) && !overwrite) {
      message("Skipping ", row$sample, " — already processed")
      next
    }

    sample_dir <- stage_sample(row, staging_dir)
    result <- run_visium_sample(row$sample, sample_dir)

    if (is.null(result)) next

    saveRDS(result, out_rds)
    results_summary[[row$sample]] <- data.frame(
      sample_id = row$sample,
      n_spots   = ncol(result$object),
      n_genes   = nrow(result$object),
      n_clusters = length(unique(Idents(result$object)))
    )

    # free memory between samples
    rm(result)
    gc(verbose = FALSE)
  }

  do.call(rbind, results_summary)
}

# ------------------------------------------------------------
# 6. RUN
# ------------------------------------------------------------
summary_table <- run_pipeline_batch(extracted_dir, staging_dir, output_dir, qc_dir)
print(summary_table)

# Save the summary table itself so you have a record of spot/gene/cluster
# counts per sample without re-running everything
write.csv(summary_table,
          file.path(output_dir, "pipeline_summary.csv"),
          row.names = FALSE)

# ============================================================
# NEXT STAGES (not implemented here — see STAR Methods in paper)
# ============================================================
# 1. RCTD deconvolution (spacexr package):
#      Requires a matched scRNA-seq reference per cancer type
#      (Table S2/S3 in the paper). This replaces the marker-based
#      cluster annotation above with proper cell-type-fraction
#      estimates per spot.
#
# 2. Compartment calling (Malignant/Stromal/Benign):
#      Requires inferCNV run per sample using the top-1%
#      non-parenchymal spots as reference, followed by the
#      CNV-denoising + Otsu-threshold logic described in the
#      Method Details "Identification of tissue compartments".
#
# 3. LCP characterization:
#      K-means clustering (k=28-32) on Z-scored cell-type
#      composition matrices, done separately for Malignant and
#      Stromal compartments, across ALL samples pooled together
#      (not per-sample) — this is a cross-sample step.
#
# 4. Niche definition, LR interactions, NCDEGs, survival:
#      Each depends on having steps 1-3 done across the full
#      sample set first.
#
# 5. Sample-level exclusion (per paper's actual criteria):
#      The paper's real QC gate is at the SAMPLE level, not the
#      spot level: exclude samples with spot number<=500, OR
#      median UMI<=500, OR median mito%>=20% (computed on that
#      sample's own spots). If qc_summary_by_sample.csv shows a
#      sample already failing these on raw data, that's a sample
#      to drop entirely (as the paper drops 48 of 421), not one to
#      rescue by loosening thresholds further.
# ============================================================
