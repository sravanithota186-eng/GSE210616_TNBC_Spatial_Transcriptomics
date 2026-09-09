# GSE210616 TNBC Spatial Transcriptomics Reproducibility

Reproducibility package for the spatial transcriptomic analysis of triple-negative breast cancer (TNBC) using GEO dataset GSE210616.

## Dataset

- GEO accession: GSE210616
- Platform: 10x Genomics Visium
- Species: Homo sapiens
- Samples: 43 spatial sections
- Patients: 22
- Cancer type: Triple-negative breast cancer (TNBC)

## Analysis workflow

RCTD deconvolution
→ tissue compartment identification
→ LCP analysis
→ spatial niche analysis
→ ligand-receptor analysis
→ tumor-cell C-SIDE NCDEG analysis
→ NCDEG meta-analysis
→ enrichment analysis
→ LR-NCDEG integration
→ reproducibility QC

## Repository contents

- scripts/ — author reference scripts and TNBC analysis scripts
- esults/final_tables/ — final analysis tables
- esults/qc/ — final reproducibility and QC audit results
- docs/ — documentation
- metadata/ — dataset and analysis metadata

## Important reproducibility note

This repository contains final summarized results and analysis scripts. Raw sequencing data, large intermediate RDS objects, Visium images, and other large input files are intentionally excluded.

The analysis is a TNBC-focused adaptation of the published Pan-Cancer SpatialOmics workflow. Therefore, TNBC-specific niche definitions and some implementation details should not be interpreted as an exact reproduction of the Pan-Cancer cohort-level results.

## Final QC

The final reproducibility audit passed all critical and full audit checks:

- Critical checks: 23/23 PASS
- Full audit checks: 51/51 PASS

## Dataset source

GEO: GSE210616
