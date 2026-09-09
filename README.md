# GSE210616 TNBC Spatial Transcriptomics Reproducibility

## Project Overview

This repository contains a **TNBC-focused computational reproduction and methodological adaptation** of the spatial transcriptomics workflow described by:

**Li J, Lin P, Wang H, et al.**  
**Pan-cancer analysis of spatial transcriptomics reveals heterogeneous tumor spatial microenvironment.**  
*Cell Reports Medicine.* 2026;7:102751.

The reference study developed a large-scale spatial transcriptomics framework for characterizing the tumor spatial microenvironment (TSME) across **373 samples from 12 cancer types**.

This project applies the methodological framework of that study specifically to the **triple-negative breast cancer (TNBC)** spatial transcriptomics cohort **GSE210616**.

The analysis integrates:

- Spatial transcriptomics quality control
- RCTD cell-type deconvolution
- Tissue compartment classification
- Local cellular program (LCP) analysis
- Malignant LCP (MLCP) analysis
- Stromal LCP (SLCP) analysis
- Spatial niche analysis
- Ligand-receptor (LR) analysis
- LR meta-analysis
- C-SIDE analysis
- Niche-related cell-type-specific differential expression (NCDEGs)
- NCDEG meta-analysis
- Gene-set enrichment analysis
- LR-NCDEG integration
- Reproducibility auditing

> **Important:** This is a TNBC-focused methodological reproduction/adaptation. It does not claim to reproduce the exact pan-cancer numerical results of the reference study.

---

# Reference Methodology

The primary methodological reference for this project is:

**Li et al., 2026 — Pan-cancer analysis of spatial transcriptomics reveals heterogeneous tumor spatial microenvironment.**

The reference study analyzed spatial transcriptomic data from multiple cancer types and developed a workflow for identifying cellular programs, spatial niches, ligand-receptor interactions, and niche-associated gene-expression changes.

The overall methodology followed in this project is:

```text
10x Visium Spatial Transcriptomics
              |
              v
        Quality Control
              |
              v
     RCTD Cell-Type Deconvolution
              |
              v
     Tissue Compartment Calling
              |
              v
       MLCP / SLCP Analysis
              |
              v
        Spatial Niches
              |
              v
   Ligand-Receptor Analysis
              |
              v
       C-SIDE / NCDEGs
              |
              v
      NCDEG Meta-analysis
              |
              v
       Enrichment Analysis
              |
              v
     LR-NCDEG Integration
              |
              v
     Reproducibility Audit

Dataset
GSE210616
GEO accession: GSE210616
Organism: Homo sapiens
Technology: 10x Genomics Visium
Cancer type: Triple-negative breast cancer (TNBC)
Patients: 22
Spatial sections: 43
Samples: GSM6433585–GSM6433627
Project: TNBC spatial transcriptomics

The dataset contains spatial transcriptomic sections from TNBC patient samples.

Large raw/intermediate files, including sequencing data, Visium images, and large RDS objects, are intentionally excluded from this repository.

Single-Cell Reference Dataset

RCTD cell-type deconvolution was performed using the breast cancer single-cell RNA-seq reference associated with:

GSE176078 — Wu et al. BRCA single-cell RNA-seq reference

A TNBC-focused reference was prepared for spatial deconvolution.

The final RCTD reference contained 11 cell types:

B
Fibroblast
Tumor cell
Endothelial cell
cDC
pDC
Macrophage
Plasma cell
CD4 T
CD8 T
NK

The TNBC-specific reference was used to estimate cell-type composition for spatial spots in GSE210616.

Project Objective

The objective of this project was to reproduce and operationalize the major computational components of the spatial transcriptomics methodology described by Li et al. (2026) on a TNBC-specific dataset.

The project emphasizes:

Methodological reproducibility
TNBC-specific adaptation
Explicit sample tracking
Stage-wise analysis
Structured result generation
Independent quality-control auditing
Reproducible computational organization
Analysis Workflow
GSE210616 TNBC Visium Dataset
              |
              v
       Spatial Data QC
              |
              v
     Spatial Coordinates
       and Spot Recovery
              |
              v
     RCTD Reference Setup
              |
              v
      RCTD Deconvolution
              |
              v
   Tissue Compartment Calling
              |
       +------+------+
       |             |
       v             v
     MLCP          SLCP
 Malignant LCP   Stromal LCP
       |             |
       +------+------+
              |
              v
      Spatial Niche Analysis
              |
       +------+------+
       |             |
       v             v
   LR Analysis     C-SIDE
       |             |
       v             v
 LR Meta-analysis  NCDEGs
       |             |
       +------+------+
              |
              v
      NCDEG Meta-analysis
              |
              v
       Enrichment Analysis
              |
              v
      LR-NCDEG Integration
              |
              v
       Reproducibility Audit
              |
              v
         Final Results
1. Spatial Transcriptomics Quality Control

The GSE210616 spatial transcriptomics sections were processed using a reproducible QC workflow.

QC included assessment of:

Number of spatial spots
UMI counts
Number of detected genes
Mitochondrial content
Ribosomal content
Spatial coordinates
Tissue-associated spots
Spatial QC visualizations

The objective was to generate consistent spatial input for downstream RCTD and spatial microenvironment analyses.

2. RCTD Cell-Type Deconvolution

Cell-type composition was estimated using RCTD from the spacexr framework.

The TNBC-focused GSE176078 single-cell reference was used for deconvolution.

Final RCTD result
Spatial sections analyzed: 43
Successful RCTD analyses: 43 / 43
Cell types: 11

RCTD generated spatial cell-type composition estimates for each Visium spot.

These estimates were subsequently used for:

Tissue compartment classification
LCP identification
Spatial niche analysis
3. Tissue Compartment Classification

The spatial tissue compartment workflow was adapted from the reference methodology.

The major compartments were:

Benign
Stromal
Malignant

The analysis incorporated cellular composition together with tumor-associated genomic information.

The workflow used:

RCTD-derived cell-type composition
inferCNV-based information
Proportion-based classification
Otsu/custom thresholding logic

The resulting compartment assignments were used for downstream MLCP and SLCP analyses.

4. Local Cellular Programs (LCPs)

The reference methodology identifies recurring local cellular programs based on spatial cell-type composition.

Because malignant and stromal spatial regions have different cellular compositions, the analysis was separated into:

Malignant compartment
        |
        v
Malignant Local Cellular Programs (MLCPs)


Stromal compartment
        |
        v
Stromal Local Cellular Programs (SLCPs)

LCPs represent recurring local combinations of cell types within the spatial tumor microenvironment.

The number and composition of TNBC LCPs are dataset-specific.

Therefore, the TNBC results should not be interpreted as an exact reproduction of the 28 MLCPs and 28 SLCPs reported in the pan-cancer reference study.

5. Spatial Niche Analysis

Spatial niches were constructed from the local neighborhood composition of LCPs.

The conceptual workflow was:

Spatial spots
     |
     v
Cell-type composition
     |
     v
LCP assignment
     |
     v
Spatial neighborhoods
     |
     v
Neighborhood LCP composition
     |
     v
Spatial niche clustering

A spatial niche represents a larger spatial microenvironment characterized by neighboring cellular programs.

The reference study identified 13 consensus niches across its pan-cancer cohort.

The current project identifies TNBC-specific spatial patterns and therefore does not force the TNBC data to reproduce those 13 pan-cancer niches.

6. Ligand-Receptor (LR) Analysis

Niche-related ligand-receptor interactions were analyzed using the conceptual framework of the reference study.

The workflow included:

CellChatDB
     |
     v
LR interaction definitions
     |
     v
Gene filtering
     |
     v
Sample-level LR testing
     |
     v
Niche comparisons
     |
     v
Rank-based meta-analysis
     |
     v
Final LR candidates
CellChat database summary
Total LR interactions: 3,233
Unique ligands: 800
Unique receptors: 780
Study-level LR analysis
Common genes across all 43 samples: 1,930
Union LR-testable genes: 2,187
Study-level LR marker records: 148,187
Final LR results
LR meta-analysis records: 10,567
Final LR candidate records: 22
Unique LR interactions: 17
Pathways represented: 12
Niche contexts represented: 11

The final candidate pathways included:

CXCL
CCL
ANGPTL
OSM
BAFF
CD40
CD45
LCK
TNF
SELPLG
IFN-II
CD86

All 22 final LR candidate records had corresponding CellChat identities.

7. C-SIDE / NCDEG Analysis

The reference methodology uses cell-type-specific differential expression to identify gene-expression changes associated with spatial niches.

In this project, C-SIDE analysis was performed for eligible MLCP and SLCP comparisons.

C-SIDE results
Eligible comparisons: 73
Completed comparisons: 73 / 73

MLCP comparisons: 14
SLCP comparisons: 59

Gene-level records: 157,713
Unique genes: 3,625

The resulting C-SIDE records were subsequently used for NCDEG meta-analysis.

8. NCDEG Meta-analysis

C-SIDE results were aggregated using the implemented meta-analysis framework.

Results
NCDEG meta-analysis records: 70,587
Unique genes: 3,625
Eligible branch × experimental-niche groups: 29

The eligible groups were:

MLCP: 8 groups
SLCP: 21 groups

This analysis was used to identify reproducible niche-associated gene-expression changes.

9. Final NCDEG Results

After applying the final filtering criteria:

Final NCDEG records: 56,512
Unique NCDEG genes: 3,597
Higher-expression records: 30,653
Lower-expression records: 25,859

The final signed-rank filtering framework included:

Higher:
signed_combined_rank >= 0.5

Lower:
signed_combined_rank <= -0.5

together with the implemented study-support criteria.

10. Replication Support

Replication support was evaluated by determining how many independent studies/comparisons supported individual NCDEG findings.

MLCP
1 study: 425 genes
2 studies: 774 genes
3 studies: 1,510 genes
SLCP
1 study:   578 genes
2 studies: 422 genes
3 studies: 405 genes
4 studies: 328 genes
5 studies: 393 genes
6+ studies: 1,371 genes

These values describe the replication structure of the implemented TNBC analysis.

11. Gene-set Enrichment Analysis

The final NCDEG results were evaluated against the curated gene-set collection recovered from the reference workflow.

The recovered collection contained approximately:

Gene-set associations: 14,991
Unique terms: 166
Unique genes: 6,512
Sources: 5
Final enrichment results
Enrichment records: 1,228
Unique enriched terms: 68

Frequently recurrent programs included:

HALLMARK_MTORC1_SIGNALING
HALLMARK_MYC_TARGETS_V1
HALLMARK_OXIDATIVE_PHOSPHORYLATION
HALLMARK_ADIPOGENESIS
HALLMARK_DNA_REPAIR
HALLMARK_E2F_TARGETS
HALLMARK_G2M_CHECKPOINT
HALLMARK_MYC_TARGETS_V2

Other recurrent programs included:

Proteasomal Degradation
Translation initiation
MYC
Cell Cycle
Respiration
Secreted programs
12. MLCP versus SLCP Enrichment

Enrichment patterns were compared between malignant and stromal branches.

MLCP-specific terms: 2
SLCP-specific terms: 10
Shared terms: 56

Term recurrence was further summarized as:

Ubiquitous terms: 21
Niche-specific terms: 5

Directionality analysis showed:

Terms represented in higher-expression results: 63
Terms represented in lower-expression results: 63
Higher-only: 5
Lower-only: 5
Both higher and lower: 58
13. LR-NCDEG Integration

The final LR and NCDEG results were integrated to characterize relationships between spatial communication signals and niche-associated transcriptional programs.

Direct gene-level overlap
Unique LR genes: 33
Unique NCDEG genes: 3,597
Direct LR-NCDEG overlap: 0

Therefore, there was no direct gene-level overlap between the final LR candidate genes and final NCDEG genes under the implemented filtering criteria.

This does not mean that the biological processes represented by LR signaling and NCDEGs are unrelated.

It means that the final filtered gene sets contained no identical genes.

14. LR versus NCDEG Niche Context

The final LR and NCDEG results were also compared by niche context.

LR candidate niche contexts: 11
NCDEG enrichment contexts: 29
Exact shared niche contexts: 0

Therefore, no exact niche-context match was observed between the final LR candidate results and NCDEG enrichment results.

Both analyses nevertheless contained higher-level:

MLCP
SLCP

branch information.

For this reason, the final integration was performed at the broader branch/evidence level rather than being interpreted as direct LR-NCDEG statistical association.

15. Branch-Level LR-NCDEG Integration
MLCP
LR records: 17
NCDEG enrichment records: 326
Integration combinations: 5,542
SLCP
LR records: 5
NCDEG enrichment records: 902
Integration combinations: 4,510
Total
LR records: 22
NCDEG enrichment records: 1,228
Integration combinations: 10,052
Important interpretation

10,052 is the number of evidence-characterization combinations generated during the integration procedure.

It is not the number of statistically significant LR-NCDEG associations.

16. Reproducibility Audit

A dedicated reproducibility audit was developed to validate the final analysis.

The audit checked:

Dataset/sample counts
RCTD completion
Cell-type reference
LR analysis
LR meta-analysis
LR filtering
C-SIDE completion
NCDEG extraction
NCDEG meta-analysis
NCDEG filtering
Enrichment results
LR-NCDEG integration
Final result organization
Final Audit
Critical checks: 23 / 23 PASS

Full audit: 51 / 51 PASS

The final audit outputs are included in the repository under:

results/qc/
Final Project Summary
Component	Final Result
Spatial sections	43
Patients	22
RCTD completion	43 / 43
RCTD cell types	11
LR study-level records	148,187
LR meta-analysis records	10,567
Final LR candidate records	22
Unique LR interactions	17
LR pathways	12
LR niche contexts	11
C-SIDE comparisons	73 / 73
C-SIDE unique genes	3,625
NCDEG meta-analysis records	70,587
Final NCDEG records	56,512
Final NCDEG unique genes	3,597
Enrichment records	1,228
Unique enrichment terms	68
Direct LR-NCDEG gene overlap	0
Exact LR-NCDEG niche-context overlap	0
Branch-level integration combinations	10,052
Critical audit	23 / 23 PASS
Full audit	51 / 51 PASS
Software Environment

The main computational environment used for the project included:

R                  4.5.1
Seurat             5.5.1
SeuratObject       5.4.0
spacexr            2.2.1
infercnv           1.24.0
CellChat           2.2.0.9001
clusterProfiler    4.2.2
Repository Structure
GSE210616_TNBC_Spatial_Transcriptomics/
|
├── README.md
├── .gitignore
|
├── metadata/
|   └── DATASET_METADATA.txt
|
├── scripts/
|   ├── author_reference/
|   └── TNBC_analysis/
|
├── results/
|   ├── final_tables/
|   |   ├── STEP64_DATASET_RECOVERY_REPORT.csv
|   |   ├── STEP64_FINAL_FIGURES_TABLES_SUMMARY.txt
|   |   ├── STEP64_TABLE1_PIPELINE_SUMMARY.csv
|   |   ├── STEP64_TABLE2_FINAL_LR_CANDIDATES.csv
|   |   ├── STEP64_TABLE3_LR_PATHWAY_CHARACTERIZATION.csv
|   |   ├── STEP64_TABLE4_CSIDE_BRANCH_SUMMARY.csv
|   |   ├── STEP64_TABLE5_NCDEG_DIRECTION_SUMMARY.csv
|   |   ├── STEP64_TABLE6_ENRICHMENT_TERM_RECURRENCE.csv
|   |   ├── STEP64_TABLE7_MLCP_SLCP_COMPARISON.csv
|   |   ├── STEP64_TABLE8_LR_NCDEG_INTEGRATION_SUMMARY.csv
|   |   └── STEP64_TABLE9_REPRODUCIBILITY_STATUS.csv
|   |
|   └── qc/
|       ├── STEP63C_CORRECTED_MASTER_AUDIT.csv
|       ├── STEP63C_CORRECTED_CRITICAL_CHECKS.csv
|       ├── STEP63C_CORRECTED_FAILED_AUDIT_CHECKS.csv
|       └── STEP63C_CORRECTED_FAILED_CRITICAL_CHECKS.csv
|
└── docs/
Final Result Files
Dataset Recovery
STEP64_DATASET_RECOVERY_REPORT.csv

Contains dataset/sample recovery and verification information.

Final Pipeline Summary
STEP64_TABLE1_PIPELINE_SUMMARY.csv

Contains the quantitative summary of the complete analysis.

Final LR Candidates
STEP64_TABLE2_FINAL_LR_CANDIDATES.csv

Contains the final niche-related ligand-receptor candidates.

LR Pathway Characterization
STEP64_TABLE3_LR_PATHWAY_CHARACTERIZATION.csv

Contains pathway-level characterization of final LR candidates.

C-SIDE Branch Summary
STEP64_TABLE4_CSIDE_BRANCH_SUMMARY.csv

Summarizes C-SIDE results for MLCP and SLCP branches.

NCDEG Direction Summary
STEP64_TABLE5_NCDEG_DIRECTION_SUMMARY.csv

Summarizes higher- and lower-expression NCDEGs.

Enrichment Term Recurrence
STEP64_TABLE6_ENRICHMENT_TERM_RECURRENCE.csv

Summarizes recurrence of enriched biological terms.

MLCP-SLCP Comparison
STEP64_TABLE7_MLCP_SLCP_COMPARISON.csv

Compares enrichment patterns between MLCP and SLCP branches.

LR-NCDEG Integration
STEP64_TABLE8_LR_NCDEG_INTEGRATION_SUMMARY.csv

Summarizes the LR-NCDEG branch-level integration.

Reproducibility Status
STEP64_TABLE9_REPRODUCIBILITY_STATUS.csv

Contains the final reproducibility status.

QC and Audit Files

The final quality-control directory contains:

STEP63C_CORRECTED_MASTER_AUDIT.csv
STEP63C_CORRECTED_CRITICAL_CHECKS.csv
STEP63C_CORRECTED_FAILED_AUDIT_CHECKS.csv
STEP63C_CORRECTED_FAILED_CRITICAL_CHECKS.csv

These files document the final validation and audit status of the computational workflow.

Reproduction Strategy

The analysis is organized as a sequential workflow:

1. Dataset recovery and metadata
2. Spatial QC
3. Spatial coordinate recovery
4. RCTD reference preparation
5. RCTD cell-type deconvolution
6. Tissue compartment classification
7. MLCP / SLCP analysis
8. Spatial niche analysis
9. Ligand-receptor analysis
10. LR meta-analysis
11. C-SIDE analysis
12. NCDEG extraction
13. NCDEG meta-analysis
14. NCDEG filtering
15. Gene-set enrichment
16. MLCP versus SLCP comparison
17. LR-NCDEG integration
18. Final result generation
19. Reproducibility audit

Scripts should be executed in the documented order because downstream stages depend on outputs from earlier stages.

Reproducibility Principles

This project follows the following reproducibility principles:

Stage-wise analysis

The workflow is divided into independent computational stages rather than one monolithic script.

Incremental result saving

Intermediate and final results are saved so that individual stages can be independently inspected.

Explicit sample tracking

All 43 GSE210616 spatial sections are explicitly tracked.

Documented cell-type reference

RCTD uses a documented TNBC-focused reference derived from GSE176078.

Structured result tables

Major analysis stages generate structured CSV outputs.

Independent auditing

A dedicated audit framework verifies expected counts and relationships between analysis stages.

Version tracking

Major R and package versions are documented.

GitHub organization

Scripts, metadata, documentation, final tables, and QC results are separated into dedicated directories.

Interpretation of the Reproduction

The primary goal of this project is methodological reproducibility, not numerical replication of the complete pan-cancer study.

The reference study analyzed:

373 samples
12 cancer types
56 LCPs
13 recurrent spatial niches

The present project analyzes:

43 spatial sections
22 patients
1 cancer type
TNBC

Therefore, differences in:

Number of LCPs
Number of spatial niches
LR candidates
NCDEG counts
Enrichment terms
Spatial patterns

are expected because the datasets and cohort structures are different.

The appropriate interpretation is that the analytical framework and computational methodology of the reference study were reproduced and adapted for the TNBC GSE210616 cohort.

Scientific Scope

The project evaluates the TNBC tumor spatial microenvironment at multiple biological levels:

Cellular composition
        |
        v
Tissue compartments
        |
        v
Local cellular programs
        |
        v
Spatial niches
        |
        v
Cell-cell communication
        |
        v
Niche-associated gene expression
        |
        v
Biological pathways
        |
        v
Integrated spatial evidence

This provides a computational framework for investigating how cellular composition, spatial organization, intercellular signaling, and niche-associated transcriptional states interact within TNBC tissue.

Important Limitations
TNBC-only cohort

The current analysis is restricted to GSE210616 and therefore does not reproduce the complete pan-cancer cohort.

Dataset-specific LCPs and niches

LCP and niche definitions depend on the dataset, spatial structure, clustering, and available cellular composition.

Therefore, TNBC-derived LCPs and niches should not automatically be interpreted as identical to the pan-cancer LCPs and niches reported by Li et al.

Direct LR-NCDEG overlap

The final analysis identified:

0 direct LR-NCDEG gene overlap

Therefore, the integration should not be interpreted as direct molecular interaction evidence.

Exact niche-context overlap

The final analysis identified:

0 exact shared LR-NCDEG niche contexts

Consequently, integration was interpreted at the broader MLCP/SLCP branch level.

Experimental validation

This repository contains computational analyses and does not provide independent experimental validation of biological mechanisms.

Project Status
Completed
[✓] GSE210616 dataset processing
[✓] Spatial QC
[✓] Spatial coordinate recovery
[✓] TNBC single-cell reference preparation
[✓] RCTD deconvolution
[✓] 43 / 43 spatial sections analyzed
[✓] Tissue compartment analysis
[✓] MLCP / SLCP analysis
[✓] Spatial niche analysis
[✓] Ligand-receptor analysis
[✓] LR meta-analysis
[✓] C-SIDE analysis
[✓] NCDEG extraction
[✓] NCDEG meta-analysis
[✓] NCDEG filtering
[✓] Gene-set enrichment
[✓] MLCP / SLCP comparison
[✓] LR-NCDEG integration
[✓] Final result tables
[✓] Reproducibility audit
[✓] GitHub repository organization
Final validation
Critical audit: 23 / 23 PASS
Full audit:     51 / 51 PASS
Relationship to the Reference Paper

The reference paper provides the methodological framework.

This project provides a TNBC-specific implementation and adaptation of that framework.

Li et al. 2026
Pan-cancer spatial transcriptomics methodology
                |
                | methodological reproduction
                v
       GSE210616 TNBC cohort
                |
                +--> RCTD
                |
                +--> Tissue compartments
                |
                +--> MLCP / SLCP
                |
                +--> Spatial niches
                |
                +--> LR interactions
                |
                +--> C-SIDE / NCDEGs
                |
                +--> Enrichment
                |
                +--> LR-NCDEG integration
                |
                +--> Reproducibility audit

The final statistics and biological findings reported in this repository are specific to the GSE210616 TNBC cohort and are not the original pan-cancer results reported by Li et al.

Suggested Academic Project Description

TNBC-focused computational reproduction and methodological adaptation of a published spatial transcriptomics framework using GSE210616, integrating RCTD-based cell-type deconvolution, tissue compartment classification, local cellular programs, spatial niche analysis, ligand-receptor signaling, C-SIDE/NCDEG analysis, gene-set enrichment, LR-NCDEG integration, and reproducibility auditing.

Keywords
Bioinformatics
Spatial Transcriptomics
Triple-Negative Breast Cancer
TNBC
GSE210616
GSE176078
10x Visium
RCTD
spacexr
Seurat
CellChat
inferCNV
C-SIDE
NCDEG
Local Cellular Programs
MLCP
SLCP
Spatial Niches
Ligand-Receptor
Tumor Microenvironment
Tumor Spatial Microenvironment
Single-Cell RNA-seq
Transcriptomics
Cancer Bioinformatics
Computational Biology
Reproducible Research
Reference

Li J, Lin P, Wang H, Tang Z, Yan X, Chen X, Yuan J, Chen W, Li H.

Pan-cancer analysis of spatial transcriptomics reveals heterogeneous tumor spatial microenvironment.

Cell Reports Medicine. 2026;7:102751.

DOI: 10.1016/j.xcrm.2026.102751

Author

Sravani Thota

M.Sc. Bioinformatics

Spatial Transcriptomics | Cancer Bioinformatics | Computational Biology

Repository Purpose

This repository provides a transparent, structured, and auditable computational implementation of a published spatial transcriptomics methodology adapted to a TNBC cohort.

The project is intended to support:

Computational reproducibility
Methodological learning
Spatial transcriptomics research
TNBC tumor microenvironment analysis
Bioinformatics portfolio development
Future extension to additional cancer cohorts
