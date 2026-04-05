# PWMR-ATTR: Cross-Phenotypic Concordance Filtering for Proteome-Wide Mendelian Randomization

## Overview

This repository contains the analytical code and curated derived results for:

> **Cross-Phenotypic Concordance Filtering for Proteome-Wide Mendelian Randomization: Empirical Calibration of Power Requirements Across Cardiometabolic and Rare-Disease Trajectories**

## Repository Structure

```
scripts/
  01_PWMR_pipeline.R                  # Main PWMR pipeline (2,094 proteins × 3 endpoints)
  02_replication_stage.R               # FinnGen R12 replication lookup
  03_triangulation_PWMR.R             # Cross-phenotypic concordance filtering
  07_sensitivity_analysis.R            # Steiger directionality, heterogeneity
  08_reverse_mr.R                      # Reverse MR analysis
  13b_locus_level_permutation.R        # Locus-level permutation testing
  15_negative_control.R                # T2D + asthma negative control analysis
  positive_control_ObesityT2DHF.R      # BMI → T2D → HF positive control trajectory
  downsampling_experiment.R            # SE-inflation power-degradation simulation
  R12_full_PWMR.R                      # Full PWMR against FinnGen R12
  alternative_permutation_sensitivity.R # Alternative anchor permutation schemes
  14_coloc_bag3_hf.R                   # BAG3–HF colocalization (coloc.abf)
  coloc_pon1_ukbppp.R                  # PON1 UKB-PPP cross-platform validation
  ukbppp_sensitivity.R                 # UKB-PPP Olink cross-platform replication
  ukbppp_mr_supplement.R               # UKB-PPP supplementary MR
  posctrl_sensitivity.R                # Positive control sensitivity analyses
  palette.R                            # Unified macaron color palette
  fix_figure_overlaps.R                # Figure generation (final versions)
  generate_missing_figures.R           # Concordance/permutation figure generation
  export_supp_table_v4.R               # Supplementary table formatting

results/
  (Derived MR results, permutation outputs, and intermediate data files)

```

## Data Sources

- **Plasma pQTL instruments**: Sun et al. (2018) via OpenGWAS (prot-a-*)
- **CTS GWAS**: UK Biobank via OpenGWAS (ukb-d-G6_CARPTU)
- **Systemic amyloidosis GWAS**: FinnGen R5 (E4_AMYLOIDOSIS, 226 cases)
- **Heart failure GWAS**: HERMES Consortium (ebi-a-GCST009541)
- **FinnGen R12 replication**: E4_AMYLNAS (573 cases)
- **HERMES 2024 (colocalization)**: European HFall meta-analysis
- **UKB-PPP Olink**: Synapse syn51364943 (European discovery pGWAS)
- **BMI**: UK Biobank Neale lab (ieu-b-40)
- **T2D**: DIAGRAM consortium (ebi-a-GCST006867)
- **Asthma**: UK Biobank/Neale lab (ukb-d-J10_ASTHMA)

## Requirements

- R 4.5.3
- TwoSampleMR (v0.7.0)
- ieugwasr (with JWT authentication)
- coloc
- ggplot2, dplyr, tidyr, patchwork, ggrepel, scales

## License

[To be specified]
