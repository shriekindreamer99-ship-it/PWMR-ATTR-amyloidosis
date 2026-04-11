# PWMR-ATTR-amyloidosis

Analytical code for:

**"Cross-Phenotypic Concordance Filtering for Proteome-Wide Mendelian Randomization: Empirical Calibration of Power Requirements Across Cardiometabolic and Rare-Disease Trajectories"**

Jincheng Xing, Jiaxu Shen

## Overview

This repository contains R scripts for the Phenotypic Trajectory Mendelian Randomization (PT-MR) framework applied to ATTR amyloidosis and cardiometabolic positive control trajectories.

## Repository Structure

```
├── scripts/
│   ├── 01_PWMR_pipeline.R              # Core PWMR pipeline (Sun et al. 2018 pQTLs → systemic amyloidosis)
│   ├── 02_replication_stage.R           # FinnGen R12 replication lookup
│   ├── 03_triangulation_PWMR.R         # Cross-phenotype concordance filtering (CTS × Amyloidosis × HF)
│   ├── 07_sensitivity_analysis.R        # Steiger, heterogeneity, pleiotropy, multi-method MR
│   ├── 08_reverse_mr.R                 # Reverse MR (disease → protein)
│   ├── 13b_locus_level_permutation.R   # Locus-level permutation test (ATTR trajectory)
│   ├── 13c_locus_permutation_3loci.R   # Corrected permutation (positive control, CA5A/rs662 LD merged)
│   ├── 14_coloc_bag3_hf.R             # BAG3-HF colocalization (PP.H4 = 0.992)
│   ├── 15_negative_control.R           # Negative control analysis (T2D + asthma)
│   ├── 30_multi_trajectory_calibration.R # Multi-trajectory calibration (4 comparison trajectories)
│   ├── 38_power_degradation_3loci.R    # SE-inflation power simulation (corrected 3 loci)
│   ├── positive_control_ObesityT2DHF.R # Positive control trajectory (BMI → T2D → HF)
│   ├── coloc_pon1_ukbppp.R            # PON1 cross-platform replication (UKB-PPP Olink)
│   ├── ukbppp_mr_supplement.R          # UKB-PPP supplementary MR analyses
│   ├── export_supp_table_v4.R          # Supplementary Table S1 generation
│   └── ...
├── results/                            # Key intermediate CSV outputs
└── README.md
```

## Data Sources

| Dataset | Source | Access |
|---------|--------|--------|
| Plasma pQTLs (2,094 proteins) | Sun et al. 2018 | OpenGWAS (prot-a-*) |
| CTS GWAS | UK Biobank | OpenGWAS |
| Systemic amyloidosis | FinnGen R5 (E4_AMYLOIDOSIS) | finngen.fi |
| Heart failure | HERMES Consortium | CVDKP |
| T2D (negative control) | DIAGRAM | OpenGWAS (ebi-a-GCST006867) |
| Asthma (negative control) | UK Biobank/Neale lab | OpenGWAS (ukb-d-J10_ASTHMA) |
| Cross-platform replication | UKB-PPP (Olink) | Synapse (syn51364943) |
| CAD | CARDIoGRAMplusC4D | OpenGWAS (ebi-a-GCST003116) |
| CKD | CKDGen | OpenGWAS (ieu-a-1102) |
| LDL cholesterol | GLGC | OpenGWAS (ieu-a-300) |
| SBP | ICBP | OpenGWAS (ieu-b-38) |
| Stroke | MEGASTROKE | OpenGWAS (ebi-a-GCST006906) |

## Key Dependencies

```r
install.packages(c("dplyr", "ggplot2", "openxlsx"))
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("MRCIEU/ieugwasr")
install.packages("coloc")
```

## Reproducibility Notes

- All analyses use publicly available summary-level GWAS data
- OpenGWAS API access requires `ieugwasr::get_opengwas_jwt()` for authentication
- UKB-PPP data requires Synapse access approval
- Random seeds are reported for all permutation analyses (see manuscript Methods)

## Citation

If you use this code, please cite the manuscript (submitted to *Human Genomics*).

## License

MIT
