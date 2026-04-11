# ============================================================
# 30_multi_trajectory_calibration.R
# Multi-trajectory calibration analysis for PT-MR framework
#
# Runs 4 additional well-powered trajectories to contextualize
# the positive control permutation P = 0.059 with a specificity
# gradient. All use Sun et al. 2018 (2,094 proteins).
#
# Trajectories:
#   BMI → T2D → CKD  (CKDGen; Wuttke et al. 2019)
#   BMI → CAD → HF   (CARDIoGRAMplusC4D; Nikpay et al. 2015)
#   LDL → CAD → HF   (GLGC; Willer et al. 2013)
#   BMI → SBP → Stroke (ICBP; Evangelou et al. 2018 + MEGASTROKE)
#
# Manuscript ref: Methods "Multi-Trajectory Calibration"
#                 Results "Multi-trajectory calibration..."
# Output: results/30_multi_trajectory_calibration.csv
# ============================================================

library(TwoSampleMR)
library(dplyr)

# ---- 0. Common setup ----
# Load pre-extracted BMI instruments (shared across trajectories)
bmi <- read.csv("results/PosCtrl_BMI_all_results.csv", stringsAsFactors = FALSE)
t2d <- read.csv("results/19_mr_t2d_results.csv", stringsAsFactors = FALSE)
hf  <- read.csv("results/10_mr_hf_results.csv", stringsAsFactors = FALSE)
all_exposures <- readRDS("results/01_extraction_progress.rds")

# Pre-extracted protein instruments
protein_ids <- unique(all_exposures$id.exposure)

clean <- function(df, suffix) {
  df %>%
    mutate(id.exposure = as.character(id.exposure)) %>%
    select(id.exposure, b, se, pval) %>%
    rename_with(~ paste0(., "_", suffix), c(b, se, pval))
}

get_lead_snps <- function(ids) {
  all_exposures %>%
    filter(id.exposure %in% ids) %>%
    group_by(id.exposure) %>%
    arrange(pval.exposure) %>%
    slice(1) %>%
    ungroup() %>%
    select(id.exposure, lead_snp = SNP) %>%
    mutate(locus_snp = ifelse(lead_snp == "rs1157745", "rs662", lead_snp))
}

run_concordance_permutation <- function(tri, n_perm = 10000, seed = 42) {
  leads <- get_lead_snps(tri$id.exposure)
  tri <- tri %>% left_join(leads, by = "id.exposure")
  
  tri <- tri %>% mutate(
    concordant = (pval_n1 < 0.05) & (pval_n2 < 0.05) & (pval_n3 < 0.05) &
      (sign(b_n1) == sign(b_n2)) & (sign(b_n1) == sign(b_n3))
  )
  
  obs_loci <- length(unique(tri$locus_snp[tri$concordant]))
  obs_prot <- sum(tri$concordant)
  
  if (obs_loci == 0) {
    return(data.frame(n_proteins = nrow(tri), concordant_proteins = 0,
                      concordant_loci = 0, empirical_p = 1.0))
  }
  
  set.seed(seed)
  null_counts <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    pp1 <- sample(tri$pval_n1); pb1 <- sample(tri$b_n1)
    pp3 <- sample(tri$pval_n3); pb3 <- sample(tri$b_n3)
    ps <- (pp1 < 0.05) & (tri$pval_n2 < 0.05) & (pp3 < 0.05)
    pd <- (sign(pb1) == sign(tri$b_n2)) & (sign(pb1) == sign(pb3))
    pc <- ps & pd
    null_counts[i] <- if (any(pc)) length(unique(tri$locus_snp[pc])) else 0
  }
  
  data.frame(n_proteins = nrow(tri), concordant_proteins = obs_prot,
             concordant_loci = obs_loci, empirical_p = mean(null_counts >= obs_loci))
}

# ---- 1. Extract outcome GWAS data ----
# Outcome IDs
outcomes <- list(
  CKD    = "ieu-a-1102",       # CKDGen consortium
  CAD    = "ebi-a-GCST003116", # CARDIoGRAMplusC4D
  LDL    = "ieu-a-300",        # GLGC
  SBP    = "ieu-b-38",         # ICBP
  Stroke = "ebi-a-GCST006906"  # MEGASTROKE
)

# Extract MR results for each outcome (batch processing)
extract_outcome_mr <- function(outcome_id, batch_size = 50) {
  all_results <- data.frame()
  batches <- split(protein_ids, ceiling(seq_along(protein_ids) / batch_size))
  
  for (j in seq_along(batches)) {
    cat(sprintf("  Batch %d/%d\n", j, length(batches)))
    tryCatch({
      exp_dat <- mv_extract_exposures(batches[[j]])
      out_dat <- extract_outcome_data(exp_dat$SNP, outcome_id)
      harm <- harmonise_data(exp_dat, out_dat)
      res <- mr(harm, method_list = c("mr_ivw", "mr_wald_ratio"))
      all_results <- rbind(all_results, res)
      Sys.sleep(3)
    }, error = function(e) {
      cat(sprintf("  Batch %d failed: %s\n", j, e$message))
    })
  }
  return(all_results)
}

cat("Extracting CKD results...\n")
mr_ckd <- extract_outcome_mr(outcomes$CKD)
write.csv(mr_ckd, "results/30_mr_ckd_results.csv", row.names = FALSE)

cat("Extracting CAD results...\n")
mr_cad <- extract_outcome_mr(outcomes$CAD)
write.csv(mr_cad, "results/30_mr_cad_results.csv", row.names = FALSE)

cat("Extracting LDL results...\n")
mr_ldl <- extract_outcome_mr(outcomes$LDL)  # Note: LDL as first node, not protein→LDL
write.csv(mr_ldl, "results/30_mr_ldl_results.csv", row.names = FALSE)

cat("Extracting SBP results...\n")
mr_sbp <- extract_outcome_mr(outcomes$SBP)
write.csv(mr_sbp, "results/30_mr_sbp_results.csv", row.names = FALSE)

cat("Extracting Stroke results...\n")
mr_stroke <- extract_outcome_mr(outcomes$Stroke)
write.csv(mr_stroke, "results/30_mr_stroke_results.csv", row.names = FALSE)

# ---- 2. Run concordance + permutation for each trajectory ----
results <- data.frame()

# Trajectory 1: BMI → T2D → CKD
cat("\n=== BMI → T2D → CKD ===\n")
tri1 <- clean(bmi, "n1") %>%
  inner_join(clean(t2d, "n2"), by = "id.exposure") %>%
  inner_join(clean(mr_ckd, "n3"), by = "id.exposure")
r1 <- run_concordance_permutation(tri1)
r1$trajectory <- "BMI → T2D → CKD"
results <- rbind(results, r1)

# Trajectory 2: BMI → CAD → HF
cat("\n=== BMI → CAD → HF ===\n")
tri2 <- clean(bmi, "n1") %>%
  inner_join(clean(mr_cad, "n2"), by = "id.exposure") %>%
  inner_join(clean(hf, "n3"), by = "id.exposure")
r2 <- run_concordance_permutation(tri2)
r2$trajectory <- "BMI → CAD → HF"
results <- rbind(results, r2)

# Trajectory 3: LDL → CAD → HF
cat("\n=== LDL → CAD → HF ===\n")
tri3 <- clean(mr_ldl, "n1") %>%
  inner_join(clean(mr_cad, "n2"), by = "id.exposure") %>%
  inner_join(clean(hf, "n3"), by = "id.exposure")
r3 <- run_concordance_permutation(tri3)
r3$trajectory <- "LDL → CAD → HF"
results <- rbind(results, r3)

# Trajectory 4: BMI → SBP → Stroke
cat("\n=== BMI → SBP → Stroke ===\n")
tri4 <- clean(bmi, "n1") %>%
  inner_join(clean(mr_sbp, "n2"), by = "id.exposure") %>%
  inner_join(clean(mr_stroke, "n3"), by = "id.exposure")
r4 <- run_concordance_permutation(tri4)
r4$trajectory <- "BMI → SBP → Stroke"
results <- rbind(results, r4)

# ---- 3. Save results ----
write.csv(results, "results/30_multi_trajectory_calibration.csv", row.names = FALSE)

cat("\n=== Multi-Trajectory Calibration Results ===\n")
print(results[, c("trajectory", "n_proteins", "concordant_loci", "empirical_p")])
