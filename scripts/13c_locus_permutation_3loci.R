# ============================================================
# 13c_locus_permutation_3loci.R
# Corrected locus-level permutation test for positive control
# 
# Correction: CA5A instrument rs1157745 (chr7:94,940,235) is
# 2.8 kb from rs662 (chr7:94,937,434) with r² = 1.0 in EUR.
# These are collapsed into a single locus (PON1/rs662).
# Result: 5 proteins / 3 independent loci (was 4).
#
# Manuscript ref: Results, "Positive control trajectory..."
# Output: results/13c_locus_permutation_3loci.csv
#         results/13c_null_counts_3loci.rds
# ============================================================

library(dplyr)

# ---- 1. Load positive control MR results ----
bmi <- read.csv("results/PosCtrl_BMI_all_results.csv", stringsAsFactors = FALSE)
t2d <- read.csv("results/19_mr_t2d_results.csv", stringsAsFactors = FALSE)
hf  <- read.csv("results/10_mr_hf_results.csv", stringsAsFactors = FALSE)

# Load instrument data for lead SNP identification
all_exposures <- readRDS("results/01_extraction_progress.rds")

# ---- 2. Build three-endpoint intersection ----
clean <- function(df, suffix) {
  df %>%
    mutate(id.exposure = as.character(id.exposure)) %>%
    select(id.exposure, b, se, pval) %>%
    rename_with(~ paste0(., "_", suffix), c(b, se, pval))
}

tri <- clean(bmi, "bmi") %>%
  inner_join(clean(t2d, "t2d"), by = "id.exposure") %>%
  inner_join(clean(hf, "hf"), by = "id.exposure")

cat("Proteins with results across all 3 endpoints:", nrow(tri), "\n")

# ---- 3. Identify lead SNPs ----
lead_snps <- all_exposures %>%
  filter(id.exposure %in% tri$id.exposure) %>%
  group_by(id.exposure) %>%
  arrange(pval.exposure) %>%
  slice(1) %>%
  ungroup() %>%
  select(id.exposure, lead_snp = SNP)

tri <- tri %>% left_join(lead_snps, by = "id.exposure")

# ---- 4. CRITICAL CORRECTION: Collapse CA5A into PON1 locus ----
# rs1157745 (CA5A) and rs662 (TNS4/STX10) are in perfect LD (r² = 1.0)
# Both map to the PON1 region at chr7q21.3 (2.8 kb apart)
tri <- tri %>%
  mutate(locus_snp = ifelse(lead_snp == "rs1157745", "rs662", lead_snp))

cat("Unique locus SNPs after LD correction:", 
    length(unique(tri$locus_snp)), "\n")

# ---- 5. Identify concordant proteins ----
tri <- tri %>%
  mutate(
    sig_bmi = pval_bmi < 0.05,
    sig_t2d = pval_t2d < 0.05,
    sig_hf  = pval_hf  < 0.05,
    dir_consistent = (sign(b_bmi) == sign(b_t2d)) & 
                     (sign(b_bmi) == sign(b_hf)),
    concordant = sig_bmi & sig_t2d & sig_hf & dir_consistent
  )

obs_proteins <- sum(tri$concordant)
obs_loci <- length(unique(tri$locus_snp[tri$concordant]))
cat("Concordant proteins:", obs_proteins, "\n")
cat("Independent loci:", obs_loci, "\n")

# ---- 6. Locus-level permutation test (10,000 iterations) ----
n_perm <- 10000
seeds <- c(42, 123, 456, 789, 2024)

results_all <- data.frame()

for (s in seeds) {
  set.seed(s)
  null_counts <- numeric(n_perm)
  
  for (i in seq_len(n_perm)) {
    # Shuffle BMI and HF P-values and betas; keep T2D fixed (anchor)
    perm_p_bmi <- sample(tri$pval_bmi)
    perm_b_bmi <- sample(tri$b_bmi)
    perm_p_hf  <- sample(tri$pval_hf)
    perm_b_hf  <- sample(tri$b_hf)
    
    perm_sig <- (perm_p_bmi < 0.05) & (tri$pval_t2d < 0.05) & (perm_p_hf < 0.05)
    perm_dir <- (sign(perm_b_bmi) == sign(tri$b_t2d)) & 
                (sign(perm_b_bmi) == sign(perm_b_hf)) &
                (sign(tri$b_t2d) == sign(perm_b_hf))
    perm_conc <- perm_sig & perm_dir
    
    if (any(perm_conc)) {
      null_counts[i] <- length(unique(tri$locus_snp[perm_conc]))
    } else {
      null_counts[i] <- 0
    }
  }
  
  emp_p <- mean(null_counts >= obs_loci)
  
  results_all <- rbind(results_all, data.frame(
    seed = s,
    n_perm = n_perm,
    observed_loci = obs_loci,
    null_mean = mean(null_counts),
    null_median = median(null_counts),
    null_max = max(null_counts),
    ge_observed = sum(null_counts >= obs_loci),
    empirical_p = emp_p
  ))
  
  cat(sprintf("Seed %d: P = %.4f\n", s, emp_p))
  
  # Save null distribution for first seed (used for Fig 2C histogram)
  if (s == seeds[1]) {
    saveRDS(null_counts, "results/13c_null_counts_3loci.rds")
  }
}

# ---- 7. Save results ----
write.csv(results_all, "results/13c_locus_permutation_3loci.csv", row.names = FALSE)

cat("\n=== Summary ===\n")
cat("Observed concordant loci:", obs_loci, "\n")
cat("Empirical P range:", sprintf("%.3f-%.3f", min(results_all$empirical_p), 
                                   max(results_all$empirical_p)), "\n")
cat("Concordant proteins:", paste(tri$id.exposure[tri$concordant], collapse = ", "), "\n")
