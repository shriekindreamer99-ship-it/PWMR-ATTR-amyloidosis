# ============================================================
# 38_power_degradation_3loci.R
# Power-degradation simulation with corrected locus count
#
# Simulates progressive sample-size reduction in the positive
# control trajectory by inflating standard errors. Uses the
# corrected locus collapsing (CA5A/rs1157745 merged with
# rs662 at PON1 locus, r² = 1.0).
#
# Key finding: concordant hits drop from 5 proteins / 3 loci
# at full power to 0 at 50% power (Neff ≈ 90,000).
#
# Manuscript ref: Results "Power-degradation simulation..."
#                 Fig 3A, 3B
# Output: results/38_power_degradation_3loci.csv
# ============================================================

library(dplyr)

# ---- 1. Load data ----
bmi <- read.csv("results/PosCtrl_BMI_all_results.csv", stringsAsFactors = FALSE)
t2d <- read.csv("results/19_mr_t2d_results.csv", stringsAsFactors = FALSE)
hf  <- read.csv("results/10_mr_hf_results.csv", stringsAsFactors = FALSE)
all_exposures <- readRDS("results/01_extraction_progress.rds")

clean <- function(df, suffix) {
  df %>%
    mutate(id.exposure = as.character(id.exposure),
           protein_name = sub(" \\|\\|.*", "", 
                              ifelse("protein_name" %in% names(.), protein_name, exposure))) %>%
    select(id.exposure, protein_name, b, se, pval) %>%
    rename_with(~ paste0(., "_", suffix), c(b, se, pval))
}

tri <- clean(bmi, "bmi") %>%
  inner_join(clean(t2d, "t2d") %>% select(-protein_name), by = "id.exposure") %>%
  inner_join(clean(hf, "hf") %>% select(-protein_name), by = "id.exposure")

cat("Three-node intersection:", nrow(tri), "proteins\n")

# ---- 2. Lead SNP + corrected locus collapsing ----
lead_snps <- all_exposures %>%
  filter(id.exposure %in% tri$id.exposure) %>%
  group_by(id.exposure) %>%
  arrange(pval.exposure) %>%
  slice(1) %>%
  ungroup() %>%
  select(id.exposure, lead_snp = SNP)

tri <- tri %>% left_join(lead_snps, by = "id.exposure")

# CRITICAL: rs1157745 (CA5A) → rs662 (PON1 locus)
tri <- tri %>%
  mutate(locus_snp = ifelse(lead_snp == "rs1157745", "rs662", lead_snp))

# ---- 3. Full-power verification ----
tri_full <- tri %>%
  mutate(
    concordant = (pval_bmi < 0.05) & (pval_t2d < 0.05) & (pval_hf < 0.05) &
      (sign(b_bmi) == sign(b_t2d)) & (sign(b_bmi) == sign(b_hf))
  )

cat("Full power: ", sum(tri_full$concordant), " proteins, ",
    length(unique(tri_full$locus_snp[tri_full$concordant])), " loci\n")

# ---- 4. SE inflation simulation ----
fractions <- c(1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005)
neff_hf_full <- 180076  # 4 * 47309 * (47309+309910) / (2*47309+309910+309910)... approx
n_perm <- 1000
set.seed(42)

results <- data.frame()

for (f in fractions) {
  cat(sprintf("\nFraction = %.3f (Neff_HF ≈ %d)...\n", f, round(neff_hf_full * f)))
  
  # Inflate SEs for T2D and HF (simulating smaller GWAS)
  tri_sim <- tri %>%
    mutate(
      se_t2d_infl = se_t2d / sqrt(f),
      se_hf_infl  = se_hf  / sqrt(f),
      pval_t2d_infl = 2 * pnorm(-abs(b_t2d / se_t2d_infl)),
      pval_hf_infl  = 2 * pnorm(-abs(b_hf  / se_hf_infl)),
      concordant = (pval_bmi < 0.05) & (pval_t2d_infl < 0.05) & (pval_hf_infl < 0.05) &
        (sign(b_bmi) == sign(b_t2d)) & (sign(b_bmi) == sign(b_hf))
    )
  
  n_prot <- sum(tri_sim$concordant)
  n_loci <- length(unique(tri_sim$locus_snp[tri_sim$concordant]))
  
  # Permutation test (only if there are concordant loci)
  if (n_loci > 0) {
    null_counts <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
      pp1 <- sample(tri_sim$pval_bmi);  pb1 <- sample(tri_sim$b_bmi)
      pp3 <- sample(tri_sim$pval_hf_infl); pb3 <- sample(tri_sim$b_hf)
      ps <- (pp1 < 0.05) & (tri_sim$pval_t2d_infl < 0.05) & (pp3 < 0.05)
      pd <- (sign(pb1) == sign(tri_sim$b_t2d)) & 
            (sign(pb1) == sign(pb3)) & 
            (sign(tri_sim$b_t2d) == sign(pb3))
      pc <- ps & pd
      null_counts[i] <- if (any(pc)) length(unique(tri_sim$locus_snp[pc])) else 0
    }
    emp_p <- mean(null_counts >= n_loci)
  } else {
    emp_p <- NA
  }
  
  cat(sprintf("  Proteins: %d, Loci: %d, P: %s\n",
              n_prot, n_loci, ifelse(is.na(emp_p), "NA", sprintf("%.4f", emp_p))))
  
  results <- rbind(results, data.frame(
    fraction = f,
    neff_hf = round(neff_hf_full * f),
    n_concordant_proteins = n_prot,
    n_concordant_loci = n_loci,
    empirical_p = emp_p
  ))
}

# ---- 5. Save ----
write.csv(results, "results/38_power_degradation_3loci.csv", row.names = FALSE)

cat("\n=== Power Degradation Results (3 loci) ===\n")
print(results)
cat("\nKey: Concordant signals vanish below Neff ≈ 90,000\n")
