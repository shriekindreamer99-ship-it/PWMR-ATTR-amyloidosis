#!/usr/bin/env Rscript
# ============================================================
# Sensitivity Analyses for Positive Control Concordant Hits
# 
# Checks: shared IVs, locus collapsing, permutation, Steiger,
#         reverse MR, F-statistics, heterogeneity
#
# Run: Rscript posctrl_sensitivity.R > sensitivity_log.txt 2>&1 &
# ============================================================

cat("=== Positive Control Sensitivity Analyses ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

library(TwoSampleMR)

# The 5 concordant proteins
hits <- data.frame(
  id = c("prot-a-332", "prot-a-3070", "prot-a-2883", "prot-a-1677", "prot-a-2231"),
  name = c("Carbonic anhydrase 5A", "Tensin-4", "Syntaxin-10", 
           "Importin subunit alpha-1", "PDGFRL"),
  stringsAsFactors = FALSE
)

# ============================================================
# 1. SHARED INSTRUMENT CHECK & LOCUS COLLAPSING
# ============================================================
cat("=== 1. Shared Instrument Check ===\n\n")

all_instruments <- data.frame()
for (i in 1:nrow(hits)) {
  exp <- extract_instruments(hits$id[i], p1 = 5e-8, clump = TRUE,
                              r2 = 0.001, kb = 10000)
  if (!is.null(exp) && nrow(exp) > 0) {
    exp$protein_name <- hits$name[i]
    all_instruments <- rbind(all_instruments, exp)
    cat(sprintf("  %s: %d IVs — lead SNP(s): %s\n", 
                hits$name[i], nrow(exp), paste(exp$SNP, collapse = ", ")))
  }
}

# Check for shared SNPs
snp_table <- table(all_instruments$SNP)
shared_snps <- names(snp_table[snp_table > 1])

if (length(shared_snps) > 0) {
  cat(sprintf("\n  *** SHARED IVs DETECTED: %d SNP(s) shared across proteins ***\n", length(shared_snps)))
  for (s in shared_snps) {
    proteins_sharing <- unique(all_instruments$protein_name[all_instruments$SNP == s])
    cat(sprintf("    %s → shared by: %s\n", s, paste(proteins_sharing, collapse = ", ")))
  }
  
  # Locus collapsing
  cat("\n  Locus-level interpretation: proteins sharing IVs = single locus unit\n")
  
  # Group by lead SNP
  locus_groups <- split(all_instruments, all_instruments$SNP)
  n_independent_loci <- length(unique(all_instruments$SNP))
  cat(sprintf("  %d proteins → %d independent locus-level signals\n\n", 
              nrow(hits), n_independent_loci))
} else {
  cat("\n  No shared IVs. All 5 signals are independently instrumented.\n\n")
}

# Save instrument details
write.csv(all_instruments[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                               "beta.exposure", "se.exposure", "pval.exposure", 
                               "id.exposure", "protein_name")],
          "PosCtrl_instruments.csv", row.names = FALSE)

# ============================================================
# 2. F-STATISTICS
# ============================================================
cat("=== 2. F-statistics ===\n\n")

for (i in 1:nrow(all_instruments)) {
  r <- all_instruments[i, ]
  N <- 3301  # Sun et al. sample size
  beta <- r$beta.exposure
  se <- r$se.exposure
  F_stat <- beta^2 / se^2
  R2 <- F_stat / (F_stat + N - 2)
  cat(sprintf("  %s | %s: F=%.1f, R²=%.4f %s\n",
              r$protein_name, r$SNP, F_stat, R2,
              ifelse(F_stat > 10, "[OK]", "[WEAK]")))
}

# ============================================================
# 3. STEIGER DIRECTIONALITY
# ============================================================
cat("\n=== 3. Steiger Directionality ===\n\n")

bmi_id <- "ieu-b-40"

for (i in 1:nrow(hits)) {
  tryCatch({
    exp <- extract_instruments(hits$id[i], p1 = 5e-8, clump = TRUE,
                                r2 = 0.001, kb = 10000)
    if (is.null(exp) || nrow(exp) == 0) next
    
    out <- extract_outcome_data(exp$SNP, bmi_id)
    if (is.null(out) || nrow(out) == 0) next
    
    harm <- harmonise_data(exp, out, action = 2)
    if (is.null(harm) || nrow(harm) == 0) next
    
    steiger <- directionality_test(harm)
    cat(sprintf("  %s: correct_direction=%s, steiger_pval=%.2e\n",
                hits$name[i], steiger$correct_causal_direction, steiger$steiger_pval))
    
  }, error = function(e) {
    cat(sprintf("  %s: Steiger failed (%s)\n", hits$name[i], e$message))
  })
  
  Sys.sleep(1)
}

# ============================================================
# 4. REVERSE MR (BMI → protein)
# ============================================================
cat("\n=== 4. Reverse MR (BMI → each protein) ===\n\n")

tryCatch({
  bmi_instruments <- extract_instruments(bmi_id, p1 = 5e-8, clump = TRUE,
                                          r2 = 0.001, kb = 10000)
  cat(sprintf("  BMI instruments: %d SNPs\n\n", nrow(bmi_instruments)))
  
  for (i in 1:nrow(hits)) {
    tryCatch({
      out <- extract_outcome_data(bmi_instruments$SNP, hits$id[i])
      if (is.null(out) || nrow(out) == 0) {
        cat(sprintf("  %s: no outcome data\n", hits$name[i]))
        next
      }
      
      harm <- harmonise_data(bmi_instruments, out, action = 2)
      mr_res <- mr(harm, method_list = c("mr_ivw"))
      
      cat(sprintf("  BMI → %s: beta=%.4f, P=%.2e %s\n",
                  hits$name[i], mr_res$b, mr_res$pval,
                  ifelse(mr_res$pval < 0.05, "[REVERSE SIGNAL]", "[no reverse]")))
      
    }, error = function(e) {
      cat(sprintf("  %s: reverse MR failed (%s)\n", hits$name[i], e$message))
    })
    Sys.sleep(1)
  }
}, error = function(e) {
  cat(sprintf("  BMI instrument extraction failed: %s\n", e$message))
})

# ============================================================
# 5. PERMUTATION TEST (locus-level, same as ATTR analysis)
# ============================================================
cat("\n=== 5. Locus-level Permutation Test ===\n\n")

# Load full BMI primary results
bmi_primary <- read.csv("PosCtrl_BMI_all_results.csv", stringsAsFactors = FALSE)
bmi_primary <- bmi_primary[bmi_primary$method %in% c("Inverse variance weighted", "Wald ratio"), ]

# Load S1 for T2D and HF
library(openxlsx)
s1_paths <- c("Supplementary_Table_S1_fixed.xlsx",
              "~/Desktop/PWMR/Supplementary_Table_S1_fixed.xlsx")
s1 <- NULL
for (p in s1_paths) {
  if (file.exists(p)) {
    s1 <- read.xlsx(p, startRow = 2)
    break
  }
}

if (!is.null(s1)) {
  # Merge BMI with S1 T2D/HF
  merged <- merge(
    bmi_primary[, c("id.exposure", "b", "pval")],
    s1[, c("OpenGWAS.ID", "Beta.(T2D)", "P.(T2D)", "Beta.(HF)", "P.(HF)")],
    by.x = "id.exposure", by.y = "OpenGWAS.ID"
  )
  names(merged) <- c("id", "bmi_beta", "bmi_p", "t2d_beta", "t2d_p", "hf_beta", "hf_p")
  merged$t2d_beta <- as.numeric(merged$t2d_beta)
  merged$t2d_p <- as.numeric(merged$t2d_p)
  merged$hf_beta <- as.numeric(merged$hf_beta)
  merged$hf_p <- as.numeric(merged$hf_p)
  merged <- merged[complete.cases(merged), ]
  
  cat(sprintf("  Proteins with all 3 endpoints: %d\n", nrow(merged)))
  
  # Count observed concordant hits
  count_concordant <- function(df) {
    n <- 0
    snps_seen <- c()
    for (j in 1:nrow(df)) {
      if (df$bmi_p[j] < 0.05 & df$t2d_p[j] < 0.05 & df$hf_p[j] < 0.05) {
        signs <- sign(c(df$bmi_beta[j], df$t2d_beta[j], df$hf_beta[j]))
        signs <- signs[signs != 0]
        if (length(unique(signs)) == 1) {
          n <- n + 1
        }
      }
    }
    return(n)
  }
  
  observed <- count_concordant(merged)
  cat(sprintf("  Observed concordant hits: %d\n", observed))
  
  # Permutation: shuffle BMI and T2D labels (hold HF fixed as anchoring phenotype)
  # Actually for positive control: hold BMI fixed (the "discovery" endpoint), shuffle T2D and HF
  set.seed(42)
  n_perm <- 10000
  null_dist <- numeric(n_perm)
  
  cat(sprintf("  Running %d permutations (shuffling T2D + HF labels)...\n", n_perm))
  
  for (p in 1:n_perm) {
    perm <- merged
    perm$t2d_p <- sample(perm$t2d_p)
    perm$t2d_beta <- sample(perm$t2d_beta)
    perm$hf_p <- sample(perm$hf_p)
    perm$hf_beta <- sample(perm$hf_beta)
    null_dist[p] <- count_concordant(perm)
  }
  
  emp_p <- mean(null_dist >= observed)
  cat(sprintf("  Null mean: %.3f, Null max: %d\n", mean(null_dist), max(null_dist)))
  cat(sprintf("  Observed: %d, Empirical P = %.4f\n", observed, emp_p))
  
  # Repeat with different seeds
  cat("\n  Stability across seeds:\n")
  for (seed in c(11, 101, 2024, 31415)) {
    set.seed(seed)
    nd <- numeric(n_perm)
    for (p in 1:n_perm) {
      perm <- merged
      perm$t2d_p <- sample(perm$t2d_p)
      perm$t2d_beta <- sample(perm$t2d_beta)
      perm$hf_p <- sample(perm$hf_p)
      perm$hf_beta <- sample(perm$hf_beta)
      nd[p] <- count_concordant(perm)
    }
    cat(sprintf("    Seed %5d: empirical P = %.4f\n", seed, mean(nd >= observed)))
  }
  
} else {
  cat("  WARNING: S1 not loaded. Skipping permutation.\n")
}

# ============================================================
# 6. NEGATIVE CONTROL CHECK
# ============================================================
cat("\n=== 6. Negative Control: Do hits associate with CTS or Asthma? ===\n\n")

if (!is.null(s1)) {
  for (i in 1:nrow(hits)) {
    s1_row <- s1[s1$`OpenGWAS.ID` == hits$id[i], ]
    if (nrow(s1_row) > 0) {
      cts_p <- as.numeric(s1_row$`P.(CTS)`[1])
      asthma_p <- as.numeric(s1_row$`P.(Asthma)`[1])
      cat(sprintf("  %s: CTS P=%.4f %s | Asthma P=%.4f %s\n",
                  hits$name[i], 
                  ifelse(is.na(cts_p), NA, cts_p),
                  ifelse(!is.na(cts_p) & cts_p < 0.05, "[SIG]", ""),
                  ifelse(is.na(asthma_p), NA, asthma_p),
                  ifelse(!is.na(asthma_p) & asthma_p < 0.05, "[SIG]", "")))
    }
  }
}

# ============================================================
# SUMMARY
# ============================================================
cat("\n")
cat("========================================\n")
cat("=== SENSITIVITY SUMMARY ===\n")
cat("========================================\n\n")
cat("Check the log above for:\n")
cat("  1. Shared IVs → how many independent loci?\n")
cat("  2. F-statistics → all >10?\n")
cat("  3. Steiger → correct direction?\n")
cat("  4. Reverse MR → any reverse signals?\n")
cat("  5. Permutation → empirical P stable?\n")
cat("  6. Negative control → phenotype-specific?\n")

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("=== Done ===\n")
