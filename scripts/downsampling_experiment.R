#!/usr/bin/env Rscript
# ============================================================
# Downsampling Experiment: Power-gradient analysis of PT-MR
#
# Logic: Take the positive control (BMI→T2D→HF) and artificially
# inflate the SE of the weakest endpoint (HF) to simulate
# progressively smaller sample sizes. Show that concordant hits
# disappear predictably as power drops toward ATTR levels.
#
# This does NOT require re-downloading data or API calls.
# It uses the existing BMI + S1 (T2D/HF) results and simulates
# reduced power by scaling standard errors.
#
# Runtime: ~30 minutes (permutation at each power level)
# ============================================================

cat("=== Downsampling Experiment ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

library(openxlsx)

# === Load data ===
bmi_all <- read.csv("PosCtrl_BMI_all_results.csv", stringsAsFactors = FALSE)
bmi_pri <- bmi_all[bmi_all$method %in% c("Inverse variance weighted", "Wald ratio"), ]

s1 <- read.xlsx("投稿包/Supplementary_Table_S1_fixed.xlsx", startRow = 2)

merged <- merge(
  bmi_pri[, c("id.exposure", "b", "se", "pval", "nsnp")],
  s1[, c("OpenGWAS.ID", "Beta.(T2D)", "SE.(T2D)", "P.(T2D)", 
         "Beta.(HF)", "SE.(HF)", "P.(HF)")],
  by.x = "id.exposure", by.y = "OpenGWAS.ID"
)

names(merged) <- c("id", "bmi_b", "bmi_se", "bmi_p", "nsnp",
                    "t2d_b", "t2d_se", "t2d_p", "hf_b", "hf_se", "hf_p")

# Convert to numeric
for (col in c("t2d_b","t2d_se","t2d_p","hf_b","hf_se","hf_p")) {
  merged[[col]] <- as.numeric(merged[[col]])
}
merged <- merged[complete.cases(merged), ]
cat("Proteins with all 3 endpoints:", nrow(merged), "\n\n")

# === Concordance counter ===
count_concordant <- function(df) {
  sum(sapply(1:nrow(df), function(j) {
    if (df$bmi_p[j] < 0.05 & df$t2d_p[j] < 0.05 & df$hf_p[j] < 0.05) {
      s <- sign(c(df$bmi_b[j], df$t2d_b[j], df$hf_b[j]))
      s <- s[s != 0]
      return(length(unique(s)) == 1)
    }
    FALSE
  }))
}

# === Downsampling simulation ===
# The key idea: if we had fewer cases in the outcome GWAS,
# the SE would be larger by factor sqrt(N_original / N_reduced).
# We simulate this by scaling SE and recalculating P-values.
#
# HF HERMES: N_eff ≈ 180,076
# Amyloidosis R5: N_eff ≈ 903
# We create a gradient between these.

# Also downsample T2D to simulate a "rare disease" version
# T2D DIAGRAM: N ≈ 898,130 (very large)

# Power levels to simulate (fraction of original N_eff)
# 1.0 = full power, down to ATTR-equivalent
hf_neff_original <- 180076
amy_neff <- 903

# Simulate reducing BOTH T2D and HF (keeping BMI fixed as the "discovery" endpoint)
# This mimics: "what if the trajectory endpoints were increasingly underpowered?"

fractions <- c(1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005)
# 0.005 × 180076 ≈ 900 → approximately ATTR amyloidosis level

cat("=== Power Gradient Results ===\n\n")
cat(sprintf("%-10s %10s %10s %8s %10s %10s\n", 
            "Fraction", "~N_eff_HF", "~N_eff_T2D", "Hits", "Perm_P", "Null_mean"))
cat(paste(rep("-", 65), collapse=""), "\n")

results <- data.frame()

for (frac in fractions) {
  # Scale SEs
  sim <- merged
  scale_factor <- 1 / sqrt(frac)  # SE scales as 1/sqrt(N)
  
  sim$hf_se_new <- sim$hf_se * scale_factor
  sim$hf_p <- 2 * pnorm(-abs(sim$hf_b / sim$hf_se_new))  # recalculate P
  
  sim$t2d_se_new <- sim$t2d_se * scale_factor
  sim$t2d_p <- 2 * pnorm(-abs(sim$t2d_b / sim$t2d_se_new))  # recalculate P
  
  # Count concordant hits
  obs <- count_concordant(sim)
  
  # Quick permutation (1000 iterations for speed)
  set.seed(42)
  null <- numeric(1000)
  for (p in 1:1000) {
    pm <- sim
    pm$t2d_p <- sample(pm$t2d_p); pm$t2d_b <- sample(pm$t2d_b)
    pm$hf_p <- sample(pm$hf_p); pm$hf_b <- sample(pm$hf_b)
    null[p] <- count_concordant(pm)
  }
  emp_p <- mean(null >= obs)
  
  approx_neff_hf <- round(hf_neff_original * frac)
  approx_neff_t2d <- round(898130 * frac)
  
  cat(sprintf("%-10s %10d %10d %8d %10.4f %10.3f\n",
              sprintf("%.1f%%", frac*100), approx_neff_hf, approx_neff_t2d,
              obs, emp_p, mean(null)))
  
  results <- rbind(results, data.frame(
    fraction = frac,
    neff_hf = approx_neff_hf,
    neff_t2d = approx_neff_t2d,
    concordant_hits = obs,
    empirical_p = emp_p,
    null_mean = mean(null),
    null_max = max(null)
  ))
}

# Save results
write.csv(results, "Downsampling_results.csv", row.names = FALSE)

cat("\n=== Key Comparisons ===\n")
cat(sprintf("  Full power (100%%):     %d concordant hits, P=%.4f\n", 
            results$concordant_hits[1], results$empirical_p[1]))

# Find the row closest to ATTR N_eff
attr_row <- which.min(abs(results$neff_hf - amy_neff))
cat(sprintf("  ATTR-equivalent (~%.1f%%): %d concordant hits, P=%.4f  (N_eff≈%d)\n",
            results$fraction[attr_row]*100, results$concordant_hits[attr_row],
            results$empirical_p[attr_row], results$neff_hf[attr_row]))

cat("\n=== Interpretation ===\n")
cat("If concordant hits drop to 0 at ATTR-equivalent power levels,\n")
cat("this directly demonstrates that the ATTR null result is a predictable\n")
cat("consequence of insufficient power, not a framework failure.\n")

# === Also: single-endpoint benchmark ===
cat("\n\n=== BONUS: Single-endpoint vs PT-MR comparison ===\n\n")

# How many proteins are Bonferroni-significant for BMI alone?
bonf_bmi <- sum(bmi_pri$pval < 0.05/nrow(bmi_pri))

# How many of those also pass T2D and HF (concordance)?
bonf_bmi_ids <- bmi_pri$id.exposure[bmi_pri$pval < 0.05/nrow(bmi_pri)]
bonf_in_merged <- merged[merged$id %in% bonf_bmi_ids, ]
bonf_concordant <- 0
if (nrow(bonf_in_merged) > 0) {
  for (j in 1:nrow(bonf_in_merged)) {
    r <- bonf_in_merged[j, ]
    if (r$t2d_p < 0.05 & r$hf_p < 0.05) {
      s <- sign(c(r$bmi_b, r$t2d_b, r$hf_b))
      s <- s[s != 0]
      if (length(unique(s)) == 1) bonf_concordant <- bonf_concordant + 1
    }
  }
}

cat(sprintf("Single-endpoint (BMI Bonferroni): %d significant proteins\n", bonf_bmi))
cat(sprintf("  Of those, %d also pass T2D+HF concordance\n", bonf_concordant))
cat(sprintf("PT-MR concordance filter: %d hits (from %d nominal BMI)\n",
            count_concordant(merged), sum(merged$bmi_p < 0.05)))
cat("\nThis shows PT-MR captures signals that single-endpoint Bonferroni misses\n")
cat("(or conversely, that concordance adds stringency beyond single-endpoint).\n")

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("=== Done ===\n")
