#!/usr/bin/env Rscript
# ============================================================
# UKB-PPP Supplementary MR: PON1 and CA5A
# Uses UKB-PPP Olink pQTLs as instruments
# Tests against non-UKB outcomes to avoid sample overlap
#
# No-overlap outcomes:
#   - Amyloidosis: FinnGen R12 (local file)
#   - T2D: DIAGRAM (OpenGWAS, ebi-a-GCST006867)
#   - HF: HERMES (OpenGWAS, ebi-a-GCST009541)
# Overlap outcomes (reported as sensitivity only):
#   - BMI: UKB Neale (ieu-b-40) — flagged as overlapping
# ============================================================

cat("=== UKB-PPP Supplementary MR ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

library(data.table)
library(ieugwasr)

# === Step 1: Extract UKB-PPP instruments ===
cat("Step 1: Extracting UKB-PPP instruments...\n\n")

# PON1: rs662 (GRCh37: chr7:94937446)
# From our earlier lookup: beta=+0.588, SE=0.00828, N=33985
pon1_inst <- data.frame(
  SNP = "rs662",
  chr = 7,
  pos37 = 94937446,
  effect_allele = "T",  # from UKB-PPP: ALLELE0=T, ALLELE1=C, beta is for ALLELE1
  other_allele = "C",
  beta.exposure = 0.588403,
  se.exposure = 0.00827716,
  pval.exposure = 10^(-1099.29),
  eaf.exposure = 0.999473,  # A1FREQ — this is freq of ALLELE1
  exposure = "PON1 (UKB-PPP Olink)",
  id.exposure = "ukbppp-PON1"
)

# Wait — A1FREQ=0.999473 means ALLELE1 (C) is nearly fixed
# The beta is for ALLELE1 (C allele), effect allele should be C
pon1_inst$effect_allele <- "C"
pon1_inst$other_allele <- "T"
pon1_inst$eaf.exposure <- 0.999473
# But MAF = 1-0.999 = 0.0005?? That's extremely rare
# Let me check: A1FREQ is for ALLELE1. If A1FREQ=0.999, then ALLELE1=C is major
# BETA=+0.588 means C allele increases PON1
# This seems like the variant is nearly monomorphic... 
# Actually, looking again: GENPOS=95308134 is GRCh38, ID has 94937446 (GRCh37)
# A1FREQ=0.287377 — wait, let me recheck

# From the actual output:
# 7 95308134 7:94937446:T:C:imp:v1 T C 0.287377 0.999473 33985 ADD 0.588403 0.00827716 5053.44 1099.29 NA
# ALLELE0=T, ALLELE1=C, A1FREQ=0.287377
# So C allele frequency = 0.287, T allele frequency = 0.713
# BETA=+0.588 is per C allele (the 192R arginine allele)

pon1_inst$effect_allele <- "C"
pon1_inst$other_allele <- "T" 
pon1_inst$eaf.exposure <- 0.287377
pon1_inst$beta.exposure <- 0.588403
pon1_inst$se.exposure <- 0.00827716

cat("  PON1 instrument: rs662 (C allele)\n")
cat(sprintf("    beta = %+.4f, SE = %.6f, EAF = %.3f, F = %.0f\n",
            pon1_inst$beta.exposure, pon1_inst$se.exposure, 
            pon1_inst$eaf.exposure,
            (pon1_inst$beta.exposure/pon1_inst$se.exposure)^2))

# CA5A: need to identify lead cis-pQTL from chr16
# From earlier output, top hit at pos=87893616 (GRCh38), LOG10P=7400
# We need GRCh37 position — extract from ID column
cat("\n  Loading CA5A to extract lead instrument...\n")
ca5a_file <- Sys.glob("ukbppp/CA5A*/discovery_chr16*")
if (length(ca5a_file) == 0) ca5a_file <- Sys.glob("~/Downloads/ukbppp/CA5A*/discovery_chr16*")

ca5a <- fread(cmd = paste("gzcat", shQuote(ca5a_file[1])))
ca5a[, pos37 := as.numeric(sapply(strsplit(ID, ":"), "[", 2))]
ca5a_top <- ca5a[which.max(LOG10P)]

cat(sprintf("  CA5A top hit: GRCh38=%d, GRCh37=%d, LOG10P=%.0f\n",
            ca5a_top$GENPOS, ca5a_top$pos37, ca5a_top$LOG10P))
cat(sprintf("    ID: %s\n", ca5a_top$ID))
cat(sprintf("    beta = %+.4f, SE = %.6f, A1FREQ = %.3f\n",
            ca5a_top$BETA, ca5a_top$SE, ca5a_top$A1FREQ))

# For CA5A, we need to find the rsID
# The ID format is "16:pos37:A:B:imp:v1"
# We'll use chr:pos to query OpenGWAS
ca5a_chrpos37 <- paste0("16:", ca5a_top$pos37)

ca5a_inst <- data.frame(
  SNP = ca5a_chrpos37,  # will need rsID for outcome lookup
  chr = 16,
  pos37 = ca5a_top$pos37,
  effect_allele = ca5a_top$ALLELE1,
  other_allele = ca5a_top$ALLELE0,
  beta.exposure = ca5a_top$BETA,
  se.exposure = ca5a_top$SE,
  pval.exposure = 10^(-ca5a_top$LOG10P),
  eaf.exposure = ca5a_top$A1FREQ,
  exposure = "CA5A (UKB-PPP Olink)",
  id.exposure = "ukbppp-CA5A"
)

cat(sprintf("  CA5A instrument: chr16:%d\n", ca5a_top$pos37))
cat(sprintf("    F = %.0f\n", (ca5a_inst$beta.exposure/ca5a_inst$se.exposure)^2))

# === Step 2: Look up instruments in outcomes ===
cat("\nStep 2: Looking up instruments in outcome GWAS...\n\n")

# 2a: FinnGen R12 Amyloidosis (local file, NO overlap)
cat("--- FinnGen R12 Amyloidosis (no overlap) ---\n")
r12_file <- Sys.glob("~/Desktop/PWMR/finngen_R12_E4_AMYLNAS.gz")
if (length(r12_file) > 0) {
  r12 <- fread(cmd = paste("gzcat", shQuote(r12_file[1])))
  
  # rs662 in R12
  r12_rs662 <- r12[rsids == "rs662"]
  if (nrow(r12_rs662) > 0) {
    cat(sprintf("  PON1→Amy: rs662 found. beta=%+.4f, SE=%.4f, P=%.4f\n",
                r12_rs662$beta[1], r12_rs662$sebeta[1], r12_rs662$pval[1]))
    
    # Wald ratio MR
    mr_beta <- r12_rs662$beta[1] / pon1_inst$beta.exposure
    mr_se <- r12_rs662$sebeta[1] / abs(pon1_inst$beta.exposure)
    mr_p <- 2 * pnorm(-abs(mr_beta / mr_se))
    cat(sprintf("  MR (Wald): beta=%+.4f, SE=%.4f, P=%.4f\n", mr_beta, mr_se, mr_p))
  }
  
  # CA5A lead SNP in R12 — search by position
  # Need to match by chr:pos
  r12_ca5a <- r12[`#chrom` == 16 & pos == ca5a_top$pos37]
  if (nrow(r12_ca5a) > 0) {
    cat(sprintf("\n  CA5A→Amy: found. beta=%+.4f, SE=%.4f, P=%.4f\n",
                r12_ca5a$beta[1], r12_ca5a$sebeta[1], r12_ca5a$pval[1]))
    mr_beta2 <- r12_ca5a$beta[1] / ca5a_inst$beta.exposure
    mr_se2 <- r12_ca5a$sebeta[1] / abs(ca5a_inst$beta.exposure)
    mr_p2 <- 2 * pnorm(-abs(mr_beta2 / mr_se2))
    cat(sprintf("  MR (Wald): beta=%+.4f, SE=%.4f, P=%.4f\n", mr_beta2, mr_se2, mr_p2))
  } else {
    cat("  CA5A lead SNP not found in R12 (position mismatch or not genotyped)\n")
  }
} else {
  cat("  R12 file not found\n")
}

# 2b: T2D DIAGRAM (OpenGWAS, minimal UKB overlap in older version)
cat("\n--- T2D DIAGRAM ---\n")
tryCatch({
  t2d_rs662 <- associations("rs662", "ebi-a-GCST006867")
  if (nrow(t2d_rs662) > 0) {
    cat(sprintf("  PON1→T2D: rs662 beta=%+.4f, SE=%.4f, P=%.2e\n",
                as.numeric(t2d_rs662$beta[1]), as.numeric(t2d_rs662$se[1]), as.numeric(t2d_rs662$pval[1])))
    mr_b <- as.numeric(t2d_rs662$beta[1]) / pon1_inst$beta.exposure
    mr_s <- as.numeric(t2d_rs662$se[1]) / abs(pon1_inst$beta.exposure)
    mr_p <- 2*pnorm(-abs(mr_b/mr_s))
    cat(sprintf("  MR (Wald): beta=%+.4f, SE=%.4f, P=%.4f\n", mr_b, mr_s, mr_p))
  }
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))

# 2c: HF HERMES (OpenGWAS)
cat("\n--- HF HERMES ---\n")
tryCatch({
  hf_rs662 <- associations("rs662", "ebi-a-GCST009541")
  if (nrow(hf_rs662) > 0) {
    cat(sprintf("  PON1→HF: rs662 beta=%+.4f, SE=%.4f, P=%.2e\n",
                as.numeric(hf_rs662$beta[1]), as.numeric(hf_rs662$se[1]), as.numeric(hf_rs662$pval[1])))
    mr_b <- as.numeric(hf_rs662$beta[1]) / pon1_inst$beta.exposure
    mr_s <- as.numeric(hf_rs662$se[1]) / abs(pon1_inst$beta.exposure)
    mr_p <- 2*pnorm(-abs(mr_b/mr_s))
    cat(sprintf("  MR (Wald): beta=%+.4f, SE=%.4f, P=%.4f\n", mr_b, mr_s, mr_p))
  }
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))

# 2d: BMI (overlap — reported as sensitivity)
cat("\n--- BMI (UKB, OVERLAPPING — sensitivity only) ---\n")
tryCatch({
  bmi_rs662 <- associations("rs662", "ieu-b-40")
  if (nrow(bmi_rs662) > 0) {
    cat(sprintf("  PON1→BMI: rs662 beta=%+.6f, SE=%.6f, P=%.2e\n",
                as.numeric(bmi_rs662$beta[1]), as.numeric(bmi_rs662$se[1]), as.numeric(bmi_rs662$pval[1])))
    mr_b <- as.numeric(bmi_rs662$beta[1]) / pon1_inst$beta.exposure
    mr_s <- as.numeric(bmi_rs662$se[1]) / abs(pon1_inst$beta.exposure)
    mr_p <- 2*pnorm(-abs(mr_b/mr_s))
    cat(sprintf("  MR (Wald): beta=%+.6f, SE=%.6f, P=%.4f\n", mr_b, mr_s, mr_p))
    cat("  ⚠️ NOTE: Sample overlap (UKB-PPP + UKB BMI). Interpret with caution.\n")
  }
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))

# === Step 3: Compare Sun et al. vs UKB-PPP instruments ===
cat("\n\n=== COMPARISON: Sun et al. vs UKB-PPP instrument strength ===\n")
cat(sprintf("%-20s %10s %10s %10s\n", "Platform", "Beta", "SE", "F-stat"))
cat(sprintf("%-20s %10.4f %10.6f %10.0f\n", "Sun (TNS4 trans)", 0.2358, 0.2358/sqrt(74.6), 74.6))
cat(sprintf("%-20s %10.4f %10.6f %10.0f\n", "Sun (STX10 trans)", 0.5296, 0.5296/sqrt(419.8), 419.8))
cat(sprintf("%-20s %10.4f %10.6f %10.0f\n", "UKB-PPP (PON1 cis)", 0.5884, 0.00828, 5053))

cat("\nFinished:", as.character(Sys.time()), "\n")
