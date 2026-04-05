#!/usr/bin/env Rscript
# ============================================================
# Colocalization: UKB-PPP PON1 pQTL vs BMI / T2D / HF
# at the rs662 locus (chr7, ±500kb)
#
# Uses UKB-PPP Olink PON1 data as the pQTL dataset
# Pulls outcome GWAS from OpenGWAS API
#
# Runtime: ~1-2 hours (API-dependent)
# Run: Rscript coloc_pon1_ukbppp.R > coloc_log.txt 2>&1 &
# ============================================================

cat("=== PON1 Colocalization (UKB-PPP Olink) ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

library(coloc)
library(data.table)
library(ieugwasr)

# === Step 1: Load UKB-PPP PON1 chr7 data ===
cat("Step 1: Loading UKB-PPP PON1 chr7 data...\n")

pon1_file <- Sys.glob("~/Downloads/ukbppp/PON1_P27169_OID30704_v1_Inflammation_II/discovery_chr7*")
if (length(pon1_file) == 0) {
  # Try alternative paths
  pon1_file <- Sys.glob("ukbppp/PON1_P27169_OID30704_v1_Inflammation_II/discovery_chr7*")
}
if (length(pon1_file) == 0) stop("Cannot find PON1 chr7 file")
cat("  File:", pon1_file[1], "\n")

pon1 <- fread(cmd = paste("gzcat", shQuote(pon1_file[1])))
cat("  Loaded:", nrow(pon1), "variants\n")
cat("  Columns:", paste(names(pon1), collapse=", "), "\n")

# rs662 position (GRCh38: 95308134, GRCh37: 94937446)
# UKB-PPP uses GRCh38
rs662_pos <- 95308134
window <- 500000  # ±500kb

# Filter to rs662 region
pon1_region <- pon1[GENPOS >= (rs662_pos - window) & GENPOS <= (rs662_pos + window)]
cat("  Variants in ±500kb window:", nrow(pon1_region), "\n")

# Convert LOG10P to P
pon1_region[, pval := 10^(-LOG10P)]

# Create coloc-format dataset
# Need: beta, varbeta, snp, position, type, N, MAF
pon1_region[, varbeta := SE^2]
pon1_region[, MAF := ifelse(A1FREQ < 0.5, A1FREQ, 1 - A1FREQ)]
pon1_region[, snp := paste0("chr7:", GENPOS)]

cat("  PON1 pQTL at rs662 position:", "\n")
rs662_row <- pon1_region[GENPOS == rs662_pos]
if (nrow(rs662_row) > 0) {
  cat(sprintf("    Beta = %+.4f, SE = %.4f, P = 10^-%.0f, N = %d\n",
              rs662_row$BETA[1], rs662_row$SE[1], rs662_row$LOG10P[1], rs662_row$N[1]))
}

# === Step 2: Get outcome GWAS data from OpenGWAS ===
cat("\nStep 2: Fetching outcome GWAS data from OpenGWAS...\n")

# Outcomes to test
outcomes <- list(
  list(id = "ieu-b-40", name = "BMI", type = "quant", ncase = NA, N = 461460),
  list(id = "ebi-a-GCST006867", name = "T2D (DIAGRAM)", type = "cc", ncase = 74124, N = 898130),
  list(id = "ebi-a-GCST009541", name = "HF (HERMES)", type = "cc", ncase = 47309, N = 977323)
)

# Function to run coloc for one outcome
run_coloc <- function(outcome_info, pqtl_data) {
  cat(sprintf("\n--- Colocalization: PON1 vs %s ---\n", outcome_info$name))
  
  tryCatch({
    # Get outcome data for the region
    # Use associations() with chr:pos format
    # First, get rsIDs from the region using tophits or associations
    out_data <- associations(
      variants = paste0("7:", rs662_pos - window, "-", rs662_pos + window),
      id = outcome_info$id
    )
    
    if (is.null(out_data) || nrow(out_data) == 0) {
      # Try alternative: query by position range
      cat("  Trying alternative query method...\n")
      # Use specific SNPs from our pQTL data
      # Get a sample of positions
      sample_pos <- pqtl_data$GENPOS[seq(1, nrow(pqtl_data), by = max(1, nrow(pqtl_data) %/% 500))]
      sample_chrpos <- paste0("7:", sample_pos)
      
      out_data <- associations(
        variants = sample_chrpos[1:min(100, length(sample_chrpos))],
        id = outcome_info$id,
        proxies = 0
      )
    }
    
    if (is.null(out_data) || nrow(out_data) == 0) {
      cat("  WARNING: No outcome data retrieved\n")
      return(NULL)
    }
    
    cat(sprintf("  Outcome variants retrieved: %d\n", nrow(out_data)))
    
    # Create position-based matching
    out_data$pos <- as.numeric(out_data$position)
    out_data$snp <- paste0("chr7:", out_data$pos)
    
    # Merge on position
    merged_snps <- intersect(pqtl_data$snp, out_data$snp)
    cat(sprintf("  Overlapping variants: %d\n", length(merged_snps)))
    
    if (length(merged_snps) < 50) {
      cat("  WARNING: Too few overlapping variants for coloc\n")
      return(NULL)
    }
    
    pqtl_sub <- pqtl_data[snp %in% merged_snps]
    out_sub <- out_data[out_data$snp %in% merged_snps, ]
    
    # Remove duplicates
    pqtl_sub <- pqtl_sub[!duplicated(snp)]
    out_sub <- out_sub[!duplicated(out_sub$snp), ]
    
    # Align
    pqtl_sub <- pqtl_sub[order(snp)]
    out_sub <- out_sub[order(out_sub$snp), ]
    
    # Build coloc datasets
    dataset1 <- list(
      beta = pqtl_sub$BETA,
      varbeta = pqtl_sub$varbeta,
      snp = pqtl_sub$snp,
      position = pqtl_sub$GENPOS,
      type = "quant",
      N = pqtl_sub$N[1],
      MAF = pqtl_sub$MAF
    )
    
    if (outcome_info$type == "quant") {
      dataset2 <- list(
        beta = as.numeric(out_sub$beta),
        varbeta = as.numeric(out_sub$se)^2,
        snp = out_sub$snp,
        position = out_sub$pos,
        type = "quant",
        N = outcome_info$N,
        MAF = as.numeric(out_sub$eaf)
      )
    } else {
      dataset2 <- list(
        beta = as.numeric(out_sub$beta),
        varbeta = as.numeric(out_sub$se)^2,
        snp = out_sub$snp,
        position = out_sub$pos,
        type = "cc",
        N = outcome_info$N,
        s = outcome_info$ncase / outcome_info$N,
        MAF = as.numeric(out_sub$eaf)
      )
    }
    
    # Fix MAF (ensure 0 < MAF < 0.5)
    dataset1$MAF[is.na(dataset1$MAF) | dataset1$MAF <= 0] <- 0.01
    dataset1$MAF[dataset1$MAF >= 0.5] <- 0.49
    dataset2$MAF[is.na(dataset2$MAF) | dataset2$MAF <= 0] <- 0.01
    dataset2$MAF[dataset2$MAF >= 0.5] <- 0.49
    
    # Run coloc
    result <- coloc.abf(dataset1, dataset2)
    
    cat(sprintf("  PP.H0 = %.4f (no association)\n", result$summary["PP.H0.abf"]))
    cat(sprintf("  PP.H1 = %.4f (pQTL only)\n", result$summary["PP.H1.abf"]))
    cat(sprintf("  PP.H2 = %.4f (outcome only)\n", result$summary["PP.H2.abf"]))
    cat(sprintf("  PP.H3 = %.4f (both, different SNPs)\n", result$summary["PP.H3.abf"]))
    cat(sprintf("  PP.H4 = %.4f (shared causal variant)\n", result$summary["PP.H4.abf"]))
    
    if (result$summary["PP.H4.abf"] > 0.8) {
      cat("  >>> STRONG colocalization evidence <<<\n")
    } else if (result$summary["PP.H4.abf"] > 0.5) {
      cat("  >>> Moderate colocalization evidence <<<\n")
    } else {
      cat("  >>> Weak/no colocalization <<<\n")
    }
    
    return(result$summary)
    
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    return(NULL)
  })
}

# === Step 3: Run coloc for each outcome ===
cat("\nStep 3: Running colocalization analyses...\n")

results_list <- list()
for (out in outcomes) {
  res <- run_coloc(out, pon1_region)
  if (!is.null(res)) {
    results_list[[out$name]] <- res
  }
  Sys.sleep(2)
}

# === Step 4: Also do CA5A ===
cat("\n\n=== CA5A Colocalization ===\n")

ca5a_file <- Sys.glob("~/Downloads/ukbppp/CA5A_P35218_OID20075_v1_Cardiometabolic/discovery_chr16*")
if (length(ca5a_file) == 0) {
  ca5a_file <- Sys.glob("ukbppp/CA5A_P35218_OID20075_v1_Cardiometabolic/discovery_chr16*")
}

if (length(ca5a_file) > 0) {
  cat("Loading CA5A chr16 data...\n")
  ca5a <- fread(cmd = paste("gzcat", shQuote(ca5a_file[1])))
  
  # Find the top hit position for CA5A
  ca5a_top <- ca5a[which.max(LOG10P)]
  ca5a_pos <- ca5a_top$GENPOS
  cat(sprintf("  CA5A top hit: pos=%d, LOG10P=%.1f\n", ca5a_pos, ca5a_top$LOG10P))
  
  # Filter to region
  ca5a_region <- ca5a[GENPOS >= (ca5a_pos - window) & GENPOS <= (ca5a_pos + window)]
  ca5a_region[, pval := 10^(-LOG10P)]
  ca5a_region[, varbeta := SE^2]
  ca5a_region[, MAF := ifelse(A1FREQ < 0.5, A1FREQ, 1 - A1FREQ)]
  ca5a_region[, snp := paste0("chr16:", GENPOS)]
  
  cat(sprintf("  Variants in ±500kb window: %d\n", nrow(ca5a_region)))
  
  for (out in outcomes) {
    out_name_ca5a <- paste0("CA5A_vs_", out$name)
    res <- run_coloc(out, ca5a_region)
    if (!is.null(res)) {
      results_list[[out_name_ca5a]] <- res
    }
    Sys.sleep(2)
  }
}

# === Summary ===
cat("\n\n========================================\n")
cat("=== COLOCALIZATION SUMMARY ===\n")
cat("========================================\n\n")

cat(sprintf("%-30s %8s %8s %8s %8s %8s\n", "Comparison", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"))
cat(paste(rep("-", 80), collapse = ""), "\n")

for (name in names(results_list)) {
  r <- results_list[[name]]
  cat(sprintf("%-30s %8.4f %8.4f %8.4f %8.4f %8.4f %s\n",
              name,
              r["PP.H0.abf"], r["PP.H1.abf"], r["PP.H2.abf"], r["PP.H3.abf"], r["PP.H4.abf"],
              ifelse(r["PP.H4.abf"] > 0.8, "***", ifelse(r["PP.H4.abf"] > 0.5, "**", ""))))
}

# Save results
if (length(results_list) > 0) {
  results_df <- do.call(rbind, lapply(names(results_list), function(n) {
    data.frame(comparison = n, t(results_list[[n]]))
  }))
  write.csv(results_df, "coloc_ukbppp_results.csv", row.names = FALSE)
  cat("\nSaved: coloc_ukbppp_results.csv\n")
}

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("=== Done ===\n")
