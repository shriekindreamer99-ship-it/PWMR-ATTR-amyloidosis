#!/usr/bin/env Rscript
# ============================================================
# Overnight Job: Full PWMR replication using FinnGen R12
# Expected runtime: 4-8 hours (API-dependent)
# 
# BEFORE RUNNING:
# 1. Make sure finngen_R12_E4_AMYLNAS.gz is in ~/
# 2. Make sure you have TwoSampleMR and ieugwasr installed
# 3. Run in terminal: Rscript R12_full_PWMR.R > R12_log.txt 2>&1 &
#    (the & makes it run in background so you can close terminal)
# ============================================================

cat("=== R12 Full PWMR Replication ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

library(TwoSampleMR)
library(data.table)

# === Step 1: Load R12 amyloidosis summary stats ===
cat("Step 1: Loading R12 summary stats...\n")
r12_file <- "~/finngen_R12_E4_AMYLNAS.gz"
if (!file.exists(r12_file)) stop("R12 file not found at ", r12_file)

r12 <- fread(cmd = paste("gzcat", r12_file), sep = "\t")
cat("  Loaded:", nrow(r12), "variants\n")
cat("  Columns:", paste(names(r12), collapse = ", "), "\n")

# Standardize column names for TwoSampleMR
# FinnGen R12 format: #chrom, pos, ref, alt, rsids, nearest_genes, pval, mlogp, beta, sebeta, af_alt, af_alt_cases, af_alt_controls
setnames(r12, c("#chrom"), c("chrom"), skip_absent = TRUE)
r12[, SNP := rsids]
r12[, effect_allele.outcome := alt]
r12[, other_allele.outcome := ref]
r12[, beta.outcome := beta]
r12[, se.outcome := sebeta]
r12[, pval.outcome := pval]
r12[, eaf.outcome := af_alt]
r12[, outcome := "R12_amyloidosis"]

# Index by rsid for fast lookup
setkey(r12, SNP)
cat("  Indexed by rsid. Ready.\n\n")

# === Step 2: Get list of all Sun et al. protein IDs ===
cat("Step 2: Getting protein list...\n")

# All Sun et al. 2018 proteins are prot-a-1 through prot-a-2994 (not all exist)
# We'll query available ones. TwoSampleMR can list them.
ao <- available_outcomes()
sun_proteins <- ao[grepl("^prot-a-", ao$id), ]
cat("  Found", nrow(sun_proteins), "prot-a-* outcomes in OpenGWAS\n")

# We need them as EXPOSURES, not outcomes
# The protein IDs for exposure extraction
protein_ids <- sun_proteins$id
cat("  Will test", length(protein_ids), "proteins\n\n")

# === Step 3: Loop through all proteins ===
cat("Step 3: Running MR for each protein against R12 amyloidosis...\n")
cat("  This will take several hours. Progress logged below.\n\n")

results_list <- list()
errors <- character()

for (i in seq_along(protein_ids)) {
  pid <- protein_ids[i]
  
  if (i %% 50 == 0 || i == 1) {
    cat(sprintf("  [%s] Processing %d/%d: %s\n", 
                format(Sys.time(), "%H:%M"), i, length(protein_ids), pid))
  }
  
  tryCatch({
    # Get instruments for this protein
    exp <- extract_instruments(pid, p1 = 5e-8, clump = TRUE, 
                                r2 = 0.001, kb = 10000)
    
    if (is.null(exp) || nrow(exp) == 0) next
    
    # Look up these SNPs in R12 data
    snps_needed <- exp$SNP
    r12_subset <- r12[SNP %in% snps_needed]
    
    if (nrow(r12_subset) == 0) next
    
    # Format as outcome for TwoSampleMR
    outcome_dat <- data.frame(
      SNP = r12_subset$SNP,
      beta.outcome = r12_subset$beta.outcome,
      se.outcome = r12_subset$se.outcome,
      pval.outcome = r12_subset$pval.outcome,
      effect_allele.outcome = r12_subset$effect_allele.outcome,
      other_allele.outcome = r12_subset$other_allele.outcome,
      eaf.outcome = r12_subset$eaf.outcome,
      outcome = "R12_E4_AMYLNAS",
      id.outcome = "R12_amyloidosis",
      stringsAsFactors = FALSE
    )
    
    # Harmonise
    harm <- tryCatch(
      harmonise_data(exp, outcome_dat, action = 2),
      error = function(e) NULL
    )
    
    if (is.null(harm) || nrow(harm) == 0) next
    
    # Run MR
    mr_res <- tryCatch(mr(harm), error = function(e) NULL)
    
    if (!is.null(mr_res) && nrow(mr_res) > 0) {
      mr_res$protein_id <- pid
      mr_res$protein_name <- sun_proteins$trait[sun_proteins$id == pid]
      results_list[[length(results_list) + 1]] <- mr_res
    }
    
  }, error = function(e) {
    errors <<- c(errors, paste(pid, ":", e$message))
  })
  
  # Be nice to the API
  if (i %% 10 == 0) Sys.sleep(1)
}

# === Step 4: Compile results ===
cat("\n\nStep 4: Compiling results...\n")

if (length(results_list) > 0) {
  all_results <- do.call(rbind, results_list)
  
  # Save full results
  write.csv(all_results, "R12_PWMR_all_results.csv", row.names = FALSE)
  cat("  Saved: R12_PWMR_all_results.csv (", nrow(all_results), "rows)\n")
  
  # Filter to IVW/Wald ratio only
  primary <- all_results[all_results$method %in% c("Inverse variance weighted", "Wald ratio"), ]
  
  # Find nominal hits (P < 0.05)
  nominal <- primary[primary$pval < 0.05, ]
  nominal <- nominal[order(nominal$pval), ]
  write.csv(nominal, "R12_PWMR_nominal_hits.csv", row.names = FALSE)
  cat("  Saved: R12_PWMR_nominal_hits.csv (", nrow(nominal), "nominal hits)\n")
  
  # === Step 5: Check concordance with CTS and HF ===
  # This requires your original CTS and HF results
  # For now, flag the ITIH1-locus proteins specifically
  itih1_proteins <- c("prot-a-2722", "prot-a-1599", "prot-a-575", "prot-a-225", "prot-a-1586")
  itih1_results <- primary[primary$id.exposure %in% itih1_proteins, ]
  
  cat("\n=== Key protein results in R12 ===\n")
  if (nrow(itih1_results) > 0) {
    for (j in 1:nrow(itih1_results)) {
      cat(sprintf("  %s (%s): beta=%.4f, P=%.4f\n",
                  itih1_results$protein_name[j],
                  itih1_results$id.exposure[j],
                  itih1_results$b[j],
                  itih1_results$pval[j]))
    }
  }
  
  cat("\n=== Summary ===\n")
  cat(sprintf("  Proteins tested: %d\n", length(unique(primary$id.exposure))))
  cat(sprintf("  Nominal P < 0.05: %d\n", nrow(nominal)))
  cat(sprintf("  Bonferroni P < %.2e: %d\n", 
              0.05/length(unique(primary$id.exposure)),
              sum(primary$pval < 0.05/length(unique(primary$id.exposure)))))
  
  # Show top 20
  cat("\n=== Top 20 by P-value ===\n")
  top20 <- head(nominal, 20)
  for (j in 1:nrow(top20)) {
    cat(sprintf("  %2d. %s (beta=%.3f, P=%.2e, nSNP=%d)\n",
                j, top20$protein_name[j], top20$b[j], top20$pval[j], top20$nsnp[j]))
  }
  
} else {
  cat("  WARNING: No results collected. Check errors.\n")
}

if (length(errors) > 0) {
  writeLines(errors, "R12_PWMR_errors.txt")
  cat("\n  Errors:", length(errors), "(see R12_PWMR_errors.txt)\n")
}

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("=== Done ===\n")
