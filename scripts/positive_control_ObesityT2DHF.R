#!/usr/bin/env Rscript
# ============================================================
# Positive Control: Obesity → T2D → HF trajectory
# 
# KEY INSIGHT: T2D and HF results ALREADY EXIST in your S1.
# Only BMI/obesity needs to be run fresh.
# Then we cross-reference for concordance.
#
# Runtime: ~2-4 hours
# Run: Rscript positive_control_ObesityT2DHF.R > pos_ctrl_log.txt 2>&1 &
# ============================================================

cat("=== Positive Control: Obesity → T2D → HF ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

library(TwoSampleMR)

# === Step 1: Run PWMR for all proteins against BMI ===
cat("Step 1: Getting protein list...\n")
ao <- available_outcomes()
sun_proteins <- ao[grepl("^prot-a-", ao$id), ]
protein_ids <- sun_proteins$id
cat("  Proteins to test:", length(protein_ids), "\n\n")

# BMI outcome (UKB Neale lab, ~360K samples - massively powered)
bmi_id <- "ieu-b-40"

cat("Step 2: Running MR for each protein against BMI (ieu-b-40)...\n")
cat("  This is the only new analysis needed.\n")
cat("  T2D and HF results will come from your existing S1 data.\n\n")

results_list <- list()
errors <- character()

for (i in seq_along(protein_ids)) {
  pid <- protein_ids[i]
  
  if (i %% 100 == 0 || i == 1) {
    cat(sprintf("  [%s] Processing %d/%d: %s\n",
                format(Sys.time(), "%H:%M"), i, length(protein_ids), pid))
  }
  
  tryCatch({
    exp <- extract_instruments(pid, p1 = 5e-8, clump = TRUE,
                                r2 = 0.001, kb = 10000)
    if (is.null(exp) || nrow(exp) == 0) next
    
    out <- extract_outcome_data(exp$SNP, bmi_id)
    if (is.null(out) || nrow(out) == 0) next
    
    harm <- harmonise_data(exp, out, action = 2)
    if (is.null(harm) || nrow(harm) == 0) next
    
    mr_res <- mr(harm)
    if (!is.null(mr_res) && nrow(mr_res) > 0) {
      mr_res$protein_id <- pid
      mr_res$protein_name <- sun_proteins$trait[sun_proteins$id == pid]
      results_list[[length(results_list) + 1]] <- mr_res
    }
  }, error = function(e) {
    errors <<- c(errors, paste(pid, ":", e$message))
  })
  
  if (i %% 10 == 0) Sys.sleep(0.5)
}

# === Step 2: Compile BMI results ===
cat("\nStep 3: Compiling BMI results...\n")

if (length(results_list) == 0) {
  cat("ERROR: No results. Check API access.\n")
  quit(status = 1)
}

all_bmi <- do.call(rbind, results_list)
write.csv(all_bmi, "PosCtrl_BMI_all_results.csv", row.names = FALSE)

# Filter to IVW/Wald
primary_bmi <- all_bmi[all_bmi$method %in% c("Inverse variance weighted", "Wald ratio"), ]
cat("  BMI: Tested", length(unique(primary_bmi$id.exposure)), "proteins\n")
cat("  BMI nominal (P<0.05):", sum(primary_bmi$pval < 0.05), "\n\n")

# === Step 3: Load T2D and HF from S1 ===
cat("Step 4: Loading T2D and HF from S1...\n")

library(openxlsx)
# Try to load S1 - adjust path as needed
s1_paths <- c(
  "Supplementary_Table_S1_fixed.xlsx",
  "~/Desktop/PWMR/Supplementary_Table_S1_fixed.xlsx",
  "/Users/jc/Desktop/PWMR/Supplementary_Table_S1_fixed.xlsx"
)

s1 <- NULL
for (p in s1_paths) {
  if (file.exists(p)) {
    s1 <- read.xlsx(p, startRow = 2)
    cat("  Loaded S1 from:", p, "\n")
    break
  }
}

if (is.null(s1)) {
  cat("  WARNING: S1 not found. Saving BMI results only.\n")
  cat("  You'll need to manually cross-reference with S1.\n")
  write.csv(primary_bmi, "PosCtrl_BMI_primary.csv", row.names = FALSE)
  cat("\nFinished:", as.character(Sys.time()), "\n")
  quit(status = 0)
}

# === Step 4: Triple concordance filter ===
cat("\nStep 5: Applying concordance filter (BMI + T2D + HF)...\n")

concordant <- data.frame()
near_miss <- data.frame()

for (j in 1:nrow(primary_bmi)) {
  pid <- primary_bmi$id.exposure[j]
  bmi_beta <- primary_bmi$b[j]
  bmi_p <- primary_bmi$pval[j]
  
  # Find in S1
  s1_row <- s1[s1$`OpenGWAS.ID` == pid, ]
  if (nrow(s1_row) == 0) next
  
  t2d_beta <- as.numeric(s1_row$`Beta.(T2D)`[1])
  t2d_p <- as.numeric(s1_row$`P.(T2D)`[1])
  hf_beta <- as.numeric(s1_row$`Beta.(HF)`[1])
  hf_p <- as.numeric(s1_row$`P.(HF)`[1])
  
  if (is.na(t2d_p) || is.na(hf_p)) next
  
  # Direction check
  signs <- sign(c(bmi_beta, t2d_beta, hf_beta))
  same_dir <- length(unique(signs[signs != 0])) == 1
  
  all_sig <- bmi_p < 0.05 & t2d_p < 0.05 & hf_p < 0.05
  two_sig <- sum(c(bmi_p < 0.05, t2d_p < 0.05, hf_p < 0.05)) >= 2
  
  row_data <- data.frame(
    protein_id = pid,
    protein_name = primary_bmi$protein_name[j],
    bmi_beta = bmi_beta, bmi_p = bmi_p,
    t2d_beta = t2d_beta, t2d_p = t2d_p,
    hf_beta = hf_beta, hf_p = hf_p,
    same_direction = same_dir,
    nsnp_bmi = primary_bmi$nsnp[j],
    stringsAsFactors = FALSE
  )
  
  if (all_sig && same_dir) {
    concordant <- rbind(concordant, row_data)
  } else if (two_sig && same_dir) {
    near_miss <- rbind(near_miss, row_data)
  }
}

# === Step 5: Report ===
cat("\n========================================\n")
cat("=== POSITIVE CONTROL RESULTS ===\n")
cat("========================================\n\n")

cat(sprintf("Trajectory: BMI (UKB, ~360K) → T2D (DIAGRAM, 898K) → HF (HERMES, 977K)\n"))
cat(sprintf("Proteins tested: %d\n", length(unique(primary_bmi$id.exposure))))
cat(sprintf("BMI nominal hits: %d\n\n", sum(primary_bmi$pval < 0.05)))

if (nrow(concordant) > 0) {
  cat(sprintf("*** %d CONCORDANT TRIPLE HITS ***\n\n", nrow(concordant)))
  concordant <- concordant[order(concordant$bmi_p), ]
  write.csv(concordant, "PosCtrl_CONCORDANT_HITS.csv", row.names = FALSE)
  
  for (k in 1:min(20, nrow(concordant))) {
    r <- concordant[k, ]
    dir <- ifelse(r$same_direction & r$bmi_beta > 0, "risk(+)", "protect(-)")
    cat(sprintf("  %2d. %s\n      BMI: b=%+.3f P=%.2e | T2D: b=%+.3f P=%.2e | HF: b=%+.3f P=%.2e [%s]\n",
                k, substr(r$protein_name, 1, 50),
                r$bmi_beta, r$bmi_p, r$t2d_beta, r$t2d_p, r$hf_beta, r$hf_p, dir))
  }
  cat(sprintf("\n  Total concordant hits: %d (saved to PosCtrl_CONCORDANT_HITS.csv)\n", nrow(concordant)))
} else {
  cat("No concordant triple hits found.\n")
}

if (nrow(near_miss) > 0) {
  cat(sprintf("\nNear misses (2/3 significant + same direction): %d\n", nrow(near_miss)))
  write.csv(near_miss, "PosCtrl_NEAR_MISSES.csv", row.names = FALSE)
}

if (length(errors) > 0) {
  writeLines(errors, "PosCtrl_errors.txt")
  cat(sprintf("\nErrors: %d (see PosCtrl_errors.txt)\n", length(errors)))
}

cat("\nFinished:", as.character(Sys.time()), "\n")
cat("=== Done ===\n")
