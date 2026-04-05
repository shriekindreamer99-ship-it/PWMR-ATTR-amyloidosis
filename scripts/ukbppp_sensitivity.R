#!/usr/bin/env Rscript
# ============================================================
# UKB-PPP Cross-Platform Sensitivity Analysis
# 
# Check whether the 5 positive-control concordant proteins
# have pQTLs in UKB-PPP (Olink) at the same lead SNPs.
#
# Runtime: ~5-10 minutes (API queries only, no full PWMR)
# ============================================================

cat("=== UKB-PPP Cross-Platform Sensitivity ===\n")
cat("Started:", as.character(Sys.time()), "\n\n")

library(ieugwasr)

# The 5 concordant proteins and their lead SNPs
hits <- data.frame(
  protein = c("CA5A", "TNS4", "STX10", "KPNA2", "PDGFRL"),
  gene = c("CA5A", "TNS4", "STX10", "KPNA2", "PDGFRL"),
  lead_snp = c("rs1157745", "rs662", "rs662", "rs11867412", "rs4683234"),
  sun_id = c("prot-a-332", "prot-a-3070", "prot-a-2883", "prot-a-1677", "prot-a-2231"),
  sun_beta = c(0.034, 0.031, 0.014, 0.033, -0.004),
  stringsAsFactors = FALSE
)

cat("Looking up UKB-PPP Olink pQTL IDs...\n")

# Search for UKB-PPP protein IDs
ao <- available_outcomes()
ukb_ppp <- ao[grepl("UKB-PPP|ukb-ppp|prot-b-|ebi-a-GCST90", ao$id, ignore.case = TRUE), ]

# Also check for Olink-based pQTLs
olink <- ao[grepl("olink|Olink", ao$trait, ignore.case = TRUE) | 
            grepl("UKB-PPP", ao$consortium, ignore.case = TRUE), ]

cat(sprintf("  Found %d UKB-PPP / Olink entries in OpenGWAS\n", nrow(olink)))

# Search for each protein
cat("\n=== Searching for protein matches ===\n\n")

for (i in 1:nrow(hits)) {
  gene <- hits$gene[i]
  
  # Search by gene name in all outcomes
  matches <- ao[grepl(gene, ao$trait, ignore.case = TRUE) & 
                grepl("prot|olink|UKB-PPP|Olink", ao$id, ignore.case = TRUE), ]
  
  # Also try broader search
  if (nrow(matches) == 0) {
    matches <- ao[grepl(gene, ao$trait, ignore.case = TRUE) & 
                  grepl("prot", ao$id, ignore.case = TRUE), ]
  }
  
  cat(sprintf("  %s: %d matches in OpenGWAS\n", gene, nrow(matches)))
  if (nrow(matches) > 0) {
    for (j in 1:min(5, nrow(matches))) {
      cat(sprintf("    %s | %s | N=%s | %s\n", 
                  matches$id[j], 
                  substr(matches$trait[j], 1, 50),
                  matches$sample_size[j],
                  matches$consortium[j]))
    }
  }
}

# Direct SNP lookup approach: query each lead SNP in all pQTL datasets
cat("\n=== Direct SNP-level lookup ===\n")
cat("Querying lead SNPs in association databases...\n\n")

for (i in 1:nrow(hits)) {
  snp <- hits$lead_snp[i]
  gene <- hits$gene[i]
  
  cat(sprintf("--- %s (lead SNP: %s) ---\n", gene, snp))
  
  tryCatch({
    # PheWAS-style lookup: what traits is this SNP associated with?
    phewas <- phewas(snp, pval = 5e-6)
    
    if (!is.null(phewas) && nrow(phewas) > 0) {
      # Filter for protein-related
      prot_hits <- phewas[grepl("prot|Olink|protein|UKB-PPP", 
                                phewas$id, ignore.case = TRUE), ]
      
      if (nrow(prot_hits) > 0) {
        prot_hits <- prot_hits[order(prot_hits$p), ]
        cat(sprintf("  Found %d pQTL associations (P < 5e-6):\n", nrow(prot_hits)))
        for (j in 1:min(10, nrow(prot_hits))) {
          cat(sprintf("    %s | beta=%+.4f | P=%.2e | %s\n",
                      prot_hits$id[j], prot_hits$beta[j], prot_hits$p[j],
                      substr(prot_hits$trait[j], 1, 50)))
        }
      } else {
        cat("  No pQTL associations found for this SNP\n")
      }
      
      # Also show top non-pQTL hits for context
      non_prot <- phewas[!grepl("prot|Olink", phewas$id, ignore.case = TRUE), ]
      if (nrow(non_prot) > 0) {
        non_prot <- non_prot[order(non_prot$p), ]
        cat(sprintf("  Top non-pQTL associations (%d total):\n", nrow(non_prot)))
        for (j in 1:min(3, nrow(non_prot))) {
          cat(sprintf("    %s | P=%.2e | %s\n",
                      non_prot$id[j], non_prot$p[j],
                      substr(non_prot$trait[j], 1, 60)))
        }
      }
    } else {
      cat("  No associations found at P < 5e-6\n")
    }
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
  })
  
  cat("\n")
  Sys.sleep(2)
}

# Also try: look up Sun et al pQTL effect and check if same SNP
# has consistent direction in any UKB-PPP pQTL
cat("\n=== Summary ===\n")
cat("For each protein, compare:\n")
cat("  1. Sun et al. 2018 pQTL effect at lead SNP\n")
cat("  2. UKB-PPP / Olink pQTL effect at same SNP (if available)\n")
cat("  3. Direction consistency\n")
cat("\nResults will be added to Supplementary Table S5.\n")

cat("\nFinished:", as.character(Sys.time()), "\n")
