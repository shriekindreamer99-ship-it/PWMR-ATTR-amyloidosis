# ============================================================
# alternative_permutation_sensitivity.R
# 
# 目的：测试permutation结果对anchor phenotype选择的敏感性
# 
# 当前方案（primary）：amyloidosis固定，CTS+HF shuffled
# 替代方案A：CTS固定，amyloidosis+HF shuffled
# 替代方案B：HF固定，CTS+amyloidosis shuffled
# 替代方案C：全部shuffle（无anchor）
#
# 如果所有方案的empirical P一致 → 框架对anchor选择不敏感
# ============================================================

library(dplyr)
set.seed(42)

cat("=== Alternative Permutation Sensitivity Analysis ===\n\n")

# --- 读取三表型MR结果 ---
# 请确认路径：你的三表型合并结果文件
# 预期列：id.exposure, protein, 
#          p_cts, beta_cts, p_amyloidosis, beta_amyloidosis, p_hf, beta_hf,
#          lead_snp (或SNP)

tri <- read.csv("results/11_triangulation_combined.csv")
cat("总蛋白数:", nrow(tri), "\n")

# --- 确认列名 ---
# 如果你的列名不同，在这里映射
# tri <- tri %>% rename(p_cts = P_CTS, p_amyloidosis = P_Amyloidosis, p_hf = P_HF, ...)

# --- 定义concordance filter函数 ---
count_concordant_loci <- function(data, p_col1, beta_col1, 
                                  p_col2, beta_col2, 
                                  p_col3, beta_col3,
                                  locus_col = "lead_snp") {
  # 筛选三表型均P<0.05且方向一致的蛋白
  concordant <- data %>%
    filter(
      .data[[p_col1]] < 0.05,
      .data[[p_col2]] < 0.05,
      .data[[p_col3]] < 0.05,
      sign(.data[[beta_col1]]) == sign(.data[[beta_col2]]),
      sign(.data[[beta_col2]]) == sign(.data[[beta_col3]])
    )
  
  # Collapse到locus level
  if (locus_col %in% colnames(data) && nrow(concordant) > 0) {
    n_loci <- length(unique(concordant[[locus_col]]))
  } else {
    n_loci <- nrow(concordant)  # fallback: protein count
  }
  
  return(n_loci)
}

# --- 观测值 ---
observed <- count_concordant_loci(
  tri, "p_cts", "beta_cts", "p_amyloidosis", "beta_amyloidosis", "p_hf", "beta_hf"
)
cat("观测到的concordant loci:", observed, "\n\n")

# --- Permutation函数 ---
run_permutation <- function(data, fixed_pheno, n_iter = 10000, seed = 42) {
  # fixed_pheno: "amyloidosis", "cts", "hf", or "none"
  set.seed(seed)
  
  null_counts <- numeric(n_iter)
  n <- nrow(data)
  
  for (i in 1:n_iter) {
    perm_data <- data
    
    if (fixed_pheno == "amyloidosis") {
      # Primary: amyloidosis固定, CTS+HF shuffled
      perm_data$p_cts <- data$p_cts[sample(n)]
      perm_data$beta_cts <- data$beta_cts[sample(n)]
      perm_data$p_hf <- data$p_hf[sample(n)]
      perm_data$beta_hf <- data$beta_hf[sample(n)]
      
    } else if (fixed_pheno == "cts") {
      # Alt A: CTS固定, amyloidosis+HF shuffled
      perm_data$p_amyloidosis <- data$p_amyloidosis[sample(n)]
      perm_data$beta_amyloidosis <- data$beta_amyloidosis[sample(n)]
      perm_data$p_hf <- data$p_hf[sample(n)]
      perm_data$beta_hf <- data$beta_hf[sample(n)]
      
    } else if (fixed_pheno == "hf") {
      # Alt B: HF固定, CTS+amyloidosis shuffled
      perm_data$p_cts <- data$p_cts[sample(n)]
      perm_data$beta_cts <- data$beta_cts[sample(n)]
      perm_data$p_amyloidosis <- data$p_amyloidosis[sample(n)]
      perm_data$beta_amyloidosis <- data$beta_amyloidosis[sample(n)]
      
    } else if (fixed_pheno == "none") {
      # Alt C: 全部shuffle
      perm_data$p_cts <- data$p_cts[sample(n)]
      perm_data$beta_cts <- data$beta_cts[sample(n)]
      perm_data$p_amyloidosis <- data$p_amyloidosis[sample(n)]
      perm_data$beta_amyloidosis <- data$beta_amyloidosis[sample(n)]
      perm_data$p_hf <- data$p_hf[sample(n)]
      perm_data$beta_hf <- data$beta_hf[sample(n)]
    }
    
    null_counts[i] <- count_concordant_loci(
      perm_data, "p_cts", "beta_cts", 
      "p_amyloidosis", "beta_amyloidosis", 
      "p_hf", "beta_hf"
    )
  }
  
  empirical_p <- mean(null_counts >= observed)
  null_mean <- mean(null_counts)
  
  return(list(
    scheme = fixed_pheno,
    observed = observed,
    null_mean = round(null_mean, 3),
    empirical_p = empirical_p,
    null_max = max(null_counts),
    null_distribution = null_counts
  ))
}

# --- 跑四种方案 ---
cat("Running 4 permutation schemes (10,000 iterations each)...\n\n")

schemes <- c("amyloidosis", "cts", "hf", "none")
scheme_labels <- c(
  "amyloidosis" = "Primary (amyloidosis fixed)",
  "cts"         = "Alt A (CTS fixed)",
  "hf"          = "Alt B (HF fixed)",
  "none"        = "Alt C (all shuffled)"
)

results <- list()
for (s in schemes) {
  cat("  Running:", scheme_labels[s], "...\n")
  results[[s]] <- run_permutation(tri, fixed_pheno = s, n_iter = 10000)
  cat("    Null mean =", results[[s]]$null_mean, 
      "| Empirical P =", results[[s]]$empirical_p, "\n")
}

# --- 汇总表格 ---
cat("\n=== SUMMARY TABLE ===\n")
cat(sprintf("%-30s  %8s  %10s  %12s  %8s\n", 
            "Scheme", "Observed", "Null mean", "Empirical P", "Null max"))
cat(paste(rep("-", 75), collapse=""), "\n")
for (s in schemes) {
  r <- results[[s]]
  cat(sprintf("%-30s  %8d  %10.3f  %12.4f  %8d\n",
              scheme_labels[s], r$observed, r$null_mean, r$empirical_p, r$null_max))
}

# --- 结论 ---
cat("\n=== INTERPRETATION ===\n")
p_values <- sapply(results, function(r) r$empirical_p)
if (max(p_values) - min(p_values) < 0.03) {
  cat("所有方案的empirical P高度一致（range:", 
      round(min(p_values), 4), "-", round(max(p_values), 4), 
      "）\n→ 结论：permutation结果对anchor phenotype选择不敏感\n")
  cat("\n写入稿件的一句话：\n")
  cat('"Sensitivity analysis using alternative anchoring phenotypes\n')
  cat('(CTS-anchored, HF-anchored, or fully shuffled) yielded consistent\n')
  cat('empirical P values (range: ', round(min(p_values), 3), '-', 
      round(max(p_values), 3), '), confirming that the permutation\n')
  cat('result is not sensitive to the choice of anchoring phenotype."\n')
} else {
  cat("方案间有差异，需要讨论。P range:", 
      round(min(p_values), 4), "-", round(max(p_values), 4), "\n")
}

cat("\nDone.\n")
