# ============================================================
# 13b_locus_level_permutation.R
# 目的：把共享lead SNP的蛋白collapse成1个locus单位后重跑permutation
# 解决问题：避免同一位点被重复计算3次而膨胀empirical P
# ============================================================

library(dplyr)
setwd("~/Desktop/PWMR")

cat("=== Locus-Level Permutation Test ===\n\n")

# --- 读取数据 ---
tri <- read.csv("results/11_triangulation_combined.csv")
cat("总蛋白数:", nrow(tri), "\n")

# --- 第一步：给每个蛋白标注lead SNP，用于locus collapse ---
# 你需要确认：你的数据里是否有每个蛋白的lead SNP信息
# 如果没有，我们用一个简化版：直接从已知信息标注

# 检查数据中是否有SNP列
if ("SNP" %in% colnames(tri)) {
  cat("数据中有SNP列，直接用于locus定义\n")
  locus_data <- tri %>%
    mutate(locus_id = SNP)  # 用lead SNP作为locus标识
} else {
  cat("数据中没有SNP列，尝试从其他结果文件获取...\n")
  
  # 尝试从提取数据中获取SNP信息
  # 如果你有 IV extraction 的结果文件，读取它
  if (file.exists("results/01_extraction_progress.rds")) {
    cat("从01_extraction_progress.rds读取SNP信息...\n")
    iv_data <- readRDS("results/01_extraction_progress.rds")
    
    # 获取每个蛋白的lead SNP（取P值最小的那个）
    lead_snps <- iv_data %>%
      group_by(id.exposure) %>%
      slice_min(pval.exposure, n = 1) %>%
      select(id.exposure, SNP) %>%
      ungroup()
    
    locus_data <- tri %>%
      left_join(lead_snps, by = "id.exposure") %>%
      mutate(locus_id = ifelse(is.na(SNP), id.exposure, SNP))
    
    cat("成功匹配到", sum(!is.na(locus_data$SNP)), "个蛋白的SNP\n")
  } else {
    cat("⚠️ 没有找到IV数据文件\n")
    cat("使用备选方案：手动标注已知共享SNP\n")
    
    # 手动标注已知的共享SNP蛋白
    shared_snp_proteins <- c("prot-a-2722", "prot-a-1599", "prot-a-575")  # SIGIRR, JAKMIP3, CLEC2L
    
    locus_data <- tri %>%
      mutate(locus_id = case_when(
        id.exposure %in% shared_snp_proteins ~ "rs1042779_locus",  # collapse成1个
        TRUE ~ id.exposure  # 其他蛋白各自独立
      ))
  }
}

# --- 第二步：Locus-level collapse ---
# 对于共享同一lead SNP的蛋白，只算1个locus hit
# 如果该locus的任何一个蛋白通过了三表型筛选，则该locus算通过

locus_summary <- locus_data %>%
  mutate(
    triple_hit = pval_ca < 0.05 & pval_cts < 0.05 & pval_hf < 0.05 &
      sign(beta_ca) == sign(beta_cts) & sign(beta_ca) == sign(beta_hf)
  ) %>%
  group_by(locus_id) %>%
  summarise(
    n_proteins = n(),
    any_triple = any(triple_hit),
    .groups = "drop"
  )

n_unique_loci <- nrow(locus_summary)
n_observed_locus <- sum(locus_summary$any_triple)

cat("\n--- Locus-level summary ---\n")
cat("总unique loci数:", n_unique_loci, "\n")
cat("观察到的locus-level triple hits:", n_observed_locus, "\n")

# 显示哪些locus通过了
passed_loci <- locus_summary %>% filter(any_triple)
cat("通过的loci:\n")
print(passed_loci)

# --- 第三步：Locus-level permutation ---
set.seed(42)
n_perm <- 10000

cat("\n开始locus-level permutation（", n_perm, "次）...\n")

null_locus_counts <- numeric(n_perm)

for (i in 1:n_perm) {
  if (i %% 1000 == 0) cat("  进度:", i, "/", n_perm, "\n")
  
  # 打乱CTS和HF的P值和beta（保持marginal分布）
  shuffled_cts_idx <- sample(1:nrow(locus_data))
  shuffled_hf_idx <- sample(1:nrow(locus_data))
  
  # 计算每个蛋白是否通过（用打乱后的数据）
  perm_triple <- locus_data$pval_ca < 0.05 & 
    locus_data$pval_cts[shuffled_cts_idx] < 0.05 & 
    locus_data$pval_hf[shuffled_hf_idx] < 0.05 &
    sign(locus_data$beta_ca) == sign(locus_data$beta_cts[shuffled_cts_idx]) & 
    sign(locus_data$beta_ca) == sign(locus_data$beta_hf[shuffled_hf_idx])
  
  # Collapse到locus level：每个locus只算1次
  perm_locus_data <- data.frame(
    locus_id = locus_data$locus_id,
    triple_hit = perm_triple
  ) %>%
    group_by(locus_id) %>%
    summarise(any_triple = any(triple_hit), .groups = "drop")
  
  null_locus_counts[i] <- sum(perm_locus_data$any_triple)
}

# --- 第四步：计算经验P值 ---
empirical_p_locus <- (sum(null_locus_counts >= n_observed_locus) + 1) / (n_perm + 1)

cat("\n=== Locus-Level Permutation Results ===\n")
cat("观察到的locus-level triple hits:", n_observed_locus, "\n")
cat("Null分布均值:", round(mean(null_locus_counts), 4), "\n")
cat("Null分布中位数:", median(null_locus_counts), "\n")
cat("Null分布最大值:", max(null_locus_counts), "\n")
cat("Null分布中 >=", n_observed_locus, "的次数:", sum(null_locus_counts >= n_observed_locus), "\n")
cat("Locus-level 经验P值:", signif(empirical_p_locus, 4), "\n\n")

# --- 对比protein-level结果 ---
cat("--- 与protein-level permutation对比 ---\n")
cat("Protein-level: 3 hits, P < 1e-4\n")
cat("Locus-level:  ", n_observed_locus, "hit(s), P =", signif(empirical_p_locus, 4), "\n")
cat("注意：locus-level才是正文应该报告的主要结果\n")
cat("protein-level可以放supplementary\n")

# --- 保存结果 ---
perm_locus_result <- data.frame(
  level = "locus",
  observed_hits = n_observed_locus,
  null_mean = round(mean(null_locus_counts), 4),
  null_median = median(null_locus_counts),
  null_max = max(null_locus_counts),
  null_ge_observed = sum(null_locus_counts >= n_observed_locus),
  empirical_p = signif(empirical_p_locus, 4),
  n_unique_units = n_unique_loci,
  n_permutations = n_perm
)

write.csv(perm_locus_result, "results/17b_locus_permutation_result.csv", row.names = FALSE)
cat("\n结果已保存到 results/17b_locus_permutation_result.csv\n")

# --- 画图 ---
pdf("figures/permutation_locus_level.pdf", width = 7, height = 5)
hist(null_locus_counts, breaks = seq(-0.5, max(null_locus_counts, n_observed_locus) + 1.5, by = 1),
     col = "grey80", border = "white",
     main = "Null distribution of locus-level triple-validated hits\n(10,000 permutations, shared-SNP proteins collapsed)",
     xlab = "Number of independent loci with triple validation",
     ylab = "Frequency")
abline(v = n_observed_locus, col = "red", lwd = 2, lty = 2)
text(n_observed_locus + 0.3, par("usr")[4] * 0.85,
     paste0("Observed = ", n_observed_locus, "\nP = ", signif(empirical_p_locus, 3)),
     col = "red", pos = 4, cex = 0.9)
dev.off()

png("figures/permutation_locus_level.png", width = 700, height = 500, res = 150)
hist(null_locus_counts, breaks = seq(-0.5, max(null_locus_counts, n_observed_locus) + 1.5, by = 1),
     col = "grey80", border = "white",
     main = "Null distribution of locus-level triple-validated hits\n(10,000 permutations, shared-SNP proteins collapsed)",
     xlab = "Number of independent loci with triple validation",
     ylab = "Frequency")
abline(v = n_observed_locus, col = "red", lwd = 2, lty = 2)
text(n_observed_locus + 0.3, par("usr")[4] * 0.85,
     paste0("Observed = ", n_observed_locus, "\nP = ", signif(empirical_p_locus, 3)),
     col = "red", pos = 4, cex = 0.9)
dev.off()

cat("图已保存到 figures/permutation_locus_level.pdf/png\n")

# --- 额外：输出BAG3在amyloidosis终点的结果 ---
cat("\n=== BAG3 在各终点的详细结果 ===\n")
bag3 <- locus_data %>% filter(id.exposure == "prot-a-225")
if (nrow(bag3) > 0) {
  cat("BAG3 (prot-a-225):\n")
  cat("  Amyloidosis: beta =", bag3$beta_ca, ", P =", bag3$pval_ca, "\n")
  cat("  CTS:         beta =", bag3$beta_cts, ", P =", bag3$pval_cts, "\n")
  cat("  HF:          beta =", bag3$beta_hf, ", P =", bag3$pval_hf, "\n")
} else {
  cat("⚠️ BAG3 (prot-a-225) 未在数据中找到\n")
  cat("请手动检查 results/04_mr_results_corrected.csv 中BAG3的amyloidosis结果\n")
}

cat("\n✅ 全部完成\n")
