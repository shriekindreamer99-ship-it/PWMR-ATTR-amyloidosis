# ============================================================
# 15_negative_control.R
# 目的：用生物学上无关联的三个疾病跑同样的框架，
#        证明不是随便什么病都能跑出三角验证信号
# 设计：CTS + Type 2 Diabetes + Asthma
#        （三个常见病，大样本，但无共享病理机制）
# 预期结果：0 个或极少三角验证蛋白 → 框架不是假阳性机器
# ============================================================

library(TwoSampleMR)
library(ieugwasr)
library(dplyr)

setwd("~/Desktop/PWMR")
cat("=== Negative Control: CTS + T2D + Asthma ===\n\n")

# --- 结局数据集 ---
# CTS: 同样用 ukb-d-G6_CARPTU (已有数据)
# T2D: ebi-a-GCST006867 (Mahajan 2018, ~75K cases, ~900K controls)
# Asthma: ukb-d-J10_ASTHMA (UK Biobank, large)

outcomes <- list(
  cts = "ukb-d-G6_CARPTU",
  t2d = "ebi-a-GCST006867",
  asthma = "ukb-d-J10_ASTHMA"
)

cat("Outcomes:\n")
cat("  CTS:", outcomes$cts, "\n")
cat("  T2D:", outcomes$t2d, "\n")
cat("  Asthma:", outcomes$asthma, "\n\n")

# --- 读取已有的 IV 数据 ---
cat("Loading IV extraction data...\n")
if (file.exists("results/01_extraction_progress.rds")) {
  iv_all <- readRDS("results/01_extraction_progress.rds")
  proteins <- unique(iv_all$id.exposure)
  cat("  Loaded", length(proteins), "proteins with IVs\n\n")
} else {
  stop("Need results/01_extraction_progress.rds from the main analysis")
}

# --- 对 T2D 和 Asthma 跑 MR ---
# CTS 结果已有 (results/09_mr_cts_results.csv)

cat("Reading existing CTS results...\n")
cts_results <- read.csv("results/09_mr_cts_results.csv")
cat("  CTS results:", nrow(cts_results), "protein-outcome pairs\n\n")

# T2D MR
cat("Running MR for T2D (this may take 30-60 minutes)...\n")
t2d_results <- data.frame()

for (i in seq_along(proteins)) {
  pid <- proteins[i]
  if (i %% 100 == 0) cat("  T2D progress:", i, "/", length(proteins), "\n")
  
  iv_sub <- iv_all %>% filter(id.exposure == pid)
  
  tryCatch({
    outcome_data <- extract_outcome_data(
      snps = iv_sub$SNP,
      outcomes = outcomes$t2d
    )
    
    if (!is.null(outcome_data) && nrow(outcome_data) > 0) {
      dat <- harmonise_data(iv_sub, outcome_data, action = 2)
      if (nrow(dat) > 0) {
        mr_res <- mr(dat, method_list = c("mr_ivw", "mr_wald_ratio"))
        if (nrow(mr_res) > 0) {
          mr_res$protein_id <- pid
          mr_res$protein_name <- iv_sub$exposure[1]
          t2d_results <- bind_rows(t2d_results, mr_res)
        }
      }
    }
  }, error = function(e) {})
}

cat("  T2D done:", nrow(t2d_results), "results\n\n")

# Asthma MR
cat("Running MR for Asthma (this may take 30-60 minutes)...\n")
asthma_results <- data.frame()

for (i in seq_along(proteins)) {
  pid <- proteins[i]
  if (i %% 100 == 0) cat("  Asthma progress:", i, "/", length(proteins), "\n")
  
  iv_sub <- iv_all %>% filter(id.exposure == pid)
  
  tryCatch({
    outcome_data <- extract_outcome_data(
      snps = iv_sub$SNP,
      outcomes = outcomes$asthma
    )
    
    if (!is.null(outcome_data) && nrow(outcome_data) > 0) {
      dat <- harmonise_data(iv_sub, outcome_data, action = 2)
      if (nrow(dat) > 0) {
        mr_res <- mr(dat, method_list = c("mr_ivw", "mr_wald_ratio"))
        if (nrow(mr_res) > 0) {
          mr_res$protein_id <- pid
          mr_res$protein_name <- iv_sub$exposure[1]
          asthma_results <- bind_rows(asthma_results, mr_res)
        }
      }
    }
  }, error = function(e) {})
}

cat("  Asthma done:", nrow(asthma_results), "results\n\n")

# --- 保存中间结果 ---
write.csv(t2d_results, "results/19_mr_t2d_results.csv", row.names = FALSE)
write.csv(asthma_results, "results/20_mr_asthma_results.csv", row.names = FALSE)

# --- 三角验证（负对照版） ---
cat("=== Negative Control Triangulation ===\n\n")

# 合并三个结局
neg_tri <- cts_results %>%
  select(id.exposure, b, pval) %>%
  rename(beta_cts = b, pval_cts = pval) %>%
  inner_join(
    t2d_results %>% select(protein_id, b, pval) %>%
      rename(id.exposure = protein_id, beta_t2d = b, pval_t2d = pval),
    by = "id.exposure"
  ) %>%
  inner_join(
    asthma_results %>% select(protein_id, b, pval) %>%
      rename(id.exposure = protein_id, beta_asthma = b, pval_asthma = pval),
    by = "id.exposure"
  )

cat("Proteins with all 3 outcomes:", nrow(neg_tri), "\n")

# 筛选三角验证
neg_triple <- neg_tri %>%
  filter(
    pval_cts < 0.05 & pval_t2d < 0.05 & pval_asthma < 0.05 &
    sign(beta_cts) == sign(beta_t2d) & sign(beta_cts) == sign(beta_asthma)
  )

cat("\n========================================\n")
cat("NEGATIVE CONTROL RESULT\n")
cat("========================================\n")
cat("Triple-validated proteins (CTS + T2D + Asthma):", nrow(neg_triple), "\n")

if (nrow(neg_triple) == 0) {
  cat("\n✅ PERFECT: Zero proteins passed the negative control\n")
  cat("   triplet. This confirms the framework does NOT produce\n")
  cat("   false positives from arbitrary disease combinations.\n")
} else if (nrow(neg_triple) <= 2) {
  cat("\n⚠️ Low count:", nrow(neg_triple), "proteins.\n")
  cat("   Check if these are known pleiotropic inflammation markers.\n")
  print(neg_triple)
} else {
  cat("\n❌ Unexpectedly high:", nrow(neg_triple), "proteins.\n")
  cat("   This may weaken the specificity argument.\n")
  print(neg_triple)
}

# Expected by chance: 0.05^3 * 0.25 * N_proteins
n_tested <- nrow(neg_tri)
expected_by_chance <- n_tested * 0.05^3 * 0.25
cat("\nExpected by chance alone:", round(expected_by_chance, 2), "\n")
cat("Observed:", nrow(neg_triple), "\n")

# --- 保存 ---
write.csv(neg_tri, "results/21_negative_control_triangulation.csv", row.names = FALSE)

result_summary <- data.frame(
  control_triplet = "CTS + T2D + Asthma",
  n_proteins_tested = n_tested,
  n_triple_validated = nrow(neg_triple),
  expected_by_chance = round(expected_by_chance, 2)
)
write.csv(result_summary, "results/22_negative_control_summary.csv", row.names = FALSE)

cat("\nResults saved to results/21-22_negative_control_*.csv\n")
cat("=== Done ===\n")
