# ============================================================
#  敏感性分析 (Sensitivity Analysis)
#  
#  对三角验证 + 双重验证的蛋白做：
#  1. 多方法 MR（IVW, MR-Egger, Weighted Median, Weighted Mode）
#  2. 异质性检验 (Cochran's Q)
#  3. 水平多效性检验 (MR-Egger Intercept)
#  4. 留一法 (Leave-one-out) + 画图
#  5. Steiger 方向性检验
#  
#  注意：只有 ≥3 个 SNP 的蛋白才能跑完整敏感性分析
#  1 个 SNP 的蛋白（Wald ratio）无法做这些检验，
#  这在论文里如实说明即可
# ============================================================

library(TwoSampleMR)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/PWMR")
dir.create("figures/sensitivity", showWarnings = FALSE)

cat("============================================\n")
cat(" 敏感性分析 (Sensitivity Analysis)\n")
cat("============================================\n\n")

# ---- 目标蛋白（三角 + 双重验证） ----
target_ids <- c(
  # 三角验证 (3/3)
  "prot-a-1599",  # JAKMIP3
  "prot-a-2722",  # SIGIRR
  "prot-a-575",   # CLEC2L
  # 双重验证 (2/3) 
  "prot-a-1003",  # Erythrocyte band 7
  "prot-a-3280",  # Zinc finger protein 843
  "prot-a-2798",  # Sorting nexin-1
  "prot-a-183",   # Set1/Ash2 histone methyltransferase
  "prot-a-2951",  # Tensin-2
  "prot-a-2999",  # Transmembrane protease serine 11D
  "prot-a-873",   # Dedicator of cytokinesis 9
  "prot-a-3168"   # Ubiquitin carboxyl-terminal hydrolase 25
)

# ---- 3 个结局 ----
outcomes <- c(
  "finn-b-E4_AMYLOIDOSIS",  # CA
  "ukb-d-G6_CARPTU",        # CTS
  "ebi-a-GCST009541"        # HF
)
outcome_labels <- c("Amyloidosis", "CTS", "HeartFailure")

# ---- 读取已有的 IV 数据 ----
cat("读取 IV 数据...\n")
all_exp <- readRDS("results/01_extraction_progress.rds")
target_exp <- all_exp %>% filter(id.exposure %in% target_ids)
cat("目标蛋白 IV 总数:", nrow(target_exp), "\n")
cat("覆盖蛋白数:", length(unique(target_exp$id.exposure)), "\n\n")

# 统计每个蛋白的 SNP 数
snp_counts <- target_exp %>% 
  group_by(id.exposure) %>% 
  summarise(
    protein = first(sub(" \\|\\|.*", "", exposure)),
    n_snp = n()
  ) %>%
  arrange(desc(n_snp))

cat("各蛋白 SNP 数量：\n")
for (i in 1:nrow(snp_counts)) {
  flag <- ifelse(snp_counts$n_snp[i] >= 3, "✓ 可做完整敏感性分析", 
          ifelse(snp_counts$n_snp[i] == 2, "△ 部分可做", "✗ 仅 Wald ratio"))
  cat(sprintf("  %s: %d SNPs  %s\n", snp_counts$protein[i], snp_counts$n_snp[i], flag))
}
cat("\n")

# ============================================================
# 对每个结局做敏感性分析
# ============================================================
all_mr_multi    <- data.frame()
all_hetero      <- data.frame()
all_pleio       <- data.frame()
all_steiger     <- data.frame()

for (oi in 1:length(outcomes)) {
  outc_id <- outcomes[oi]
  outc_label <- outcome_labels[oi]
  
  cat(sprintf("━━━ 结局: %s (%s) ━━━\n", outc_label, outc_id))
  
  # 提取结局数据
  cat("  提取结局数据...\n")
  outcome_dat <- tryCatch({
    extract_outcome_data(
      snps = unique(target_exp$SNP),
      outcomes = outc_id
    )
  }, error = function(e) {
    cat(sprintf("  ⚠️ 提取失败: %s\n", e$message))
    Sys.sleep(5)
    return(NULL)
  })
  
  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
    cat("  结局数据为空，跳过\n\n")
    next
  }
  
  # 数据对齐
  dat <- harmonise_data(exposure_dat = target_exp, outcome_dat = outcome_dat)
  cat("  对齐后数据:", nrow(dat), "行\n")
  
  # ---- 1. 多方法 MR ----
  cat("  [1] 多方法 MR...\n")
  mr_res <- mr(dat, method_list = c(
    "mr_ivw", "mr_egger_regression", 
    "mr_weighted_median", "mr_weighted_mode",
    "mr_wald_ratio"
  ))
  mr_res$outcome_label <- outc_label
  mr_res$protein_name <- sub(" \\|\\|.*", "", mr_res$exposure)
  all_mr_multi <- rbind(all_mr_multi, mr_res)
  
  # ---- 2. 异质性检验 (Cochran's Q) ----
  cat("  [2] 异质性检验...\n")
  hetero <- tryCatch({
    h <- mr_heterogeneity(dat)
    h$outcome_label <- outc_label
    h$protein_name <- sub(" \\|\\|.*", "", h$exposure)
    h
  }, error = function(e) {
    cat("    (无法计算，可能 SNP 不足)\n")
    data.frame()
  })
  if (nrow(hetero) > 0) all_hetero <- rbind(all_hetero, hetero)
  
  # ---- 3. 多效性检验 (MR-Egger Intercept) ----
  cat("  [3] 多效性检验...\n")
  pleio <- tryCatch({
    p <- mr_pleiotropy_test(dat)
    p$outcome_label <- outc_label
    p$protein_name <- sub(" \\|\\|.*", "", p$exposure)
    p
  }, error = function(e) {
    cat("    (无法计算，可能 SNP 不足)\n")
    data.frame()
  })
  if (nrow(pleio) > 0) all_pleio <- rbind(all_pleio, pleio)
  
  # ---- 4. 留一法 + 画图 ----
  cat("  [4] 留一法分析...\n")
  loo <- tryCatch({
    mr_leaveoneout(dat)
  }, error = function(e) {
    cat("    (无法计算)\n")
    NULL
  })
  
  if (!is.null(loo) && nrow(loo) > 0) {
    # 对每个有多 SNP 的蛋白画留一法图
    for (exp_id in unique(loo$id.exposure)) {
      loo_sub <- loo[loo$id.exposure == exp_id, ]
      if (nrow(loo_sub) > 2) {
        prot_nm <- sub(" \\|\\|.*", "", loo_sub$exposure[1])
        
        p_loo <- ggplot(loo_sub, aes(y = SNP, x = b)) +
          geom_point() +
          geom_errorbarh(aes(xmin = b - 1.96*se, xmax = b + 1.96*se), height = 0.2) +
          geom_vline(xintercept = 0, linetype = "dashed") +
          labs(
            title = sprintf("Leave-one-out: %s → %s", prot_nm, outc_label),
            x = "MR Effect (beta)", y = ""
          ) +
          theme_minimal(base_size = 12)
        
        fname <- sprintf("figures/sensitivity/loo_%s_%s.png", 
                         gsub("[^a-zA-Z0-9]", "", prot_nm), outc_label)
        ggsave(fname, p_loo, width = 8, height = max(4, nrow(loo_sub)*0.3), dpi = 200)
      }
    }
  }
  
  # ---- 5. 散点图 (Scatter plot) ----
  cat("  [5] 散点图...\n")
  for (exp_id in unique(dat$id.exposure)) {
    dat_sub <- dat[dat$id.exposure == exp_id, ]
    if (sum(dat_sub$mr_keep) >= 2) {
      prot_nm <- sub(" \\|\\|.*", "", dat_sub$exposure[1])
      
      p_scatter <- tryCatch({
        mr_scatter_plot(
          mr_res[mr_res$id.exposure == exp_id, ],
          dat_sub
        )[[1]] +
          labs(title = sprintf("%s → %s", prot_nm, outc_label)) +
          theme_minimal(base_size = 12)
      }, error = function(e) NULL)
      
      if (!is.null(p_scatter)) {
        fname <- sprintf("figures/sensitivity/scatter_%s_%s.png",
                         gsub("[^a-zA-Z0-9]", "", prot_nm), outc_label)
        ggsave(fname, p_scatter, width = 8, height = 6, dpi = 200)
      }
    }
  }
  
  # ---- 6. Steiger 方向性检验 ----
  cat("  [6] Steiger 方向性检验...\n")
  steiger <- tryCatch({
    s <- directionality_test(dat)
    s$outcome_label <- outc_label
    s$protein_name <- sub(" \\|\\|.*", "", s$exposure)
    s
  }, error = function(e) {
    cat("    (无法计算)\n")
    data.frame()
  })
  if (nrow(steiger) > 0) all_steiger <- rbind(all_steiger, steiger)
  
  cat("  ✅ 完成\n\n")
  Sys.sleep(2)
}

# ============================================================
# 汇总输出
# ============================================================
cat("\n╔══════════════════════════════════════════════╗\n")
cat("║        敏感性分析最终结果汇总                 ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")

# 多方法 MR 结果
cat("【1】多方法 MR 结果（方向一致性）\n")
cat("────────────────────────────────────────────\n")
mr_summary <- all_mr_multi %>%
  filter(!is.na(pval)) %>%
  select(protein_name, outcome_label, method, nsnp, b, pval) %>%
  mutate(direction = ifelse(b > 0, "+", "-"),
         sig = ifelse(pval < 0.05, "*", ""))

for (prot in unique(mr_summary$protein_name)) {
  cat(sprintf("\n  %s:\n", prot))
  sub <- mr_summary[mr_summary$protein_name == prot, ]
  for (i in 1:nrow(sub)) {
    r <- sub[i,]
    cat(sprintf("    %s | %-25s | beta=%6.3f %s | P=%.2e %s\n",
                r$outcome_label, r$method, r$b, r$direction, r$pval, r$sig))
  }
}

# 异质性
if (nrow(all_hetero) > 0) {
  cat("\n\n【2】异质性检验 (Cochran's Q)\n")
  cat("────────────────────────────────────────────\n")
  cat("  Q P > 0.05 表示无显著异质性 ✓\n\n")
  for (i in 1:nrow(all_hetero)) {
    r <- all_hetero[i,]
    verdict <- ifelse(r$Q_pval > 0.05, "✓ 无异质性", "⚠️ 有异质性")
    cat(sprintf("  %s → %s | %s | Q=%.2f, P=%.3f  %s\n",
                r$protein_name, r$outcome_label, r$method, r$Q, r$Q_pval, verdict))
  }
}

# 多效性
if (nrow(all_pleio) > 0) {
  cat("\n\n【3】MR-Egger 截距检验（水平多效性）\n")
  cat("────────────────────────────────────────────\n")
  cat("  P > 0.05 表示无显著多效性 ✓\n\n")
  for (i in 1:nrow(all_pleio)) {
    r <- all_pleio[i,]
    verdict <- ifelse(r$pval > 0.05, "✓ 无多效性", "⚠️ 有多效性")
    cat(sprintf("  %s → %s | intercept=%.4f, P=%.3f  %s\n",
                r$protein_name, r$outcome_label, r$egger_intercept, r$pval, verdict))
  }
}

# Steiger
if (nrow(all_steiger) > 0) {
  cat("\n\n【4】Steiger 方向性检验\n")
  cat("────────────────────────────────────────────\n")
  cat("  correct_causal_direction = TRUE 表示方向正确 ✓\n\n")
  for (i in 1:nrow(all_steiger)) {
    r <- all_steiger[i,]
    verdict <- ifelse(r$correct_causal_direction, "✓ 方向正确", "⚠️ 可能反向因果")
    cat(sprintf("  %s → %s | snp_r2.exposure=%.4f, snp_r2.outcome=%.6f | %s\n",
                r$protein_name, r$outcome_label, r$snp_r2.exposure, r$snp_r2.outcome, verdict))
  }
}

# ============================================================
# 保存
# ============================================================
write.csv(all_mr_multi, "results/13_sensitivity_mr_multi.csv", row.names = FALSE)
if (nrow(all_hetero) > 0) write.csv(all_hetero, "results/13_sensitivity_heterogeneity.csv", row.names = FALSE)
if (nrow(all_pleio) > 0)  write.csv(all_pleio, "results/13_sensitivity_pleiotropy.csv", row.names = FALSE)
if (nrow(all_steiger) > 0) write.csv(all_steiger, "results/13_sensitivity_steiger.csv", row.names = FALSE)

cat("\n\n============================================\n")
cat(" ✅ 敏感性分析全部完成！\n")
cat(" 结果文件:\n")
cat("   results/13_sensitivity_*.csv\n")
cat("   figures/sensitivity/*.png\n")
cat(" \n")
cat(" 📝 论文数据部分杀青！\n")
cat("    下一步：把所有 CSV 结果发给 Claude\n")
cat("    我来帮你搭建论文框架\n")
cat("============================================\n")
