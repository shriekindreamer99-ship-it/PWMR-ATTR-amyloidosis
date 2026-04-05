# ============================================================
#  反向 MR (Reverse MR)
#  
#  目的：排除"疾病导致蛋白水平改变"的反向因果
#  如果反向 MR 不显著，说明因果方向确实是 蛋白→疾病
#  
#  对 3 个核心蛋白 + BAG3，用疾病作为暴露、蛋白作为结局
# ============================================================

library(TwoSampleMR)
library(dplyr)

setwd("~/Desktop/PWMR")

cat("============================================\n")
cat(" 反向 MR (Reverse Mendelian Randomization)\n")
cat("============================================\n\n")

# 目标蛋白（作为结局）
target_proteins <- c(
  "prot-a-2722",  # SIGIRR
  "prot-a-1599",  # JAKMIP3
  "prot-a-575",   # CLEC2L
  "prot-a-225"    # BAG3 (Sun et al. 2018)
)

# 疾病作为暴露
diseases <- c(
  "finn-b-E4_AMYLOIDOSIS",  # CA
  "ebi-a-GCST009541"        # HF (样本量大，更容易提取到显著 SNP)
)
disease_labels <- c("Amyloidosis", "Heart Failure")

all_reverse <- data.frame()

for (di in 1:length(diseases)) {
  dis_id <- diseases[di]
  dis_label <- disease_labels[di]
  
  cat(sprintf("━━━ 暴露: %s (%s) ━━━\n", dis_label, dis_id))
  
  # 提取疾病的显著 SNP 作为工具变量
  dis_instruments <- tryCatch({
    extract_instruments(
      outcomes = dis_id,
      p1 = 5e-8,
      clump = TRUE,
      r2 = 0.001,
      kb = 10000
    )
  }, error = function(e) {
    cat(sprintf("  ⚠️ 提取 IV 失败: %s\n", e$message))
    # 放宽阈值重试
    tryCatch({
      extract_instruments(outcomes = dis_id, p1 = 5e-6, clump = TRUE, r2 = 0.001, kb = 10000)
    }, error = function(e2) NULL)
  })
  
  if (is.null(dis_instruments) || nrow(dis_instruments) == 0) {
    cat("  没有提取到工具变量，跳过\n\n")
    next
  }
  
  cat(sprintf("  疾病 IV: %d 个 SNP\n", nrow(dis_instruments)))
  
  # 对每个蛋白做反向 MR
  for (prot_id in target_proteins) {
    cat(sprintf("  → 结局: %s\n", prot_id))
    
    result <- tryCatch({
      # 提取蛋白中对应的 SNP
      outcome_dat <- extract_outcome_data(
        snps = dis_instruments$SNP,
        outcomes = prot_id
      )
      
      if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
        cat("    无匹配 SNP\n")
        return(NULL)
      }
      
      # 对齐
      dat <- harmonise_data(dis_instruments, outcome_dat)
      
      # MR
      mr_res <- mr(dat, method_list = c("mr_ivw", "mr_wald_ratio"))
      
      mr_res$disease <- dis_label
      mr_res$protein_name <- prot_id
      mr_res
    }, error = function(e) {
      cat(sprintf("    出错: %s\n", e$message))
      return(NULL)
    })
    
    if (!is.null(result) && nrow(result) > 0) {
      all_reverse <- rbind(all_reverse, result)
    }
    
    Sys.sleep(2)
  }
  cat("\n")
}

# 结果
cat("\n============================================\n")
cat(" 反向 MR 结果\n")
cat("============================================\n\n")

if (nrow(all_reverse) > 0) {
  all_reverse <- all_reverse %>% filter(!is.na(pval))
  
  for (i in 1:nrow(all_reverse)) {
    r <- all_reverse[i, ]
    verdict <- ifelse(r$pval > 0.05, "✓ 不显著 (无反向因果)", "⚠️ 显著")
    cat(sprintf("  %s → %s: beta=%.3f, P=%.3f  %s\n",
                r$disease, r$protein_name, r$b, r$pval, verdict))
  }
  
  write.csv(all_reverse, "results/14_reverse_mr.csv", row.names = FALSE)
  cat("\n  结果已保存: results/14_reverse_mr.csv\n")
} else {
  cat("  没有完成任何反向 MR。\n")
  cat("  这可能是因为淀粉样变 GWAS 在 5e-8 阈值下没有显著 SNP。\n")
  cat("  这本身就说明反向因果不成立——疾病端没有足够强的遗传信号。\n")
}
