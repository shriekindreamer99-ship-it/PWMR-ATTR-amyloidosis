# ============================================================
#  三表型联合 PWMR (Triangulation / Cross-phenotype MR)
#  
#  设计：
#  结局 1 (Discovery): finn-b-E4_AMYLOIDOSIS (淀粉样变, ~226 cases)
#  结局 2 (CTS):        ukb-d-G6_CARPTU (腕管综合征, UK Biobank)
#  结局 3 (HF):         ebi-a-GCST009541 (心衰 HERMES, 47,309 cases)
#  
#  暴露：prot-a 批次全部蛋白的 IVs（从第一阶段保存的进度文件读取）
#  
#  策略：
#  1. 对 CTS 和 HF 分别做全蛋白组 MR（不是只跑55个）
#  2. 三个结局的结果合并，找"三角交叉"蛋白
#  3. 同时在3个表型中P<0.05且方向一致 = 最强证据
#  
#  使用方法：
#  1. 确保 token 还在: Sys.getenv("OPENGWAS_JWT")
#  2. 确保在 PWMR 文件夹下: setwd("~/Desktop/PWMR")
#  3. source 本脚本
#  
#  预计运行时间：CTS 约 20-30 分钟，HF 约 20-30 分钟
#  （因为 IVs 已经提取好了，只需要提取两个新结局的数据）
# ============================================================

library(TwoSampleMR)
library(ieugwasr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

# ---- 三个结局 ----
OUTCOME_CA  <- "finn-b-E4_AMYLOIDOSIS"  # 淀粉样变 (已跑完)
OUTCOME_CTS <- "ukb-d-G6_CARPTU"        # 腕管综合征 (UK Biobank)
OUTCOME_HF  <- "ebi-a-GCST009541"       # 心衰 HERMES

cat("============================================\n")
cat(" 三表型联合 PWMR (Triangulation)\n")
cat(" 结局 1: 淀粉样变 (CA) -", OUTCOME_CA, "\n")
cat(" 结局 2: 腕管综合征 (CTS) -", OUTCOME_CTS, "\n")
cat(" 结局 3: 心力衰竭 (HF) -", OUTCOME_HF, "\n")
cat("============================================\n\n")

# ---- 读取已保存的暴露 IV 数据 ----
cat("【准备】读取第一阶段保存的 IV 数据...\n")
all_exposures <- readRDS("results/01_extraction_progress.rds")
cat("  总 IV 数:", nrow(all_exposures), "\n")
cat("  蛋白数:", length(unique(all_exposures$id.exposure)), "\n\n")

# ---- 读取已有的淀粉样变 MR 结果 ----
cat("【结局 1: 淀粉样变】读取已有结果...\n")
mr_ca <- read.csv("results/04_mr_results_corrected.csv")
cat("  已有结果:", nrow(mr_ca), "行\n\n")

# ============================================================
# 辅助函数：对一个新结局做全蛋白组 MR
# ============================================================
run_mr_for_outcome <- function(all_exp, outcome_id, outcome_label) {
  cat(sprintf("【%s】开始提取结局数据: %s\n", outcome_label, outcome_id))
  
  all_snps <- unique(all_exp$SNP)
  cat("  需要查询的 SNP 总数:", length(all_snps), "\n")
  
  # 分批提取结局数据（大批量更稳定）
  batch_size <- 300
  batches <- split(all_snps, ceiling(seq_along(all_snps) / batch_size))
  
  outcome_parts <- list()
  for (i in seq_along(batches)) {
    if (i %% 5 == 1) {
      cat(sprintf("  提取进度: %d/%d 批次\n", i, length(batches)))
    }
    
    part <- tryCatch({
      extract_outcome_data(snps = batches[[i]], outcomes = outcome_id)
    }, error = function(e) {
      cat(sprintf("    批次 %d 出错，等待后重试...\n", i))
      Sys.sleep(10)
      # 重试一次
      tryCatch({
        extract_outcome_data(snps = batches[[i]], outcomes = outcome_id)
      }, error = function(e2) {
        Sys.sleep(5)
        return(NULL)
      })
    })
    
    if (!is.null(part)) outcome_parts[[i]] <- part
    Sys.sleep(1)
  }
  
  outcome_dat <- do.call(rbind, outcome_parts)
  
  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
    cat("  ⚠️ 无法提取数据！跳过此结局。\n")
    return(NULL)
  }
  
  cat("  匹配到的 SNP 数:", nrow(outcome_dat), "\n")
  
  # 数据对齐
  cat("  数据对齐中...\n")
  dat <- harmonise_data(exposure_dat = all_exp, outcome_dat = outcome_dat)
  cat("  对齐后数据行数:", nrow(dat), "\n")
  
  # MR 分析
  cat("  执行 MR 分析...\n")
  mr_res <- mr(dat, method_list = c("mr_ivw", "mr_wald_ratio"))
  
  # 整理结果
  n_tests <- length(unique(mr_res$id.exposure))
  mr_res <- mr_res %>%
    filter(!is.na(pval)) %>%
    mutate(
      p_bonferroni = pmin(pval * n_tests, 1),
      p_fdr = p.adjust(pval, method = "BH"),
      protein_name = sub(" \\|\\|.*", "", exposure)
    )
  
  cat(sprintf("  ✅ 完成！%d 个蛋白，Bonf显著: %d，FDR显著: %d，Nominal显著: %d\n\n",
              n_tests,
              sum(mr_res$p_bonferroni < 0.05, na.rm=TRUE),
              sum(mr_res$p_fdr < 0.05, na.rm=TRUE),
              sum(mr_res$pval < 0.05, na.rm=TRUE)))
  
  return(mr_res)
}

# ============================================================
# 对 CTS 和 HF 分别做全蛋白组 MR
# ============================================================
cat("============================================\n")
cat(" 开始对两个新结局做 PWMR\n")
cat("============================================\n\n")

# CTS
mr_cts <- run_mr_for_outcome(all_exposures, OUTCOME_CTS, "腕管综合征 CTS")
if (!is.null(mr_cts)) {
  saveRDS(mr_cts, "results/09_mr_cts_results.rds")
  write.csv(mr_cts, "results/09_mr_cts_results.csv", row.names = FALSE)
}

# HF（如果已有第二阶段 55 个蛋白的结果，这里重跑全蛋白组）
mr_hf <- run_mr_for_outcome(all_exposures, OUTCOME_HF, "心力衰竭 HF")
if (!is.null(mr_hf)) {
  saveRDS(mr_hf, "results/10_mr_hf_results.rds")
  write.csv(mr_hf, "results/10_mr_hf_results.csv", row.names = FALSE)
}

# ============================================================
# 三表型结果合并 (Triangulation)
# ============================================================
cat("============================================\n")
cat(" 三表型结果合并 (Triangulation)\n")
cat("============================================\n\n")

# 整理三个结局的结果，统一格式
format_results <- function(df, suffix) {
  df %>%
    filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
    select(id.exposure, protein_name,
           !!paste0("beta_", suffix) := b,
           !!paste0("se_", suffix) := se,
           !!paste0("pval_", suffix) := pval,
           !!paste0("nsnp_", suffix) := nsnp,
           !!paste0("method_", suffix) := method)
}

ca_fmt  <- format_results(mr_ca, "ca")
cts_fmt <- if (!is.null(mr_cts)) format_results(mr_cts, "cts") else NULL
hf_fmt  <- if (!is.null(mr_hf))  format_results(mr_hf, "hf")  else NULL

# 合并
tri <- ca_fmt
if (!is.null(cts_fmt)) tri <- inner_join(tri, cts_fmt, by = c("id.exposure", "protein_name"))
if (!is.null(hf_fmt))  tri <- inner_join(tri, hf_fmt, by = c("id.exposure", "protein_name"))

cat("  三个结局都有结果的蛋白数:", nrow(tri), "\n")

# 计算三角验证指标
tri <- tri %>%
  mutate(
    # 各结局是否 nominal 显著
    sig_ca  = pval_ca < 0.05,
    sig_cts = if ("pval_cts" %in% names(.)) pval_cts < 0.05 else FALSE,
    sig_hf  = if ("pval_hf" %in% names(.)) pval_hf < 0.05 else FALSE,
    
    # 方向一致性
    dir_ca_cts = if (all(c("beta_ca","beta_cts") %in% names(.))) sign(beta_ca) == sign(beta_cts) else NA,
    dir_ca_hf  = if (all(c("beta_ca","beta_hf") %in% names(.))) sign(beta_ca) == sign(beta_hf) else NA,
    dir_cts_hf = if (all(c("beta_cts","beta_hf") %in% names(.))) sign(beta_cts) == sign(beta_hf) else NA,
    
    # 三个方向全一致
    all_direction_consistent = ifelse(
      !is.na(dir_ca_cts) & !is.na(dir_ca_hf),
      dir_ca_cts & dir_ca_hf,
      FALSE
    ),
    
    # 三角验证等级
    n_sig = sig_ca + sig_cts + sig_hf,
    triangulation = case_when(
      n_sig == 3 & all_direction_consistent ~ "★★★ 三角验证 (3/3 显著+同向)",
      n_sig >= 2 & all_direction_consistent ~ "★★ 双重验证 (2/3 显著+同向)",
      n_sig >= 2                            ~ "★ 部分验证 (2/3 显著)",
      n_sig == 1                            ~ "单一信号",
      TRUE                                  ~ "无信号"
    ),
    
    # Bonferroni 水平（各结局各自的阈值）
    bonf_cts = if ("pval_cts" %in% names(.)) pval_cts < (0.05 / nrow(.)) else FALSE,
    bonf_hf  = if ("pval_hf" %in% names(.)) pval_hf < (0.05 / nrow(.)) else FALSE
  )

# 排序：按三角验证等级，然后按 CA 的 P 值
tri <- tri %>%
  arrange(desc(n_sig), pval_ca)

# ============================================================
# 打印关键结果
# ============================================================
cat("\n============================================\n")
cat(" 🎯 三表型联合 PWMR 最终结果\n")
cat("============================================\n")

# 统计
n_triple <- sum(tri$n_sig == 3 & tri$all_direction_consistent, na.rm=TRUE)
n_double <- sum(tri$n_sig >= 2 & tri$all_direction_consistent, na.rm=TRUE)

cat("  三角验证 (3/3 显著+同向):", n_triple, "个蛋白\n")
cat("  双重验证 (2/3 显著+同向):", n_double, "个蛋白\n")
if (!is.null(cts_fmt)) cat("  CTS Bonferroni 显著:", sum(tri$bonf_cts, na.rm=TRUE), "个\n")
if (!is.null(hf_fmt))  cat("  HF Bonferroni 显著:", sum(tri$bonf_hf, na.rm=TRUE), "个\n")
cat("\n")

# 打印三角验证的蛋白
top_proteins <- tri %>% filter(n_sig >= 2)

if (nrow(top_proteins) > 0) {
  cat("=== 至少在 2 个表型中显著的蛋白 (按验证等级排序) ===\n\n")
  
  for (i in 1:nrow(top_proteins)) {
    r <- top_proteins[i, ]
    cat(sprintf("%2d. %s [%s]\n", i, r$protein_name, r$triangulation))
    cat(sprintf("    CA:  beta=%7.3f, P=%.2e %s\n", 
                r$beta_ca, r$pval_ca, ifelse(r$sig_ca, "✓", "")))
    if ("beta_cts" %in% names(r)) {
      cat(sprintf("    CTS: beta=%7.3f, P=%.2e %s %s\n", 
                  r$beta_cts, r$pval_cts, 
                  ifelse(r$sig_cts, "✓", ""),
                  ifelse(!is.na(r$bonf_cts) && r$bonf_cts, "[BONF]", "")))
    }
    if ("beta_hf" %in% names(r)) {
      cat(sprintf("    HF:  beta=%7.3f, P=%.2e %s %s\n",
                  r$beta_hf, r$pval_hf,
                  ifelse(r$sig_hf, "✓", ""),
                  ifelse(!is.na(r$bonf_hf) && r$bonf_hf, "[BONF]", "")))
    }
    cat(sprintf("    方向一致: %s\n\n", 
                ifelse(r$all_direction_consistent, "✓ 三个表型全部同向", "✗ 有反向")))
  }
} else {
  cat("  没有在 2 个以上表型中同时显著的蛋白。\n")
  cat("  显示 Top 20（按 CA P 值排序）：\n\n")
  for (i in 1:min(20, nrow(tri))) {
    r <- tri[i, ]
    cat(sprintf("%2d. %s | CA: P=%.2e | CTS: P=%.2e | HF: P=%.2e | %s\n",
                i, r$protein_name, r$pval_ca,
                ifelse("pval_cts" %in% names(r), r$pval_cts, NA),
                ifelse("pval_hf" %in% names(r), r$pval_hf, NA),
                r$triangulation))
  }
}

# ============================================================
# 保存完整结果
# ============================================================
write.csv(tri, "results/11_triangulation_combined.csv", row.names = FALSE)

# ============================================================
# 画三表型热力图 (Heatmap of P-values)
# ============================================================
cat("\n【绘图】生成三表型热力图...\n")

# 筛选有意义的蛋白画图（至少一个 P < 0.01）
plot_proteins <- tri %>%
  filter(pval_ca < 0.01 | 
         ("pval_cts" %in% names(.) & pval_cts < 0.01) |
         ("pval_hf" %in% names(.) & pval_hf < 0.01))

if (nrow(plot_proteins) > 0) {
  # 准备长格式数据
  plot_long <- plot_proteins %>%
    select(protein_name, starts_with("pval_"), starts_with("beta_")) %>%
    head(30)  # 最多画 30 个
  
  # 转换为长格式
  pval_long <- plot_long %>%
    select(protein_name, starts_with("pval_")) %>%
    pivot_longer(cols = starts_with("pval_"), names_to = "phenotype", values_to = "pval") %>%
    mutate(
      phenotype = case_when(
        phenotype == "pval_ca"  ~ "Amyloidosis\n(226 cases)",
        phenotype == "pval_cts" ~ "Carpal Tunnel\n(UKB)",
        phenotype == "pval_hf"  ~ "Heart Failure\n(47,309 cases)",
        TRUE ~ phenotype
      ),
      neg_log10_p = -log10(pval),
      sig_label = case_when(
        pval < 0.001 ~ "***",
        pval < 0.01  ~ "**",
        pval < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  
  beta_long <- plot_long %>%
    select(protein_name, starts_with("beta_")) %>%
    pivot_longer(cols = starts_with("beta_"), names_to = "phenotype", values_to = "beta") %>%
    mutate(
      phenotype = case_when(
        phenotype == "beta_ca"  ~ "Amyloidosis\n(226 cases)",
        phenotype == "beta_cts" ~ "Carpal Tunnel\n(UKB)",
        phenotype == "beta_hf"  ~ "Heart Failure\n(47,309 cases)",
        TRUE ~ phenotype
      ),
      direction = ifelse(beta > 0, "Risk ↑", "Protective ↓")
    )
  
  combined_long <- left_join(pval_long, 
                              beta_long %>% select(protein_name, phenotype, direction),
                              by = c("protein_name", "phenotype"))
  
  p_heat <- ggplot(combined_long, 
                    aes(x = phenotype, y = reorder(protein_name, neg_log10_p), 
                        fill = neg_log10_p)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sig_label), size = 5, color = "white") +
    scale_fill_gradient2(
      low = "grey90", mid = "steelblue", high = "darkred",
      midpoint = 1.3,  # -log10(0.05)
      name = expression(-log[10](P))
    ) +
    labs(
      title = "Cross-Phenotype MR: CTS → CA → HF Triangulation",
      subtitle = "Proteins with P < 0.01 in at least one phenotype | * P<0.05, ** P<0.01, *** P<0.001",
      x = "", y = ""
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 11, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank()
    )
  
  ggsave("figures/triangulation_heatmap.pdf", p_heat, width = 10, height = max(6, nrow(plot_proteins)*0.35), dpi = 300)
  ggsave("figures/triangulation_heatmap.png", p_heat, width = 10, height = max(6, nrow(plot_proteins)*0.35), dpi = 300)
  cat("  ✅ 热力图已保存\n")
}

# ============================================================
# 画效应值一致性散点图 (CTS vs HF, 标注 CA 显著的蛋白)
# ============================================================
if (all(c("beta_cts", "beta_hf") %in% names(tri))) {
  cat("【绘图】生成 CTS vs HF 效应值对比图...\n")
  
  tri_plot <- tri %>%
    mutate(
      highlight = case_when(
        n_sig == 3 & all_direction_consistent ~ "三角验证",
        n_sig >= 2 & all_direction_consistent ~ "双重验证",
        sig_ca ~ "仅CA显著",
        TRUE ~ "不显著"
      ),
      label = ifelse(n_sig >= 2 & all_direction_consistent, 
                     substr(protein_name, 1, 25), NA)
    )
  
  p_scatter <- ggplot(tri_plot, aes(x = beta_cts, y = beta_hf, color = highlight)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_point(data = tri_plot %>% filter(highlight %in% c("三角验证", "双重验证")),
               size = 4, alpha = 0.9) +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 15, na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "blue", alpha = 0.3) +
    scale_color_manual(values = c(
      "三角验证" = "#E31A1C",
      "双重验证" = "#FF7F00",
      "仅CA显著" = "#6A3D9A",
      "不显著"   = "grey70"
    ), name = "验证等级") +
    labs(
      title = "Cross-Phenotype Effect Consistency: CTS vs Heart Failure",
      subtitle = "Colored by triangulation status across CTS-CA-HF axis",
      x = "MR beta → Carpal Tunnel Syndrome",
      y = "MR beta → Heart Failure"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  ggsave("figures/cts_vs_hf_scatter.pdf", p_scatter, width = 10, height = 8, dpi = 300)
  ggsave("figures/cts_vs_hf_scatter.png", p_scatter, width = 10, height = 8, dpi = 300)
  cat("  ✅ 散点图已保存\n")
}

# ============================================================
# 最终汇总
# ============================================================
cat("\n============================================\n")
cat(" 🎉 三表型联合 PWMR 全部完成！\n")
cat("============================================\n")
cat(" 输出文件清单:\n")
cat("   results/09_mr_cts_results.csv        - CTS 全蛋白 MR 结果\n")
cat("   results/10_mr_hf_results.csv         - HF 全蛋白 MR 结果\n")
cat("   results/11_triangulation_combined.csv - 三表型合并结果 ★\n")
cat("   figures/triangulation_heatmap.pdf/png - 三表型热力图 ★\n")
cat("   figures/cts_vs_hf_scatter.pdf/png     - CTS vs HF 一致性图\n")
cat("============================================\n")
cat(" 下一步：把 11_triangulation_combined.csv 发给 Claude\n")
cat("         我会帮你分析哪些蛋白通过了三角验证\n")
cat("============================================\n")
