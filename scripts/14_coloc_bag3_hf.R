# ============================================================
# 14_coloc_bag3_hf.R
# 目的：BAG3 cis-pQTL (rs2234962) vs Heart Failure 的共定位分析
# BAG3是cis-pQTL，coloc最适合这种情况
# ============================================================

# --- 安装 coloc（如果没有的话）---
if (!require("coloc", quietly = TRUE)) {
  install.packages("remotes")
  remotes::install_github("chr1swallace/coloc")
}

library(coloc)
library(ieugwasr)
library(dplyr)

setwd("~/Desktop/PWMR")

cat("=== BAG3 Colocalization Analysis ===\n\n")

# --- 参数设置 ---
# BAG3 基因位置（GRCh37/hg19，因为OpenGWAS用的是hg19）
bag3_chr <- 10
bag3_pos_center <- 121429633  # rs2234962 位置 (GRCh37)
window <- 500000  # ±500kb

region_start <- bag3_pos_center - window
region_end <- bag3_pos_center + window

cat("BAG3 区域: chr", bag3_chr, ":", region_start, "-", region_end, " (GRCh37)\n")

# --- BAG3 pQTL exposure ID ---
# 你需要确认这个ID是否正确
# 从你之前的数据中查找BAG3的exposure ID
bag3_pqtl_id <- "prot-a-225"  # BAG3在Sun 2018中的ID

# --- HF outcome ID ---
hf_id <- "ebi-a-GCST009541"  # HERMES heart failure

# --- 提取 pQTL 区域数据 ---
cat("\n--- 步骤1：提取BAG3 pQTL区域数据 ---\n")

pqtl_data <- tryCatch({
  associations(
    variants = paste0(bag3_chr, ":", region_start, "-", region_end),
    id = bag3_pqtl_id
  )
}, error = function(e) {
  cat("❌ pQTL数据提取失败:", e$message, "\n")
  cat("可能原因：OpenGWAS API不稳定，或者prot-a-225不支持区域查询\n")
  cat("\n尝试备选方案：用tophits提取...\n")
  
  # 备选：提取该蛋白的所有显著SNP
  tryCatch({
    tophits(id = bag3_pqtl_id)
  }, error = function(e2) {
    cat("❌ 备选方案也失败:", e2$message, "\n")
    return(NULL)
  })
})

if (is.null(pqtl_data) || nrow(pqtl_data) == 0) {
  cat("\n⚠️ 无法从API提取BAG3 pQTL区域数据\n")
  cat("这在Sun 2018 (prot-a) 数据中很常见——该批次不支持区域查询\n")
  cat("coloc需要整个区域的SNP统计量，而不仅仅是lead SNP\n\n")
  cat("=== 结论：Coloc无法完成 ===\n")
  cat("建议：从Methods中删除coloc，在Limitations中说明原因：\n")
  cat("  'Colocalization analysis was not performed as the Sun et al.\n")
  cat("   pQTL dataset does not provide regional summary statistics\n")
  cat("   through the OpenGWAS API, precluding locus-level comparison.'\n")
  
  # 保存失败记录
  coloc_result <- data.frame(
    analysis = "BAG3_vs_HF",
    status = "FAILED",
    reason = "pQTL regional data not available via API",
    recommendation = "Remove coloc from Methods, explain in Limitations"
  )
  write.csv(coloc_result, "results/18_coloc_result.csv", row.names = FALSE)
  cat("\n记录已保存到 results/18_coloc_result.csv\n")
  
  stop("Coloc无法完成，请参考上述建议修改论文")
}

cat("pQTL区域获得", nrow(pqtl_data), "个SNP\n")

# --- 提取 HF 同区域数据 ---
cat("\n--- 步骤2：提取HF GWAS同区域数据 ---\n")

hf_data <- tryCatch({
  associations(
    variants = paste0(bag3_chr, ":", region_start, "-", region_end),
    id = hf_id
  )
}, error = function(e) {
  cat("❌ HF数据提取失败:", e$message, "\n")
  return(NULL)
})

if (is.null(hf_data) || nrow(hf_data) == 0) {
  cat("⚠️ HF区域无数据\n")
  
  coloc_result <- data.frame(
    analysis = "BAG3_vs_HF",
    status = "FAILED",
    reason = "HF regional data not available",
    recommendation = "Remove coloc from Methods"
  )
  write.csv(coloc_result, "results/18_coloc_result.csv", row.names = FALSE)
  stop("HF区域数据获取失败")
}

cat("HF区域获得", nrow(hf_data), "个SNP\n")

# --- 取交集 ---
cat("\n--- 步骤3：匹配共有SNP ---\n")

common_snps <- intersect(pqtl_data$rsid, hf_data$rsid)
cat("共有SNP:", length(common_snps), "\n")

if (length(common_snps) < 50) {
  cat("⚠️ 共有SNP < 50，coloc结果可能不可靠\n")
}

if (length(common_snps) == 0) {
  cat("❌ 无共有SNP\n")
  coloc_result <- data.frame(
    analysis = "BAG3_vs_HF",
    status = "FAILED",
    reason = "No overlapping SNPs between pQTL and HF datasets",
    recommendation = "Remove coloc from Methods"
  )
  write.csv(coloc_result, "results/18_coloc_result.csv", row.names = FALSE)
  stop("无共有SNP")
}

# --- 准备 coloc 输入 ---
cat("\n--- 步骤4：运行coloc ---\n")

pqtl_sub <- pqtl_data %>% 
  filter(rsid %in% common_snps) %>%
  arrange(rsid) %>%
  distinct(rsid, .keep_all = TRUE)

hf_sub <- hf_data %>% 
  filter(rsid %in% common_snps) %>%
  arrange(rsid) %>%
  distinct(rsid, .keep_all = TRUE)

# 确保顺序一致
common_ordered <- intersect(pqtl_sub$rsid, hf_sub$rsid)
pqtl_sub <- pqtl_sub %>% filter(rsid %in% common_ordered) %>% arrange(rsid)
hf_sub <- hf_sub %>% filter(rsid %in% common_ordered) %>% arrange(rsid)

cat("最终用于coloc的SNP数:", nrow(pqtl_sub), "\n")

# 获取样本量
n_pqtl <- pqtl_sub$n[1]
n_hf <- hf_sub$n[1]
cat("pQTL样本量:", n_pqtl, "\n")
cat("HF样本量:", n_hf, "\n")

# HF case比例
hf_case_prop <- 47309 / (47309 + 930014)  # ≈ 0.0484

# 运行 coloc
result <- tryCatch({
  coloc.abf(
    dataset1 = list(
      beta = pqtl_sub$beta,
      varbeta = pqtl_sub$se^2,
      type = "quant",
      N = n_pqtl,
      snp = pqtl_sub$rsid,
      sdY = 1  # 标准化蛋白水平
    ),
    dataset2 = list(
      beta = hf_sub$beta,
      varbeta = hf_sub$se^2,
      type = "cc",
      N = n_hf,
      s = hf_case_prop,
      snp = hf_sub$rsid
    )
  )
}, error = function(e) {
  cat("❌ coloc运行出错:", e$message, "\n")
  return(NULL)
})

if (is.null(result)) {
  coloc_result <- data.frame(
    analysis = "BAG3_vs_HF",
    status = "FAILED",
    reason = "coloc.abf() error",
    recommendation = "Remove coloc from Methods"
  )
  write.csv(coloc_result, "results/18_coloc_result.csv", row.names = FALSE)
  stop("coloc运行失败")
}

# --- 输出结果 ---
cat("\n========================================\n")
cat("   BAG3 pQTL vs Heart Failure: Coloc Results\n")
cat("========================================\n\n")
cat("PP.H0 (无关联):          ", round(result$summary["PP.H0.abf"], 4), "\n")
cat("PP.H1 (仅pQTL有信号):   ", round(result$summary["PP.H1.abf"], 4), "\n")
cat("PP.H2 (仅HF有信号):     ", round(result$summary["PP.H2.abf"], 4), "\n")
cat("PP.H3 (两个独立信号):   ", round(result$summary["PP.H3.abf"], 4), "\n")
cat("PP.H4 (共享因果变异):   ", round(result$summary["PP.H4.abf"], 4), "\n")
cat("\n")

pp_h4 <- result$summary["PP.H4.abf"]
if (pp_h4 > 0.75) {
  cat("✅ PP.H4 > 0.75 → 强共定位证据！BAG3 pQTL和HF信号共享同一因果变异\n")
  cat("这可以加入正文作为重要验证\n")
} else if (pp_h4 > 0.5) {
  cat("⚠️ PP.H4 在 0.5-0.75 → 中等共定位证据\n")
  cat("可以在正文提及，但需要谨慎措辞\n")
} else {
  cat("❌ PP.H4 < 0.5 → 不支持共定位\n")
  cat("建议从正文删除coloc，或在补充材料中报告阴性结果\n")
}

# --- 保存 ---
coloc_result <- data.frame(
  analysis = "BAG3_vs_HF",
  status = "SUCCESS",
  n_snps = nrow(pqtl_sub),
  PP_H0 = round(result$summary["PP.H0.abf"], 4),
  PP_H1 = round(result$summary["PP.H1.abf"], 4),
  PP_H2 = round(result$summary["PP.H2.abf"], 4),
  PP_H3 = round(result$summary["PP.H3.abf"], 4),
  PP_H4 = round(result$summary["PP.H4.abf"], 4),
  interpretation = ifelse(pp_h4 > 0.75, "Strong colocalization",
                          ifelse(pp_h4 > 0.5, "Moderate colocalization",
                                 "No colocalization"))
)

write.csv(coloc_result, "results/18_coloc_result.csv", row.names = FALSE)
cat("\n结果已保存到 results/18_coloc_result.csv\n")

# --- 如果成功，画区域关联图 ---
if (pp_h4 > 0.3) {
  cat("\n--- 生成区域关联图 ---\n")
  
  # 合并数据用于画图
  plot_data <- data.frame(
    pos = pqtl_sub$position,
    pqtl_logp = -log10(pqtl_sub$p),
    hf_logp = -log10(hf_sub$p)
  )
  
  pdf("figures/coloc_bag3_hf_region.pdf", width = 10, height = 6)
  par(mfrow = c(2, 1), mar = c(2, 4, 2, 1))
  
  plot(plot_data$pos, plot_data$pqtl_logp, pch = 16, cex = 0.6, col = "steelblue",
       xlab = "", ylab = expression(-log[10](P)),
       main = paste0("BAG3 pQTL (chr10:", region_start, "-", region_end, ")"))
  abline(h = -log10(5e-8), lty = 2, col = "red")
  
  par(mar = c(4, 4, 1, 1))
  plot(plot_data$pos, plot_data$hf_logp, pch = 16, cex = 0.6, col = "coral",
       xlab = paste0("Chromosome ", bag3_chr, " position (GRCh37)"),
       ylab = expression(-log[10](P)),
       main = paste0("Heart Failure GWAS — PP.H4 = ", round(pp_h4, 3)))
  abline(h = -log10(5e-8), lty = 2, col = "red")
  
  dev.off()
  cat("图已保存到 figures/coloc_bag3_hf_region.pdf\n")
}

cat("\n✅ Coloc分析完成\n")
