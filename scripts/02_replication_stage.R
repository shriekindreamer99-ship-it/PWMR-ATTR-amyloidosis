# ============================================================
#  第二阶段 MR 验证 (Replication Stage)
#  
#  设计思路：
#  - Discovery: finn-b-E4_AMYLOIDOSIS (淀粉样变, 226 cases)
#    → 筛出 55 个 nominal P < 0.05 的候选蛋白
#  - Replication: ebi-a-GCST009541 (HERMES 心衰, 47,309 cases)
#    → 在大样本心衰队列中验证这 55 个蛋白
#  
#  校正阈值：0.05 / 55 = 9.09e-4
#  
#  使用方法：先确保 token 还在，然后直接 source 这个脚本
#  Sys.setenv(OPENGWAS_JWT = "你的token")
# ============================================================

library(TwoSampleMR)
library(ieugwasr)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ---- 设置 ----
# HERMES 心衰 GWAS: 47,309 cases / 930,014 controls (欧洲人群)
# 这是目前 OpenGWAS 上最大的心衰队列之一
REPLICATION_OUTCOME <- "ebi-a-GCST009541"

# Bonferroni 校正阈值 (55 个蛋白)
N_CANDIDATES <- 55
REPL_BONF_THRESH <- 0.05 / N_CANDIDATES  # = 9.09e-4

cat("============================================\n")
cat(" 第二阶段 MR 验证 (Replication)\n")
cat(" Discovery: finn-b-E4_AMYLOIDOSIS (淀粉样变)\n")
cat(" Replication:", REPLICATION_OUTCOME, "(心衰 HERMES)\n")
cat(" 候选蛋白数:", N_CANDIDATES, "\n")
cat(" 校正阈值: P <", formatC(REPL_BONF_THRESH, format="e", digits=2), "\n")
cat("============================================\n\n")

# ---- 55 个候选蛋白 ID (Discovery P < 0.05) ----
candidate_ids <- c(
  "prot-a-2601",  # Ribonucleoside-diphosphate reductase large subunit
  "prot-a-268",   # BPI fold-containing family B member 1
  "prot-a-250",   # Myc box-dependent-interacting protein 1
  "prot-a-2120",  # Netrin-G1
  "prot-a-1513",  # Interleukin-23 receptor
  "prot-a-2319",  # GDP-fucose protein O-fucosyltransferase 2
  "prot-a-1229",  # Glycine N-methyltransferase
  "prot-a-170",   # Arrestin domain-containing protein 3
  "prot-a-330",   # Carbonic anhydrase 3
  "prot-a-3199",  # Vascular endothelial growth factor C
  "prot-a-2361",  # Prolargin
  "prot-a-132",   # Apolipoprotein E (isoform E2)
  "prot-a-1329",  # Probable E3 ubiquitin-protein ligase HERC1
  "prot-a-2281",  # Polycystin-2
  "prot-a-3280",  # Zinc/RING finger protein 4
  "prot-a-161",   # ADP-ribosylation factor-like protein 1
  "prot-a-1586",  # Inter-alpha-trypsin inhibitor heavy chain H1
  "prot-a-2050",  # Nidogen-2
  "prot-a-971",   # Endoplasmic reticulum aminopeptidase 1
  "prot-a-1357",  # Non-histone chromosomal protein HMG-14
  "prot-a-2898",  # SUN domain-containing protein 3
  "prot-a-1792",  # Leucine-rich repeat neuronal protein 1
  "prot-a-2242",  # PDZK1-interacting protein 1
  "prot-a-1643",  # Killer cell immunoglobulin-like receptor 2DL5A
  "prot-a-2983",  # Tight junction protein ZO-1
  "prot-a-366",   # Calpastatin
  "prot-a-1503",  # Interleukin-1 receptor-like 2
  "prot-a-207",   # Axin-2
  "prot-a-2380",  # Activated Protein C
  "prot-a-1738",  # Leukocyte immunoglobulin-like receptor subfamily A member 4
  "prot-a-1273",  # Gremlin-1
  "prot-a-2267",  # PILR alpha-associated neural protein
  "prot-a-3175",  # Cytochrome b-c1 complex subunit 7
  "prot-a-1834",  # ER mannosyl-oligosaccharide 1,2-alpha-mannosidase
  "prot-a-1599",  # Janus kinase and microtubule-interacting protein 3
  "prot-a-3073",  # TOM1-like protein 1
  "prot-a-661",   # Cysteine-rich secretory protein 2
  "prot-a-1347",  # HLA class II histocompatibility antigen, DQ alpha 2 chain
  "prot-a-2722",  # Single Ig IL-1-related receptor
  "prot-a-575",   # C-type lectin domain family 2 member L
  "prot-a-3036",  # Tumor necrosis factor receptor superfamily member 11A
  "prot-a-1043",  # Protein FAM19A5
  "prot-a-2314",  # Bifunctional polynucleotide phosphatase/kinase
  "prot-a-1150",  # Furin
  "prot-a-3242",  # DNA repair protein XRCC4
  "prot-a-2398",  # PH and SEC7 domain-containing protein 1
  "prot-a-2875",  # Stromal interaction molecule 1
  "prot-a-2930",  # Tubulin-specific chaperone A
  "prot-a-1378",  # Heparan-sulfate 6-O-sulfotransferase 1
  "prot-a-2824",  # Kunitz-type protease inhibitor 2
  "prot-a-57",    # Protein argonaute-1
  "prot-a-2556",  # Ribonuclease 4
  "prot-a-2618",  # Protein S100-A4
  "prot-a-735",   # CUB and zona pellucida-like domain-containing protein 1
  "prot-a-762"    # Dystroglycan
)

cat("【第 1 步】提取 55 个候选蛋白的工具变量...\n")

# 读取第一阶段保存的 IV 数据
if (file.exists("results/01_extraction_progress.rds")) {
  all_exp <- readRDS("results/01_extraction_progress.rds")
  cat("  从保存的进度文件中读取 IV 数据\n")
} else {
  stop("找不到第一阶段的数据文件，请确认 results/ 文件夹存在")
}

# 只保留 55 个候选蛋白的 IV
candidate_exp <- all_exp %>% filter(id.exposure %in% candidate_ids)
cat("  候选蛋白的 IV 数量:", nrow(candidate_exp), "\n")
cat("  覆盖的候选蛋白数:", length(unique(candidate_exp$id.exposure)), "\n\n")

# ============================================================
# 第 2 步：从心衰 GWAS 提取结局数据
# ============================================================
cat("【第 2 步】从心衰 GWAS 提取结局数据...\n")
cat("  结局:", REPLICATION_OUTCOME, "\n")

outcome_dat <- tryCatch({
  extract_outcome_data(
    snps     = unique(candidate_exp$SNP),
    outcomes = REPLICATION_OUTCOME
  )
}, error = function(e) {
  cat("  一次性提取出错，改为分批...\n")
  snp_list <- unique(candidate_exp$SNP)
  batch_size <- 100
  batches <- split(snp_list, ceiling(seq_along(snp_list) / batch_size))
  
  outcome_parts <- list()
  for (i in seq_along(batches)) {
    cat(sprintf("    批次 %d/%d\n", i, length(batches)))
    part <- tryCatch({
      extract_outcome_data(snps = batches[[i]], outcomes = REPLICATION_OUTCOME)
    }, error = function(e2) {
      Sys.sleep(5)
      return(NULL)
    })
    if (!is.null(part)) outcome_parts[[i]] <- part
    Sys.sleep(2)
  }
  do.call(rbind, outcome_parts)
})

if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
  stop("无法从心衰 GWAS 提取数据，请检查网络")
}

cat("  匹配到的 SNP 数:", nrow(outcome_dat), "\n\n")

# ============================================================
# 第 3 步：数据对齐 + MR 分析
# ============================================================
cat("【第 3 步】数据对齐 + MR 分析...\n")

dat <- harmonise_data(
  exposure_dat = candidate_exp,
  outcome_dat  = outcome_dat
)

mr_repl <- mr(dat, method_list = c("mr_ivw", "mr_wald_ratio"))

cat("  MR 分析完成，结果行数:", nrow(mr_repl), "\n\n")

# ============================================================
# 第 4 步：多重校正 + 与 Discovery 结果合并
# ============================================================
cat("【第 4 步】多重校正 + 结果合并...\n")

# 读取 Discovery 阶段结果
discovery <- read.csv("results/04_mr_results_corrected.csv") %>%
  filter(id.exposure %in% candidate_ids) %>%
  select(id.exposure, protein_name, 
         disc_method = method, disc_nsnp = nsnp,
         disc_beta = b, disc_se = se, disc_pval = pval)

# 整理 Replication 结果
replication <- mr_repl %>%
  filter(!is.na(pval)) %>%
  mutate(
    protein_name = sub(" \\|\\|.*", "", exposure),
    repl_bonf = pmin(pval * N_CANDIDATES, 1),
    repl_fdr  = p.adjust(pval, method = "BH")
  ) %>%
  select(id.exposure, protein_name,
         repl_method = method, repl_nsnp = nsnp,
         repl_beta = b, repl_se = se, repl_pval = pval,
         repl_bonf, repl_fdr)

# 合并两阶段结果
combined <- inner_join(discovery, replication, by = c("id.exposure", "protein_name"))

# 标记验证结果
combined <- combined %>%
  mutate(
    replicated_bonf = repl_bonf < 0.05,
    replicated_fdr  = repl_fdr < 0.05,
    replicated_nominal = repl_pval < 0.05,
    # 方向一致性检查：两阶段 beta 同号
    direction_consistent = sign(disc_beta) == sign(repl_beta)
  )

# 排序
combined <- combined %>% arrange(repl_pval)

cat("\n============================================\n")
cat(" 两阶段 MR 验证结果汇总\n")
cat("============================================\n")
cat("  候选蛋白总数:", nrow(combined), "\n")
cat("  Replication Bonferroni 显著 (P <", formatC(REPL_BONF_THRESH, format="e", digits=2), "):", 
    sum(combined$replicated_bonf, na.rm=TRUE), "个\n")
cat("  Replication FDR < 0.05:", sum(combined$replicated_fdr, na.rm=TRUE), "个\n")
cat("  Replication nominal P < 0.05:", sum(combined$replicated_nominal, na.rm=TRUE), "个\n")
cat("  方向一致 (beta 同号):", sum(combined$direction_consistent, na.rm=TRUE), "/", 
    nrow(combined), "\n")
cat("============================================\n\n")

# 打印 Top 结果
cat("=== 按 Replication P 值排序的完整结果 ===\n\n")
for (i in 1:nrow(combined)) {
  r <- combined[i, ]
  flag <- ""
  if (!is.na(r$replicated_bonf) && r$replicated_bonf) flag <- " ★★★ BONFERRONI"
  else if (!is.na(r$replicated_fdr) && r$replicated_fdr) flag <- " ★★ FDR"
  else if (!is.na(r$replicated_nominal) && r$replicated_nominal) flag <- " ★ NOMINAL"
  
  dir_mark <- ifelse(!is.na(r$direction_consistent) && r$direction_consistent, "✓同向", "✗反向")
  
  cat(sprintf("%2d. %s%s\n", i, r$protein_name, flag))
  cat(sprintf("    Discovery:    beta=%.3f, P=%.2e\n", r$disc_beta, r$disc_pval))
  cat(sprintf("    Replication:  beta=%.3f, P=%.2e (Bonf=%.4f) [%s]\n\n",
              r$repl_beta, r$repl_pval, r$repl_bonf, dir_mark))
}

# ============================================================
# 第 5 步：保存结果 + 画图
# ============================================================
write.csv(combined, "results/08_two_stage_combined.csv", row.names = FALSE)

# 画两阶段对比图 (Miami plot style)
plot_data <- combined %>%
  mutate(
    neg_log10_disc = -log10(disc_pval),
    neg_log10_repl = -log10(repl_pval),
    status = case_when(
      replicated_bonf ~ "Bonferroni 验证",
      replicated_fdr  ~ "FDR 验证", 
      replicated_nominal ~ "Nominal 验证",
      TRUE ~ "未验证"
    ),
    label = ifelse(repl_pval < 0.05, substr(protein_name, 1, 25), NA)
  )

p <- ggplot(plot_data, aes(x = disc_beta, y = repl_beta, color = status)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "blue", alpha = 0.5) +
  scale_color_manual(values = c(
    "Bonferroni 验证" = "#E31A1C",
    "FDR 验证"        = "#FF7F00",
    "Nominal 验证"    = "#33A02C",
    "未验证"          = "grey60"
  ), name = "验证状态") +
  labs(
    title = "Two-Stage MR: Discovery (Amyloidosis) vs Replication (Heart Failure)",
    subtitle = paste0("N candidates = ", N_CANDIDATES, 
                      " | Replication: ", REPLICATION_OUTCOME,
                      " | Bonf threshold: P < ", formatC(REPL_BONF_THRESH, format="e", digits=1)),
    x = "Discovery beta (Amyloidosis)",
    y = "Replication beta (Heart Failure)",
    caption = "Diagonal line = perfect agreement\nGreen/Orange/Red = validated proteins"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey40", size = 11),
    legend.position = "bottom"
  )

ggsave("figures/two_stage_scatter.pdf", p, width = 10, height = 8, dpi = 300)
ggsave("figures/two_stage_scatter.png", p, width = 10, height = 8, dpi = 300)

# 火山图 (Replication 阶段)
p2 <- ggplot(plot_data, aes(x = repl_beta, y = neg_log10_repl, color = status)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, na.rm = TRUE) +
  geom_hline(yintercept = -log10(REPL_BONF_THRESH), linetype = "dashed", color = "#E31A1C") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#33A02C") +
  scale_color_manual(values = c(
    "Bonferroni 验证" = "#E31A1C",
    "FDR 验证"        = "#FF7F00",
    "Nominal 验证"    = "#33A02C",
    "未验证"          = "grey60"
  ), name = "验证状态") +
  labs(
    title = "Replication MR: 55 Candidate Proteins → Heart Failure",
    subtitle = paste0("Outcome: ", REPLICATION_OUTCOME, " (HERMES, 47,309 HF cases)"),
    x = "MR Effect Size (beta)",
    y = expression(-log[10](P)),
    caption = paste0("Red dashed = Bonferroni (", formatC(REPL_BONF_THRESH, format="e", digits=1), 
                     ")\nGreen dashed = nominal P = 0.05")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 15),
    legend.position = "bottom"
  )

ggsave("figures/replication_volcano.pdf", p2, width = 10, height = 8, dpi = 300)
ggsave("figures/replication_volcano.png", p2, width = 10, height = 8, dpi = 300)

cat("\n✅ 所有结果已保存！\n")
cat("   results/08_two_stage_combined.csv  - 两阶段合并结果\n")
cat("   figures/two_stage_scatter.pdf/png  - 两阶段效应值对比图\n")
cat("   figures/replication_volcano.pdf/png - Replication 火山图\n")
cat("============================================\n")
