# ============================================================
#  全蛋白质组学孟德尔随机化 (PWMR) Pipeline
#  目标：识别与淀粉样变有因果关系的循环蛋白质
#  结局 (Outcome): finn-b-E4_AMYLOIDOSIS (FinnGen 淀粉样变)
#  暴露 (Exposure): prot-a 批次 (血液循环蛋白 pQTL)
#
#  使用方法：
#  1. 先运行 00_install_packages.R 安装所有包
#  2. 在 RStudio 里打开这个文件
#  3. 全选 → 点 "Run" 或按 Cmd+Shift+Enter
#
#  预计运行时间：2-6 小时（取决于网速，因为要从网上逐个下载数据）
#  如果中途断了，不用怕，脚本有自动保存进度的功能
# ============================================================

# ---- 加载所有需要的包 ----
library(TwoSampleMR)
library(ieugwasr)
library(coloc)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

# ---- 全局设置 ----
OUTCOME_ID  <- "finn-b-E4_AMYLOIDOSIS"   # 结局：淀粉样变
PVAL_THRESH <- 5e-8     # P 值阈值：只要最显著的 SNP
CLUMP_R2    <- 0.001    # LD clumping：非常严格，确保 SNP 之间互相独立
CLUMP_KB    <- 10000    # clumping 窗口：10Mb
BATCH_ID    <- "prot-a" # 蛋白质组学数据批次

# ---- 创建输出文件夹 ----
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

cat("============================================\n")
cat(" PWMR Pipeline 启动\n")
cat(" 结局 ID:", OUTCOME_ID, "\n")
cat(" 暴露批次:", BATCH_ID, "\n")
cat("============================================\n\n")

# ============================================================
# 第 1 步：获取所有 prot-a 蛋白的列表
# ============================================================
cat("【第 1 步】正在从 OpenGWAS 获取蛋白列表...\n")

ao <- available_outcomes()
prot_ids <- ao$id[grepl(paste0("^", BATCH_ID, "-"), ao$id)]

cat("  找到", length(prot_ids), "个蛋白质\n\n")

# ============================================================
# 第 2 步：逐个提取每个蛋白的工具变量 (IVs)
#          这是最耗时的步骤，有防断线机制
# ============================================================
cat("【第 2 步】开始逐个提取蛋白的工具变量 (IVs)...\n")
cat("  这一步需要较长时间，请耐心等待。\n")
cat("  如果某个蛋白报错，会自动跳过，不会崩溃。\n\n")

# 检查是否有之前跑了一半的进度文件
progress_file <- "results/01_extraction_progress.rds"
if (file.exists(progress_file)) {
  cat("  发现之前的进度文件，从断点继续...\n")
  all_exposures <- readRDS(progress_file)
  done_ids <- unique(all_exposures$id.exposure)
  prot_ids_todo <- setdiff(prot_ids, done_ids)
  cat("  已完成:", length(done_ids), "个, 剩余:", length(prot_ids_todo), "个\n")
} else {
  all_exposures <- data.frame()
  prot_ids_todo <- prot_ids
}

# 设置计数器
counter <- 0
total <- length(prot_ids_todo)
fail_count <- 0

for (pid in prot_ids_todo) {
  counter <- counter + 1

  # 每处理 10 个蛋白打印一次进度
  if (counter %% 10 == 0 || counter == 1) {
    cat(sprintf("  进度: %d/%d (失败跳过: %d)\n", counter, total, fail_count))
  }

  # 核心：tryCatch 防断线
  result <- tryCatch({
    # 提取该蛋白的显著 SNP 并做 clumping
    exp_dat <- extract_instruments(
      outcomes = pid,
      p1       = PVAL_THRESH,
      clump    = TRUE,
      r2       = CLUMP_R2,
      kb       = CLUMP_KB
    )

    # 如果这个蛋白没有显著 SNP，跳过
    if (is.null(exp_dat) || nrow(exp_dat) == 0) {
      return(NULL)
    }

    exp_dat
  },
  error = function(e) {
    fail_count <<- fail_count + 1
    # 网络错误时等 5 秒再继续
    Sys.sleep(5)
    return(NULL)
  })

  # 把成功的结果存起来
  if (!is.null(result)) {
    all_exposures <- rbind(all_exposures, result)
  }

  # 每处理 50 个蛋白自动保存一次进度（防止全部白跑）
  if (counter %% 50 == 0) {
    saveRDS(all_exposures, progress_file)
    cat("  [自动保存进度]\n")
  }

  # 每次请求之间暂停 1 秒，避免被服务器封 IP
  Sys.sleep(1)
}

# 最终保存
saveRDS(all_exposures, progress_file)

n_proteins_with_iv <- length(unique(all_exposures$id.exposure))
n_total_snps <- nrow(all_exposures)
cat("\n  ✅ IV 提取完成！\n")
cat("  有工具变量的蛋白数:", n_proteins_with_iv, "\n")
cat("  总 SNP 数:", n_total_snps, "\n\n")

# ============================================================
# 第 3 步：提取结局数据 + 数据对齐 + MR 分析
# ============================================================
cat("【第 3 步】提取结局数据并执行 MR 分析...\n")

# 提取结局中对应的 SNP
outcome_dat <- tryCatch({
  extract_outcome_data(
    snps     = unique(all_exposures$SNP),
    outcomes = OUTCOME_ID
  )
}, error = function(e) {
  cat("  ⚠️ 提取结局数据出错:", e$message, "\n")
  cat("  尝试分批提取...\n")

  # 如果一次提取太多 SNP 报错，就分批来
  snp_list <- unique(all_exposures$SNP)
  batch_size <- 200
  batches <- split(snp_list, ceiling(seq_along(snp_list) / batch_size))

  outcome_parts <- list()
  for (i in seq_along(batches)) {
    cat(sprintf("    批次 %d/%d\n", i, length(batches)))
    part <- tryCatch({
      extract_outcome_data(snps = batches[[i]], outcomes = OUTCOME_ID)
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
  stop("❌ 无法从结局数据中提取到任何 SNP，请检查网络连接后重试。")
}

cat("  结局中匹配到的 SNP 数:", nrow(outcome_dat), "\n")

# 数据对齐（Harmonise）
cat("  正在对齐暴露和结局数据...\n")
dat <- harmonise_data(
  exposure_dat = all_exposures,
  outcome_dat  = outcome_dat
)

cat("  对齐后的数据行数:", nrow(dat), "\n")

# 执行 MR 分析（批量，对每个蛋白分别做）
cat("  正在执行 MR 分析...\n")
mr_results <- mr(dat, method_list = c("mr_ivw", "mr_wald_ratio"))

cat("  ✅ MR 分析完成！共", nrow(mr_results), "条结果\n\n")

# 保存中间结果
saveRDS(dat, "results/02_harmonised_data.rds")
saveRDS(mr_results, "results/03_mr_results_raw.rds")

# ============================================================
# 第 4 步：多重检验校正（Bonferroni + FDR）
# ============================================================
cat("【第 4 步】多重检验校正...\n")

# 只保留有 P 值的结果
mr_results <- mr_results %>% filter(!is.na(pval))

# 计算校正后的 P 值
n_tests <- length(unique(mr_results$id.exposure))  # 检验次数 = 蛋白数
cat("  检验次数 (蛋白数):", n_tests, "\n")

mr_results <- mr_results %>%
  mutate(
    # Bonferroni 校正：P 值乘以检验次数（最严格）
    p_bonferroni = pmin(pval * n_tests, 1),
    # FDR 校正：允许一定比例的假阳性（稍微宽松）
    p_fdr = p.adjust(pval, method = "BH"),
    # 提取蛋白名称（从 exposure 列中提取）
    protein_name = sub(" \\|\\|.*", "", exposure)
  )

# 筛选阳性蛋白
bonf_sig <- mr_results %>% filter(p_bonferroni < 0.05)
fdr_sig  <- mr_results %>% filter(p_fdr < 0.05)

cat("  Bonferroni 校正后显著的蛋白:", nrow(bonf_sig), "个\n")
cat("  FDR 校正后显著的蛋白:", nrow(fdr_sig), "个\n")

# 保存结果
write.csv(mr_results, "results/04_mr_results_corrected.csv", row.names = FALSE)
write.csv(bonf_sig, "results/05_bonferroni_significant.csv", row.names = FALSE)
write.csv(fdr_sig, "results/06_fdr_significant.csv", row.names = FALSE)

cat("  ✅ 结果已保存到 results/ 文件夹\n\n")

# ============================================================
# 第 5 步：共定位分析 (Colocalization)
#          只对 FDR 显著的蛋白做（省时间）
# ============================================================
cat("【第 5 步】共定位分析 (Colocalization)...\n")

if (nrow(fdr_sig) == 0) {
  cat("  ⚠️ 没有 FDR 显著的蛋白，跳过共定位分析。\n")
  cat("  这可能是因为样本量不够大（704 例），不代表你的分析有错。\n")
  cat("  你可以考虑：\n")
  cat("    1. 查看 FDR < 0.2 的蛋白作为提示性结果\n")
  cat("    2. 换一个更大的结局队列重新跑\n\n")

  # 即使没有 Bonferroni 显著的，也看看 Top 10
  top10 <- mr_results %>%
    arrange(pval) %>%
    head(10)
  cat("  Top 10 最小 P 值的蛋白：\n")
  print(top10 %>% select(protein_name, method, nsnp, b, pval, p_fdr, p_bonferroni))

  coloc_results <- data.frame()

} else {
  # 对每个显著蛋白做共定位
  sig_exposure_ids <- unique(fdr_sig$id.exposure)
  coloc_results <- data.frame()

  for (exp_id in sig_exposure_ids) {
    prot_name <- fdr_sig$protein_name[fdr_sig$id.exposure == exp_id][1]
    cat("  正在分析:", prot_name, "\n")

    result <- tryCatch({
      # 获取该蛋白的所有关联 SNP（不限 P 值）
      exp_full <- associations(id = exp_id, proxies = FALSE)
      out_full <- associations(id = OUTCOME_ID, proxies = FALSE)

      # 找两个数据集共有的 SNP
      common_snps <- intersect(exp_full$rsid, out_full$rsid)

      if (length(common_snps) < 100) {
        cat("    共有 SNP 太少 (", length(common_snps), ")，跳过\n")
        return(NULL)
      }

      exp_sub <- exp_full %>% filter(rsid %in% common_snps) %>%
        arrange(rsid) %>% distinct(rsid, .keep_all = TRUE)
      out_sub <- out_full %>% filter(rsid %in% common_snps) %>%
        arrange(rsid) %>% distinct(rsid, .keep_all = TRUE)

      # 确保顺序一致
      exp_sub <- exp_sub[match(sort(common_snps), exp_sub$rsid), ]
      out_sub <- out_sub[match(sort(common_snps), out_sub$rsid), ]

      # 运行 coloc
      coloc_res <- coloc.abf(
        dataset1 = list(
          beta    = exp_sub$beta,
          varbeta = exp_sub$se^2,
          type    = "quant",
          N       = max(exp_sub$n, na.rm = TRUE),
          snp     = exp_sub$rsid
        ),
        dataset2 = list(
          beta    = out_sub$beta,
          varbeta = out_sub$se^2,
          type    = "cc",
          N       = max(out_sub$n, na.rm = TRUE),
          s       = 704 / max(out_sub$n, na.rm = TRUE),
          snp     = out_sub$rsid
        )
      )

      # 提取关键指标 PP.H4 = 共享同一因果变异的概率
      data.frame(
        protein      = prot_name,
        exposure_id  = exp_id,
        PP.H0        = coloc_res$summary["PP.H0.abf"],
        PP.H1        = coloc_res$summary["PP.H1.abf"],
        PP.H2        = coloc_res$summary["PP.H2.abf"],
        PP.H3        = coloc_res$summary["PP.H3.abf"],
        PP.H4        = coloc_res$summary["PP.H4.abf"],
        n_snps       = length(common_snps)
      )
    },
    error = function(e) {
      cat("    ⚠️ 出错跳过:", e$message, "\n")
      return(NULL)
    })

    if (!is.null(result)) {
      coloc_results <- rbind(coloc_results, result)
    }

    Sys.sleep(3)  # 避免 API 过载
  }

  if (nrow(coloc_results) > 0) {
    cat("\n  共定位结果 (PP.H4 > 0.75 表示强证据):\n")
    print(coloc_results %>% select(protein, PP.H4, n_snps))
    write.csv(coloc_results, "results/07_colocalization.csv", row.names = FALSE)
  }
}

cat("  ✅ 共定位分析完成！\n\n")

# ============================================================
# 第 6 步：画火山图 (Volcano Plot)
# ============================================================
cat("【第 6 步】绘制火山图...\n")

# 准备画图数据
plot_data <- mr_results %>%
  filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
  mutate(
    neg_log10_p = -log10(pval),
    significant = case_when(
      p_bonferroni < 0.05 ~ "Bonferroni 显著",
      p_fdr < 0.05        ~ "FDR 显著",
      TRUE                ~ "不显著"
    ),
    # 如果蛋白名太长，截取前 20 个字符
    label = ifelse(p_fdr < 0.05, substr(protein_name, 1, 25), NA)
  )

# Bonferroni 和 FDR 的阈值线
bonf_line <- -log10(0.05 / n_tests)
fdr_candidates <- plot_data$pval[plot_data$p_fdr < 0.05]
fdr_line <- if (length(fdr_candidates) > 0) -log10(max(fdr_candidates)) else NA

# 画图
p <- ggplot(plot_data, aes(x = b, y = neg_log10_p, color = significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Bonferroni 显著" = "#E31A1C",
      "FDR 显著"        = "#FF7F00",
      "不显著"          = "grey60"
    ),
    name = "显著性"
  ) +
  geom_hline(yintercept = bonf_line, linetype = "dashed", color = "#E31A1C",
             linewidth = 0.5) +
  {if (!is.na(fdr_line)) geom_hline(yintercept = fdr_line, linetype = "dashed",
                                     color = "#FF7F00", linewidth = 0.5)} +
  geom_text_repel(
    aes(label = label),
    size = 3, max.overlaps = 20, na.rm = TRUE,
    box.padding = 0.5, point.padding = 0.3
  ) +
  labs(
    title    = "Proteome-wide MR: Circulating Proteins → Amyloidosis",
    subtitle = paste0("Outcome: ", OUTCOME_ID, " | Exposure batch: ", BATCH_ID),
    x        = "MR Effect Size (beta)",
    y        = expression(-log[10](P)),
    caption  = paste0(
      "Red dashed = Bonferroni (", formatC(0.05/n_tests, format="e", digits=1), ")\n",
      "Orange dashed = FDR 0.05 threshold\n",
      "N proteins tested = ", n_tests
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey40"),
    legend.position = "bottom"
  )

# 保存图片（高清）
ggsave("figures/volcano_plot.pdf", p, width = 12, height = 8, dpi = 300)
ggsave("figures/volcano_plot.png", p, width = 12, height = 8, dpi = 300)

cat("  ✅ 火山图已保存到 figures/ 文件夹\n\n")

# ============================================================
# 最终汇总报告
# ============================================================
cat("============================================\n")
cat(" 🎉 Pipeline 运行完毕！\n")
cat("============================================\n")
cat(" 结局:", OUTCOME_ID, "\n")
cat(" 暴露批次:", BATCH_ID, "\n")
cat(" 测试蛋白数:", n_tests, "\n")
cat(" Bonferroni 显著:", nrow(bonf_sig), "个\n")
cat(" FDR 显著:", nrow(fdr_sig), "个\n")
if (nrow(coloc_results) > 0) {
  coloc_strong <- sum(coloc_results$PP.H4 > 0.75, na.rm = TRUE)
  cat(" 共定位 PP.H4 > 0.75:", coloc_strong, "个\n")
}
cat("\n 输出文件清单:\n")
cat("   results/04_mr_results_corrected.csv  - MR 全部结果\n")
cat("   results/05_bonferroni_significant.csv - Bonferroni 显著蛋白\n")
cat("   results/06_fdr_significant.csv        - FDR 显著蛋白\n")
if (file.exists("results/07_colocalization.csv")) {
  cat("   results/07_colocalization.csv         - 共定位结果\n")
}
cat("   figures/volcano_plot.pdf              - 火山图 (PDF)\n")
cat("   figures/volcano_plot.png              - 火山图 (PNG)\n")
cat("============================================\n")
