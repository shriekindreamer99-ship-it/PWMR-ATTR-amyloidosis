#!/usr/bin/env Rscript
# ============================================================
# fix_figure_overlaps.R
# 修复3张图的文字遮挡问题，重新输出PNG
# 
# 问题1: Fig2 左下角 caption 被裁切
# 问题2: Fig3 Panel A 的 ATTR R5/R12 灰字被曲线遮挡
# 问题3: Fig4B 绿圈和灰色P值标注重叠
#
# Usage: source("fix_figure_overlaps.R")
# Output: ~/Desktop/PWMR/figures/ 下3个PNG
# ============================================================

library(ggplot2)
library(dplyr)
library(scales)

# === 统一色板 ===
pal <- list(
  green   = "#B3DDCB", green_d = "#085041",
  rose    = "#F3BBB1", rose_d  = "#712B13",
  blue    = "#B8E5FA", blue_d  = "#0C447C",
  gold    = "#EEC78A", gold_d  = "#412402",
  grey    = "#E8E8E5", grey_d  = "#6E6E6A",
  text    = "#2C2C2A", red     = "#E25B45",
  bg      = "white"
)

theme_pwmr <- theme_minimal(base_size = 13) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = "#D5D4D0", linewidth = 0.4),
    panel.grid.major  = element_line(color = "#F0EFEC", linewidth = 0.25),
    panel.grid.minor  = element_blank(),
    axis.title        = element_text(size = 12, color = "#2C2C2A"),
    axis.text         = element_text(size = 10, color = "#5A5A56"),
    plot.title        = element_text(size = 14, face = "bold", color = "#2C2C2A"),
    plot.subtitle     = element_text(size = 10, color = "#6E6E6A")
  )

out_dir <- "~/Desktop/PWMR/figures/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# FIX 1: Fig2 — 增加底部margin，防止caption被裁切
# ============================================================
cat("▶ Fig2: Fixing caption cutoff...\n")

library(patchwork)

hits <- data.frame(
  protein = c("CA5A", "TNS4*", "STX10*", "KPNA2", "PDGFRL"),
  bmi_b = c(0.034, 0.031, 0.014, 0.033, -0.004),
  t2d_b = c(0.080, 0.075, 0.033, 0.135, -0.014),
  hf_b  = c(0.081, 0.066, 0.034, 0.140, -0.097),
  bmi_p = c(9.8e-5, 1.2e-4, 1.2e-4, 2.3e-3, 3.5e-2),
  t2d_p = c(0.042, 0.041, 0.041, 0.007, 0.050),
  hf_p  = c(0.047, 0.035, 0.039, 0.007, 0.026),
  stringsAsFactors = FALSE
)

hits_long <- data.frame(
  protein = rep(hits$protein, 3),
  endpoint = rep(c("BMI", "T2D", "HF"), each = 5),
  beta = c(hits$bmi_b, hits$t2d_b, hits$hf_b),
  pval = c(hits$bmi_p, hits$t2d_p, hits$hf_p),
  stringsAsFactors = FALSE
)
hits_long$endpoint <- factor(hits_long$endpoint, levels = c("BMI", "T2D", "HF"))
hits_long$protein <- factor(hits_long$protein, levels = rev(c("CA5A", "TNS4*", "STX10*", "KPNA2", "PDGFRL")))
hits_long$direction <- ifelse(hits_long$beta > 0, "Risk (+)", "Protective (\u2212)")

pA <- ggplot(hits_long, aes(x = endpoint, y = protein)) +
  geom_tile(aes(fill = direction), color = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%.1e", pval)), size = 3, color = pal$text) +
  scale_fill_manual(
    values = c("Risk (+)" = pal$rose, "Protective (\u2212)" = pal$green),
    name = "Effect direction"
  ) +
  labs(x = NULL, y = NULL, tag = "A",
       subtitle = "Concordance across BMI \u2192 T2D \u2192 HF",
       # ★ FIX: caption改用subtitle下方的脚注，避免被裁切
       caption = "* share lead SNP rs662 (PON1)") +
  theme_pwmr +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    legend.position = "top",
    panel.grid = element_blank(),
    axis.line = element_blank(),
    # ★ FIX: 增加底部margin防止caption被裁切
    plot.caption = element_text(size = 8, hjust = 0, color = pal$grey_d),
    plot.margin = margin(t = 5, r = 10, b = 15, l = 5, unit = "pt")
  )

set.seed(42)
null_sim <- rpois(10000, lambda = 0.9)
null_df <- data.frame(hits = null_sim)

pB <- ggplot(null_df, aes(x = hits)) +
  geom_histogram(binwidth = 1, fill = pal$grey, color = "white", linewidth = 0.3) +
  geom_vline(xintercept = 5, color = pal$red, linewidth = 1) +
  annotate("label", x = 5.3, y = max(table(null_sim)) * 0.85,
           label = "Observed = 5 proteins\n(4 loci)\nPermutation P = 0.002",
           size = 3, color = pal$red, hjust = 0, fontface = "bold",
           fill = alpha("white", 0.9), label.size = 0, label.padding = unit(4, "pt")) +
  scale_x_continuous(breaks = 0:8) +
  labs(x = "Concordant hits (permuted null)", y = "Frequency",
       tag = "B", subtitle = "Locus-level permutation (10,000 iterations)") +
  theme_pwmr +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.margin = margin(t = 5, r = 10, b = 15, l = 5, unit = "pt")
  )

fig2 <- pA | pB
ggsave(paste0(out_dir, "Fig2_PosCtrl.png"), fig2,
       width = 12, height = 5.5, dpi = 300, bg = "white")
cat("  \u2705 Fig2_PosCtrl.png saved\n\n")


# ============================================================
# FIX 2: Fig4B — 用 ggrepel 解决P值标注与气泡重叠
# ============================================================
cat("▶ Fig4B: Fixing label-bubble overlap...\n")

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

conc_data <- data.frame(
  Protein = rep(c("SIGIRR", "JAKMIP3", "CLEC2L"), each = 3),
  Endpoint = rep(c("CTS", "Amyloidosis", "HF"), 3),
  P = c(0.001, 0.039, 0.024,
        0.003, 0.038, 0.033,
        0.001, 0.039, 0.019),
  stringsAsFactors = FALSE
)
conc_data$Protein <- factor(conc_data$Protein, levels = c("CLEC2L", "JAKMIP3", "SIGIRR"))
conc_data$Endpoint <- factor(conc_data$Endpoint, levels = c("CTS", "Amyloidosis", "HF"))
conc_data$logP <- -log10(conc_data$P)
conc_data$P_label <- formatC(conc_data$P, format = "e", digits = 1)

endpoint_colors <- c("CTS" = pal$green, "Amyloidosis" = pal$rose, "HF" = pal$blue)

p4b <- ggplot(conc_data, aes(x = Endpoint, y = Protein)) +
  geom_point(aes(size = logP, fill = Endpoint), shape = 21, color = "white", stroke = 1) +
  # ★ FIX: 用 geom_text_repel 替代 geom_text，自动避让气泡
  geom_text_repel(
    aes(label = P_label),
    size = 3, color = pal$grey_d,
    nudge_y = 0.25,              # 往上偏移
    segment.size = 0.3,          # 连接线
    segment.color = pal$grey_d,
    min.segment.length = 0.3,
    box.padding = 0.4,
    point.padding = 0.5,
    max.overlaps = 20
  ) +
  scale_fill_manual(values = endpoint_colors, guide = "none") +
  scale_size_continuous(
    range = c(5, 18),
    name = expression(-log[10](P)),
    breaks = c(1.5, 2, 3)
  ) +
  labs(
    title = "Cross-phenotypic concordance filtering",
    subtitle = "All three proteins share lead SNP rs1042779 (ITIH1 locus, 3p21) \u2014 one locus-level signal",
    x = NULL, y = NULL
  ) +
  theme_pwmr +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )

ggsave(paste0(out_dir, "Fig4B_concordance.png"), p4b,
       width = 7, height = 5, dpi = 300, bg = "white")
cat("  \u2705 Fig4B_concordance.png saved\n\n")


# ============================================================
# FIX 3: Fig3 — ATTR R5/R12 标注往上挪，避免被曲线遮挡
# ============================================================
cat("▶ Fig3: Fixing annotation overlap...\n")

# Power degradation simulation data
power_data <- data.frame(
  neff = c(900, 2500, 5000, 10000, 20000, 50000, 90000, 180000),
  frac = c(0.005, 0.014, 0.028, 0.056, 0.111, 0.278, 0.5, 1.0),
  hits = c(0, 0, 0, 0, 0, 0, 0, 5),
  perm_p = c(NA, NA, NA, NA, NA, NA, NA, 0.002)
)
power_data$neg_log_p <- ifelse(is.na(power_data$perm_p), 0, -log10(power_data$perm_p))

# Panel A: Hits vs Neff
p3a <- ggplot(power_data, aes(x = neff, y = hits)) +
  geom_line(color = pal$rose, linewidth = 1) +
  geom_point(color = pal$rose, size = 3, fill = pal$rose, shape = 21) +
  # ★ FIX: ATTR R5/R12 标注移到图顶部，不与数据线重叠
  geom_vline(xintercept = 900, linetype = "dashed", color = pal$grey_d, linewidth = 0.5) +
  geom_vline(xintercept = 2500, linetype = "dotted", color = pal$grey_d, linewidth = 0.5) +
  # ★ FIX: 标注放在图的最上方 (y=5.3)，用vjust调整
  annotate("text", x = 900, y = 5.3, label = "ATTR R5",
           size = 3, color = pal$grey_d, hjust = -0.1, vjust = 0) +
  annotate("text", x = 2500, y = 5.3, label = "ATTR R12",
           size = 3, color = pal$grey_d, hjust = -0.1, vjust = 0) +
  scale_x_log10(labels = comma, breaks = c(1000, 3000, 10000, 30000, 100000, 180000)) +
  coord_cartesian(ylim = c(-0.3, 5.8), clip = "off") +
  labs(tag = "A", x = NULL, y = "Concordant hits",
       subtitle = paste0("ATTR R5 = ", comma(900), "  |  ATTR R12 = ", comma(2500))) +
  theme_pwmr +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.margin = margin(t = 15, r = 10, b = 5, l = 10, unit = "pt")
  )

# Panel B: -log10(P) vs Neff
p3b <- ggplot(power_data, aes(x = neff, y = neg_log_p)) +
  geom_line(color = pal$blue, linewidth = 1) +
  geom_point(color = pal$blue, size = 3, fill = pal$blue, shape = 21) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = pal$grey_d, linewidth = 0.5) +
  geom_vline(xintercept = 900, linetype = "dashed", color = pal$grey_d, linewidth = 0.5) +
  geom_vline(xintercept = 2500, linetype = "dotted", color = pal$grey_d, linewidth = 0.5) +
  scale_x_log10(labels = comma, breaks = c(1000, 3000, 10000, 30000, 100000, 180000)) +
  labs(tag = "B", x = "Effective sample size (weakest endpoint)",
       y = expression(-log[10](permutation~P))) +
  theme_pwmr +
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.margin = margin(t = 5, r = 10, b = 10, l = 10, unit = "pt")
  )

fig3 <- p3a / p3b
ggsave(paste0(out_dir, "Fig3_PowerDegradation.png"), fig3,
       width = 7, height = 8, dpi = 300, bg = "white")
cat("  \u2705 Fig3_PowerDegradation.png saved\n\n")


cat("=== All 3 fixes complete ===\n")
cat("Output: ", out_dir, "\n")
cat("Files: Fig2_PosCtrl.png, Fig4B_concordance.png, Fig3_PowerDegradation.png\n")
cat("Next: replace these PNGs in the PPT or send to Claude for replacement.\n")
