# ============================================================
# Supplementary Table S1: Full 2094-protein MR results
# PWMR ATTR project — export script v3
# 使用方法：
#   1. install.packages("openxlsx")   # 第一次需要安装
#   2. setwd("~/Desktop/PWMR")
#   3. source("export_supp_table_v3.R")
# 输出: results/Supplementary_Table_S1_Full_MR_Results.xlsx
# ============================================================

library(dplyr)
library(openxlsx)

# ---- 1. 文件路径（根据实际文件名）--------------------------

cts_file  <- "results/09_mr_cts_results.csv"
ca_file   <- "results/04_mr_results_corrected.csv"   # systemic amyloidosis
hf_file   <- "results/10_mr_hf_results.csv"
t2d_file  <- "results/19_mr_t2d_results.csv"
asth_file <- "results/20_mr_asthma_results.csv"

cat("使用文件:\n")
for (f in c(cts_file, ca_file, hf_file, t2d_file, asth_file)) {
  flag <- if (file.exists(f)) "OK" else "MISSING!"
  cat(" ", flag, f, "\n")
}

# ---- 2. 读取数据 --------------------------------------------

cts  <- read.csv(cts_file,  stringsAsFactors = FALSE)
ca   <- read.csv(ca_file,   stringsAsFactors = FALSE)
hf   <- read.csv(hf_file,   stringsAsFactors = FALSE)
t2d  <- read.csv(t2d_file,  stringsAsFactors = FALSE)
asth <- read.csv(asth_file, stringsAsFactors = FALSE)

# 确保id.exposure列是字符型
for (df_name in c("cts","ca","hf","t2d","asth")) {
  df <- get(df_name)
  df$id.exposure <- as.character(df$id.exposure)
  assign(df_name, df)
}

cat("\n数据行数 — CTS:", nrow(cts), "| CA:", nrow(ca),
    "| HF:", nrow(hf), "| T2D:", nrow(t2d), "| Asthma:", nrow(asth), "\n\n")

# ---- 3. 清洗函数 --------------------------------------------

clean_mr <- function(df, suffix) {
  if (!"protein_name" %in% names(df)) df$protein_name <- df$exposure
  df %>%
    mutate(
      id.exposure  = as.character(id.exposure),
      protein_name = gsub(" \\|\\| id:.*", "", protein_name),
      b    = round(b,    4),
      se   = round(se,   4),
      pval = signif(pval, 3)
    ) %>%
    select(id.exposure, protein_name, nsnp, b, se, pval) %>%
    rename_with(~ paste0(., "_", suffix), c(nsnp, b, se, pval))
}

cts_c  <- clean_mr(cts,  "CTS")
ca_c   <- clean_mr(ca,   "CA")
hf_c   <- clean_mr(hf,   "HF")
t2d_c  <- clean_mr(t2d,  "T2D")
asth_c <- clean_mr(asth, "Asthma")

# ---- 4. 合并 ------------------------------------------------

merged <- cts_c %>%
  full_join(ca_c   %>% select(-protein_name), by = "id.exposure") %>%
  full_join(hf_c   %>% select(-protein_name), by = "id.exposure") %>%
  full_join(t2d_c  %>% select(-protein_name), by = "id.exposure") %>%
  full_join(asth_c %>% select(-protein_name), by = "id.exposure")

# ---- 5. 标注 ------------------------------------------------

merged <- merged %>%
  mutate(
    triple_primary = (!is.na(pval_CTS) & pval_CTS < 0.05) &
                     (!is.na(pval_CA)  & pval_CA  < 0.05) &
                     (!is.na(pval_HF)  & pval_HF  < 0.05),
    dir_concordant = ifelse(
      triple_primary,
      (sign(b_CTS) == sign(b_CA)) & (sign(b_CA) == sign(b_HF)),
      NA
    ),
    neg_ctrl_pass  = (is.na(pval_T2D)    | pval_T2D    > 0.05) &
                     (is.na(pval_Asthma) | pval_Asthma > 0.05),
    annotation = case_when(
      id.exposure == "prot-a-2722" ~ "ITIH1-locus — SIGIRR (primary finding)",
      id.exposure == "prot-a-1599" ~ "ITIH1-locus — JAKMIP3 (primary finding)",
      id.exposure == "prot-a-575"  ~ "ITIH1-locus — CLEC2L (primary finding)",
      id.exposure == "prot-a-225"  ~ "BAG3 — orthogonal HF signal (cis-pQTL)",
      triple_primary & isTRUE(dir_concordant) ~ "Triple-validated, direction concordant",
      triple_primary                           ~ "Triple-validated, direction discordant",
      !is.na(pval_CTS) & pval_CTS < 0.05      ~ "CTS nominal only",
      !is.na(pval_CA)  & pval_CA  < 0.05      ~ "CA nominal only",
      !is.na(pval_HF)  & pval_HF  < 0.05      ~ "HF nominal only",
      TRUE                                     ~ "Not significant"
    )
  ) %>%
  arrange(
    desc(id.exposure %in% c("prot-a-2722","prot-a-1599","prot-a-575")),
    desc(triple_primary & isTRUE(dir_concordant)),
    desc(triple_primary),
    replace(pval_CTS, is.na(pval_CTS), 1)
  ) %>%
  select(
    `Protein ID`                   = id.exposure,
    `Protein Name`                 = protein_name,
    `N SNPs (CTS)`                 = nsnp_CTS,
    `Beta (CTS)`                   = b_CTS,
    `SE (CTS)`                     = se_CTS,
    `P (CTS)`                      = pval_CTS,
    `Beta (CA)`                    = b_CA,
    `SE (CA)`                      = se_CA,
    `P (CA)`                       = pval_CA,
    `Beta (HF)`                    = b_HF,
    `SE (HF)`                      = se_HF,
    `P (HF)`                       = pval_HF,
    `Beta (T2D)`                   = b_T2D,
    `SE (T2D)`                     = se_T2D,
    `P (T2D)`                      = pval_T2D,
    `Beta (Asthma)`                = b_Asthma,
    `SE (Asthma)`                  = se_Asthma,
    `P (Asthma)`                   = pval_Asthma,
    `Triple validated (CTS+CA+HF)` = triple_primary,
    `Direction concordant`         = dir_concordant,
    `Neg ctrl PASS (T2D+Asthma)`   = neg_ctrl_pass,
    `Annotation`                   = annotation
  )

cat("总蛋白数:", nrow(merged), "\n")
cat("三角验证:", sum(merged$`Triple validated (CTS+CA+HF)`, na.rm=TRUE), "个\n\n")

# ---- 6. 构建Excel -----------------------------------------

wb <- createWorkbook()

hdr <- createStyle(fontName="Arial", fontSize=10, fontColour="#FFFFFF", fgFill="#2C5282",
                   halign="CENTER", valign="CENTER", textDecoration="bold",
                   wrapText=TRUE, border="Bottom", borderStyle="medium", borderColour="#FFFFFF")
body_c <- createStyle(fontName="Arial", fontSize=9, halign="CENTER")
body_l <- createStyle(fontName="Arial", fontSize=9, halign="LEFT")
hit_c  <- createStyle(fontName="Arial", fontSize=9, halign="CENTER", fgFill="#EBF8FF")
hit_l  <- createStyle(fontName="Arial", fontSize=9, halign="LEFT",   fgFill="#EBF8FF")
alt_c  <- createStyle(fontName="Arial", fontSize=9, halign="CENTER", fgFill="#F7FAFC")
alt_l  <- createStyle(fontName="Arial", fontSize=9, halign="LEFT",   fgFill="#F7FAFC")
sig_c  <- createStyle(fontName="Arial", fontSize=9, halign="CENTER",
                       fontColour="#1A6B3C", textDecoration="bold")
top_c  <- createStyle(fontName="Arial", fontSize=9, halign="CENTER", fgFill="#FFF3CD")
top_l  <- createStyle(fontName="Arial", fontSize=9, halign="LEFT",   fgFill="#FFF3CD")
title_s <- createStyle(fontName="Arial", fontSize=10, fontColour="#2C5282",
                        textDecoration="bold", wrapText=TRUE)

col_w   <- c(14,35,8,8,7,9,8,7,9,8,7,9,8,7,9,8,7,9,10,10,10,28)
p_cols  <- c(6,9,12,15,18)
top_ids <- c("prot-a-2722","prot-a-1599","prot-a-575","prot-a-225")

apply_row_styles <- function(sheet, df) {
  n <- nrow(df)
  for (i in seq_len(n)) {
    r      <- i + 2
    is_top <- df$`Protein ID`[i] %in% top_ids
    is_hit <- isTRUE(df$`Triple validated (CTS+CA+HF)`[i]) & !is_top
    is_alt <- (i %% 2 == 0) & !is_top & !is_hit
    fc <- if (is_top) top_c else if (is_hit) hit_c else if (is_alt) alt_c else body_c
    fl <- if (is_top) top_l else if (is_hit) hit_l else if (is_alt) alt_l else body_l
    addStyle(wb, sheet, fl, rows=r, cols=1:2,  gridExpand=TRUE)
    addStyle(wb, sheet, fc, rows=r, cols=3:22, gridExpand=TRUE)
    for (pc in p_cols) {
      v <- df[[pc]][i]
      if (!is.na(v) && v < 0.05)
        addStyle(wb, sheet, sig_c, rows=r, cols=pc, stack=TRUE)
    }
  }
}

# Sheet 1: 全部蛋白
s1 <- "Full MR Results (2094 proteins)"
addWorksheet(wb, s1)
writeData(wb, s1,
  data.frame(x=paste0(
    "Supplementary Table S1: MR results for all 2,094 circulating proteins across the ATTR-enriched clinical axis and negative controls. ",
    "CTS=carpal tunnel syndrome (UKB, 12,312 cases); CA=systemic amyloidosis (FinnGen E85, 226 cases); ",
    "HF=heart failure (HERMES, 47,309 cases); T2D/Asthma=negative controls. ",
    "Green bold = P<0.05. Yellow = ITIH1-locus signals & BAG3. Blue = other triple-validated."
  )),
  startRow=1, startCol=1, colNames=FALSE)
addStyle(wb, s1, title_s, rows=1, cols=1)
mergeCells(wb, s1, cols=1:22, rows=1)
setRowHeights(wb, s1, rows=1, heights=55)
writeData(wb, s1, merged, startRow=2, headerStyle=hdr, withFilter=TRUE)
setRowHeights(wb, s1, rows=2, heights=45)
setColWidths(wb, s1, cols=1:22, widths=col_w)
freezePane(wb, s1, firstActiveRow=3, firstActiveCol=3)
apply_row_styles(s1, merged)

# Sheet 2: 主要信号
s2 <- "Primary Hits Summary"
addWorksheet(wb, s2)
hits <- merged %>% filter(`Protein ID` %in% top_ids | `Triple validated (CTS+CA+HF)` == TRUE)
writeData(wb, s2,
  data.frame(x="Primary signals: ITIH1-locus proteins, BAG3, and all triple-validated proteins"),
  startRow=1, startCol=1, colNames=FALSE)
addStyle(wb, s2, title_s, rows=1, cols=1)
mergeCells(wb, s2, cols=1:22, rows=1)
setRowHeights(wb, s2, rows=1, heights=35)
writeData(wb, s2, hits, startRow=2, headerStyle=hdr, withFilter=FALSE)
setRowHeights(wb, s2, rows=2, heights=45)
setColWidths(wb, s2, cols=1:22, widths=col_w)
apply_row_styles(s2, hits)

# ---- 7. 保存 ------------------------------------------------

out <- "results/Supplementary_Table_S1_Full_MR_Results.xlsx"
saveWorkbook(wb, out, overwrite=TRUE)
cat("已保存至:", file.path(getwd(), out), "\n")
cat("Sheet 1:", nrow(merged), "个蛋白 | Sheet 2:", nrow(hits), "个主要信号\n")
