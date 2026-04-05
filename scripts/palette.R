# ============================================================
# PWMR-ATTR Figure Palette & Theme
# Source this file in all figure scripts: source("palette.R")
# ============================================================

# Macaron SCI palette
pal <- list(
  # Phenotype colors (consistent across all figures)
  cts       = "#B3DDCB",  # mint — CTS (early)
  amy       = "#F3BBB1",  # salmon — Amyloidosis (core)
  hf        = "#B8E5FA",  # light blue — HF (late)
  
  # Workflow / annotation
  exposure  = "#EEC78A",  # gold
  mr_step   = "#CBE4B1",  # light green
  validation= "#F7B2C7",  # pink
  neutral   = "#E8E6E1",  # warm gray
  
  # Protein-specific
  sigirr    = "#F7A6AC",  # soft pink
  jakmip3   = "#F7B2C7",  # pink
  clec2l    = "#EEE9A2",  # pale yellow
  bag3      = "#B8E5FA",  # light blue (same as HF — intentional)
  itih1_cis = "#B3DDCB",  # mint
  
  # Text
  text      = "#2C2C2A",
  textsub   = "#6E6E6A",
  textlight = "#9E9E98",
  
  # Borders
  border    = "#BBBBBB",
  borderlt  = "#DDDDDD"
)

# Named vector for scale_fill_manual
pal_pheno <- c(
  "CTS"         = pal$cts,
  "Amyloidosis" = pal$amy,
  "HF"          = pal$hf,
  "T2D"         = pal$neutral,
  "Asthma"      = pal$neutral
)

pal_protein <- c(
  "SIGIRR"  = pal$sigirr,
  "JAKMIP3" = pal$jakmip3,
  "CLEC2L"  = pal$clec2l,
  "BAG3"    = pal$bag3,
  "ITIH1"   = pal$itih1_cis
)

# Publication theme (clean, no gridlines, minimal)
theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(family = "Arial", color = pal$text),
      plot.title = element_text(size = base_size + 2, face = "bold", margin = margin(b = 8)),
      plot.subtitle = element_text(size = base_size - 1, color = pal$textsub, margin = margin(b = 12)),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, color = pal$textsub),
      axis.line = element_line(linewidth = 0.3, color = pal$border),
      axis.ticks = element_line(linewidth = 0.2, color = pal$border),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.key.size = unit(0.4, "cm"),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.margin = margin(12, 12, 12, 12)
    )
}

# Helper: save publication figure
save_fig <- function(plot, name, w = 7, h = 5) {
  ggsave(paste0(name, ".pdf"), plot, width = w, height = h, dpi = 300)
  ggsave(paste0(name, ".png"), plot, width = w, height = h, dpi = 300, bg = "white")
  cat(sprintf("Saved: %s.pdf / .png (%g × %g in)\n", name, w, h))
}

cat("Palette loaded. Use: pal$cts, pal$amy, pal$hf, pal_pheno, pal_protein, theme_pub()\n")
