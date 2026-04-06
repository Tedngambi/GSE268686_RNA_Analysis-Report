# ============================================================
#  GSE288686 · QuantSeq RNA-seq Analysis
#  Bc. Teddy Ng'ambi | Masaryk University, Brno
#
#  CRISPR/Cas9 Epigenetic Regulator KO in Arabidopsis thaliana
#    Set A: PRDM9 KO (StPRDM9g4) vs No-guide Control (StPRDM9noG)
#    Set E: SDG2  KO (StSDG2g4)  vs No-guide Control (StSDG2noG)
#    Cross-group: rdr6 epigenetic mutant (present in both sets)
#
#  Figures produced (matching the HTML interactive report):
#    Fig  1 – Library size barplots
#    Fig  2 – log2(CPM+0.5) density distributions
#    Fig  3 – Pearson correlation heatmaps
#    Fig  4 – PCA biplots
#    Fig  5 – Volcano plots  (KO vs control)
#    Fig  6 – MA plots       (KO vs control)
#    Fig  7 – DEG count bar charts across threshold grid
#    Fig  8 – Z-score expression heatmaps (top 40 DEGs)
#    Fig  9 – Individual gene boxplots (top 6 DEGs, Set E)
#    Fig 10 – Volcano plots  (rdr6 vs CRISPR KO, both sets)
#
#  Statistical method: Welch t-test on log2(CPM+0.5) values,
#  Benjamini-Hochberg FDR correction per comparison.
#  Primary thresholds: padj < 0.1 AND |log2FC| > 0.5.
#
#  Required packages:
#    install.packages(c("tidyverse","readxl","ggrepel",
#                       "pheatmap","RColorBrewer","patchwork",
#                       "ggtext","scales"))
#    BiocManager::install("edgeR")   # for CPM / calcNormFactors
# ============================================================

library(tidyverse)
library(readxl)
library(edgeR)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(scales)
library(patchwork)

# ── Output directory ─────────────────────────────────────────
dir.create("figures", showWarnings = FALSE)

# ============================================================
# SECTION 1 · DATA IMPORT & SAMPLE METADATA
# ============================================================

setA_raw <- read_excel("GSE288686_QuantSeq_SetA_Gene_Counts_Table.xlsx") |>
  column_to_rownames("Gene_ID") |>
  as.matrix()

setE_raw <- read_excel("GSE288686_QuantSeq_SetE_Gene_Counts_Table.xlsx") |>
  column_to_rownames("Gene_ID") |>
  as.matrix()

# Assign clean column names (order matches the Excel files)
colnames(setA_raw) <- c(
  "PRDM9g4_r1", "PRDM9g4_r2", "PRDM9g4_r3",
  "PRDM9noG_r1", "PRDM9noG_r2", "PRDM9noG_r3",
  "rdr6A_r1",  "rdr6A_r2",  "rdr6A_r3"
)
colnames(setE_raw) <- c(
  "SDG2g4_r1", "SDG2g4_r2", "SDG2g4_r3",
  "rdr6E_r1",  "rdr6E_r2",  "rdr6E_r3",
  "SDG2noG_r1", "SDG2noG_r2", "SDG2noG_r3"
)

# Sample metadata data frames
meta_A <- tibble(
  sample = colnames(setA_raw),
  group  = c(rep("PRDM9 KO", 3), rep("PRDM9 Ctrl", 3), rep("rdr6", 3)),
  set    = "Set A"
)
meta_E <- tibble(
  sample = colnames(setE_raw),
  group  = c(rep("SDG2 KO", 3), rep("rdr6", 3), rep("SDG2 Ctrl", 3)),
  set    = "Set E"
)

# Colour palettes (consistent across all figures)
pal_A <- c("PRDM9 KO" = "#C0392B", "PRDM9 Ctrl" = "#2471A3", "rdr6" = "#E67E22")
pal_E <- c("SDG2 KO"  = "#C0392B", "SDG2 Ctrl"  = "#2471A3", "rdr6" = "#E67E22")

cat("Raw counts loaded.\n",
    "  Set A:", nrow(setA_raw), "genes ×", ncol(setA_raw), "samples\n",
    "  Set E:", nrow(setE_raw), "genes ×", ncol(setE_raw), "samples\n")

# ============================================================
# SECTION 2 · CPM NORMALISATION & FILTERING
# ============================================================

# TMM normalisation via edgeR, then log2(CPM + 0.5)
normalise_counts <- function(raw) {
  dge <- DGEList(counts = round(raw)) |>
    calcNormFactors(method = "TMM")
  cpm(dge, log = TRUE, prior.count = 0.5)
}

log_A_all <- normalise_counts(setA_raw)
log_E_all <- normalise_counts(setE_raw)

# Filter: keep genes with CPM > 1 in ≥ 2 samples
dge_A <- DGEList(counts = round(setA_raw)) |> calcNormFactors("TMM")
dge_E <- DGEList(counts = round(setE_raw)) |> calcNormFactors("TMM")

keepA <- rowSums(cpm(dge_A) > 1) >= 2
keepE <- rowSums(cpm(dge_E) > 1) >= 2

log_A <- normalise_counts(setA_raw[keepA, ])
log_E <- normalise_counts(setE_raw[keepE, ])

lib_A <- colSums(setA_raw) / 1e6   # library sizes in millions
lib_E <- colSums(setE_raw) / 1e6

cat("After CPM > 1 filter:\n",
    "  Set A:", nrow(log_A), "genes (was", nrow(log_A_all), ")\n",
    "  Set E:", nrow(log_E), "genes (was", nrow(log_E_all), ")\n")

# ============================================================
# SECTION 3 · DIFFERENTIAL EXPRESSION (Welch t-test + BH)
# ============================================================

# Run pairwise Welch t-test for all genes, return tidy data frame
# Add gene column to DE table
de_E2 <- de_E |>
  tibble::rownames_to_column("gene")

# Top 6 DEGs
top6_genes <- de_E2 |>
  filter(sig != "Not significant") |>
  arrange(padj) |>
  slice_head(n = 6) |>
  pull(gene)

# Build tidy data frame
box_df <- as_tibble(log_E[top6_genes, , drop = FALSE],
                    rownames = "gene") |>
  pivot_longer(-gene, names_to = "sample", values_to = "log2cpm") |>
  left_join(meta_E, by = "sample") |>
  left_join(select(de_E2, gene, log2FC, padj), by = "gene") |>
  mutate(
    gene  = factor(gene, levels = top6_genes),
    group = factor(group, levels = c("SDG2 Ctrl", "SDG2 KO", "rdr6")),
    label = paste0(
      gene, "\n",
      "log\u2082FC = ", round(log2FC, 2),
      " | padj = ", formatC(padj, format = "e", digits = 2)
    )
  )

fig9 <- ggplot(box_df, aes(x = group, y = log2cpm, fill = group, colour = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.35, linewidth = 0.7) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.9) +
  scale_fill_manual(values = pal_E) +
  scale_colour_manual(values = pal_E) +
  facet_wrap(~label, nrow = 2, scales = "free_y") +
  labs(
    title    = "Fig 9 · Set E: Individual Gene Expression — Top 6 DEGs",
    subtitle = "Points = biological replicates (n = 3 per group)",
    x        = NULL,
    y        = "log\u2082(CPM + 0.5)",
    fill     = NULL,
    colour   = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x    = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    strip.text      = element_text(size = 8)
  )
# Set A: PRDM9 KO vs PRDM9noG control
de_A <- run_de(log_A,
               ko_cols   = c("PRDM9g4_r1",  "PRDM9g4_r2",  "PRDM9g4_r3"),
               ctrl_cols = c("PRDM9noG_r1", "PRDM9noG_r2", "PRDM9noG_r3"))

# Set E: SDG2 KO vs SDG2noG control
de_E <- run_de(log_E,
               ko_cols   = c("SDG2g4_r1",  "SDG2g4_r2",  "SDG2g4_r3"),
               ctrl_cols = c("SDG2noG_r1", "SDG2noG_r2", "SDG2noG_r3"))

# Classify DEGs at primary threshold (padj < 0.1, |log2FC| > 0.5)
classify_degs <- function(df, lfc_cut = 0.5, p_cut = 0.1) {
  df |> mutate(
    sig = case_when(
      log2FC >=  lfc_cut & padj < p_cut ~ "Up-regulated",
      log2FC <= -lfc_cut & padj < p_cut ~ "Down-regulated",
      TRUE                              ~ "Not significant"
    )
  )
}


run_de <- function(log_mat, ko_cols, ctrl_cols) {
  
  results <- apply(log_mat, 1, function(g) {
    
    g_ko   <- g[ko_cols]
    g_ctrl <- g[ctrl_cols]
    
    # Skip genes with zero variance
    if (sd(g_ko) == 0 && sd(g_ctrl) == 0) {
      return(c(log2FC = NA, pval = NA))
    }
    
    log2FC <- mean(g_ko) - mean(g_ctrl)
    
    pval <- tryCatch(
      t.test(g_ko, g_ctrl, var.equal = FALSE)$p.value,
      error = function(e) NA
    )
    
    c(log2FC = log2FC, pval = pval)
  })
  
  results <- as.data.frame(t(results))
  results$padj <- p.adjust(results$pval, method = "BH")
  results
}

de_A <- run_de(
  log_A,
  ko_cols   = c("PRDM9g4_r1", "PRDM9g4_r2", "PRDM9g4_r3"),
  ctrl_cols = c("PRDM9noG_r1", "PRDM9noG_r2", "PRDM9noG_r3")
)

de_E <- run_de(
  log_E,
  ko_cols   = c("SDG2g4_r1", "SDG2g4_r2", "SDG2g4_r3"),
  ctrl_cols = c("SDG2noG_r1", "SDG2noG_r2", "SDG2noG_r3")
)


de_A <- classify_degs(de_A)
de_E <- classify_degs(de_E)
cat("\nSet A (PRDM9 KO vs Ctrl):\n"); print(table(de_A$sig))
cat("\nSet E (SDG2 KO vs Ctrl):\n");  print(table(de_E$sig))

# rdr6 comparisons for Figure 10
de_rdr6_A <- run_de(log_A,
                    ko_cols   = c("rdr6A_r1",   "rdr6A_r2",   "rdr6A_r3"),
                    ctrl_cols = c("PRDM9g4_r1", "PRDM9g4_r2", "PRDM9g4_r3"))
de_rdr6_A <- classify_degs(de_rdr6_A)

de_rdr6_E <- run_de(log_E,
                    ko_cols   = c("rdr6E_r1",  "rdr6E_r2",  "rdr6E_r3"),
                    ctrl_cols = c("SDG2g4_r1", "SDG2g4_r2", "SDG2g4_r3"))
de_rdr6_E <- classify_degs(de_rdr6_E)

# ============================================================
# FIGURE 1 · LIBRARY SIZE BARPLOTS
# ============================================================

plot_libsize <- function(lib_sizes, meta, pal, title) {
  df <- tibble(
    sample     = names(lib_sizes),
    reads_M    = lib_sizes,
    group      = meta$group[match(names(lib_sizes), meta$sample)]
  ) |>
    mutate(sample = fct_inorder(sample))
  
  mean_reads <- mean(df$reads_M)
  
  ggplot(df, aes(x = sample, y = reads_M, fill = group)) +
    geom_col(colour = "white", linewidth = 0.3) +
    geom_hline(yintercept = mean_reads, linetype = "dashed",
               colour = "#555", linewidth = 0.7) +
    annotate("text", x = Inf, y = mean_reads,
             label = paste0("mean = ", round(mean_reads, 2), "M"),
             hjust = 1.05, vjust = -0.4, size = 3, colour = "#555") +
    scale_fill_manual(values = pal) +
    scale_y_continuous(labels = label_comma(suffix = " M"),
                       expand = expansion(mult = c(0, 0.08))) +
    labs(title = title,
         subtitle = "Dashed line = group mean library size",
         x = NULL, y = "Total reads (millions)", fill = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "top",
          panel.grid.major.x = element_blank())
}

fig1_A <- plot_libsize(lib_A, meta_A, pal_A,
                       "Fig 1A · Set A: Library Sizes — PRDM9 KO experiment")
fig1_E <- plot_libsize(lib_E, meta_E, pal_E,
                       "Fig 1B · Set E: Library Sizes — SDG2 KO experiment")

fig1 <- fig1_A / fig1_E
ggsave("figures/fig01_library_sizes.pdf", fig1, width = 10, height = 9, dpi = 150)
ggsave("figures/fig01_library_sizes.png", fig1, width = 10, height = 9, dpi = 150)
cat("Figure 1 saved.\n")

# ============================================================
# FIGURE 2 · EXPRESSION DENSITY DISTRIBUTIONS
# ============================================================

plot_density <- function(log_mat, meta, pal, title) {
  as_tibble(log_mat, rownames = "Gene_ID") |>
    pivot_longer(-Gene_ID, names_to = "sample", values_to = "log2cpm") |>
    left_join(meta, by = "sample") |>
    mutate(sample = fct_inorder(sample)) |>
    ggplot(aes(x = log2cpm, colour = group, group = sample)) +
    geom_density(linewidth = 0.7, alpha = 0.8) +
    scale_colour_manual(values = pal) +
    labs(title = title,
         subtitle = "One curve per sample; colour = group",
         x = "log\u2082(CPM + 0.5)", y = "Density", colour = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")
}

fig2_A <- plot_density(log_A, meta_A, pal_A,
                       "Fig 2A · Set A: Normalised Expression Density")
fig2_E <- plot_density(log_E, meta_E, pal_E,
                       "Fig 2B · Set E: Normalised Expression Density")

fig2 <- fig2_A | fig2_E
ggsave("figures/fig02_density_distributions.pdf", fig2, width = 12, height = 5, dpi = 150)
ggsave("figures/fig02_density_distributions.png", fig2, width = 12, height = 5, dpi = 150)
cat("Figure 2 saved.\n")

# ============================================================
# FIGURE 3 · PEARSON CORRELATION HEATMAPS
# ============================================================

plot_cor_heatmap <- function(log_mat, meta, title, filename) {
  # Reorder columns: KO | Ctrl | rdr6 for visual clarity
  cor_mat <- cor(log_mat, method = "pearson")
  
  ann_col <- data.frame(
    Group = meta$group[match(colnames(cor_mat), meta$sample)],
    row.names = colnames(cor_mat)
  )
  
  # Green-yellow palette matching dashboard
  pal_cor <- colorRampPalette(c("#FFFDE7", "#A5D6A7", "#1B5E20"))(100)
  
  pheatmap(cor_mat,
           color            = pal_cor,
           breaks           = seq(0.85, 1.0, length.out = 101),
           annotation_col   = ann_col,
           annotation_row   = ann_col,
           display_numbers  = TRUE,
           number_format    = "%.3f",
           fontsize_number  = 7,
           fontsize         = 9,
           border_color     = "white",
           main             = title,
           filename         = filename)
}



# PNG versions
plot_cor_heatmap(log_A, meta_A,
                 "Fig 3A · Set A: Pearson Correlation Heatmap",
                 "figures/fig03A_correlation_setA.png")
plot_cor_heatmap(log_E, meta_E,
                 "Fig 3B · Set E: Pearson Correlation Heatmap",
                 "figures/fig03B_correlation_setE.png")
cat("Figure 3 saved.\n")

# ============================================================
# FIGURE 4 · PCA BIPLOTS
# ============================================================

plot_pca <- function(log_mat, meta, pal, title) {
  pca     <- prcomp(t(log_mat), scale. = TRUE)
  var_exp <- summary(pca)$importance[2, 1:2] * 100
  
  df <- as_tibble(pca$x[, 1:2], rownames = "sample") |>
    left_join(meta, by = "sample")
  
  ggplot(df, aes(x = PC1, y = PC2, colour = group, label = sample)) +
    geom_point(size = 5, alpha = 0.9) +
    geom_text_repel(size = 3, show.legend = FALSE,
                    box.padding = 0.4, max.overlaps = 15) +
    scale_colour_manual(values = pal) +
    labs(title    = title,
         x        = paste0("PC1 (", round(var_exp[1], 1), "% variance)"),
         y        = paste0("PC2 (", round(var_exp[2], 1), "% variance)"),
         colour   = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")
}

fig4_A <- plot_pca(log_A, meta_A, pal_A,
                   "Fig 4A · Set A: PCA — PRDM9 experiment")
fig4_E <- plot_pca(log_E, meta_E, pal_E,
                   "Fig 4B · Set E: PCA — SDG2 experiment")

fig4 <- fig4_A | fig4_E
ggsave("figures/fig04_pca.pdf", fig4, width = 12, height = 5.5, dpi = 150)
ggsave("figures/fig04_pca.png", fig4, width = 12, height = 5.5, dpi = 150)
cat("Figure 4 saved.\n")

# ============================================================
# FIGURE 5 · VOLCANO PLOTS  (KO vs Control)
# ============================================================

PADJ_CUT <- 0.1
LFC_CUT  <- 0.5

sig_pal <- c(
  "Up-regulated"    = "#C0392B",
  "Down-regulated"  = "#2471A3",
  "Not significant" = "#CCCCCC"
)

plot_volcano <- function(de_df, title, label_n = 12) {
  top_genes <- de_df |>
    filter(sig != "Not significant") |>
    arrange(padj) |>
    slice_head(n = label_n)
  
  ggplot(de_df, aes(x = log2FC, y = -log10(pmax(padj, 1e-10)),
                    colour = sig)) +
    geom_point(alpha = 0.5, size = 1.4) +
    geom_point(data = top_genes, size = 2.8) +
    geom_text_repel(data = top_genes, aes(label = Gene_ID),
                    size = 2.8, show.legend = FALSE,
                    box.padding = 0.35, max.overlaps = 14) +
    geom_hline(yintercept = -log10(PADJ_CUT),
               linetype = "dashed", colour = "#555", linewidth = 0.7) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT),
               linetype = "dotted", colour = "#999", linewidth = 0.6) +
    scale_colour_manual(values = sig_pal) +
    annotate("text", x = Inf, y = -log10(PADJ_CUT),
             label = paste0("padj = ", PADJ_CUT),
             hjust = 1.05, vjust = -0.4, size = 3, colour = "#555") +
    labs(title    = title,
         subtitle = paste0("Thresholds: padj < ", PADJ_CUT,
                           " | |log\u2082FC| > ", LFC_CUT,
                           "   [Up: ", sum(de_df$sig == "Up-regulated"),
                           " | Down: ", sum(de_df$sig == "Down-regulated"), "]"),
         x        = "log\u2082 Fold Change (KO vs Control)",
         y        = "-log\u2081\u2080(adjusted p-value)",
         colour   = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")
}

plot_volcano <- function(de_df, title, label_n = 12) {
  
  de_df <- de_df |>
    tibble::rownames_to_column("gene")
  
  top_genes <- de_df |>
    dplyr::filter(sig != "Not significant") |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = label_n)
  
  ggplot(de_df, aes(x = log2FC,
                    y = -log10(pmax(padj, 1e-10)),
                    colour = sig)) +
    geom_point(alpha = 0.5, size = 1.4) +
    geom_point(data = top_genes, size = 2.8) +
    ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 2.8,
      show.legend = FALSE,
      box.padding = 0.35,
      max.overlaps = 14
    ) +
    geom_hline(yintercept = -log10(PADJ_CUT),
               linetype = "dashed", colour = "#555", linewidth = 0.7) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT),
               linetype = "dotted", colour = "#999", linewidth = 0.6) +
    scale_colour_manual(values = sig_pal) +
    annotate("text", x = Inf, y = -log10(PADJ_CUT),
             label = paste0("padj = ", PADJ_CUT),
             hjust = 1.05, vjust = -0.4, size = 3, colour = "#555") +
    labs(
      title    = title,
      subtitle = paste0(
        "Thresholds: padj < ", PADJ_CUT,
        " | |log₂FC| > ", LFC_CUT,
        "   [Up: ", sum(de_df$sig == "Up-regulated"),
        " | Down: ", sum(de_df$sig == "Down-regulated"), "]"
      ),
      x = "log₂ Fold Change (KO vs Control)",
      y = "-log₁₀(adjusted p-value)",
      colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")
}


fig5_A <- plot_volcano(de_A, "Fig 5A · Set A: Volcano — PRDM9 KO vs Control")
fig5_E <- plot_volcano(de_E, "Fig 5B · Set E: Volcano — SDG2 KO vs Control")

fig5 <- fig5_A | fig5_E
ggsave("figures/fig05_volcano_KO_vs_ctrl.pdf", fig5, width = 13, height = 6, dpi = 150)
ggsave("figures/fig05_volcano_KO_vs_ctrl.png", fig5, width = 13, height = 6, dpi = 150)
cat("Figure 5 saved.\n")

# ============================================================
# FIGURE 6 · MA PLOTS  (KO vs Control)
# ============================================================
de_A$mean_expr <- rowMeans(log_A)
de_E$mean_expr <- rowMeans(log_E)
plot_ma <- function(de_df, title, label_n = 10) {
  
  de_df <- de_df |>
    tibble::rownames_to_column("gene")
  
  top_genes <- de_df |>
    dplyr::filter(sig != "Not significant") |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = label_n)
  
  ggplot(de_df, aes(x = mean_expr, y = log2FC, colour = sig)) +
    geom_point(alpha = 0.45, size = 1.2) +
    geom_point(data = top_genes, size = 2.5) +
    ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 2.6,
      show.legend = FALSE,
      box.padding = 0.35,
      max.overlaps = 12
    ) +
    geom_hline(yintercept = 0, colour = "#555", linewidth = 0.8) +
    geom_hline(
      yintercept = c(-LFC_CUT, LFC_CUT),
      linetype = "dashed",
      colour = "#aaa",
      linewidth = 0.5
    ) +
    scale_colour_manual(values = sig_pal) +
    labs(
      title    = title,
      subtitle = "M = log₂FC  |  A = mean log₂(CPM+0.5)",
      x        = "A — mean log₂(CPM + 0.5)",
      y        = "M — log₂ Fold Change (KO vs Control)",
      colour   = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")
}

fig6_A <- plot_ma(de_A, "Fig 6A · Set A: MA Plot — PRDM9 KO vs Control")
fig6_E <- plot_ma(de_E, "Fig 6B · Set E: MA Plot — SDG2 KO vs Control")

fig6 <- fig6_A | fig6_E
ggsave("figures/fig06_MA_plots.pdf", fig6, width = 13, height = 6, dpi = 150)
ggsave("figures/fig06_MA_plots.png", fig6, width = 13, height = 6, dpi = 150)
cat("Figure 6 saved.\n")

# ============================================================
# FIGURE 7 · DEG COUNTS ACROSS THRESHOLD GRID
# ============================================================

# 4 threshold combinations: 2 padj × 2 |log2FC| cuts
threshold_grid <- expand.grid(
  lfc_cut  = c(0.5, 1.0),
  padj_cut = c(0.1, 0.05)
)

count_degs <- function(de_df, lfc_cut, padj_cut) {
  up   <- sum(de_df$log2FC >=  lfc_cut & de_df$padj < padj_cut, na.rm = TRUE)
  down <- sum(de_df$log2FC <= -lfc_cut & de_df$padj < padj_cut, na.rm = TRUE)
  c(up = up, down = -down)  # negative so down bars go below zero
}

thresh_df <- threshold_grid |>
  as_tibble() |>
  mutate(
    label   = paste0("padj<", padj_cut, "\n|LFC|>", lfc_cut),
    setA_up = mapply(function(l, p) count_degs(de_A, l, p)["up"],   lfc_cut, padj_cut),
    setA_dn = mapply(function(l, p) count_degs(de_A, l, p)["down"], lfc_cut, padj_cut),
    setE_up = mapply(function(l, p) count_degs(de_E, l, p)["up"],   lfc_cut, padj_cut),
    setE_dn = mapply(function(l, p) count_degs(de_E, l, p)["down"], lfc_cut, padj_cut)
  )

# Tidy to long format for ggplot
thresh_long <- thresh_df |>
  pivot_longer(cols = starts_with("set"),
               names_to = c("Set", "Direction"),
               names_pattern = "set([AE])_([a-z]+)") |>
  mutate(
    Set       = ifelse(Set == "A", "Set A — PRDM9 KO", "Set E — SDG2 KO"),
    Direction = ifelse(Direction == "up", "Up-regulated", "Down-regulated"),
    label     = factor(label, levels = unique(label))
  )

fig7 <- ggplot(thresh_long, aes(x = label, y = value,
                                fill = Direction, alpha = Set)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_hline(yintercept = 0, colour = "#333") +
  scale_fill_manual(values = c("Up-regulated" = "#C0392B",
                               "Down-regulated" = "#2471A3")) +
  scale_alpha_manual(values = c("Set A — PRDM9 KO" = 0.55,
                                "Set E — SDG2 KO"  = 1.00)) +
  scale_y_continuous(labels = abs) +
  labs(title    = "Fig 7 · DEG Counts Across Threshold Combinations",
       subtitle = "Faded bars = Set A (PRDM9, no significant DEGs); solid = Set E (SDG2)",
       x        = "Threshold combination",
       y        = "Number of DEGs",
       fill     = NULL, alpha = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

ggsave("figures/fig07_DEG_threshold_grid.pdf", fig7, width = 9, height = 5.5, dpi = 150)
ggsave("figures/fig07_DEG_threshold_grid.png", fig7, width = 9, height = 5.5, dpi = 150)
cat("Figure 7 saved.\n")

# ============================================================
# FIGURE 8 · Z-SCORE EXPRESSION HEATMAPS  (top 40 DEGs)
# ============================================================
plot_deg_heatmap <- function(de_df, log_mat, meta, n_top = 40,
                             title, filename) {
  
  # Add gene column
  de_df <- de_df |>
    tibble::rownames_to_column("gene")
  
  # Top significant genes
  top_genes <- de_df |>
    dplyr::filter(sig != "Not significant") |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = n_top) |>
    dplyr::pull(gene)
  
  # If no significant genes, fallback to smallest p-values
  if (length(top_genes) == 0) {
    top_genes <- de_df |>
      dplyr::arrange(pval) |>
      dplyr::slice_head(n = n_top) |>
      dplyr::pull(gene)
  }
  
  # Extract matrix
  mat  <- log_mat[top_genes, , drop = FALSE]
  zmat <- t(scale(t(mat)))   # row-wise z-score
  
  # Column annotation
  ann_col <- data.frame(
    Group = meta$group[match(colnames(zmat), meta$sample)],
    row.names = colnames(zmat)
  )
  
  # Gap positions at group boundaries
  grp_vec  <- ann_col$Group
  gaps_col <- which(diff(as.integer(factor(grp_vec))) != 0)
  
  pheatmap(
    zmat,
    color          = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    breaks         = seq(-2.5, 2.5, length.out = 101),
    annotation_col = ann_col,
    cluster_cols   = FALSE,
    gaps_col       = gaps_col,
    fontsize_row   = 6.5,
    fontsize       = 9,
    border_color   = NA,
    main           = title,
    filename       = filename
  )
}

plot_deg_heatmap(de_A, log_A, meta_A,
                 title    = "Fig 8A · Set A: Z-score Heatmap — Top 40 candidates (PRDM9)",
                 filename = "figures/fig08A_heatmap_setA.pdf")
plot_deg_heatmap(de_E, log_E, meta_E,
                 title    = "Fig 8B · Set E: Z-score Heatmap — Top 40 DEGs (SDG2 KO)",
                 filename = "figures/fig08B_heatmap_setE.pdf")

plot_deg_heatmap(de_A, log_A, meta_A,
                 title    = "Fig 8A · Set A: Z-score Heatmap — Top 40 candidates (PRDM9)",
                 filename = "figures/fig08A_heatmap_setA.png")
plot_deg_heatmap(de_E, log_E, meta_E,
                 title    = "Fig 8B · Set E: Z-score Heatmap — Top 40 DEGs (SDG2 KO)",
                 filename = "figures/fig08B_heatmap_setE.png")
cat("Figure 8 saved.\n")

# ============================================================
# FIGURE 9 · INDIVIDUAL GENE BOXPLOTS  (top 6 DEGs, Set E only)
# ============================================================

# Add gene column to DE table
de_E2 <- de_E |>
  tibble::rownames_to_column("gene")

# Top 6 DEGs
top6_genes <- de_E2 |>
  filter(sig != "Not significant") |>
  arrange(padj) |>
  slice_head(n = 6) |>
  pull(gene)

# Build tidy data frame
box_df <- as_tibble(log_E[top6_genes, , drop = FALSE],
                    rownames = "gene") |>
  pivot_longer(-gene, names_to = "sample", values_to = "log2cpm") |>
  left_join(meta_E, by = "sample") |>
  left_join(select(de_E2, gene, log2FC, padj), by = "gene") |>
  mutate(
    gene  = factor(gene, levels = top6_genes),
    group = factor(group, levels = c("SDG2 Ctrl", "SDG2 KO", "rdr6")),
    label = paste0(
      gene, "\n",
      "log\u2082FC = ", round(log2FC, 2),
      " | padj = ", formatC(padj, format = "e", digits = 2)
    )
  )

fig9 <- ggplot(box_df, aes(x = group, y = log2cpm, fill = group, colour = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.35, linewidth = 0.7) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.9) +
  scale_fill_manual(values = pal_E) +
  scale_colour_manual(values = pal_E) +
  facet_wrap(~label, nrow = 2, scales = "free_y") +
  labs(
    title    = "Fig 9 · Set E: Individual Gene Expression — Top 6 DEGs",
    subtitle = "Points = biological replicates (n = 3 per group)",
    x        = NULL,
    y        = "log\u2082(CPM + 0.5)",
    fill     = NULL,
    colour   = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x    = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    strip.text      = element_text(size = 8)
  )

ggsave("figures/fig09_gene_boxplots_setE.pdf", fig9, width = 12, height = 7, dpi = 150)
ggsave("figures/fig09_gene_boxplots_setE.png", fig9, width = 12, height = 7, dpi = 150)
cat("Figure 9 saved.\n")

# ============================================================
# FIGURE 10 · VOLCANO PLOTS  (rdr6 vs CRISPR KO)
# ============================================================

fig10_A <- plot_volcano(
  de_rdr6_A,
  title   = "Fig 10A · Set A: Volcano — rdr6 vs PRDM9g4 KO",
  label_n = 10
)
fig10_E <- plot_volcano(
  de_rdr6_E,
  title   = "Fig 10B · Set E: Volcano — rdr6 vs SDG2g4 KO",
  label_n = 10
)

# Relabel axes to reflect rdr6 as "numerator"
relabel_rdr6 <- function(p) {
  p + labs(x = "log\u2082 Fold Change (rdr6 vs CRISPR KO)")
}

fig10 <- relabel_rdr6(fig10_A) | relabel_rdr6(fig10_E)
ggsave("figures/fig10_volcano_rdr6_vs_KO.pdf", fig10, width = 13, height = 6, dpi = 150)
ggsave("figures/fig10_volcano_rdr6_vs_KO.png", fig10, width = 13, height = 6, dpi = 150)
cat("Figure 10 saved.\n")

# ============================================================
# SAVE DE RESULTS TO CSV
# ============================================================

write_csv(arrange(de_A, padj), "DEG_SetA_PRDM9_results.csv")
write_csv(arrange(de_E, padj), "DEG_SetE_SDG2_results.csv")
write_csv(arrange(de_rdr6_A, padj), "DEG_SetA_rdr6vsPRDM9_results.csv")
write_csv(arrange(de_rdr6_E, padj), "DEG_SetE_rdr6vsSDG2_results.csv")

cat("\n================================================================\n")
cat("All 10 figures saved to figures/  (PDF + PNG)\n")
cat("DE results saved to CSV.\n\n")
cat("Set A (PRDM9 KO vs Ctrl): 0 padj-significant DEGs at padj<0.1\n")
cat("Set E (SDG2 KO vs Ctrl): ", sum(de_E$sig != "Not significant"),
    "DEGs (padj<0.1, |log2FC|>0.5)\n")
cat("  Up-regulated:   ", sum(de_E$sig == "Up-regulated"),   "\n")
cat("  Down-regulated: ", sum(de_E$sig == "Down-regulated"), "\n")
cat("================================================================\n")
