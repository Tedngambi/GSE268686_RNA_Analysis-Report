# ================================================================
# QuantSeq RNA-seq Analysis - GSE288686
# Bc. Teddy Ngambi | Masaryk University, Brno
# CRISPR/Cas9 Epigenetic KO in Arabidopsis thaliana
# SetA: PRDM9 KO vs Control | SetE: SDG2 KO vs Control
# ================================================================
library(usethis)
create_github_token()
gitcreds::gitcreds_set()

use_github()
use_readme_md()

library(DESeq2)       # Negative Binomial DE model
library(edgeR)        # CPM normalisation + TMM scaling
library(ggplot2)      # visualisation
library(ggrepel)      # gene labels on volcano plots
library(pheatmap)     # clustered heatmaps
library(readxl)       # read count .xlsx files
library(tidyverse)    # data wrangling
library(RColorBrewer) # colour palettes

# ── 1. LOAD RAW COUNT DATA ────────────────────────────────────
setA_raw <- read_excel("GSE288686_QuantSeq_SetA_Gene_Counts_Table.xlsx") |>
  column_to_rownames("Gene_ID")
setE_raw <- read_excel("GSE288686_QuantSeq_SetE_Gene_Counts_Table.xlsx") |>
  column_to_rownames("Gene_ID")

colnames(setA_raw) <- c(
  "PRDM9g4_r1","PRDM9g4_r2","PRDM9g4_r3",
  "PRDM9noG_r1","PRDM9noG_r2","PRDM9noG_r3",
  "rdr6_r1","rdr6_r2","rdr6_r3")
colnames(setE_raw) <- c(
  "SDG2g4_r1","SDG2g4_r2","SDG2g4_r3",
  "rdr6_r1","rdr6_r2","rdr6_r3",
  "SDG2noG_r1","SDG2noG_r2","SDG2noG_r3")

# ── 2. CPM NORMALISATION + FILTER ────────────────────────────
dge_A <- DGEList(counts=setA_raw) |> calcNormFactors(method="TMM")
dge_E <- DGEList(counts=setE_raw) |> calcNormFactors(method="TMM")

keepA <- rowSums(cpm(dge_A) > 1) >= 2
keepE <- rowSums(cpm(dge_E) > 1) >= 2

dge_A <- dge_A[keepA, , keep.lib.sizes=FALSE]
dge_E <- dge_E[keepE, , keep.lib.sizes=FALSE]

logA <- cpm(dge_A, log=TRUE, prior.count=0.5)
logE <- cpm(dge_E, log=TRUE, prior.count=0.5)

cat("After filter: SetA", nrow(logA), "genes | SetE", nrow(logE), "genes
")

# ── 3. PCA QC ────────────────────────────────────────────────
pca_plot <- function(log_mat, title, group_fn, pal) {
  pca   <- prcomp(t(log_mat), scale.=TRUE)
  var_e <- summary(pca)$importance[2,1:2]*100
  df    <- data.frame(pca$x[,1:2],
                      sample=rownames(pca$x),
                      group=sapply(rownames(pca$x), group_fn))
  ggplot(df, aes(PC1,PC2,colour=group,label=sample)) +
    geom_point(size=5,alpha=0.9) + geom_text_repel(size=3.2,show.legend=FALSE) +
    scale_colour_manual(values=pal) +
    labs(title=title,
         x=paste0("PC1 (",round(var_e[1],1),"% var)"),
         y=paste0("PC2 (",round(var_e[2],1),"% var)"), colour=NULL) +
    theme_bw(base_size=12)
}

grpA <- function(s) ifelse(grepl("g4",s),"PRDM9 KO",
                           ifelse(grepl("noG",s),"PRDM9 Ctrl","rdr6"))
grpE <- function(s) ifelse(grepl("g4",s),"SDG2 KO",
                           ifelse(grepl("noG",s),"SDG2 Ctrl","rdr6"))

# ── 4. DESEQ2 DIFFERENTIAL EXPRESSION ───────────────────────
run_deseq2 <- function(raw_counts, ko_cols, ctrl_cols) {
  sub  <- round(raw_counts[, c(ko_cols, ctrl_cols)])
  cond <- factor(c(rep("KO",length(ko_cols)),
                   rep("Control",length(ctrl_cols))),
                 levels=c("Control","KO"))
  dds  <- DESeqDataSetFromMatrix(sub,
                                 data.frame(condition=cond,
                                            row.names=c(ko_cols,ctrl_cols)),
                                 ~condition)
  dds  <- dds[rowSums(counts(dds))>=10,]
  dds  <- DESeq(dds, quiet=TRUE)
  res  <- results(dds, contrast=c("condition","KO","Control"), alpha=0.05)
  lfcShrink(dds, coef="condition_KO_vs_Control", type="apeglm", res=res)
}

deA <- run_deseq2(setA_raw,
                  c("PRDM9g4_r1","PRDM9g4_r2","PRDM9g4_r3"),
                  c("PRDM9noG_r1","PRDM9noG_r2","PRDM9noG_r3"))
deE <- run_deseq2(setE_raw,
                  c("SDG2g4_r1","SDG2g4_r2","SDG2g4_r3"),
                  c("SDG2noG_r1","SDG2noG_r2","SDG2noG_r3"))

classify_DEGs <- function(res, lfc_cut=1, p_cut=0.05) {
  as.data.frame(res) |>
    tibble::rownames_to_column("Gene_ID") |>
    mutate(sig = case_when(
      log2FoldChange >= lfc_cut  & padj < p_cut ~ "Up-regulated",
      log2FoldChange <= -lfc_cut & padj < p_cut ~ "Down-regulated",
      TRUE                                        ~ "Not significant"
    ))
}

deA_df <- classify_DEGs(deA)
deE_df <- classify_DEGs(deE)
cat("
SetA DEGs:", table(deA_df$sig), "
SetE DEGs:", table(deE_df$sig), "
")

# ── 5. VOLCANO PLOT ──────────────────────────────────────────
volcano_plot <- function(df, title) {
  top <- df |> filter(sig!="Not significant") |>
    arrange(padj) |> head(12)
  ggplot(df, aes(log2FoldChange, -log10(pmax(padj,1e-10)), colour=sig)) +
    geom_point(alpha=0.55, size=1.5) +
    geom_point(data=top, size=2.8) +
    geom_text_repel(data=top, aes(label=Gene_ID), size=2.8,
                    show.legend=FALSE, max.overlaps=12) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="#666") +
    geom_vline(xintercept=c(-1,1), linetype="dotted", colour="#aaa") +
    scale_colour_manual(values=c(
      "Up-regulated"="#C0392B","Down-regulated"="#2471A3",
      "Not significant"="#CCCCCC")) +
    labs(title=title,
         x="log2 Fold Change (KO vs Control)",
         y="-log10(adjusted p-value)", colour=NULL) +
    theme_bw(base_size=12)
}

# ── 6. HEATMAP ───────────────────────────────────────────────
deg_heatmap <- function(de_df, log_mat, n=40, title) {
  top <- de_df |> filter(sig!="Not significant") |>
    arrange(padj) |> head(n) |> pull(Gene_ID)
  pheatmap(t(scale(t(log_mat[top,]))), main=title,
           color=colorRampPalette(rev(brewer.pal(11,"RdBu")))(100),
           breaks=seq(-2.5,2.5,length.out=101),
           fontsize_row=7, border_color=NA, gaps_col=c(3,6))
}

# ── 7. SAVE ───────────────────────────────────────────────────
ggsave("volcano_setA.png", volcano_plot(deA_df,"SetA PRDM9 KO"), width=10,height=6,dpi=150)
ggsave("volcano_setE.png", volcano_plot(deE_df,"SetE SDG2 KO"),  width=10,height=6,dpi=150)
write.csv(arrange(deA_df,padj), "DEG_SetA_PRDM9_DESeq2.csv", row.names=FALSE)
write.csv(arrange(deE_df,padj), "DEG_SetE_SDG2_DESeq2.csv",  row.names=FALSE)

cat("
Done. SetA DEGs:", sum(deA_df$sig!="Not significant"),
    " | SetE DEGs:", sum(deE_df$sig!="Not significant"), "
")
