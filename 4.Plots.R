# ============================================================
# 1. QQ PLOT
# ============================================================

png(file.path(outDir, "QQ_plot.png"), width=800, height=800)

qqt(fit2_age$t, df=fit2_age$df.total,
    main = "QQ Plot (limma moderated t)")

dev.off()

pvals <- all_DMPs$P.Value
pvals <- pvals[!is.na(pvals) & pvals > 0]

png(file.path(outDir, "QQ_plot_pvalues.png"), width=800, height=800)

expected <- -log10(ppoints(length(pvals)))
observed <- -log10(sort(pvals))

plot(expected, observed,
     pch=16, cex=0.6,
     xlab="Expected -log10(p)",
     ylab="Observed -log10(p)",
     main="QQ Plot (p-values)")

abline(0,1,col="red", lwd=2)

dev.off()

lambda <- median(qchisq(1 - pvals, 1)) / qchisq(0.5, 1)
lambda

# ============================================================
# 2. REGULATION INTERPRETATION PLOT
# Predictive Activating vs Predictive Silencing DMPs
# ============================================================

activated_genes <- sig_DMPs_anno %>%
  filter(Regulation == "Promoter Hypomethylation") %>%
  arrange(adj.P.Val) %>%
  mutate(Gene = sapply(strsplit(UCSC_RefGene_Name, ";"),
                        `[`, 1)) %>%
  filter(!is.na(Gene) & Gene != "") %>%
  distinct(Gene, .keep_all=TRUE) %>%
  slice_head(n=20)

silenced_genes <- sig_DMPs_anno %>%
  filter(Regulation == "Promoter Hypermethylation") %>%
  arrange(adj.P.Val) %>%
  mutate(Gene = sapply(strsplit(UCSC_RefGene_Name, ";"),
                        `[`, 1)) %>%
  filter(!is.na(Gene) & Gene != "") %>%
  distinct(Gene, .keep_all=TRUE) %>%
  slice_head(n=20)

plot_reg <- bind_rows(
  mutate(activated_genes, Regulation="Promoter Hypomethylation"),
  mutate(silenced_genes,  Regulation="Promoter Hypermethylation")
)

if (nrow(plot_reg) > 0) {

  plot_reg$neglog10p <- -log10(plot_reg$adj.P.Val)
  plot_reg$Gene      <- factor(
    plot_reg$Gene,
    levels=rev(unique(plot_reg$Gene))
  )

  p_reg_scatter <- ggplot(
    plot_reg,
    aes(x=logFC, y=Gene,
        color=Regulation, size=neglog10p)
  ) +
    geom_point(alpha=0.85) +
    geom_vline(xintercept=0, linewidth=0.3,
               color="grey60", linetype="solid") +
    scale_color_manual(
      values=c("Promoter Hypomethylation"="#2166AC",
               "Promoter Hypermethylation"="#B2182B"),
      name=NULL
    ) +
    scale_size_continuous(
      name=expression(-log[10]~"(FDR)"),
      range=c(3,9)
    ) +
    facet_wrap(~Regulation, scales="free_y") +
    labs(
      x = "M-value difference (logFC)",
      y = NULL
    ) +
    theme_classic(base_size=12) +
    theme(
      axis.title   = element_text(face="bold"),
      axis.text.y  = element_text(color="black", size=9,
                                   face="italic"),
      axis.text.x  = element_text(color="black"),
      strip.text   = element_text(face="bold", size=11),
      plot.title   = element_text(face="bold", size=11),
      legend.title = element_text(face="bold"),
      plot.margin  = margin(10,15,10,10)
    )

  ggsave(file.path(outDir, "Regulation_TopGenes_Scatter.png"),
         p_reg_scatter, width=12, height=10, dpi=300)
  ggsave(file.path(outDir, "Regulation_TopGenes_Scatter.pdf"),
         p_reg_scatter, width=12, height=10)

  cat("Regulation scatter plot saved\n")
}

cat("\n=== Regulation analysis complete ===\n")
cat("Files saved:\n")
cat("  Regulation_TopGenes_Scatter.png/pdf\n")
cat("  Promoter_Hypomethylatio_Genes.csv\n")
cat("  Promoter_Hypermethylation_Genes.csv\n")

# =========================================================
# 3. HEATMAP — TOP 500 DMPs 
# =========================================================

# GET SIGNIFICANT DMPs WITH EFFECT SIZE FILTER
# Use both statistical AND biological significance thresholds

# Primary filter: FDR + meaningful effect size
sig_DMPs <- subset(all_DMPs,
                   adj.P.Val < 0.1 &
                   abs(logFC) > 0.2)

cat("DMPs passing FDR<0.1 and |logFC|>0.2:", nrow(sig_DMPs), "\n")

# Sort by adjusted p-value (most significant first)
sig_DMPs <- sig_DMPs[order(sig_DMPs$adj.P.Val), ]

# SELECT TOP N PROBES

n_probes <- min(500, nrow(sig_DMPs))
top_cpgs <- rownames(sig_DMPs)[1:n_probes]
cat("Probes in heatmap:", n_probes, "\n")

# BUILD BETA MATRIX WITH CORRECT SAMPLE ORDER

beta_cols      <- colnames(beta)
is_control     <- targets$Group[match(beta_cols,
                                       colnames(beta))] == "Control"

control_cols   <- beta_cols[is_control]
puv_cols       <- beta_cols[!is_control]

# Order: Controls left, PUV right
ordered_cols   <- c(control_cols, puv_cols)

# Subset beta matrix
heatmap_mat    <- beta[top_cpgs, ordered_cols]

cat("Heatmap matrix:", nrow(heatmap_mat), "probes x",
    ncol(heatmap_mat), "samples\n")

# SCALE FIRST, THEN CLUSTER

heatmap_scaled <- t(scale(t(heatmap_mat)))

heatmap_scaled[heatmap_scaled >  4] <-  4
heatmap_scaled[heatmap_scaled < -4] <- -4

dist_rows  <- dist(heatmap_scaled, method = "euclidean")
hc_rows    <- hclust(dist_rows, method = "complete")

# ANNOTATION

annotation_col <- data.frame(
  Group = factor(
    c(rep("Control", length(control_cols)),
      rep("PUV",     length(puv_cols))),
    levels = c("Control", "PUV")
  ),
  row.names = ordered_cols
)

stopifnot(all(rownames(annotation_col) == colnames(heatmap_mat)))
cat("Annotation check passed\n")

ann_colors <- list(
  Group = c(Control = "#2166AC",   # blue
            PUV     = "#D6604D")   # red
)

# COLOR PALETTE

heatmap_colors <- colorRampPalette(
  c("#2166AC",   # deep blue (hypo)
    "#F7F7F7",   # white (mid)
    "#D6604D")   # red (hyper)
)(100)

# PLOT

pheatmap(
  heatmap_scaled,

  cluster_rows     = hc_rows,
  cluster_cols     = FALSE,
  
  treeheight_row   = 30,
  treeheight_col   = 0,    

  color            = heatmap_colors,
  breaks           = seq(-4, 4, length.out = 101),

  show_rownames    = FALSE,
  show_colnames    = FALSE,

  annotation_col   = annotation_col,
  annotation_colors = ann_colors,

  border_color     = NA,
  fontsize         = 5,
  annotation_legend = TRUE,
  
  filename = file.path(
    outDir,
    "Heatmap_top500_DMPs.png"
  ),

  width  = 6.69,
  height = 8.5,
  dpi = 300
)
cat("PNG Heatmap saved\n")

pdf(file.path(outDir, "Heatmap_top500_DMPs.pdf"),
    width  = 6.69,
    height = 8.5)

pheatmap(
  heatmap_scaled,
  cluster_rows      = hc_rows,
  cluster_cols      = FALSE,
  treeheight_row    = 30,
  treeheight_col    = 0,
  color             = heatmap_colors,
  breaks            = seq(-4, 4, length.out = 101),
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  annotation_col    = annotation_col,
  annotation_colors = ann_colors,
  border_color      = NA,
  fontsize          = 5,
  annotation_legend = TRUE,
  main = paste0("Top ", n_probes,
                " Significant DMPs (FDR < 0.1, |logFC| > 0.2)\n",
                "PUV vs Control")
)

dev.off()
cat("PDF heatmap saved\n")

# =========================================================
# 4. VOLCANO PLOT
# =========================================================

# CATEGORY ASSIGNMENT ON FULL DATASET

all_DMPs$Category <- "Not significant"

all_DMPs$Category[
  all_DMPs$adj.P.Val < 0.1 &
  all_DMPs$logFC >  0.2
] <- "Hypermethylated"

all_DMPs$Category[
  all_DMPs$adj.P.Val < 0.1 &
  all_DMPs$logFC < -0.2
] <- "Hypomethylated"

# Count each category
cat("=== Volcano category counts ===\n")
print(table(all_DMPs$Category))

# ANNOTATION FOR LABELS

ann      <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
ann_sub  <- ann[, c("Name", "UCSC_RefGene_Name", "chr", "pos")]

ann_sub$Gene <- sapply(
  strsplit(ann_sub$UCSC_RefGene_Name, ";"),
  function(x) if (length(x) == 0 || x[1] == "") NA else x[1]
)

all_DMPs$Gene <- ann_sub$Gene[match(all_DMPs$CpG, ann_sub$Name)]
all_DMPs$chr  <- ann_sub$chr[match(all_DMPs$CpG, ann_sub$Name)]

# SELECT PROBES TO LABEL

label_hyper <- subset(all_DMPs,
                       Category == "Hypermethylated" &
                       !is.na(Gene) & Gene != "")
label_hyper <- label_hyper[order(label_hyper$adj.P.Val), ][1:min(8, nrow(label_hyper)), ]

label_hypo  <- subset(all_DMPs,
                       Category == "Hypomethylated" &
                       !is.na(Gene) & Gene != "")
label_hypo  <- label_hypo[order(label_hypo$adj.P.Val), ][1:min(8, nrow(label_hypo)), ]

label_df    <- rbind(label_hyper, label_hypo)

cat("\nProbes to be labelled:\n")
print(label_df[, c("CpG", "Gene", "delta_beta", "logFC", "adj.P.Val", "Category")])

# AXIS LIMITS FROM DATA

x_max <- max(abs(all_DMPs$logFC), na.rm = TRUE)
x_lim <- c(-ceiling(x_max * 10) / 10 - 0.2,
             ceiling(x_max * 10) / 10 + 0.2)

y_max <- max(-log10(all_DMPs$adj.P.Val[
  is.finite(-log10(all_DMPs$adj.P.Val))
]), na.rm = TRUE)
y_lim <- c(0, ceiling(y_max) + 0.3)

cat("\nx-axis range:", round(x_lim, 3), "\n")
cat("y-axis range:", round(y_lim, 3), "\n")

# PLOT

cat_colors <- c(
  "Hypermethylated" = "#B2182B",   # dark red
  "Hypomethylated"  = "#2166AC",   # dark blue
  "Not significant" = "grey75"
)

cat_sizes <- c(
  "Hypermethylated" = 1.0,
  "Hypomethylated"  = 1.0,
  "Not significant" = 0.4         
)

cat_alpha <- c(
  "Hypermethylated" = 0.8,
  "Hypomethylated"  = 0.8,
  "Not significant" = 0.3         
)

all_DMPs$plot_order <- ifelse(all_DMPs$Category == "Not significant", 1, 2)
plot_data_ordered   <- all_DMPs[order(all_DMPs$plot_order), ]

xlim_use <- 1

p_volcano <- ggplot(
  plot_data_ordered,
  aes(x     = logFC,
      y     = -log10(adj.P.Val),
      color = Category,
      size  = Category,
      alpha = Category)
) +

  geom_point(stroke = 0) +

  geom_hline(yintercept = -log10(0.1),
             linetype = "dashed",
             linewidth = 0.5,
             color = "grey40") +

  geom_vline(xintercept = c(-0.2, 0.2),
             linetype = "dashed",
             linewidth = 0.5,
             color = "grey40") +

  geom_label_repel(
    data         = label_df,
    aes(label    = Gene),
    size         = 2.8,
    fontface     = "italic",
    color        = "black",
    fill         = alpha("white", 0.85),
    label.size   = 0.25,
    box.padding  = 0.4,
    point.padding = 0.3,
    max.overlaps = 20,
    segment.size = 0.3,
    segment.color = "grey50",
    show.legend  = FALSE
  ) +

  scale_color_manual(
    values = cat_colors,
    name   = NULL,
    guide  = guide_legend(override.aes = list(size = 3, alpha = 1))
  ) +
  scale_size_manual(values = cat_sizes, guide = "none") +
  scale_alpha_manual(values = cat_alpha, guide = "none") +

  scale_x_continuous(
  name = "M-value difference (logFC)",

  limits = c(-xlim_use, xlim_use),

  breaks = seq(-1, 1, by = 0.5),

  labels = scales::number_format(accuracy = 0.1)
 ) +
  scale_y_continuous(
    name   = expression(-log[10]~"(FDR-adjusted"~italic(p)~")"),
    limits = y_lim,
    expand = expansion(mult = c(0, 0.1))
  ) +

  annotate("text",
           x = x_lim[2] * 0.85,
           y = y_lim[2] * 0.97,
           label = paste0("Hyper: ",
                          sum(all_DMPs$Category == "Hypermethylated")),
           color = "#B2182B", size = 3.5, fontface = "bold", hjust = 1) +

  annotate("text",
           x = x_lim[1] * 0.85,
           y = y_lim[2] * 0.97,
           label = paste0("Hypo: ",
                          sum(all_DMPs$Category == "Hypomethylated")),
           color = "#2166AC", size = 3.5, fontface = "bold", hjust = 0) +

  theme_classic(base_size = 13) +
  theme(
    axis.title      = element_text(face = "bold", size = 13),
    axis.text       = element_text(color = "black", size = 11),
    axis.line       = element_line(color = "black", linewidth = 0.5),
    axis.ticks      = element_line(color = "black"),
    legend.position  = "top",
    legend.text      = element_text(size = 10),
    legend.key.size  = unit(0.4, "cm"),
    legend.spacing.x = unit(0.3, "cm"),
    panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10)
  )

ggsave(
  file.path(outDir, "Volcano_final.png"),
  p_volcano,
  width  = 6.69,
  height = 6.5,
  dpi    = 300
)

ggsave(
  file.path(outDir, "Volcano_final.pdf"),
  p_volcano,
  width  = 6.69,
  height = 6.5
)

cat("Volcano plot saved as PNG and PDF\n")

# =========================================================
# 5. HIERARCHICAL CLUSTERING DENDROGRAM
# =========================================================

sample_map <- c(
  "MD9/24"  = "S22", "MD23/24" = "S28", "MD31/24" = "S33",
  "MD33/24" = "S35", "MD35/24" = "S37", "MD25/24" = "S30",
  "MD32/24" = "S34", "MD34/24" = "S36", "MD38/24" = "S38",
  "MD42/24" = "S39", "MD1/24"  = "S21", "MD17/24" = "S25",
  "MD21/24" = "S27", "MD44/24" = "S40", "MD16/24" = "S24",
  "MD19/24" = "S26", "MD27/24" = "S31", "MD10/24" = "S23",
  "MD24/24" = "S29", "MD28/24" = "S32", "MC1/24"  = "S1",
  "MC8/24"  = "S2",  "MC11/24" = "S3", "MC13/24" = "S4",
  "MC14/24" = "S5",  "MC23/24" = "S6", "MC26/24" = "S7",
  "MC30/24" = "S8",  "MC32/24" = "S9", "MC34/24" = "S10",
  "MC41/24" = "S11", "MC42/24" = "S12", "MC43/24" = "S13",
  "MC47/24" = "S14", "MC48/24" = "S15", "MC50/24" = "S16",
  "MC53/24" = "S17", "MC56/24" = "S18", "MC60/24" = "S19"
)

targets$Plot_ID <- sample_map[targets$Sample_Name]

# top 50k variable probes 

top_var  <- order(rowVars(beta), decreasing = TRUE)[1:50000]
beta_use <- beta[top_var, ]
colnames(beta_use) <- targets$Plot_ID

cat("Probes used for clustering:", nrow(beta_use), "\n")
cat("Samples:", ncol(beta_use), "\n")

# CLUSTERING

cat("Computing distance matrix... (may take 1-2 min with all probes)\n")
dist_mat <- dist(t(beta_use), method = "euclidean")

cat("Clustering...\n")
hc <- hclust(dist_mat, method = "complete")

# Convert to dendrogram object
dend <- as.dendrogram(hc)

# COLORING
# S1-S19 = Control = blue
# S21-S40 = PUV = red

get_color <- function(label) {
  num <- as.numeric(sub("S", "", label))
  ifelse(num <= 20, "#2166AC", "#D6604D")
}

leaf_colors <- get_color(labels(dend))
labels_colors(dend) <- leaf_colors
labels_cex(dend)    <- 0.7

# PLOT — clean, publication style

png(file.path(outDir, "Dendrogram_top_probes_final.png"),
    width  = 2007,
    height = 1677,
    res    = 300)

par(mar = c(9, 5, 4, 2),
    bg  = "white")

plot(dend,
     ylab    = "Height (Euclidean distance)",
     main    = "",
     sub     = "",
     leaflab = "perpendicular",
     nodePar = list(pch = NA),
     axes    = TRUE)

legend("topright",
       legend   = c("Control", "PUV"),
       text.col = c("#2166AC", "#D6604D"),
       bty      = "n",
       cex      = 0.82,
       text.font = 2)

dev.off()

pdf(file.path(outDir, "Dendrogram_top_probes_final.pdf"),
    width  = 6.69,
    height = 8.5)

par(mar = c(8, 5, 4, 2))

plot(dend,
     ylab    = "Height (Euclidean distance)",
     main    = "",
     leaflab = "perpendicular",
     nodePar = list(pch = NA))

legend("topright",
       legend   = c("Control", "PUV"),
       text.col = c("#2166AC", "#D6604D"),
       bty      = "n",
       cex      = 0.82,
       text.font = 2)

dev.off()

cat("Final dendrogram saved as PNG and PDF\n")

# =========================================================
# 6. PCA PLOT
# =========================================================

# PCA ON BETA VALUES
# Top 20k variable probes, no scaling needed

probe_var  <- rowVars(beta)
top_probes <- order(probe_var, decreasing = TRUE)[1:20000]
beta_pca   <- beta[top_probes, ]

pca_beta <- prcomp(
  t(beta_pca),
  center = TRUE,
  scale. = FALSE    # beta already 0-1, scaling not needed
)

pct_beta <- round(100 * pca_beta$sdev^2 / sum(pca_beta$sdev^2), 1)

cat("=== Beta PCA variance explained ===\n")
for (i in 1:5) cat(sprintf("  PC%d: %.1f%%\n", i, pct_beta[i]))

# BUILD DATAFRAMES

pca_beta_df <- data.frame(
  PC1       = pca_beta$x[, 1],
  PC2       = pca_beta$x[, 2],
  Group     = targets$Group,
  Sample_ID = targets$Plot_ID,
  Age       = targets$Age
)

cat("\n=== PC1 correlation analysis (beta PCA) ===\n")

# Age correlation
r_age  <- cor(pca_beta_df$PC1, pca_beta_df$Age,
               use = "complete.obs")
cat(sprintf("PC1 vs Age:   r = %.3f (p = %.4f)\n",
            r_age,
            cor.test(pca_beta_df$PC1, pca_beta_df$Age)$p.value))

# Group correlation
r_grp  <- cor(pca_beta_df$PC1, as.numeric(pca_beta_df$Group),
               use = "complete.obs")
cat(sprintf("PC1 vs Group: r = %.3f (p = %.4f)\n",
            r_grp,
            cor.test(pca_beta_df$PC1,
                     as.numeric(pca_beta_df$Group))$p.value))

cat("\n=== PC2 correlation analysis (beta PCA) ===\n")
r_age2 <- cor(pca_beta_df$PC2, pca_beta_df$Age,
               use = "complete.obs")
cat(sprintf("PC2 vs Age:   r = %.3f\n", r_age2))
r_grp2 <- cor(pca_beta_df$PC2, as.numeric(pca_beta_df$Group),
               use = "complete.obs")
cat(sprintf("PC2 vs Group: r = %.3f\n", r_grp2))

# PLOT FUNCTION 

plot_pca <- function(df, pct) {

  ggplot(df, aes(x = PC1, y = PC2, color = Group)) +

    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey65") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey65") +

    geom_point(size = 4, alpha = 0.88, stroke = 0) +

    scale_color_manual(
      values = c("Control" = "#2166AC", "PUV" = "#D6604D"),
      name   = NULL,
      guide  = guide_legend(
        override.aes = list(size = 5, alpha = 1)
      )
    ) +

    labs(
      x        = paste0("PC1 (", pct[1], "%)"),
      y        = paste0("PC2 (", pct[2], "%)")
    ) +

    theme_bw(base_size = 13) +
    theme(
      axis.title       = element_text(face = "bold", size = 13),
      axis.text        = element_text(color = "black", size = 11),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black",
                                       fill = NA, linewidth = 0.6),
      legend.position  = "right",
      legend.text      = element_text(size = 10, face = "bold"),
      legend.key.size  = unit(0.55, "cm"),
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 7, color = "grey40"),
      plot.margin      = margin(10, 15, 10, 10)
    )
}

# PLOTS

p_beta <- plot_pca(
  df            = pca_beta_df,
  pct           = pct_beta
)

ggsave(file.path(outDir, "PCA_beta_values.png"),
       p_beta, width = 6.69, height = 8.5, dpi = 300)
ggsave(file.path(outDir, "PCA_beta_values.pdf"),
       p_beta, width = 6.69, height = 8.5)

cat("\nPCA plots saved\n")   

# ============================================================
# 7. VIOLIN PLOT — Distribution of methylation changes across genomic features
# ============================================================

plot_df <- sig_DMPs_anno

# PRIORITY-BASED ASSIGNMENT

plot_df$Genomic_Feature <- "Other"

plot_df$Genomic_Feature[
  grepl("TSS200",
         plot_df$UCSC_RefGene_Group)
] <- "TSS200"

plot_df$Genomic_Feature[
  grepl("TSS1500",
         plot_df$UCSC_RefGene_Group) &
  plot_df$Genomic_Feature == "Other"
] <- "TSS1500"

plot_df$Genomic_Feature[
  grepl("5'UTR|5UTR",
         plot_df$UCSC_RefGene_Group) &
  plot_df$Genomic_Feature == "Other"
] <- "5'UTR"

plot_df$Genomic_Feature[
  grepl("1stExon|exon_1",
         plot_df$UCSC_RefGene_Group) &
  plot_df$Genomic_Feature == "Other"
] <- "1st Exon"

plot_df$Genomic_Feature[
  grepl("exon_[2-9]|exon_[1-9][0-9]|Body",
         plot_df$UCSC_RefGene_Group) &
  plot_df$Genomic_Feature == "Other"
] <- "Body"

plot_df$Genomic_Feature[
  grepl("3'UTR|3UTR",
         plot_df$UCSC_RefGene_Group) &
  plot_df$Genomic_Feature == "Other"
] <- "3'UTR"

plot_df <- plot_df[
  plot_df$Genomic_Feature != "Other",
]

# ORDER FEATURES

plot_df$Genomic_Feature <- factor(
  plot_df$Genomic_Feature,
  levels = c(
    "TSS200",
    "TSS1500",
    "5'UTR",
    "1st Exon",
    "Body",
    "3'UTR"
  )
)

# PLOT

p_violin <- ggplot(
  plot_df,
  aes(
    x = Genomic_Feature,
    y = logFC,
    fill = Genomic_Feature
  )
) +

geom_violin(
  trim = TRUE,
  adjust = 0.8,
  color = "grey25",
  linewidth = 0.3
) +

coord_cartesian(ylim = c(-1, 1)) +

scale_y_continuous(
  breaks = c(-1, 0, 1)
) +

scale_fill_manual(
  values = c(
    "1st Exon" = "#B22222",   # dark red
    "3'UTR"    = "#808000",   # olive
    "5'UTR"    = "#006400",   # dark green
    "Body"     = "#008B8B",   # dark cyan
    "TSS1500"  = "#1E5AA8",   # deep blue
    "TSS200"   = "#8B3FBF"    # dark purple
  )
)+

theme_minimal(base_size = 11) +

theme(
  legend.position = "none",

  axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ),

  panel.grid.minor = element_blank()
) +

labs(
  x = NULL,
  y = "Log Fold Change (M-value difference)"
)

ggsave(
  file.path(outDir, "Violin_Genomic_Features.png"),
  p_violin,
  width = 6.69,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(outDir, "Violin_Genomic_Features.pdf"),
  p_violin,
  width = 6.69,
  height = 5.5
)

# ============================================================
# 8. PROPORTION OF HYPER/HYPO DMPs ACROSS CpG REGIONS
# ============================================================

plot_df <- sig_DMPs_anno

plot_df$Methylation_Status <- ifelse(
  plot_df$logFC > 0,
  "Hypermethylated",
  "Hypomethylated"
)

plot_df$CpG_Region <- plot_df$Relation_Primary

# Optional cleanup
plot_df$CpG_Region <- recode(
  plot_df$CpG_Region,

  "Island" = "Island",

  "N_Shore" = "N Shore",
  "S_Shore" = "S Shore",

  "N_Shelf" = "N Shelf",
  "S_Shelf" = "S Shelf",

  "OpenSea" = "Open Sea"
)

plot_df <- plot_df[!is.na(plot_df$CpG_Region), ]

# CALCULATE PROPORTIONS

prop_df <- plot_df %>%
  group_by(CpG_Region, Methylation_Status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CpG_Region) %>%
  mutate(Proportion = n / sum(n))

# ORDER CpG REGIONS

prop_df$CpG_Region <- factor(
  prop_df$CpG_Region,
  levels = c(
    "Island",
    "N Shore",
    "S Shore",
    "N Shelf",
    "S Shelf",
    "Open Sea"
  )
)

# PLOT

p_cpg <- ggplot(
  prop_df,
  aes(
    x = CpG_Region,
    y = Proportion,
    fill = Methylation_Status
  )
) +

geom_bar(
  stat = "identity",
  width = 0.7,
  color = "black",
  linewidth = 0.2
) +

scale_fill_manual(
  values = c(
    "Hypomethylated" = "#4C72B0",
    "Hypermethylated" = "#D62728"
  )
) +

scale_y_continuous(
  limits = c(0, 1),
  expand = c(0, 0)
) +

theme_minimal(base_size = 10) +

theme(
  axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ),
  legend.position = "bottom",
  legend.direction = "horizontal",

  legend.title = element_blank(),
  
  panel.grid.minor = element_blank()
) +

labs(
  x = NULL,
  y = "Proportion of DMPs"
)

ggsave(
  file.path(outDir, "CpG_region_proportions.png"),
  p_cpg,
  width = 3.35,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(outDir, "CpG_region_proportions.pdf"),
  p_cpg,
  width = 3.35,
  height = 5.5
)

# ============================================================
# 9. CHROMOSOMAL DISTRIBUTION OF DMPs
# ============================================================

plot_df <- sig_DMPs_anno

plot_df <- plot_df[
  !is.na(plot_df$chr) &
  !is.na(plot_df$pos),
]

plot_df$chr <- gsub("chr", "", plot_df$chr)

plot_df <- plot_df[
  plot_df$chr %in% c(as.character(1:22), "X", "Y"),
]

plot_df$chr <- factor(
  plot_df$chr,
  levels = c(as.character(1:22))
)

plot_df$Methylation_Status <- ifelse(
  plot_df$logFC > 0,
  "Hypermethylated",
  "Hypomethylated"
)

# CREATE CHROMOSOME LENGTH TABLE
# hg38 approximate chromosome sizes

chr_lengths <- data.frame(
  chr = factor(
    as.character(1:22),
    levels = rev(as.character(1:22))
  ),

  length = c(
    248956422, 242193529, 198295559, 190214555,
    181538259, 170805979, 159345973, 145138636,
    138394717, 133797422, 135086622, 133275309,
    114364328, 107043718, 101991189, 90338345,
    83257441, 80373285, 58617616, 64444167,
    46709983, 50818468
  )
)

base_df <- chr_lengths

# PLOT

p_chr <- ggplot() +

geom_segment(
  data = base_df,

  aes(
    x = 0,
    xend = length,
    y = chr,
    yend = chr
  ),

  linewidth = 3,
  color = "grey65"
) +

geom_point(
  data = plot_df,

  aes(
    x = pos,
    y = chr,
    color = Methylation_Status
  ),

  size = 0.7,
  alpha = 0.85
) +

scale_color_manual(
  values = c(
    "Hypomethylated" = "blue",
    "Hypermethylated" = "red"
  )
) +

scale_x_continuous(
  labels = scales::label_number(scale_cut = scales::cut_short_scale())
) +

theme_minimal(base_size = 11) +

theme(
  panel.grid.major.y = element_blank(),
  panel.grid.minor = element_blank(),

  legend.title = element_blank(),

  axis.text.y = element_text(face = "bold")
) +

labs(
  x = "Genomic Position",
  y = "Chromosome"
)

ggsave(
  file.path(outDir, "Chromosomal_distribution_DMPs.png"),
  p_chr,
  width = 6.69,
  height = 5.0,
  dpi = 300
)

ggsave(
  file.path(outDir, "Chromosomal_distribution_DMPs.pdf"),
  p_chr,
  width = 6.69,
  height = 5.0
)
