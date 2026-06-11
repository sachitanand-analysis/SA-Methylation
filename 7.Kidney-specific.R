# ============================================================
# 1. KIDNEY-SPECIFIC DMPs 
# ============================================================

kidney_genes <- scan(
  "kidney_genes.txt", what=character()
)

cat("\n=== Kidney-specific DMP analysis ===\n")

kidney_hits <- sig_DMPs_anno %>%
  separate_rows(UCSC_RefGene_Name, sep=";") %>%
  mutate(UCSC_RefGene_Name = trimws(UCSC_RefGene_Name)) %>%
  filter(UCSC_RefGene_Name %in% kidney_genes)

# CpG-level object (use for ALL plots and summaries)
kidney_hits_cpg <- kidney_hits %>%
  distinct(CpG, .keep_all = TRUE)
  
kidney_hits_cpg$Direction <- ifelse(kidney_hits_cpg$logFC > 0,
                                 "Hypermethylated",
                                 "Hypomethylated")

# Top kidney DMPs
top_kidney <- kidney_hits_cpg %>%
  arrange(adj.P.Val) %>%
  distinct(UCSC_RefGene_Name, .keep_all=TRUE) %>%
  slice_head(n=100) %>%
  select(UCSC_RefGene_Name, Name, chr, pos, logFC,
         adj.P.Val, Direction, Relation_to_Island,
         Relation_Primary, UCSC_RefGene_Group)

write.csv(top_kidney,
          file.path(outDir,"Top_Kidney_DMPs.csv"),
          row.names=FALSE)

all_unique_kidney <- kidney_hits_cpg %>%
  arrange(adj.P.Val) %>%
  distinct(UCSC_RefGene_Name, .keep_all=TRUE)

write.csv(all_unique_kidney,
          file.path(outDir,"All_Unique_Kidney_Genes.csv"),
          row.names=FALSE)

write.csv(kidney_hits_cpg,
          file.path(outDir,"All_Kidney_DMPs.csv"),
          row.names=FALSE)

# ============================================================
# 2. KIDNEY-SPECIFIC HYPER / HYPO METHYLATED GENES
# ============================================================

# Hyper-methylated kidney genes

kidney_hyper <- kidney_hits_cpg %>%
  
  filter(Direction == "Hypermethylated") %>%
  
  arrange(adj.P.Val) %>%
  
  distinct(UCSC_RefGene_Name, .keep_all = TRUE)

# Hypo-methylated kidney genes

kidney_hypo <- kidney_hits_cpg %>%
  
  filter(Direction == "Hypomethylated") %>%
  
  arrange(adj.P.Val) %>%
  
  distinct(UCSC_RefGene_Name, .keep_all = TRUE)

# SAVE

write.csv(
  kidney_hyper,
  file.path(outDir,
  "Kidney_Genes_Hypermethylated.csv"),
  row.names = FALSE
)

write.csv(
  kidney_hypo,
  file.path(outDir,
  "Kidney_Genes_Hypomethylated.csv"),
  row.names = FALSE
)
# ============================================================
# 3. KIDNEY SUMMARY STATISTICS 
# ============================================================

cat("Total kidney-associated DMPs:", nrow(kidney_hits_cpg), "\n")
cat("Hypermethylated kidney DMPs:",
    sum(kidney_hits_cpg$Direction=="Hypermethylated"), "\n")
cat("Hypomethylated kidney DMPs:",
    sum(kidney_hits_cpg$Direction=="Hypomethylated"), "\n")

cat("\nCpG island context of kidney DMPs:\n")
print(table(kidney_hits_cpg$Relation_Primary))

cat("\nKidney DMPs by gene feature group:\n")
print(table(kidney_hits_cpg$UCSC_RefGene_Group))

# Direction by context
cat("\nKidney DMPs — direction by context:\n")
cat("\n===== Kidney-specific gene summary =====\n")

cat("Unique kidney genes:",
    length(unique(kidney_hits_cpg$UCSC_RefGene_Name)), "\n")

cat("Hypermethylated kidney genes:",
    nrow(kidney_hyper), "\n")

cat("Hypomethylated kidney genes:",
    nrow(kidney_hypo), "\n")
 
table(kidney_hits_cpg$Direction)

# Save kidney summaries
kidney_summary <- data.frame(
  Metric  = c("Total kidney DMPs",
               "Unique genes",
               "Hypermethylated",
               "Hypomethylated"),
  Count   = c(nrow(kidney_hits_cpg),
               length(unique(kidney_hits_cpg$UCSC_RefGene_Name)),
               sum(kidney_hits_cpg$Direction=="Hypermethylated"),
               sum(kidney_hits_cpg$Direction=="Hypomethylated"))
)
write.csv(kidney_summary,
          file.path(outDir,"Kidney_DMP_Summary.csv"),
          row.names=FALSE)

# ============================================================
# 4. KIDNEY HEATMAPS
# ============================================================

targets$Array_ID <- basename(targets$Basename)

valid_cpgs  <- intersect(top_kidney$Name, rownames(beta))
heatmap_mat <- beta[valid_cpgs, ]
selected2   <- top_kidney[match(valid_cpgs, top_kidney$Name), ]
row_labels  <- paste0(selected2$UCSC_RefGene_Name, "_", selected2$Name)
rownames(heatmap_mat) <- row_labels

sample_order <- c(
  targets$Array_ID[targets$Group == "Control"],
  targets$Array_ID[targets$Group == "PUV"]
)
heatmap_mat   <- heatmap_mat[, sample_order]
annotation_col <- data.frame(Group = targets$Group[
  match(sample_order, targets$Array_ID)])
rownames(annotation_col) <- sample_order

ann_colors <- list(Group = c(Control="#4575B4", PUV="#D73027"))

if (length(dev.list()) > 0) graphics.off()
png(file.path(outDir,"Kidney_Top_Heatmap.png"),
    width=2007, height=2550, res=300)
tryCatch({
  pheatmap(heatmap_mat, scale="row",
           cluster_cols=FALSE, cluster_rows=TRUE,
           clustering_method="ward.D2",
           annotation_col=annotation_col,
           annotation_colors=ann_colors,
           show_rownames=TRUE, show_colnames=FALSE,
           fontsize_row=8, border_color=NA
           main="Top Kidney DMPs")
}, finally=dev.off())

pdf(file.path(outDir,"Kidney_Top_Heatmap.pdf"),
    width=6.69, height=8.5)
tryCatch({
  pheatmap(heatmap_mat, scale="row",
           cluster_cols=FALSE, cluster_rows=TRUE,
           clustering_method="ward.D2",
           annotation_col=annotation_col,
           annotation_colors=ann_colors,
           show_rownames=TRUE, show_colnames=FALSE,
           fontsize_row=8, border_color=NA,
           main="Top Kidney DMPs")
}, finally=dev.off())

cat("Kidney heatmaps saved\n")

# ============================================================
# 5. Kidney Specific VOLCANO 
# ============================================================

kidney_hits_cpg$Category <- "Not significant"

kidney_hits_cpg$Category[
  kidney_hits_cpg$adj.P.Val < 0.1 &
  kidney_hits_cpg$logFC >  0.2
] <- "Hypermethylated"

kidney_hits_cpg$Category[
  kidney_hits_cpg$adj.P.Val < 0.1 &
  kidney_hits_cpg$logFC < -0.2
] <- "Hypomethylated"

cat("=== Volcano category counts ===\n")
print(table(kidney_hits_cpg$Category))

ann      <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
ann_sub  <- ann[, c("Name", "UCSC_RefGene_Name", "chr", "pos")]

ann_sub$Gene <- sapply(
  strsplit(ann_sub$UCSC_RefGene_Name, ";"),
  function(x) if (length(x) == 0 || x[1] == "") NA else x[1]
)

kidney_hits_cpg$Gene <- ann_sub$Gene[match(kidney_hits_cpg$CpG, ann_sub$Name)]
kidney_hits_cpg$chr  <- ann_sub$chr[match(kidney_hits_cpg$CpG, ann_sub$Name)]


label_hyper <- subset(kidney_hits_cpg,
                       Category == "Hypermethylated" &
                       !is.na(Gene) & Gene != "")
label_hyper <- label_hyper[
  order(label_hyper$adj.P.Val),
]

label_hyper <- label_hyper[
  !duplicated(label_hyper$Gene),
]

label_hyper <- label_hyper[
  1:min(8, nrow(label_hyper)),
]

label_hypo  <- subset(kidney_hits_cpg,
                       Category == "Hypomethylated" &
                       !is.na(Gene) & Gene != "")
label_hypo <- label_hypo[
  order(label_hypo$adj.P.Val),
]

label_hypo <- label_hypo[
  !duplicated(label_hypo$Gene),
]

label_hypo <- label_hypo[
  1:min(8, nrow(label_hypo)),
]

label_df    <- rbind(label_hyper, label_hypo)

cat("\nProbes to be labelled:\n")
print(label_df[, c("CpG", "Gene", "delta_beta", "logFC", "adj.P.Val", "Category")])

x_max <- max(abs(kidney_hits_cpg$logFC), na.rm = TRUE)
x_lim <- c(-ceiling(x_max * 10) / 10 - 0.2,
             ceiling(x_max * 10) / 10 + 0.2)

y_max <- max(-log10(kidney_hits_cpg$adj.P.Val[
  is.finite(-log10(kidney_hits_cpg$adj.P.Val))
]), na.rm = TRUE)
y_lim <- c(0, ceiling(y_max) + 0.3)

cat("\nx-axis range:", round(x_lim, 3), "\n")
cat("y-axis range:", round(y_lim, 3), "\n")

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

kidney_hits_cpg$plot_order <- ifelse(kidney_hits_cpg$Category == "Not significant", 1, 2)
plot_data_ordered   <- kidney_hits_cpg[order(kidney_hits_cpg$plot_order), ]

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
                          sum(kidney_hits_cpg$Category == "Hypermethylated")),
           color = "#B2182B", size = 3.5, fontface = "bold", hjust = 1) +

  annotate("text",
           x = x_lim[1] * 0.85,
           y = y_lim[2] * 0.97,
           label = paste0("Hypo: ",
                          sum(kidney_hits_cpg$Category == "Hypomethylated")),
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
  file.path(outDir, "Volcano_final_kidney.png"),
  p_volcano,
  width  = 6.69,
  height = 6.5,
  dpi    = 300
)

ggsave(
  file.path(outDir, "Volcano_final_kidney.pdf"),
  p_volcano,
  width  = 6.69,
  height = 6.5
)

cat("Volcano plot saved as PNG and PDF\n")

# ============================================================
# 6. VIOLIN PLOT
# ============================================================

plot_df <- kidney_hits_cpg

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
  title = "Kidney-specific DMP Distribution",

  x = NULL,

  y = "Log Fold Change (M-value difference)"
)

ggsave(
  file.path(outDir, "Violin_Genomic_Features_Kidney.png"),
  p_violin,
  width = 6.69,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(outDir, "Violin_Genomic_Features_Kidney.pdf"),
  p_violin,
  width = 6.69,
  height = 5.5
)

# ============================================================
# 7. PROPORTION OF HYPER/HYPO DMPs ACROSS CpG REGIONS
# ============================================================

plot_df <- kidney_hits_cpg

plot_df$Methylation_Status <- ifelse(
  plot_df$logFC > 0,
  "Hypermethylated",
  "Hypomethylated"
)

plot_df$CpG_Region <- plot_df$Relation_Primary

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

prop_df <- plot_df %>%
  group_by(CpG_Region, Methylation_Status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CpG_Region) %>%
  mutate(Proportion = n / sum(n))

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
  
 plot.title = element_text(
  size = 8,
  face = "bold",
  hjust = 0.5
 ),

  panel.grid.minor = element_blank()
) +

labs(
  title = "Kidney-specific CpG region distribution",
  x = NULL,

  y = "Proportion of DMPs"
)

ggsave(
  file.path(outDir, "CpG_region_proportions_Kidney.png"),
  p_cpg,
  width = 3.35,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(outDir, "CpG_region_proportions_Kidney.pdf"),
  p_cpg,
  width = 3.35,
  height = 5.5
) 
# ============================================================
# 8. CHROMOSOMAL DISTRIBUTION OF KIDNEY DMPs
# ============================================================

plot_df <- kidney_hits_cpg

plot_df <- plot_df[
  !is.na(plot_df$chr),
]

plot_df$chr <- gsub("chr", "", plot_df$chr)

plot_df <- plot_df[
  plot_df$chr %in% as.character(1:22),
]

plot_df$Direction <- ifelse(
  plot_df$logFC > 0,
  "Hypermethylated",
  "Hypomethylated"
)

chr_counts <- plot_df %>%

  group_by(chr, Direction) %>%

  summarise(
    Count = n(),
    .groups = "drop"
  )

chr_counts$chr <- factor(
  chr_counts$chr,
  levels = as.character(1:22)
)

p_chr <- ggplot(

  chr_counts,

  aes(
    x = chr,
    y = Count,
    fill = Direction
  )
) +

geom_col(
  position = "stack",
  width = 0.75,
  color = "black",
  linewidth = 0.2
) +

scale_fill_manual(
  values = c(
    "Hypomethylated" = "#4C72B0",
    "Hypermethylated" = "#D62728"
  )
) +

theme_bw(base_size = 11) +

labs(
  title = "Chromosomal Distribution of Kidney-specific DMPs",
  x = "Chromosome",

  y = "Number of DMPs"
) +

theme(

  plot.title = element_text(
    face = "bold",
    size = 11,
    hjust = 0.5
  ),

  axis.title = element_text(
    face = "bold",
    size = 10
  ),

  axis.text = element_text(
    size = 9,
    color = "black"
  ),

  legend.title = element_blank(),

  legend.text = element_text(size = 9),

  panel.grid.minor = element_blank()
)

ggsave(
  file.path(
    outDir,
    "Kidney_DMP_Chromosome_Distribution.png"
  ),

  p_chr,

  width = 6.69,
  height = 3.5,
  dpi = 300
)

ggsave(
  file.path(
    outDir,
    "Kidney_DMP_Chromosome_Distribution.pdf"
  ),

  p_chr,

  width = 6.69,
  height = 3.5
)

# ============================================================
# 9. KIDNEY GO ENRICHMENT
# Universe = all significant DMP genes 
# ============================================================

# Rebuild universe from significant DMPs only
sig_gene_symbols_all <- unique(unlist(strsplit(
  sig_DMPs_anno$UCSC_RefGene_Name[
    !is.na(sig_DMPs_anno$UCSC_RefGene_Name) &
    sig_DMPs_anno$UCSC_RefGene_Name != ""
  ], ";"
)))
sig_gene_symbols_all <- trimws(
  sig_gene_symbols_all[sig_gene_symbols_all != ""]
)

sig_universe_df <- bitr(
  sig_gene_symbols_all,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
sig_universe_entrez <- unique(sig_universe_df$ENTREZID)

cat("Universe for kidney GO (all sig DMP genes):",
    length(sig_universe_entrez), "\n")
    
kidney_sig_genes <- unique(
  kidney_hits_cpg$UCSC_RefGene_Name
)

kidney_sig_genes <- kidney_sig_genes[
  kidney_sig_genes != ""
]

kidney_gene_df <- bitr(
  kidney_sig_genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

kidney_entrez <- unique(
  kidney_gene_df$ENTREZID
)

cat("Kidney genes mapped to Entrez:",
    length(kidney_entrez), "\n")    
cat("Kidney genes being tested:", length(kidney_entrez), "\n")

kidney_bp <- enrichGO(
  gene          = kidney_entrez,
  universe      = sig_universe_entrez,  
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

kidney_cc <- enrichGO(
  gene          = kidney_entrez,
  universe      = sig_universe_entrez,  
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

kidney_mf <- enrichGO(
  gene          = kidney_entrez,
  universe      = sig_universe_entrez,  
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

cat("Kidney BP terms :",
    nrow(as.data.frame(kidney_bp)), "\n")
cat("Kidney CC terms :",
    nrow(as.data.frame(kidney_cc)), "\n")
cat("Kidney MF terms :",
    nrow(as.data.frame(kidney_mf)), "\n")

write.csv(as.data.frame(kidney_bp),
          file.path(outDir, "Kidney_GO_BP.csv"),
          row.names=FALSE)
write.csv(as.data.frame(kidney_cc),
          file.path(outDir, "Kidney_GO_CC.csv"),
          row.names=FALSE)
write.csv(as.data.frame(kidney_mf),
          file.path(outDir, "Kidney_GO_MF.csv"),
          row.names=FALSE)
          
# ============================================================
# GO PLOTS
# ============================================================

plot_go <- function(go_obj, ontology_name, filename_prefix) {

  if (is.null(go_obj) || nrow(as.data.frame(go_obj)) == 0) {

    cat("No significant", ontology_name, "terms\n")
    return(NULL)
  }

  p <- dotplot(
    go_obj,
    showCategory = 20,
    font.size = 10,
    title = paste0(
      "Kidney-specific GO ",
      ontology_name,
      " Enrichment"
    )
  ) +

  theme_bw(base_size = 11) +

  theme(

    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),

    axis.title = element_text(
      face = "bold",
      size = 10
    ),

    axis.text.x = element_text(
      size = 9,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),

    legend.title = element_text(
      face = "bold",
      size = 9
    ),

    legend.text = element_text(
      size = 8
    ),

    panel.grid.minor = element_blank()
  )

  ggsave(
    file.path(
      outDir,
      paste0(filename_prefix, "_Dotplot.png")
    ),

    p,

    width = 6.69,
    height = 8.0,
    dpi = 300
  )

  ggsave(
    file.path(
      outDir,
      paste0(filename_prefix, "_Dotplot.pdf")
    ),

    p,

    width = 6.69,
    height = 8.0
  )

  cat(ontology_name, "plot saved\n")
}

plot_go(
  kidney_bp,
  "BP",
  "Kidney_GO_BP"
)

plot_go(
  kidney_cc,
  "CC",
  "Kidney_GO_CC"
)

plot_go(
  kidney_mf,
  "MF",
  "Kidney_GO_MF"
)

