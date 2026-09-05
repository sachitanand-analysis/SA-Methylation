# ============================================================
# 1. GO PATHWAYS (Clusterprofiler)
# ============================================================

sig_genes <- sig_DMPs_anno %>%
  
  separate_rows(UCSC_RefGene_Name, sep = ";") %>%
  
  mutate(
    UCSC_RefGene_Name = trimws(UCSC_RefGene_Name)
  ) %>%
  
  filter(UCSC_RefGene_Name != "") %>%
  
  pull(UCSC_RefGene_Name) %>%
  
  unique()
  
# Build background gene universe from array-tested probes
bg_gene_symbols <- unique(unlist(strsplit(
  all_DMPs_anno$UCSC_RefGene_Name[
    !is.na(all_DMPs_anno$UCSC_RefGene_Name) &
    all_DMPs_anno$UCSC_RefGene_Name != ""
  ], ";"
)))
bg_gene_symbols <- trimws(bg_gene_symbols[bg_gene_symbols != ""])

bg_gene_df <- bitr(bg_gene_symbols, fromType="SYMBOL",
                    toType="ENTREZID", OrgDb=org.Hs.eg.db)
bg_entrez   <- unique(bg_gene_df$ENTREZID)

# Significant gene list
sig_gene_symbols <- unique(unlist(strsplit(
  sig_DMPs_anno$UCSC_RefGene_Name[
    !is.na(sig_DMPs_anno$UCSC_RefGene_Name) &
    sig_DMPs_anno$UCSC_RefGene_Name != ""
  ], ";"
)))
sig_gene_symbols <- trimws(sig_gene_symbols[sig_gene_symbols != ""])

sig_gene_df <- bitr(sig_gene_symbols, fromType="SYMBOL",
                     toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez      <- unique(sig_gene_df$ENTREZID)

cat("Sig genes mapped to Entrez:", length(entrez), "\n")
cat("Background genes mapped:", length(bg_entrez), "\n")

# ============================================================
# GO BIOLOGICAL PATHWAYS
# ============================================================

goe_bp <- enrichGO(
  gene          = entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

# ============================================================
# GO CELLULAR COMPONENTS 
# ============================================================

goe_cc <- enrichGO(
  gene          = entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

# ============================================================
# GO MOLECULAR FUNCTIONS
# ============================================================

goe_mf <- enrichGO(
  gene          = entrez,
  universe      = bg_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

# ------------------------------------------------------------
# GO REDUNDANCY REDUCTION
# ------------------------------------------------------------

if (!is.null(goe_bp) && nrow(as.data.frame(goe_bp)) > 0) {
  goe_bp <- clusterProfiler::simplify(
    goe_bp, cutoff = 0.7, by = "p.adjust", select_fun = min
  )
}
if (!is.null(goe_cc) && nrow(as.data.frame(goe_cc)) > 0) {
  goe_cc <- clusterProfiler::simplify(
    goe_cc, cutoff = 0.7, by = "p.adjust", select_fun = min
  )
}
if (!is.null(goe_mf) && nrow(as.data.frame(goe_mf)) > 0) {
  goe_mf <- clusterProfiler::simplify(
    goe_mf, cutoff = 0.7, by = "p.adjust", select_fun = min
  )
}

cat("\n=== Post-simplify() GO term counts (Section 14, total DMP set) ===\n")
cat("GO BP:", nrow(as.data.frame(goe_bp)), "\n")
cat("GO CC:", nrow(as.data.frame(goe_cc)), "\n")
cat("GO MF:", nrow(as.data.frame(goe_mf)), "\n")

# ============================================================
# SAVE
# ============================================================

write.csv(
  as.data.frame(goe_bp),
  file.path(
    outDir,
  "GO_Biological_Process.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(goe_cc),
  file.path(
    outDir,
  "GO_Cellular_Component.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(goe_mf),
  file.path(
    outDir,
  "GO_Molecular_Function.csv"),
  row.names = FALSE
)

while (!is.null(dev.list())) dev.off()

# ============================================================
# PLOT
# ============================================================

p_bp <- dotplot(
  goe_bp,
  showCategory = 20,
  title = "GO Biological Process Enrichment"
)

print(p_bp)

ggsave(file.path(outDir,
  "GO_BP_Dotplot.png"),
  p_bp,
  width = 6.69,
  height = 8.0,
  dpi = 300
)

ggsave(file.path(outDir,
  "GO_BP_Dotplot.pdf"),
  p_bp,
  width = 6.69,
  height = 8.0
)

p_cc <- dotplot(
  goe_cc,
  showCategory = 20,
  title = "GO Cellular Component Enrichment"
)

print(p_cc)

ggsave(file.path(outDir,
  "GO_CC_Dotplot.png"),
  p_cc,
  width = 6.69,
  height = 8.0,
  dpi = 300
)

ggsave(file.path(outDir,
  "GO_CC_Dotplot.pdf"),
  p_cc,
  width = 6.69,
  height = 8.0
)

p_mf <- dotplot(
  goe_mf,
  showCategory = 20,
  title = "GO Molecular Function Enrichment"
)

print(p_mf)

ggsave(file.path(outDir,
  "GO_MF_Dotplot.png"),
  p_mf,
  width = 6.69,
  height = 8.0,
  dpi = 300
)

ggsave(file.path(outDir,
  "GO_MF_Dotplot.pdf"),
  p_mf,
  width = 6.69,
  height = 8.0
)

# ============================================================
# 2. KEGG (Clusterprofiler)
# ============================================================

cat("Extracting genes from significant DMPs...\n")

# All significant DMP genes (background for enrichment)
extract_genes <- function(gene_col) {
  genes <- unlist(strsplit(gene_col, ";"))
  genes <- unique(trimws(genes))
  genes <- genes[genes != "" & !is.na(genes)]
  return(genes)
}

# Significant DMP genes
sig_genes_list <- lapply(
  sig_DMPs_anno$UCSC_RefGene_Name,
  extract_genes
)
sig_genes <- unique(unlist(sig_genes_list))
sig_genes  <- sig_genes[sig_genes != ""]

cat("Unique genes in significant DMPs:", length(sig_genes), "\n")

# Background — all tested DMP genes
bg_genes_list <- lapply(
  all_DMPs_anno$UCSC_RefGene_Name,
  extract_genes
)
bg_genes <- unique(unlist(bg_genes_list))
bg_genes  <- bg_genes[bg_genes != ""]

cat("Unique genes in background:", length(bg_genes), "\n")

# Separate hyper and hypo genes
hyper_anno <- subset(sig_DMPs_anno, logFC > 0)
hypo_anno  <- subset(sig_DMPs_anno, logFC < 0)

hyper_genes <- unique(unlist(lapply(
  hyper_anno$UCSC_RefGene_Name, extract_genes)))
hypo_genes  <- unique(unlist(lapply(
  hypo_anno$UCSC_RefGene_Name,  extract_genes)))

cat("Hypermethylated genes:", length(hyper_genes), "\n")
cat("Hypomethylated genes:", length(hypo_genes), "\n")

# ============================================================
# CONVERT GENE SYMBOLS TO ENTREZ IDs
# clusterProfiler requires Entrez IDs for KEGG
# ============================================================

convert_to_entrez <- function(gene_symbols) {
  mapped <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  cat("Mapped:", nrow(mapped), "of", length(gene_symbols),
      "genes to Entrez IDs\n")
  return(mapped$ENTREZID)
}

cat("\nConverting significant genes...\n")
sig_entrez  <- convert_to_entrez(sig_genes)

cat("Converting background genes...\n")
bg_entrez   <- convert_to_entrez(bg_genes)

cat("Converting hyper genes...\n")
hyper_entrez <- convert_to_entrez(hyper_genes)

cat("Converting hypo genes...\n")
hypo_entrez  <- convert_to_entrez(hypo_genes)

# ============================================================
# FULL KEGG ENRICHMENT (all pathways)
# ============================================================

cat("\nRunning KEGG enrichment on all pathways...\n")

kegg_all <- enrichKEGG(
  gene          = sig_entrez,
  universe      = bg_entrez,
  organism      = "hsa",
  keyType       = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,    
  minGSSize     = 5,
  maxGSSize     = 2000
)

kegg_df <- as.data.frame(kegg_all)
cat("Total significant KEGG pathways:", nrow(kegg_df), "\n")

if (nrow(kegg_df) == 0) {
  cat("No significant pathways at p.adj<0.05\n")
}

write.csv(kegg_df,
          file.path(outDir, "KEGG_sig_pathways.csv"),
          row.names = FALSE)

# ============================================================
# GET KEGG PATHWAY CATEGORIES
# ============================================================

cat("\nFetching KEGG pathway categories...\n")

# Get all human KEGG pathways with their categories
tryCatch({
  pathway_list <- keggList("pathway", "hsa")
  pathway_df   <- data.frame(
    PathwayID   = sub("path:", "", names(pathway_list)),
    Description = as.character(pathway_list),
    stringsAsFactors = FALSE
  )
  cat("Total KEGG human pathways fetched:", nrow(pathway_df), "\n")
}, error = function(e) {
  cat("KEGG API not accessible. Using manual category lists.\n")
})

write.csv(pathway_df,
          file.path(outDir, "KEGG_human_pathways_catalog.csv"),
          row.names = FALSE)
          
# ============================================================
# PLOT TOP 20 KEGG PATHWAYS
# ============================================================

# Check if pathways exist

if (nrow(kegg_df) > 0) {

  kegg_top20 <- kegg_df %>%
    
    arrange(p.adjust) %>%
    
    slice_head(n = 20)

  kegg_top20$neglog10FDR <- -log10(kegg_top20$p.adjust)
  
  kegg_top20$Description <- str_wrap(
  kegg_top20$Description,
  width = 35
)

  kegg_top20$Description <- factor(
    kegg_top20$Description,
    levels = rev(kegg_top20$Description)
  )

  p_kegg <- ggplot(
    kegg_top20,
    
    aes(
      x     = neglog10FDR,
      y     = Description,
      size  = Count,
      color = neglog10FDR
    )
  ) +

    geom_point(alpha = 0.9) +

    scale_color_gradient(
      low  = "#56B1F7",
      high = "#D73027"
    ) +

    geom_vline(
      xintercept = -log10(0.05),
      linetype   = "dashed",
      color      = "grey40",
      linewidth  = 0.5
    ) +

    annotate(
      "text",
      x = -log10(0.05) + 0.1,
      y = 1,
      label = "FDR = 0.1",
      size = 3,
      color = "grey40",
      hjust = 0
    ) +

    theme_bw(base_size = 13) +

    labs(
      
      x = expression(-log[10]("Adjusted p-value")),
      
      y = "KEGG Pathway",
      
      size = "Gene Count",
      
      color = expression(-log[10]("FDR"))
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 11
      ),

      axis.text.y = element_text(size = 12)
    )

  ggsave(
    file.path(outDir,
    "KEGG_Top20_Pathways.png"),
    
    p_kegg,
    
    width = 6.69,
    height = 8.0,
    dpi = 300
  )

  ggsave(
    file.path(outDir,
    "KEGG_Top20_Pathways.pdf"),
    
    p_kegg,
    
    width = 6.69,
    height = 8.0
  )

  cat("Top 20 KEGG pathway plot saved\n")

} else {

  cat("No KEGG pathways available for plotting\n")
}          

cat("\n=== KEGG analysis complete ===\n")     

# ============================================================
# 3. REACTOME — ReactomePA 
# ============================================================

cat("\nRunning Reactome enrichment (ReactomePA)...\n")

reactome_res <- enrichPathway(
  gene          = entrez,
  universe      = bg_entrez,   # array background
  organism      = "human",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

reactome_df <- as.data.frame(reactome_res)
cat("Reactome pathways significant (FDR<0.05):",
    nrow(reactome_df), "\n")

write.csv(reactome_df,
          file.path(outDir,"Reactome_Enrichment.csv"),
          row.names=FALSE)

if (nrow(reactome_df) > 0) {
  if (length(dev.list()) > 0) graphics.off()
  p_reactome <- dotplot(reactome_res, showCategory=20,
                         title="Reactome Pathway Enrichment")
  ggsave(file.path(outDir,"Reactome_Dotplot.png"),
         p_reactome, width=6.69, height=8.0, dpi=300)
  ggsave(file.path(outDir,"Reactome_Dotplot.pdf"),
         p_reactome, width=6.69, height=8.0)
  cat("Reactome dotplot saved\n")
} else {
  cat("No significant Reactome pathways at FDR<0.05\n")
}
