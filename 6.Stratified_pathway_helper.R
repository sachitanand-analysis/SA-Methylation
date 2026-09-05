# ============================================================
# 1. STRATIFIED PATHWAY ANALYSIS
#     Hyper vs Hypo vs Promoter vs Non-Promoter loci
#     (GO Biological Process, KEGG, Reactome)
# ============================================================

# ------------------------------------------------------------
# REUSABLE ENRICHMENT FUNCTION
# ------------------------------------------------------------

run_stratified_enrichment <- function(gene_symbols, bg_entrez, label, outDir) {

  gene_symbols <- unique(trimws(gene_symbols))
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]

  cat("\n=====", label, "=====\n")
  cat("Input gene symbols:", length(gene_symbols), "\n")

  if (length(gene_symbols) < 3) {
    cat("Fewer than 3 genes — skipping enrichment for", label, "\n")
    return(invisible(NULL))
  }

  gene_map <- tryCatch(
    bitr(gene_symbols, fromType = "SYMBOL",
         toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )

  if (is.null(gene_map) || nrow(gene_map) < 3) {
    cat("Fewer than 3 genes mapped to Entrez — skipping", label, "\n")
    return(invisible(NULL))
  }

  entrez_ids <- unique(gene_map$ENTREZID)
  cat("Mapped to Entrez:", length(entrez_ids), "\n")

  # ---- GO Biological Process ----
  goe_bp <- tryCatch(
    enrichGO(gene = entrez_ids, universe = bg_entrez,
             OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = "BP", pAdjustMethod = "BH",
             qvalueCutoff = 0.1, readable = TRUE),
    error = function(e) { cat("GO BP failed:", conditionMessage(e), "\n"); NULL }
  )

  # GO redundancy reduction — collapse DAG parent/child/sibling
  # terms (semantic similarity >= 0.7) to the lowest-p.adjust
  # representative, so small gene sets don't report inflated term
  # counts purely from GO hierarchy structure.
  if (!is.null(goe_bp) && nrow(as.data.frame(goe_bp)) > 0) {
    goe_bp <- tryCatch(
      clusterProfiler::simplify(goe_bp, cutoff = 0.7,
                                 by = "p.adjust", select_fun = min),
      error = function(e) { cat("GO BP simplify() failed, using raw result:",
                                 conditionMessage(e), "\n"); goe_bp }
    )
  }

  if (!is.null(goe_bp) && nrow(as.data.frame(goe_bp)) > 0) {
    write.csv(as.data.frame(goe_bp),
              file.path(outDir, paste0("GO_BP_", label, ".csv")),
              row.names = FALSE)
    p_bp <- dotplot(goe_bp, showCategory = 20,
                     title = paste0("GO Biological Process — ", label))
    ggsave(file.path(outDir, paste0("GO_BP_Dotplot_", label, ".png")),
           p_bp, width = 6.69, height = 8.0, dpi = 600)
    ggsave(file.path(outDir, paste0("GO_BP_Dotplot_", label, ".pdf")),
           p_bp, width = 6.69, height = 8.0)
    cat("GO BP:", nrow(as.data.frame(goe_bp)), "significant terms saved\n")
  } else {
    cat("No significant GO BP terms for", label, "\n")
  }
  
  # ---- GO Cellular Component ----
  goe_cc <- tryCatch(
    enrichGO(gene = entrez_ids, universe = bg_entrez,
             OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = "CC", pAdjustMethod = "BH",
             qvalueCutoff = 0.1, readable = TRUE),
    error = function(e) { cat("GO CC failed:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(goe_cc) && nrow(as.data.frame(goe_cc)) > 0) {
    goe_cc <- tryCatch(
      clusterProfiler::simplify(goe_cc, cutoff = 0.7,
                                 by = "p.adjust", select_fun = min),
      error = function(e) { cat("GO CC simplify() failed, using raw result:",
                                 conditionMessage(e), "\n"); goe_cc }
    )
  }

  if (!is.null(goe_cc) && nrow(as.data.frame(goe_cc)) > 0) {
    write.csv(as.data.frame(goe_cc),
              file.path(outDir, paste0("GO_CC_", label, ".csv")),
              row.names = FALSE)
    p_cc <- dotplot(goe_cc, showCategory = 20,
                     title = paste0("GO Cellular Component — ", label))
    ggsave(file.path(outDir, paste0("GO_CC_Dotplot_", label, ".png")),
           p_cc, width = 6.69, height = 8.0, dpi = 300)
    ggsave(file.path(outDir, paste0("GO_CC_Dotplot_", label, ".pdf")),
           p_cc, width = 6.69, height = 8.0)
    cat("GO CC:", nrow(as.data.frame(goe_cc)), "significant terms saved\n")
  } else {
    cat("No significant GO CC terms for", label, "\n")
  }
  
  # ---- GO Molecular Functions ----
  goe_mf <- tryCatch(
    enrichGO(gene = entrez_ids, universe = bg_entrez,
             OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = "MF", pAdjustMethod = "BH",
             qvalueCutoff = 0.1, readable = TRUE),
    error = function(e) { cat("GO MF failed:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(goe_mf) && nrow(as.data.frame(goe_mf)) > 0) {
    goe_mf <- tryCatch(
      clusterProfiler::simplify(goe_mf, cutoff = 0.7,
                                 by = "p.adjust", select_fun = min),
      error = function(e) { cat("GO MF simplify() failed, using raw result:",
                                 conditionMessage(e), "\n"); goe_mf }
    )
  }

  if (!is.null(goe_mf) && nrow(as.data.frame(goe_mf)) > 0) {
    write.csv(as.data.frame(goe_mf),
              file.path(outDir, paste0("GO_MF_", label, ".csv")),
              row.names = FALSE)
    p_mf <- dotplot(goe_mf, showCategory = 20,
                     title = paste0("GO Molecular Functions — ", label))
    ggsave(file.path(outDir, paste0("GO_MF_Dotplot_", label, ".png")),
           p_mf, width = 6.69, height = 8.0, dpi = 600)
    ggsave(file.path(outDir, paste0("GO_MF_Dotplot_", label, ".pdf")),
           p_mf, width = 6.69, height = 8.0)
    cat("GO MF:", nrow(as.data.frame(goe_mf)), "significant terms saved\n")
  } else {
    cat("No significant GO MF terms for", label, "\n")
  }
  

  # ---- KEGG ----
  kegg_res <- tryCatch(
    enrichKEGG(gene = entrez_ids, universe = bg_entrez,
               organism = "hsa", keyType = "ncbi-geneid",
               pAdjustMethod = "BH", pvalueCutoff = 0.05,
               qvalueCutoff = 0.1, minGSSize = 5, maxGSSize = 2000),
    error = function(e) { cat("KEGG failed:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    kegg_df <- as.data.frame(kegg_res)
    write.csv(kegg_df,
              file.path(outDir, paste0("KEGG_", label, ".csv")),
              row.names = FALSE)
    p_kegg <- dotplot(kegg_res, showCategory = 20,
      title = paste0("KEGG Pathway Enrichment — ", label)
  )
    ggsave(file.path(outDir, paste0("KEGG_Dotplot_", label, ".png")),
      p_kegg, width = 6.69, height = 8.0, dpi = 600)

    ggsave(file.path(outDir, paste0("KEGG_Dotplot_", label, ".pdf")),
      p_kegg, width = 6.69, height = 8.0)

    cat("KEGG:", nrow(kegg_df), "significant pathways saved\n")
  } else {
    cat("No significant KEGG pathways for", label, "\n")
  }

  # ---- Reactome ----
  reactome_res <- tryCatch(
    enrichPathway(gene = entrez_ids, universe = bg_entrez,
                  organism = "human", pAdjustMethod = "BH",
                  qvalueCutoff = 0.1, readable = TRUE),
    error = function(e) { cat("Reactome failed:", conditionMessage(e), "\n"); NULL }
  )

  if (!is.null(reactome_res) && nrow(as.data.frame(reactome_res)) > 0) {
    write.csv(as.data.frame(reactome_res),
              file.path(outDir, paste0("Reactome_", label, ".csv")),
              row.names = FALSE)
    p_react <- dotplot(reactome_res, showCategory = 20,
                        title = paste0("Reactome Enrichment — ", label))
    ggsave(file.path(outDir, paste0("Reactome_Dotplot_", label, ".png")),
           p_react, width = 6.69, height = 8.0, dpi = 600)
    ggsave(file.path(outDir, paste0("Reactome_Dotplot_", label, ".pdf")),
           p_react, width = 6.69, height = 8.0)
    cat("Reactome:", nrow(as.data.frame(reactome_res)), "significant pathways saved\n")
  } else {
    cat("No significant Reactome pathways for", label, "\n")
  }

  return(list(GO_BP = goe_bp, KEGG = kegg_res, Reactome = reactome_res))
}

# ------------------------------------------------------------
# HYPER vs HYPO 
# ------------------------------------------------------------

res_hyper <- run_stratified_enrichment(hyper_genes, bg_entrez, "Hypermethylated", outDir)
res_hypo  <- run_stratified_enrichment(hypo_genes,  bg_entrez, "Hypomethylated",  outDir)

# ------------------------------------------------------------
# PROMOTER vs NON-PROMOTER
# ------------------------------------------------------------

promoter_regex <- "TSS200|TSS1500|5'UTR|1stExon|5UTR|exon_1"

sig_DMPs_anno$Loci_Context <- ifelse(
  grepl(promoter_regex, sig_DMPs_anno$UCSC_RefGene_Group),
  "Promoter", "NonPromoter"
)

cat("\nPromoter vs Non-Promoter split among significant DMPs:\n")
print(table(sig_DMPs_anno$Loci_Context))

promoter_genes <- unique(unlist(lapply(
  subset(sig_DMPs_anno, Loci_Context == "Promoter")$UCSC_RefGene_Name,
  extract_genes
)))

nonpromoter_genes <- unique(unlist(lapply(
  subset(sig_DMPs_anno, Loci_Context == "NonPromoter")$UCSC_RefGene_Name,
  extract_genes
)))

res_promoter    <- run_stratified_enrichment(promoter_genes,    bg_entrez, "Promoter",    outDir)
res_nonpromoter <- run_stratified_enrichment(nonpromoter_genes, bg_entrez, "NonPromoter", outDir)

cat("\n=== Stratified pathway analysis complete ===\n")

