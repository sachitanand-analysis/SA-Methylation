# ============================================================
# 1. DMRcate — REGION-LEVEL DIFFERENTIAL METHYLATION
# ANNOTATE CpGs FOR DMR CALLING
# ============================================================

library(DMRcate)

dmr_annotated <- cpg.annotate(
  datatype      = "array",
  object        = mvals,
  what          = "M",
  arraytype     = "EPICv2",
  epicv2Remap   = TRUE,
  epicv2Filter  = "mean",
  analysis.type = "differential",
  design        = design_age,
  coef          = "GroupPUV",
  fdr           = 0.1
)

cat("\nIndividually significant probes feeding into DMR calling:",
    sum(dmr_annotated@ranges$is.sig), "\n")

# ------------------------------------------------------------
# 2 — CALL DMRs
# lambda = 1000 (max allowed distance between CpGs in a region) and
# C = 2 are the DMRcate defaults recommended for EPIC-family arrays.
# ------------------------------------------------------------

set.seed(1234)
dmrcate_res <- dmrcate(
  dmr_annotated,
  lambda = 1000,
  C      = 2,
  pcutoff = "fdr",
  min.cpgs = 2
)

cat("Total DMRs called:", nrow(dmrcate_res@coord), "\n")

# ------------------------------------------------------------
# 3 — EXTRACT GENOMIC RANGES + ANNOTATE TO GENES
# ------------------------------------------------------------

dmr_ranges <- extractRanges(dmrcate_res, genome = "hg38")

library(GenomicRanges)

cpg_gr <- GRanges(
    seqnames = seqnames(rowRanges(mSet)),
    ranges   = ranges(rowRanges(mSet)),
    CpG      = rownames(mSet)
)

hits <- findOverlaps(
    dmr_ranges,
    cpg_gr,
    ignore.strand = TRUE
)

dmr_cpgs <- split(
    mcols(cpg_gr)$CpG[subjectHits(hits)],
    queryHits(hits)
)

dmr_cpgs <- lapply(dmr_cpgs, function(x) paste(x, collapse = ";"))

dmr_df <- as.data.frame(dmr_ranges)

dmr_df$CpGs <- NA_character_
dmr_df$CpGs[as.integer(names(dmr_cpgs))] <- unlist(dmr_cpgs)
    
dmr_df$Recovered_CpGs <- ifelse(
    is.na(dmr_df$CpGs),
    0,
    sapply(strsplit(dmr_df$CpGs, ";"), length)
)

cat("\nRecovered CpGs per DMR:\n")
print(
    table(
        Recovered = dmr_df$Recovered_CpGs,
        Expected = dmr_df$no.cpgs
    )
)

cat("DMRs extracted as GRanges:", nrow(dmr_df), "\n")


get_genes_for_cpgs <- function(cpg_string) {

  if (is.na(cpg_string) || cpg_string == "")
    return(NA_character_)

  cpgs <- strsplit(cpg_string, ";")[[1]]

  genes <- ann$UCSC_RefGene_Name[
    match(cpgs, rownames(ann))
  ]

  genes <- unique(
    trimws(
      unlist(
        strsplit(
          genes[!is.na(genes)],
          ";"
        )
      )
    )
  )

  genes <- genes[
    genes != "" &
    !is.na(genes)
  ]

  if (length(genes) == 0)
    return(NA_character_)

  paste(genes, collapse = ";")

}

dmr_df$Genes_from_manifest <- sapply(
    dmr_df$CpGs,
    get_genes_for_cpgs
)

cat("\nDMRs with NA in extractRanges() gene column but recovered from manifest:\n")
print(
    sum(
        is.na(dmr_df$overlapping.genes) &
        !is.na(dmr_df$Genes_from_manifest)
    )
)

cat("DMRs with recovered genes:",
    sum(!is.na(dmr_df$Genes_from_manifest)),
    "of",
    nrow(dmr_df),
    "\n")
    
# ------------------------------------------------------------
# 4 — SUMMARY: DIRECTION AND SIZE OF DMRs
# ------------------------------------------------------------

dmr_df$Direction <- ifelse(dmr_df$meandiff > 0, "Hyper", "Hypo")

write.csv(dmr_df,
          file.path(outDir, "DMRcate_Regions_PUV_vs_Control.csv"),
          row.names = FALSE)

cat("\n=== DMR Summary ===\n")
cat("Total DMRs:", nrow(dmr_df), "\n")
print(table(dmr_df$Direction))
cat("Median width (bp):", median(dmr_df$width), "\n")
cat("Median No.CpGs per DMR:", median(dmr_df$no.cpgs), "\n")


# ============================================================
# 5. PATHWAY ANALYSIS OF DMRs
#     Total / Hyper / Hypo / Promoter / Non-Promoter
# ============================================================

library(dplyr)
library(tidyr)

dmr_df$DMR_ID <- paste0("DMR_", seq_len(nrow(dmr_df)))

dmr_cpg_rows <- dmr_df %>%
  dplyr::select(DMR_ID, Direction, CpGs) %>%
  dplyr::filter(!is.na(CpGs) & CpGs != "") %>%
  tidyr::separate_rows(CpGs, sep = ";") %>%
  dplyr::rename(CpG = CpGs)

cat("\nDMRs with at least one recovered CpG:",
    length(unique(dmr_cpg_rows$DMR_ID)), "of", nrow(dmr_df), "\n")

ann_gene_lookup <- data.frame(
  CpG                = rownames(ann_master),
  UCSC_RefGene_Group = ann_master$UCSC_RefGene_Group,
  UCSC_RefGene_Name  = ann_master$UCSC_RefGene_Name,
  stringsAsFactors   = FALSE
)

dmr_cpg_rows <- dplyr::left_join(dmr_cpg_rows, ann_gene_lookup, by = "CpG")

dmr_cpg_rows$Loci_Context <- ifelse(
  grepl(promoter_regex, dmr_cpg_rows$UCSC_RefGene_Group),
  "Promoter", "NonPromoter"
)

cat("\nPromoter vs Non-Promoter split among recovered DMR CpGs (CpG-level — for reference only, NOT used for gene-set construction below):\n")
print(table(dmr_cpg_rows$Loci_Context))

dmr_level_context <- dmr_cpg_rows %>%
  group_by(DMR_ID) %>%
  summarise(
    DMR_Loci_Context = if (any(Loci_Context == "Promoter")) "Promoter" else "NonPromoter",
    .groups = "drop"
  )

cat("\nPromoter vs Non-Promoter split at the DMR level (this is what pathway analysis below now uses):\n")
print(table(dmr_level_context$DMR_Loci_Context))

dmr_cpg_rows <- dplyr::left_join(dmr_cpg_rows, dmr_level_context, by = "DMR_ID")


total_dmr_genes <- unique(unlist(lapply(
  dmr_cpg_rows$UCSC_RefGene_Name, extract_genes
)))

hyper_dmr_genes <- unique(unlist(lapply(
  subset(dmr_cpg_rows, Direction == "Hyper")$UCSC_RefGene_Name, extract_genes
)))

hypo_dmr_genes <- unique(unlist(lapply(
  subset(dmr_cpg_rows, Direction == "Hypo")$UCSC_RefGene_Name, extract_genes
)))

promoter_dmr_genes <- unique(unlist(lapply(
  subset(dmr_cpg_rows, DMR_Loci_Context == "Promoter")$UCSC_RefGene_Name, extract_genes
)))

nonpromoter_dmr_genes <- unique(unlist(lapply(
  subset(dmr_cpg_rows, DMR_Loci_Context == "NonPromoter")$UCSC_RefGene_Name, extract_genes
)))

cat("\nDMR-associated gene set sizes:\n")
cat("Total:        ", length(total_dmr_genes), "\n")
cat("Hyper:        ", length(hyper_dmr_genes), "\n")
cat("Hypo:         ", length(hypo_dmr_genes), "\n")
cat("Promoter:     ", length(promoter_dmr_genes), "\n")
cat("Non-Promoter: ", length(nonpromoter_dmr_genes), "\n")

# ------------------------------------------------------------
# 6. Run the same GO BP / KEGG / Reactome enrichment used
# for the DMP-level stratified analysis, on each DMR gene set
# ------------------------------------------------------------

res_dmr_total        <- run_stratified_enrichment(total_dmr_genes,        bg_entrez, "Total_DMR",        outDir)
res_dmr_hyper        <- run_stratified_enrichment(hyper_dmr_genes,        bg_entrez, "Hypermethylated_DMR", outDir)
res_dmr_hypo         <- run_stratified_enrichment(hypo_dmr_genes,         bg_entrez, "Hypomethylated_DMR",  outDir)
res_dmr_promoter     <- run_stratified_enrichment(promoter_dmr_genes,     bg_entrez, "Promoter_DMR",     outDir)
res_dmr_nonpromoter  <- run_stratified_enrichment(nonpromoter_dmr_genes,  bg_entrez, "NonPromoter_DMR",  outDir)

cat("\n=== DMR-level stratified pathway analysis complete ===\n")

# ============================================================
# 7. DMR VISUALIZATION — GENOMIC FEATURE & CpG CONTEXT
# ============================================================


library(dplyr)
library(tidyr)
library(ggplot2)

dmr_df$DMR_ID <- paste0("DMR_", seq_len(nrow(dmr_df)))


assign_genomic_feature <- function(x) {
  if (is.na(x) || x == "") return("Non-Promoter")
  parts <- unique(trimws(unlist(strsplit(x, ";"))))
  if (any(c("TSS200", "TSS1500", "5'UTR", "1stExon") %in% parts)) return("Promoter")
  if (any(c("Body", "ExonBnd", "3'UTR") %in% parts)) return("Gene Body")
  return("Non-Promoter")
}

anno_lookup <- data.frame(
  CpG                 = rownames(ann_master),
  UCSC_RefGene_Group  = ann_master$UCSC_RefGene_Group,
  Relation_Primary    = ann_master$Relation_Primary,
  stringsAsFactors    = FALSE
)

dmr_cpgs <- dmr_df %>%
  dplyr::select(DMR_ID, Direction, CpGs) %>%
  dplyr::filter(!is.na(CpGs) & CpGs != "") %>%
  tidyr::separate_rows(CpGs, sep = ";") %>%
  dplyr::rename(CpG = CpGs) %>%
  dplyr::left_join(anno_lookup, by = "CpG")

dmr_cpgs$GenomicFeature <- sapply(dmr_cpgs$UCSC_RefGene_Group, assign_genomic_feature)

dmr_cpgs$CpG_Context_Detailed <- factor(
  dmr_cpgs$Relation_Primary,
  levels = c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea")
)


feature_priority <- function(features) {
  features <- features[!is.na(features)]
  if (length(features) == 0) return(NA_character_)
  if ("Promoter" %in% features) return("Promoter")
  if ("Gene Body" %in% features) return("Gene Body")
  return("Intergenic")
}

context_priority <- function(contexts) {
  contexts <- as.character(contexts[!is.na(contexts)])
  if (length(contexts) == 0) return(NA_character_)
  if (any(contexts %in% c("Island"))) return("Island")
  if (any(contexts %in% c("N_Shore", "S_Shore"))) return("Shore")
  if (any(contexts %in% c("N_Shelf", "S_Shelf"))) return("Shelf")
  return("OpenSea")
}

dmr_feature <- dmr_cpgs %>%
  group_by(DMR_ID) %>%
  summarise(
    Direction      = dplyr::first(Direction),
    GenomicFeature = feature_priority(GenomicFeature),
    .groups = "drop"
  )

dmr_context <- dmr_cpgs %>%
  group_by(DMR_ID) %>%
  summarise(
    Direction          = dplyr::first(Direction),
    CpG_Context_Simple = context_priority(CpG_Context_Detailed),
    .groups = "drop"
  )

dmr_feature$Methylation_Status <- ifelse(dmr_feature$Direction == "Hyper", "Hypermethylated", "Hypomethylated")
dmr_context$Methylation_Status <- ifelse(dmr_context$Direction == "Hyper", "Hypermethylated", "Hypomethylated")

#--------------------------------------------------------------
# 8. PROPORTIONS (within each genomic feature / CpG region, not raw DMR counts)
#--------------------------------------------------------------

feature_counts <- dmr_feature %>%
  dplyr::filter(!is.na(GenomicFeature)) %>%
  count(GenomicFeature, Methylation_Status) %>%
  group_by(GenomicFeature) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()

feature_counts$GenomicFeature <- factor(
  feature_counts$GenomicFeature,
  levels = c("Promoter", "Gene Body", "Intergenic")
)

context_counts <- dmr_context %>%
  dplyr::filter(!is.na(CpG_Context_Simple)) %>%
  count(CpG_Context_Simple, Methylation_Status) %>%
  group_by(CpG_Context_Simple) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()

context_counts$CpG_Region <- factor(
  context_counts$CpG_Context_Simple,
  levels = c("Island", "Shore", "Shelf", "OpenSea")
)

#--------------------------------------------------------------
# 9. PLOT — Genomic feature 
#--------------------------------------------------------------

p_dmr_feature <- ggplot(
  feature_counts,
  aes(x = GenomicFeature, y = Proportion, fill = Methylation_Status)
) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Hypomethylated" = "#4C72B0", "Hypermethylated" = "#D62728")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Genomic feature distribution of DMRs", x = NULL, y = "Proportion of DMRs")

ggsave(file.path(outDir, "DMR_Genomic_Feature_proportions.png"),
       p_dmr_feature, width = 3.35, height = 5.5, dpi = 600)
ggsave(file.path(outDir, "DMR_Genomic_Feature_proportions.pdf"),
       p_dmr_feature, width = 3.35, height = 5.5)

#--------------------------------------------------------------
# 10. PLOT — CpG region context 
#--------------------------------------------------------------

p_dmr_context <- ggplot(
  context_counts,
  aes(x = CpG_Region, y = Proportion, fill = Methylation_Status)
) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Hypomethylated" = "#4C72B0", "Hypermethylated" = "#D62728")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "CpG region distribution of DMRs", x = NULL, y = "Proportion of DMRs")

ggsave(file.path(outDir, "DMR_CpG_region_proportions.png"),
       p_dmr_context, width = 3.35, height = 5.5, dpi = 600)
ggsave(file.path(outDir, "DMR_CpG_region_proportions.pdf"),
       p_dmr_context, width = 3.35, height = 5.5)

cat("\n=== DMR genomic feature / CpG region proportion plots saved ===\n")

cat("\n=== CpG Region Counts and Percentages ===\n")

cpg_region_counts <- table(sig_DMPs_anno$Relation_Primary)

cpg_region_pct <- round(
  100 * cpg_region_counts / sum(cpg_region_counts),
  2
)

cpg_region_summary <- data.frame(
  CpG_Region = names(cpg_region_counts),
  Count = as.integer(cpg_region_counts),
  Percent = cpg_region_pct,
  row.names = NULL
)

print(cpg_region_summary)

cat("\n=== Promoter vs Non-Promoter ===\n")

promoter_counts <- table(sig_DMPs_anno$Loci_Context)

promoter_pct <- round(
  100 * promoter_counts / sum(promoter_counts),
  2
)

promoter_summary <- data.frame(
  Region = names(promoter_counts),
  Count = as.integer(promoter_counts),
  Percent = promoter_pct,
  row.names = NULL
)

print(promoter_summary)

cat("\n=== Hyper/Hypo Distribution within CpG Regions ===\n")

cpg_direction <- prop.table(
  table(sig_DMPs_anno$Relation_Primary,
        sig_DMPs_anno$Direction),
  margin = 1
) * 100

print(round(cpg_direction, 2))

cat("\n=== Hyper/Hypo Distribution within Promoter Status ===\n")

promoter_direction <- prop.table(
  table(sig_DMPs_anno$Loci_Context,
        sig_DMPs_anno$Direction),
  margin = 1
) * 100

print(round(promoter_direction, 2))

round(prop.table(
  table(dmr_context$CpG_Context_Simple,
        dmr_context$Methylation_Status),
  margin = 1
) * 100, 2)
