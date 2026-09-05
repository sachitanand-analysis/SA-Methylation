# ============================================================
# 1. MANIFEST ANNOTATION — CpG island context
# ============================================================

manifest <- read.csv(
  "EPIC-8v2-0_A1.csv",
  skip = 7, stringsAsFactors = FALSE
)

raw_vals <- manifest$Relation_to_UCSC_CpG_Island

assign_primary_context <- function(x) {
  if (is.na(x) || x == "") return("OpenSea")
  parts <- unique(trimws(unlist(strsplit(x, ";"))))
  if ("Island"  %in% parts) return("Island")
  if ("N_Shore" %in% parts) return("N_Shore")
  if ("S_Shore" %in% parts) return("S_Shore")
  if ("N_Shelf" %in% parts) return("N_Shelf")
  if ("S_Shelf" %in% parts) return("S_Shelf")
  return("OpenSea")
}

assign_simple_context <- function(x) {
  if (is.na(x) || x == "") return("OpenSea")
  parts <- unique(trimws(unlist(strsplit(x, ";"))))
  if ("Island" %in% parts)                    return("Island")
  if (any(c("N_Shore","S_Shore") %in% parts)) return("Shore")
  if (any(c("N_Shelf","S_Shelf") %in% parts)) return("Shelf")
  return("OpenSea")
}

is_multi_context <- function(x) {
  if (is.na(x) || x == "") return(FALSE)
  parts <- unique(trimws(unlist(strsplit(x, ";"))))
  sum(c("Island" %in% parts,
        any(c("N_Shore","S_Shore") %in% parts),
        any(c("N_Shelf","S_Shelf") %in% parts))) > 1
}

cat("Assigning CpG island contexts from manifest...\n")
manifest$Relation_Primary <- sapply(raw_vals, assign_primary_context)
manifest$Relation_Simple  <- sapply(raw_vals, assign_simple_context)
manifest$Is_Multi_Context <- sapply(raw_vals, is_multi_context)
cat("Done.\n")

idx <- match(rownames(ann), manifest$IlmnID)
cat("Matched probes:", sum(!is.na(idx)), "\n")

ann$Relation_Primary  <- manifest$Relation_Primary[idx]
ann$Relation_Simple   <- manifest$Relation_Simple[idx]
ann$Is_Multi_Context  <- manifest$Is_Multi_Context[idx]

na_idx <- is.na(ann$Relation_Primary)
if (sum(na_idx) > 0) {
  existing <- ann$Relation_to_Island[na_idx]
  ann$Relation_Primary[na_idx] <- as.character(existing)
  ann$Relation_Simple[na_idx]  <- as.character(existing)
  ann$Is_Multi_Context[na_idx] <- FALSE
}

cat("Final annotation counts:\n")
print(table(ann$Relation_Primary))

# ============================================================
# 2. ANNOTATE DMPs
# ============================================================

ann_cols <- c("chr","pos","strand","Name",
              "UCSC_RefGene_Name","UCSC_RefGene_Group",
              "Relation_to_Island","Relation_Primary",
              "Relation_Simple","Is_Multi_Context",
              "Islands_Name","Regulatory_Feature_Group")
ann_cols    <- ann_cols[ann_cols %in% colnames(ann)]
ann_sub     <- as.data.frame(ann)[, ann_cols]

ann_master  <- ann_sub

common_cpgs   <- intersect(rownames(all_DMPs), rownames(ann_sub))
all_DMPs_anno <- cbind(all_DMPs[common_cpgs, ], ann_sub[common_cpgs, ])

common_sig    <- intersect(rownames(sig_DMPs), rownames(ann_sub))
sig_DMPs_anno <- cbind(sig_DMPs[common_sig, ], ann_sub[common_sig, ])

cat("Annotated DMPs:", nrow(all_DMPs_anno), "\n")
cat("Annotated sig DMPs:", nrow(sig_DMPs_anno), "\n")

cat("\n=== Direction label check ===\n")
cat("Unique values in sig_DMPs_anno$Direction:\n")
print(table(sig_DMPs_anno$Direction, useNA = "always"))

cat("\nUnique values in all_DMPs_anno$Direction:\n")
print(table(all_DMPs_anno$Direction, useNA = "always"))

write.table(all_DMPs_anno,
            file.path(outDir, "All_DMPs_annotated.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(sig_DMPs_anno,
            file.path(outDir, "Significant_DMPs_annotated.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ============================================================
# 3a. REGULATION INTERPRETATION
# ============================================================

# ---------- ALL DMPs ----------

all_DMPs_anno$Regulation <- "Unclear"

# Hypermethylated promoter-associated CpGs
all_DMPs_anno$Regulation[
  all_DMPs_anno$Direction == "Hyper" &
  grepl(
    "TSS200|TSS1500|5'UTR|1stExon|5UTR|exon_1",
    all_DMPs_anno$UCSC_RefGene_Group
  )
] <- "Promoter Hypermethylation"

# Hypomethylated promoter-associated CpGs
all_DMPs_anno$Regulation[
  all_DMPs_anno$Direction == "Hypo" &
  grepl(
    "TSS200|TSS1500|5'UTR|1stExon|5UTR|exon_1",
    all_DMPs_anno$UCSC_RefGene_Group
  )
] <- "Promoter Hypomethylation"

# ---------- SIGNIFICANT DMPs ----------

sig_DMPs_anno$Regulation <- "Unclear"

sig_DMPs_anno$Regulation[
  sig_DMPs_anno$Direction == "Hyper" &
  grepl(
    "TSS200|TSS1500|5'UTR|1stExon|5UTR|exon_1",
    sig_DMPs_anno$UCSC_RefGene_Group
  )
] <- "Promoter Hypermethylation"

sig_DMPs_anno$Regulation[
  sig_DMPs_anno$Direction == "Hypo" &
  grepl(
    "TSS200|TSS1500|5'UTR|1stExon|5UTR|exon_1",
    sig_DMPs_anno$UCSC_RefGene_Group
  )
] <- "Promoter Hypomethylation"

# ============================================================
# SUMMARY COUNTS
# ============================================================

reg_summary <- sig_DMPs_anno %>%
  count(Regulation, Direction) %>%
  filter(Regulation != "Unclear")

cat("\n=== Regulation interpretation summary ===\n")
print(reg_summary)
cat("Promoter Hypomethylation:",
    sum(sig_DMPs_anno$Regulation == "Promoter Hypomethylation"), "\n")
cat("Promoter Hypermethylation:",
    sum(sig_DMPs_anno$Regulation == "Promoter Hypermethylation"), "\n")
cat("Unclear:",
    sum(sig_DMPs_anno$Regulation == "Unclear"), "\n")
    
# ============================================================
# 4. HYPO & HYPER METHYL ONLY
# ============================================================

hyper_sig <- subset(sig_DMPs_anno, logFC > 0)
hypo_sig  <- subset(sig_DMPs_anno, logFC < 0)

write.table(hyper_sig, file.path(outDir, "Hyper_sig.tsv"), sep="\t", quote=FALSE)
write.table(hypo_sig,  file.path(outDir, "Hypo_sig.tsv"),  sep="\t", quote=FALSE)

hyper_all <- subset(all_DMPs_anno, logFC > 0)
hypo_all  <- subset(all_DMPs_anno, logFC < 0)

write.table(hyper_all,
            file.path(outDir, "Hyper_all.tsv"),
            sep="\t", quote=FALSE)

write.table(hypo_all,
            file.path(outDir, "Hypo_all.tsv"),
            sep="\t", quote=FALSE)
            
# ============================================================
# 5. PERCENTAGE SUMMARY STATISTICS
# ============================================================

cat("\n====================================================\n")
cat("PERCENTAGE SUMMARY OF SIGNIFICANT DMPs\n")
cat("====================================================\n")

# ============================================================
# A. HYPER vs HYPO PERCENTAGE
# ============================================================

n_total_sig <- nrow(sig_DMPs_anno)

n_hyper <- sum(sig_DMPs_anno$Direction == "Hyper", na.rm = TRUE)
n_hypo  <- sum(sig_DMPs_anno$Direction == "Hypo",  na.rm = TRUE)

pct_hyper <- round((n_hyper / n_total_sig) * 100, 2)
pct_hypo  <- round((n_hypo  / n_total_sig) * 100, 2)

cat("\n=== Hyper vs Hypomethylation ===\n")
cat("Total significant DMPs:", n_total_sig, "\n")
cat("Hypermethylated:", n_hyper, "(", pct_hyper, "% )\n")
cat("Hypomethylated :", n_hypo,  "(", pct_hypo,  "% )\n")

# ============================================================
# B. CpG ISLAND CONTEXT DISTRIBUTION
# ============================================================

cat("\n=== CpG Island Context Distribution ===\n")

cpg_counts <- table(sig_DMPs_anno$Relation_Primary)

for (ctx in names(cpg_counts)) {

  pct_ctx <- round((cpg_counts[ctx] / n_total_sig) * 100, 2)

  cat(ctx, ":",
      cpg_counts[ctx],
      "(", pct_ctx, "% )\n")
}

# ============================================================
# C. HYPER/HYPO WITHIN CpG CONTEXT
# ============================================================

cat("\n=== Hyper/Hypo Distribution within CpG Context ===\n")

cpg_direction <- prop.table(
  table(sig_DMPs_anno$Relation_Primary,
        sig_DMPs_anno$Direction),
  margin = 1
) * 100

print(round(cpg_direction, 2))
