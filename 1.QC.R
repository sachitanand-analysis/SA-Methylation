# ============================================================
# 1. LIBRARIES
# ============================================================

library(minfi)
library(limma)
library(sva)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(genefilter)
library(dendextend)
library(RColorBrewer)
library(ggrepel)
library(matrixStats)
library(KEGGREST)
library(igraph)
library(ggraph)
library(enrichR)       
library(ReactomePA)
library(scales)
library(enrichplot)
library(stringr)
library(ggalluvial)
library(readr)
library(stringr)

# --- Directories ---
baseDir  <- ""
outDir   <- "PUV_vs_Control"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. LOAD DATA
# ============================================================

targets <- read.metharray.sheet(baseDir)
targets$Group <- factor(targets$Group, levels = c("Control", "PUV"))

rgSet <- read.metharray.exp(targets = targets, extended = TRUE)
rgSet@annotation <- c(array      = "IlluminaHumanMethylationEPICv2",
                       annotation = "20a1.hg38") 
                       
# ============================================================
# 3. QUALITY CONTROL
# ============================================================

detP <- detectionP(rgSet)

failed_samples <- colMeans(detP) > 0.05
if (any(failed_samples)) {
  message("Removing failed samples: ",
          paste(colnames(rgSet)[failed_samples], collapse = ", "))
  rgSet   <- rgSet[, !failed_samples]
  targets <- targets[!failed_samples, ]
  detP    <- detP[, !failed_samples]
}

mSet_qc <- preprocessRaw(rgSet)
qc      <- getQC(mSet_qc)
png(file.path(outDir, "QC_plot.png"), width = 800, height = 600)
plotQC(qc)
dev.off()

# ============================================================
# 4. NORMALIZATION
# ============================================================

mSet_noob <- preprocessNoob(rgSet)
saveRDS(mSet_noob, file.path(outDir, "mSet_noob.rds"))

# ============================================================
# 5a. PROBE FILTERING
# ============================================================

keep_detP      <- rowSums(detP < 0.01) == ncol(detP)
mSet_detP      <- mSet_noob[keep_detP, ]
mSet_detP      <- mapToGenome(mSet_detP)
mSet_after_snp <- dropLociWithSnps(mSet_detP,
                                    snps = c("SBE", "CpG"), maf = 0.05)

ann      <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
keep_auto <- !(featureNames(mSet_after_snp) %in%
               ann$Name[ann$chr %in% c("chrX", "chrY")])
mSet     <- mSet_after_snp[keep_auto, ]

cat("Probes remaining after all filtering:", nrow(mSet), "\n")
saveRDS(mSet, file.path(outDir, "mSet_filtered.rds"))

# ============================================================
# 5b. PROBE FILTERING SUMMARY
# ============================================================

n_total     <- nrow(mSet_noob)
n_filtered  <- nrow(mSet)
n_removed   <- n_total - n_filtered
pct_removed <- round(100 * n_removed / n_total, 2)

cat("\n=== Probe Filtering Summary ===\n")
cat("Total probes before filtering:", n_total, "\n")
cat("Probes after filtering:",        n_filtered, "\n")
cat("Probes removed:",                n_removed, "\n")
cat("Percentage removed:",            pct_removed, "%\n")
cat("After detectionP filter:",       sum(keep_detP), "\n")
cat("After SNP removal:",             nrow(mSet_after_snp), "\n")
cat("After sex chr removal:",         nrow(mSet), "\n")

lost_probes <- setdiff(rownames(mSet_noob), rownames(mSet))
write.table(lost_probes, file.path(outDir, "lost_probes.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)


# Cross-reactive probes removal 

PIDSLEY_CSV <- "pidsley2024.csv"
pid_manifest <- read.csv(PIDSLEY_CSV, stringsAsFactors = FALSE,
                           row.names = 1)   

cat("pidsley2024.csv loaded. Rows:", nrow(pid_manifest), "\n")


cat("\n=== Cross-Reactive Probe Diagnostic ===\n")
cat("Total probes in manifest:           ", nrow(pid_manifest), "\n")


offtarget_table <- table(pid_manifest$Num_offtargets == 0)
cat("Probes with Num_offtargets == 0 (clean):    ",
    sum(pid_manifest$Num_offtargets == 0, na.rm = TRUE), "\n")
cat("Probes with Num_offtargets  > 0 (cross-reactive): ",
    sum(pid_manifest$Num_offtargets  > 0, na.rm = TRUE), "\n")


cat("Probes with missing position (MissingPos == Y): ",
    sum(pid_manifest$MissingPos == "Y", na.rm = TRUE), "\n")


overlap <- sum(featureNames(mSet) %in% rownames(pid_manifest))
cat("Your probes matching manifest:      ", overlap,
    "of", nrow(mSet), "\n\n")

qc_summary <- data.frame(
  Step  = character(),
  Count = numeric(),
  Note  = character(),
  stringsAsFactors = FALSE
)

# BUILD PROBE REMOVAL LISTS

cross_reactive <- rownames(pid_manifest)[
  !is.na(pid_manifest$Num_offtargets) &
  pid_manifest$Num_offtargets > 0
]
cat("Cross-reactive probes identified:   ", length(cross_reactive), "\n")

# --- Non-mapping probes 

non_mapping <- rownames(pid_manifest)[
  !is.na(pid_manifest$MissingPos) &
  pid_manifest$MissingPos == "Y"
]
cat("Non-mapping probes (chr0):           ", length(non_mapping), "\n")

# --- Combine: all probes to remove ---
probes_to_remove <- unique(c(cross_reactive, non_mapping))
cat("Total probes flagged for removal:   ", length(probes_to_remove), "\n\n")

# APPLY FILTER TO mSet

probes_before_cr <- nrow(mSet)

all_probe_names <- featureNames(mSet)

keep_cr <- !all_probe_names %in% probes_to_remove

removed_probe_names <- all_probe_names[!keep_cr]

mSet <- mSet[keep_cr, ]

cat("=== Cross-Reactive Probe Removal Summary ===\n")
cat("Probes before filter:  ", probes_before_cr, "\n")
cat("Probes removed:        ", length(removed_probe_names), "\n")
cat("  of which cross-reactive (Num_offtargets > 0):",
    sum(cross_reactive %in% removed_probe_names), "\n")
cat("  of which non-mapping  (chr0/MissingPos == Y):",
    sum(non_mapping %in% removed_probe_names), "\n")
cat("Probes remaining:      ", nrow(mSet), "\n\n")

# UPDATE QC SUMMARY TABLE

qc_summary <- rbind(
  qc_summary,
  data.frame(
    Step = c(
      "After cross-reactive/non-mapping probe removal"
    ),

    Count = c(
      nrow(mSet)
    ),

    Note = c(
      paste0(
        "Removed ",
        length(removed_probe_names),
        " probes using Peters et al. 2024 EPICv2 manifest"
      )
    ),

    stringsAsFactors = FALSE
  )
)

write.csv(qc_summary, file.path(outDir, "QC_Summary.csv"), row.names = FALSE)

# SANITY CHECK
cat("=== Sanity check: any cross-reactive probes remain? ===\n")
remaining_cr <- sum(featureNames(mSet) %in% cross_reactive)
cat("Cross-reactive probes still in mSet:", remaining_cr,
    "(should be 0)\n\n")

# Verify all remaining probes have genomic coordinates
# (no chr0 probes should remain)
cat("=== Sanity check: any chr0 probes remain? ===\n")
remaining_nm <- sum(featureNames(mSet) %in% non_mapping)
cat("Non-mapping probes still in mSet:   ", remaining_nm,
    "(should be 0)\n\n")

