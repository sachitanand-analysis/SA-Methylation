# ============================================================
# 11. missMethyl PROBE-NUMBER-ADJUSTED ENRICHMENT
#     PUV vs Control
# ============================================================

library(missMethyl)

cat("\n====================================================\n")
cat("missMethyl Enrichment Analysis — PUV vs Control\n")
cat("====================================================\n")


# ------------------------------------------------------------
# 1. DEFINE CpG SETS
# ------------------------------------------------------------

# Significant DMPs from the primary analysis
sig_cpg <- unique(sig_DMPs$CpG)

# All CpGs actually tested in the DMP analysis
all_cpg <- unique(rownames(mvals))

# Keep only significant CpGs that were actually tested
sig_cpg <- intersect(sig_cpg, all_cpg)

cat("Background CpGs:", length(all_cpg), "\n")
cat("Significant CpGs:", length(sig_cpg), "\n")


# ------------------------------------------------------------
# 2. HYPER- AND HYPOMETHYLATED CpGs
# ------------------------------------------------------------

hyper_cpg <- unique(
    sig_DMPs$CpG[sig_DMPs$Direction == "Hyper"]
)

hypo_cpg <- unique(
    sig_DMPs$CpG[sig_DMPs$Direction == "Hypo"]
)
 
hyper_cpg <- intersect(hyper_cpg, all_cpg)
hypo_cpg  <- intersect(hypo_cpg, all_cpg)

cat("Hypermethylated CpGs:", length(hyper_cpg), "\n")
cat("Hypomethylated CpGs:", length(hypo_cpg), "\n")


# ============================================================
# 3. GO ENRICHMENT — ALL SIGNIFICANT DMPs
# ============================================================

GO_all <- gometh(
    sig.cpg   = sig_cpg,
    all.cpg   = all_cpg,
    collection = "GO",
    array.type = "EPIC_V2",
    prior.prob = TRUE
)

write.table(
    GO_all,
    file = file.path(
        outDir,
        "missMethyl_GO_All_DMPs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 4. KEGG ENRICHMENT — ALL SIGNIFICANT DMPs
# ============================================================

KEGG_all <- gometh(
    sig.cpg   = sig_cpg,
    all.cpg   = all_cpg,
    collection = "KEGG",
    array.type = "EPIC_V2",
    prior.prob = TRUE
)

write.table(
    KEGG_all,
    file = file.path(
        outDir,
        "missMethyl_KEGG_All_DMPs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 5. GO — HYPERMETHYLATED DMPs
# ============================================================

GO_hyper <- gometh(
    sig.cpg   = hyper_cpg,
    all.cpg   = all_cpg,
    collection = "GO",
    array.type = "EPIC_V2",
    prior.prob = TRUE
)

write.table(
    GO_hyper,
    file = file.path(
        outDir,
        "missMethyl_GO_Hyper_DMPs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 6. GO — HYPOMETHYLATED DMPs
# ============================================================

GO_hypo <- gometh(
    sig.cpg   = hypo_cpg,
    all.cpg   = all_cpg,
    collection = "GO",
    array.type = "EPIC_V2",
    prior.prob = TRUE
)

write.table(
    GO_hypo,
    file = file.path(
        outDir,
        "missMethyl_GO_Hypo_DMPs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 7. KEGG — HYPERMETHYLATED DMPs
# ============================================================

KEGG_hyper <- gometh(
    sig.cpg   = hyper_cpg,
    all.cpg   = all_cpg,
    collection = "KEGG",
    array.type = "EPIC_V2",
    prior.prob = TRUE
)

write.table(
    KEGG_hyper,
    file = file.path(
        outDir,
        "missMethyl_KEGG_Hyper_DMPs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 8. KEGG — HYPOMETHYLATED DMPs
# ============================================================

KEGG_hypo <- gometh(
    sig.cpg   = hypo_cpg,
    all.cpg   = all_cpg,
    collection = "KEGG",
    array.type = "EPIC_V2",
    prior.prob = TRUE
)

write.table(
    KEGG_hypo,
    file = file.path(
        outDir,
        "missMethyl_KEGG_Hypo_DMPs.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 9. SAVE FDR < 0.05 RESULTS
# ============================================================

write.table(
    subset(GO_all, FDR < 0.05),
    file = file.path(
        outDir,
        "missMethyl_GO_All_DMPs_FDR05.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    subset(GO_hyper, FDR < 0.05),
    file = file.path(
        outDir,
        "missMethyl_GO_Hyper_DMPs_FDR05.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    subset(GO_hypo, FDR < 0.05),
    file = file.path(
        outDir,
        "missMethyl_GO_Hypo_DMPs_FDR05.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    subset(KEGG_all, FDR < 0.05),
    file = file.path(
        outDir,
        "missMethyl_KEGG_All_DMPs_FDR05.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    subset(KEGG_hyper, FDR < 0.05),
    file = file.path(
        outDir,
        "missMethyl_KEGG_Hyper_DMPs_FDR05.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    subset(KEGG_hypo, FDR < 0.05),
    file = file.path(
        outDir,
        "missMethyl_KEGG_Hypo_DMPs_FDR05.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 10. QUICK SUMMARY FDR 0.05
# ============================================================

cat("\n====================================================\n")
cat("missMethyl Analysis Summary\n")
cat("====================================================\n")

cat("GO — all DMPs, FDR < 0.05:",
    sum(GO_all$FDR < 0.05, na.rm = TRUE), "\n")

cat("GO — hyper DMPs, FDR < 0.05:",
    sum(GO_hyper$FDR < 0.05, na.rm = TRUE), "\n")

cat("GO — hypo DMPs, FDR < 0.05:",
    sum(GO_hypo$FDR < 0.05, na.rm = TRUE), "\n")

cat("KEGG — all DMPs, FDR < 0.05:",
    sum(KEGG_all$FDR < 0.05, na.rm = TRUE), "\n")

cat("KEGG — hyper DMPs, FDR < 0.05:",
    sum(KEGG_hyper$FDR < 0.05, na.rm = TRUE), "\n")

cat("KEGG — hypo DMPs, FDR < 0.05:",
    sum(KEGG_hypo$FDR < 0.05, na.rm = TRUE), "\n")

# ============================================================
# 11. With Plot bias
# ============================================================

png(
    filename = file.path(outDir, "PUV_vs_Control_gometh_probe_bias_GO.png"),
    width = 5400,
    height = 4200,
    res = 600
)

go_all_pb <- gometh(
    sig.cpg   = sig_cpg,
    all.cpg   = all_cpg,
    collection = "GO",
    array.type = "EPIC_V2",
    prior.prob = TRUE,
    plot.bias = TRUE
)

dev.off()

png(
    filename = file.path(outDir, "PUV_vs_Control_gometh_probe_bias_KEGG.png"),
    width = 1800,
    height = 1400,
    res = 200
)

kegg_all_pb <- gometh(
    sig.cpg   = sig_cpg,
    all.cpg   = all_cpg,
    collection = "KEGG",
    array.type = "EPIC_V2",
    prior.prob = TRUE,
    plot.bias = TRUE
)
dev.off()

