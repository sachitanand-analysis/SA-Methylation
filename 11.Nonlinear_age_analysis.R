# ============================================================
# NONLINEAR AGE SENSITIVITY ANALYSIS
# PUV vs Control — EPIC v2
#
# Primary model:
#   ~ Group + Age + SV1 + SV2 + SV3
#
# Sensitivity model:
#   ~ Group + ns(Age, df = 3) + SVs
#
# Same:
#   - 865,432 QC-retained CpGs
#   - M-value matrix
#   - FDR < 0.1
#   - |logFC| > 0.2
# ============================================================

library(limma)
library(sva)
library(splines)


cat("\n============================================\n")
cat("STEP 1: BASIC CHECKS\n")
cat("============================================\n")

cat("Number of CpGs:", nrow(mvals), "\n")
cat("Number of samples:", ncol(mvals), "\n")


cat("Age is numeric:",
    is.numeric(targets$Age), "\n")

cat("Missing age values:",
    sum(is.na(targets$Age)), "\n")


cat("\nGroup counts:\n")
print(table(targets$Group))

targets$Sample_ID <- basename(targets$Basename)

cat("\nSamples matching mvals:",
    sum(colnames(mvals) %in% targets$Sample_ID),
    "\n")

cat("Targets matching mvals:",
    sum(targets$Sample_ID %in% colnames(mvals)),
    "\n")

cat("Duplicate Sample_IDs:",
    anyDuplicated(targets$Sample_ID),
    "\n")

targets <- targets[
    match(colnames(mvals), targets$Sample_ID),
    ,
    drop = FALSE
]

rownames(targets) <- targets$Sample_ID

cat("\nSample order identical:",
    identical(colnames(mvals), rownames(targets)),
    "\n")

stopifnot(
    identical(colnames(mvals), rownames(targets))
)

targets$Group <- factor(
    targets$Group,
    levels = c("Control", "PUV")
)

cat("\nGroup levels:\n")
print(levels(targets$Group))

stopifnot(
    levels(targets$Group)[1] == "Control"
)

cat("\n============================================\n")
cat("STEP 2: NONLINEAR AGE MODEL\n")
cat("============================================\n")


age_spline <- ns(
    targets$Age,
    df = 3
)


# ============================================================
# CREATE FULL AND NULL SVA MODELS
# ============================================================

# Full model:
# Group + nonlinear age

mod_spline <- model.matrix(
    ~ Group + age_spline,
    data = targets
)


# Null model:
# nonlinear age only

mod0_spline <- model.matrix(
    ~ age_spline,
    data = targets
)


cat("\nFull SVA model columns:\n")
print(colnames(mod_spline))

cat("\nNull SVA model columns:\n")
print(colnames(mod0_spline))


# ============================================================
# RECALCULATE SVA
# ============================================================

cat("\n============================================\n")
cat("STEP 3: RECALCULATE SVA\n")
cat("============================================\n")


sva_spline <- sva(
    mvals,
    mod_spline,
    mod0_spline
)


# Number of SVs estimated
n_sv_spline <- ncol(sva_spline$sv)

cat("\nNumber of SVs estimated:",
    n_sv_spline,
    "\n")


if (n_sv_spline > 0) {

    for (i in seq_len(n_sv_spline)) {

        targets[[paste0("SV", i)]] <-
            sva_spline$sv[, i]

    }

}

cat("\n============================================\n")
cat("STEP 4: LIMMA MODEL\n")
cat("============================================\n")


design_formula <- "~ Group + age_spline"


if (n_sv_spline > 0) {

    sv_terms <- paste(
        paste0("SV", seq_len(n_sv_spline)),
        collapse = " + "
    )

    design_formula <- paste(
        design_formula,
        "+",
        sv_terms
    )

}


cat("\nDesign formula:\n")
print(design_formula)


design_age_spline <- model.matrix(
    as.formula(design_formula),
    data = targets
)


cat("\nDesign matrix columns:\n")
print(colnames(design_age_spline))


stopifnot(
    "GroupPUV" %in% colnames(design_age_spline)
)


# ============================================================
# RUN LIMMA
# ============================================================

cat("\n============================================\n")
cat("STEP 5: LIMMA\n")
cat("============================================\n")


fit_age_spline <- lmFit(
    mvals,
    design_age_spline
)


fit_age_spline <- eBayes(
    fit_age_spline
)

res_age_spline <- topTable(
    fit_age_spline,
    coef = "GroupPUV",
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
)


# Add CpG IDs
res_age_spline$CpG <- rownames(res_age_spline)


sig_DMPs_age_spline <- subset(
    res_age_spline,
    adj.P.Val < 0.1 &
    abs(logFC) > 0.2
)


# Direction
sig_DMPs_age_spline$Direction <- ifelse(
    sig_DMPs_age_spline$logFC > 0,
    "Hyper",
    "Hypo"
)


# ============================================================
# SUMMARY
# ============================================================

cat("\n============================================\n")
cat("NONLINEAR-AGE RESULTS\n")
cat("============================================\n")

cat(
    "Significant DMPs:",
    nrow(sig_DMPs_age_spline),
    "\n"
)

cat(
    "Hypermethylated:",
    sum(sig_DMPs_age_spline$Direction == "Hyper"),
    "\n"
)

cat(
    "Hypomethylated:",
    sum(sig_DMPs_age_spline$Direction == "Hypo"),
    "\n"
)


# ============================================================
# COMPARE WITH PRIMARY LINEAR-AGE ANALYSIS
# ============================================================

cat("\n============================================\n")
cat("STEP 6: COMPARISON WITH PRIMARY MODEL\n")
cat("============================================\n")


# Primary DMPs
linear_DMPs <- unique(sig_DMPs$CpG)


# Nonlinear-age DMPs
spline_DMPs <- unique(sig_DMPs_age_spline$CpG)


# Shared
common_DMPs <- intersect(
    linear_DMPs,
    spline_DMPs
)


# Union
union_DMPs <- union(
    linear_DMPs,
    spline_DMPs
)


# Linear-only
linear_only <- setdiff(
    linear_DMPs,
    spline_DMPs
)


# Spline-only
spline_only <- setdiff(
    spline_DMPs,
    linear_DMPs
)


# Retention
retention <- 100 *
    length(common_DMPs) /
    length(linear_DMPs)


# Jaccard
jaccard <- length(common_DMPs) /
    length(union_DMPs)


cat(
    "Linear-age DMPs:",
    length(linear_DMPs),
    "\n"
)

cat(
    "Nonlinear-age DMPs:",
    length(spline_DMPs),
    "\n"
)

cat(
    "Common DMPs:",
    length(common_DMPs),
    "\n"
)

cat(
    "Linear-only DMPs:",
    length(linear_only),
    "\n"
)

cat(
    "Spline-only DMPs:",
    length(spline_only),
    "\n"
)

cat(
    "Percentage of primary DMPs retained:",
    round(retention, 2),
    "%\n"
)

cat(
    "Jaccard index:",
    round(jaccard, 4),
    "\n"
)


# ============================================================
# GENOME-WIDE logFC CONCORDANCE
# ============================================================

cat("\n============================================\n")
cat("STEP 7: GENOME-WIDE CONCORDANCE\n")
cat("============================================\n")


# Both analyses should contain the same 865,432 CpGs

common_all <- intersect(
    rownames(all_DMPs),
    rownames(res_age_spline)
)


cat(
    "Common CpGs across full models:",
    length(common_all),
    "\n"
)


# Pearson
pearson_all <- cor(
    all_DMPs[common_all, "logFC"],
    res_age_spline[common_all, "logFC"],
    method = "pearson",
    use = "complete.obs"
)


# Spearman
spearman_all <- cor(
    all_DMPs[common_all, "logFC"],
    res_age_spline[common_all, "logFC"],
    method = "spearman",
    use = "complete.obs"
)


# Direction concordance
linear_logFC <- all_DMPs[
    common_all,
    "logFC"
]

spline_logFC <- res_age_spline[
    common_all,
    "logFC"
]


same_direction <- sign(linear_logFC) ==
    sign(spline_logFC)


direction_concordance <- mean(
    same_direction,
    na.rm = TRUE
) * 100


cat(
    "Genome-wide Pearson r:",
    pearson_all,
    "\n"
)

cat(
    "Genome-wide Spearman rho:",
    spearman_all,
    "\n"
)

cat(
    "Genome-wide direction concordance:",
    round(direction_concordance, 2),
    "%\n"
)


# ============================================================
# SAVE RESULTS
# ============================================================

write.csv(
    res_age_spline,
    file = file.path(
        outDir,
        "PUV_vs_Control_NonlinearAge_All_CpGs.csv"
    ),
    row.names = FALSE
)


write.csv(
    sig_DMPs_age_spline,
    file = file.path(
        outDir,
        "PUV_vs_Control_NonlinearAge_Significant_DMPs_FDR0.1.csv"
    ),
    row.names = FALSE
)


# Save shared DMP comparison
comparison_sets <- data.frame(
    CpG = union_DMPs,
    Linear_Age = union_DMPs %in% linear_DMPs,
    Nonlinear_Age = union_DMPs %in% spline_DMPs
)

write.csv(
    comparison_sets,
    file = file.path(
        outDir,
        "PUV_vs_Control_Linear_vs_NonlinearAge_DMP_Comparison.csv"
    ),
    row.names = FALSE
)


# ============================================================
# SAVE SUMMARY
# ============================================================

summary_results <- data.frame(
    Metric = c(
        "Linear-age DMPs",
        "Nonlinear-age DMPs",
        "Common DMPs",
        "Linear-only DMPs",
        "Spline-only DMPs",
        "DMP retention (%)",
        "Jaccard index",
        "Genome-wide Pearson r",
        "Genome-wide Spearman rho",
        "Genome-wide direction concordance (%)",
        "Number of nonlinear-age SVs"
    ),
    Value = c(
        length(linear_DMPs),
        length(spline_DMPs),
        length(common_DMPs),
        length(linear_only),
        length(spline_only),
        retention,
        jaccard,
        pearson_all,
        spearman_all,
        direction_concordance,
        n_sv_spline
    )
)


write.csv(
    summary_results,
    file = file.path(
        outDir,
        "PUV_vs_Control_NonlinearAge_Sensitivity_Summary.csv"
    ),
    row.names = FALSE
)


cat("\n============================================\n")
cat("ANALYSIS COMPLETE\n")
cat("============================================\n")

cat(
    "Results saved to:",
    outDir,
    "\n"
)

# ============================================================
# UpSet plot: Linear-age vs nonlinear-age DMPs
# ============================================================

library(ggplot2)
library(ComplexUpset)

# ------------------------------------------------------------
# 1. Create membership table
# ------------------------------------------------------------

linear_DMPs <- unique(sig_DMPs$CpG)
spline_DMPs <- unique(sig_DMPs_age_spline$CpG)

all_DMPs_compare <- union(
    linear_DMPs,
    spline_DMPs
)

upset_df <- data.frame(
    CpG = all_DMPs_compare,
    Linear_age = all_DMPs_compare %in% linear_DMPs,
    Nonlinear_age = all_DMPs_compare %in% spline_DMPs
)

# Check counts
cat("Linear-age DMPs:",
    sum(upset_df$Linear_age), "\n")

cat("Nonlinear-age DMPs:",
    sum(upset_df$Nonlinear_age), "\n")

cat("Shared DMPs:",
    sum(upset_df$Linear_age & upset_df$Nonlinear_age), "\n")

# ------------------------------------------------------------
# 2. UpSet plot
# ------------------------------------------------------------

p_upset <- upset(
    upset_df,
    intersect = c("Linear_age", "Nonlinear_age"),
    name = "DMP overlap",

    base_annotations = list(
        "Intersection size" =
            intersection_size(
                text = list(
                    vjust = -0.5,
                    size = 5
                )
            )
    ),

    set_sizes = upset_set_size(
        geom = geom_bar(
            width = 0.6
        )
    ),

    width_ratio = 0.25
) +

    theme_classic(base_size = 13) +

    theme(
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        strip.text = element_text(
            size = 12,
            face = "bold"
        )
    )

p_upset

ggsave(
    filename = file.path(
        outDir,
        "Figure_Nonlinear_Age_UpSet.pdf"
    ),
    plot = p_upset,
    width = 7,
    height = 5.5,
    units = "in"
)

ggsave(
    filename = file.path(
        outDir,
        "Figure_Nonlinear_Age_UpSet.png"
    ),
    plot = p_upset,
    width = 7,
    height = 5.5,
    units = "in",
    dpi = 600
)

# ============================================================
# Genome-wide logFC concordance plot
# Linear-age vs nonlinear-age model
# ============================================================

library(ggplot2)

comparison_all <- data.frame(
    CpG = common_all,
    logFC_linear = all_DMPs[common_all, "logFC"],
    logFC_spline = res_age_spline[common_all, "logFC"]
)

# Remove missing values
comparison_all <- comparison_all[
    complete.cases(comparison_all),
]

# Correlations
pearson_all <- cor(
    comparison_all$logFC_linear,
    comparison_all$logFC_spline,
    method = "pearson"
)

spearman_all <- cor(
    comparison_all$logFC_linear,
    comparison_all$logFC_spline,
    method = "spearman"
)

# Direction concordance
direction_concordance <- mean(
    sign(comparison_all$logFC_linear) ==
    sign(comparison_all$logFC_spline)
) * 100


# ============================================================
# Plot
# ============================================================

p_logFC_concordance <- ggplot(
    comparison_all,
    aes(
        x = logFC_linear,
        y = logFC_spline
    )
) +

    geom_point(
        size = 0.4,
        alpha = 0.25
    ) +

    # Identity line
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        linewidth = 0.7
    ) +

    # Statistics moved to upper-left corner
    annotate(
        "text",
        x = -4.7,
        y = 5.7,
        label = paste0(
            "Pearson r = ",
            sprintf("%.3f", pearson_all),
            "\n",
            "Spearman \u03c1 = ",
            sprintf("%.3f", spearman_all),
            "\n",
            "Direction concordance = ",
            sprintf("%.1f", direction_concordance),
            "%"
        ),
        hjust = 0,
        vjust = 1,
        size = 4.2
    ) +

    labs(
        x = "logFC: linear-age model",
        y = "logFC: nonlinear-age model"
    ) +

    theme_classic(base_size = 13) +

    theme(
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        plot.margin = margin(10, 15, 10, 10)
    )

p_logFC_concordance

ggsave(
    filename = file.path(
        outDir,
        "Figure_NonlinearAge_GenomeWide_logFC_Concordance2.png"
    ),
    plot = p_logFC_concordance,
    width = 7,
    height = 6,
    units = "in",
    dpi = 600
)
