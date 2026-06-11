# ============================================================
# 1. EXTRACT BETA AND M-VALUES
# ============================================================

beta  <- getBeta(mSet)
mvals <- getM(mSet)

saveRDS(beta,  file.path(outDir, "beta_values.rds"))
saveRDS(mvals, file.path(outDir, "m_values.rds"))

png(file.path(outDir, "Beta_density_plot.png"), width = 900, height = 600)
densityPlot(beta, sampGroups = targets$Group,
            main = "Beta Value Distribution — Post-Normalization",
            legend = TRUE)
dev.off()


# ============================================================
# 2. SVA — LEEK METHOD 
# ============================================================

targets$Age_scaled <- scale(targets$Age)[,1]

mod_age  <- model.matrix(~ Group + Age_scaled, data = targets)
mod0_age <- model.matrix(~ Age_scaled, data = targets)

# Estimate number of surrogate variables
n_sv <- num.sv(mvals, mod_age, method = "leek")
cat("Estimated surrogate variables:", n_sv, "\n")


sva_age <- sva(
  mvals,
  mod_age,
  mod0_age,
  n.sv = n_sv
)

# Remove batch effect from M-values for downstream use

mvals_corrected <- removeBatchEffect(mvals,
                                      covariates = sva_age$sv,
                                      design = mod_age)

saveRDS(mvals_corrected, file.path(outDir, "mvals_SVA_corrected.rds"))

# ============================================================
# 3. DIFFERENTIAL METHYLATION ANALYSIS
# ============================================================

# Design matrix (SVs included as covariates )
design_age <- cbind(mod_age, sva_age$sv)
colnames(design_age) <- c(colnames(mod_age),
                           paste0("SV", seq_len(ncol(sva_age$sv))))

fit_age  <- lmFit(mvals, design_age)   
fit2_age <- eBayes(fit_age)

# Extract all DMPs
all_DMPs <- topTable(fit2_age,
                      coef        = "GroupPUV",
                      number      = Inf,
                      adjust.method = "BH",
                      sort.by     = "p")
all_DMPs$CpG <- rownames(all_DMPs)

delta_beta       <- rowMeans(beta[, targets$Group == "PUV"]) -
                    rowMeans(beta[, targets$Group == "Control"])
names(delta_beta) <- rownames(beta)
all_DMPs$delta_beta <- delta_beta[match(all_DMPs$CpG, names(delta_beta))]

all_DMPs$Direction <- ifelse(all_DMPs$logFC > 0, "Hyper",
                      ifelse(all_DMPs$logFC < 0, "Hypo", "NoChange"))

write.table(all_DMPs,
            file.path(outDir, "All_DMPs_PUV_vs_Control.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sig_DMPs <- subset(all_DMPs, adj.P.Val < 0.1 & abs(logFC) > 0.2)
sig_DMPs$delta_beta <- delta_beta[match(sig_DMPs$CpG, names(delta_beta))]
sig_DMPs$Direction  <- ifelse(sig_DMPs$logFC > 0, "Hyper",
                       ifelse(sig_DMPs$logFC < 0, "Hypo", "NoChange"))

cat("Significant DMPs:", nrow(sig_DMPs), "\n")
cat("  Hypermethylated:", sum(sig_DMPs$Direction == "Hyper"), "\n")
cat("  Hypomethylated: ", sum(sig_DMPs$Direction == "Hypo"),  "\n")

write.table(sig_DMPs,
            file.path(outDir, "Significant_DMPs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
