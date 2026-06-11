library(ComplexUpset)
library(ggplot2)
library(dplyr)

# ============================================================
# Create significance categories
# ============================================================

upset_df <- data.frame(
  CpG = rownames(all_DMPs)
)

upset_df$FDR05_Hyper <- with(
  all_DMPs,
  adj.P.Val < 0.05 & logFC > 0.2
)

upset_df$FDR05_Hypo <- with(
  all_DMPs,
  adj.P.Val < 0.05 & logFC < -0.2
)

upset_df$FDR1_Hyper <- with(
  all_DMPs,
  adj.P.Val < 0.1 & logFC > 0.2
)

upset_df$FDR1_Hypo <- with(
  all_DMPs,
  adj.P.Val < 0.1 & logFC < -0.2
)

# ============================================================
# KEEP ONLY SIGNIFICANT PROBES
# ============================================================

upset_df <- subset(
  upset_df,
  FDR05_Hyper |
  FDR05_Hypo |
  FDR1_Hyper |
  FDR1_Hypo
)

# ============================================================
# Plot
# ============================================================

png(
  "Sensitivity_UpSet_SignificantOnly.png",
  width = 2007,
  height = 1677,
  res = 300
)

upset(
  upset_df,
  c(
    "FDR05_Hyper",
    "FDR05_Hypo",
    "FDR1_Hyper",
    "FDR1_Hypo"
  ),

  width_ratio = 0.18,

  base_annotations = list(
    "Intersection size" =
      intersection_size(
        text = list(size = 6)
      )
  )
)

dev.off()

