# ============================================================
# 1. CREATE BED FILE FOR HOMER 
# ============================================================

bed_df <- sig_DMPs_anno[
  
  !is.na(sig_DMPs_anno$chr) &
  !is.na(sig_DMPs_anno$pos),
  
  c("chr", "pos", "Name")
]

# Ensure UCSC chromosome format
bed_df$chr <- as.character(bed_df$chr)

bed_df$chr <- ifelse(
  grepl("^chr", bed_df$chr),
  bed_df$chr,
  paste0("chr", bed_df$chr)
)

# Create 200bp windows centered on CpGs
bed_df$start <- pmax(0, bed_df$pos - 100)

bed_df$end <- bed_df$pos + 100

bed_out <- bed_df[
  ,
  c("chr", "start", "end", "Name")
]

write.table(
  bed_out,
  
  file.path(
    outDir,
    "Significant_DMPs_200bp.bed"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

cat("BED file created\n")

# ============================================================
# 2. HOMER RESULTS - Run in Terminal
# ============================================================

findMotifsGenome.pl \
Significant_DMPs_200bp.bed \
hg38 \
HOMER_Motif_Output \
-size given \
-mask \
-p 8

# ============================================================
# 3. MOTIF BASED TF ENRICHMENT (HOMER)
# ============================================================

homer_res <- read.delim("knownResults.txt")

head(homer_res)

plot_df <- homer_res[1:30, ]

plot_df$TF <- str_extract(
  plot_df$Motif.Name,
  "^[^/]+"
)

plot_df$TF <- gsub(
  "\\(",
  " (",
  plot_df$TF
)

plot_df$TF <- trimws(plot_df$TF)

plot_df$TargetPct <- as.numeric(
  gsub(
    "%",
    "",
    plot_df$X..of.Target.Sequences.with.Motif
  )
)

plot_df$BgPct <- as.numeric(
  gsub(
    "%",
    "",
    plot_df$X..of.Background.Sequences.with.Motif
  )
)

plot_df$logP <- -log10(plot_df$P.value)

plot_df$FoldEnrichment <- (
  plot_df$TargetPct /
  plot_df$BgPct
)

plot_df$Family <- str_extract(
  plot_df$TF,
  "\\(.*?\\)"
)

plot_df$Family <- gsub(
  "[\\(\\)]",
  "",
  plot_df$Family
)

plot_df$TF <- factor(
  plot_df$TF,
  levels = rev(plot_df$TF)
)

# PLOT

p_nature <- ggplot(
  plot_df,
  
  aes(
    x = logP,
    y = TF,
    size = FoldEnrichment,
    color = logP
  )
) +

geom_point(alpha = 0.95) +

scale_color_gradient(
  low = "#6BAED6",
  high = "#CB181D"
) +

theme_bw(base_size = 16) +

theme(
  
  axis.text.y = element_text(size = 10),
  
  axis.text.x = element_text(size = 12),
  
  plot.title = element_text(
    size = 11,
    face = "bold",
    hjust = 0.5
  )
) + 
labs(
  
  x = expression(-log[10]("P-value")),
  
  y = "Transcription Factor",
  
  size = "Fold\nEnrichment",
  
  color = expression(-log[10]("P-value"))
)

print(p_nature)

ggsave(
  file.path(
    outDir,
    "HOMER.png"
  ),
  
  p_nature,
  
  width = 6.69,
  height = 8.5,
  dpi = 300
)

ggsave(
  file.path(
    outDir,
    "HOMER.pdf"
  ),

  p_nature,

  width = 6.69,
  height = 8.5
)

# ============================================================
# 4. TF ANALYSIS — ChEA/ENCODE (binding-based)
# ============================================================

cat("\nRunning TF binding site enrichment (ChEA/ENCODE)...\n")

sig_genes_tf <- sig_DMPs_anno %>%
  separate_rows(UCSC_RefGene_Name, sep=";") %>%
  mutate(UCSC_RefGene_Name = trimws(UCSC_RefGene_Name)) %>%
  filter(UCSC_RefGene_Name != "") %>%
  pull(UCSC_RefGene_Name) %>%
  unique()

cat("Genes for TF analysis:", length(sig_genes_tf), "\n")

dbs_tf <- c(
  "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
  "ChEA_2022"
)

tf_results <- enrichr(sig_genes_tf, dbs_tf)

encode_tf <- tf_results[["ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"]]
chea_tf   <- tf_results[["ChEA_2022"]]

write.csv(encode_tf,
          file.path(outDir,"ENCODE_TFBS_Enrichment.csv"),
          row.names=FALSE)
write.csv(chea_tf,
          file.path(outDir,"ChEA_TFBS_Enrichment.csv"),
          row.names=FALSE)

plot_df_tf <- encode_tf %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n=30)

plot_df_tf$logP <- -log10(plot_df_tf$Adjusted.P.value)
plot_df_tf$TF   <- gsub("_.*", "", plot_df_tf$Term)

p_tf <- ggplot(plot_df_tf,
               aes(x=logP, y=reorder(TF,logP),
                   color=Combined.Score, size=Combined.Score)) +
  geom_point() +
  scale_color_gradientn(
    colors=c("#2166AC","#92C5DE","#FFFFBF","#F4A582","#B2182B"),
    name="Combined score"
  ) +
  scale_size_continuous(range=c(3,8), name="Combined score") +
  theme_minimal(base_size=11) +
  labs(x=expression(-log[10]~"(Adjusted P-value)"),
       y="Transcription Factor"
   )

ggsave(file.path(outDir,"TF_ENCODE_ChEA_Enrichment.png"),
       p_tf, width=6.69, height=8.5, dpi=300)
ggsave(file.path(outDir,"TF_ENCODE_ChEA_Enrichment.pdf"),
       p_tf, width=6.69, height=8.5)
cat("TF ChEA/ENCODE plot saved\n")

##########################################
# 5. TOP HOMER HITS - Run in Terminal
##########################################

findMotifsGenome.pl \
Significant_DMPs_200bp.bed \
hg38 \
homer_out \
-find HOMER_Motif_Output/knownResults/known1.motif \
> ELK4_hits.txt

findMotifsGenome.pl \
Significant_DMPs_200bp.bed \
hg38 \
homer_out \
-find HOMER_Motif_Output/knownResults/known3.motif \
> ELK1_hits.txt

findMotifsGenome.pl \
Significant_DMPs_200bp.bed \
hg38 \
homer_out \
-find HOMER_Motif_Output/knownResults/known4.motif \
> FLI1_hits.txt

findMotifsGenome.pl \
Significant_DMPs_200bp.bed \
hg38 \
homer_out \
-find HOMER_Motif_Output/knownResults/known5.motif \
> ELF1_hits.txt

findMotifsGenome.pl \
Significant_DMPs_200bp.bed \
hg38 \
homer_out \
-find HOMER_Motif_Output/knownResults/known7.motif \
> ETV4_hits.txt

findMotifsGenome.pl \
Significant_DMPs_200bp.bed \
hg38 \
homer_out \
-find HOMER_Motif_Output/knownResults/known9.motif \
> GABPA_hits.txt

# ============================================================
# 6. ETS MOTIF - HOMER motif-hit based analysis
# ============================================================

# RECREATE BED DATAFRAME

bed_df <- sig_DMPs_anno %>%
  
  filter(
    !is.na(chr),
    !is.na(pos)
  ) %>%
  
  mutate(
    start = pmax(1, pos - 100),
    end   = pos + 100
  ) %>%
  
  select(
    chr,
    start,
    end,
    CpG,
    UCSC_RefGene_Name,
    logFC,
    delta_beta,
    Direction,
    Relation_Simple
  )
  
# READ HOMER MOTIF-HIT FILES

cat("\nReading HOMER motif-hit files...\n")

elk4 <- read.delim(
  "ELK4_hits.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

elk4$Motif <- "ELK4(ETS)"

elk1 <- read.delim(
  "ELK1_hits.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

elk1$Motif <- "ELK1(ETS)"

fli1 <- read.delim(
  "FLI1_hits.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

fli1$Motif <- "FLI1(ETS)"

elf1 <- read.delim(
  "ELF1_hits.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

elf1$Motif <- "ELF1(ETS)"

etv4 <- read.delim(
  "ETV4_hits.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

etv4$Motif <- "ETV4(ETS)"

gabpa <- read.delim(
  "GABPA_hits.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

gabpa$Motif <- "GABPA(ETS)"

# MERGE MOTIF-HIT TABLES

motif_hits <- bind_rows(
  elk4,
  elk1,
  fli1,
  elf1,
  etv4,
  gabpa
)

cat(
  "Total motif hits:",
  nrow(motif_hits),
  "\n"
)

# MATCH HOMER MOTIF HITS TO DMP ANNOTATIONS

cat("\nMatching motif hits to DMP genes...\n")

# Create annotation lookup table
bed_df2 <- sig_DMPs_anno %>%
  
  select(
    CpG,
    UCSC_RefGene_Name,
    logFC,
    delta_beta,
    Direction,
    Relation_Simple
  )

# Match HOMER PositionID to CpG IDs
motif_gene_df <- inner_join(
  
  motif_hits,
  
  bed_df2,
  
  by = c(
    "PositionID" = "CpG"
  )
)

cat(
  "Motif-gene overlaps:",
  nrow(motif_gene_df),
  "\n"
)

# EXPAND MULTI-GENE ANNOTATIONS

motif_gene_df <- motif_gene_df %>%
  
  separate_rows(
    UCSC_RefGene_Name,
    sep = ";"
  ) %>%
  
  mutate(
    UCSC_RefGene_Name =
      trimws(UCSC_RefGene_Name)
  ) %>%
  
  filter(
    UCSC_RefGene_Name != "",
    !is.na(UCSC_RefGene_Name)
  )

# STEP 6 — SUMMARISE MOTIF-GENE RELATIONSHIPS

sankey_df <- motif_gene_df %>%
  
  group_by(
    Motif,
    Gene = UCSC_RefGene_Name
  ) %>%
  
  summarise(
    
    N_DMPs = n(),
    
    Direction =
      names(
        sort(
          table(Direction),
          decreasing = TRUE
        )
      )[1],
    
    .groups = "drop"
  )
cat(
  "Unique motif-gene edges:",
  nrow(sankey_df),
  "\n"
)

gene_counts <- sankey_df %>%
  
  group_by(Gene) %>%
  
  summarise(
    total_DMPs = sum(N_DMPs),
    .groups = "drop"
  )

sankey_df <- left_join(
  
  sankey_df,
  
  gene_counts,
  
  by = "Gene"
)

top20_genes <- sankey_df %>%
  
  group_by(Gene) %>%
  
  summarise(
    total = sum(N_DMPs),
    .groups = "drop"
  ) %>%
  
  arrange(desc(total)) %>%
  
  slice_head(n = 20) %>%
  
  pull(Gene)

sankey_df <- sankey_df %>%
  
  filter(
    Gene %in% top20_genes
  )

cat(
  "After top-gene filtering:",
  nrow(sankey_df),
  "edges\n"
)

# Order motifs
sankey_df$Motif <- factor(
  sankey_df$Motif,
  
  levels = c(
    "ELF1(ETS)",
    "ELK1(ETS)",
    "ELK4(ETS)",
    "ETV4(ETS)",
    "FLI1(ETS)",
    "GABPA(ETS)"
  )
)

# Order genes by number of motif associations
gene_order <- sankey_df %>%
  
  group_by(Gene) %>%
  
  summarise(
    n_connections = n(),
    .groups = "drop"
  ) %>%
  
  arrange(desc(n_connections), Gene)

sankey_df$Gene <- factor(
  sankey_df$Gene,
  levels = gene_order$Gene
)

# PLOT 

motif_colors <- c(
  "ELF1(ETS)"  = "#00BFC4",
  "ELK1(ETS)"  = "#619CFF",
  "ELK4(ETS)"  = "#F8766D",
  "ETV4(ETS)"  = "#C77CFF",
  "FLI1(ETS)"  = "#7CAE00",
  "GABPA(ETS)" = "#FF9F1C"
)

p_sankey <- ggplot(
  
  sankey_df,
  
  aes(
    axis1 = Motif,
    axis2 = Gene,
    y     = N_DMPs
  )
) +

  geom_alluvium(
    
    aes(fill = Motif),
    
    alpha      = 0.35,
    
    width      = 0.03,
    
    curve_type = "sigmoid"
  ) +

  geom_stratum(
    
    fill      = "grey88",
    
    color     = "grey35",
    
    width     = 0.12,
    
    linewidth = 0.25
  ) +

  geom_text(
    
    stat = "stratum",
    
    aes(label = after_stat(stratum)),
    
    size = 3,
    
    fontface = "bold"
  ) +

  scale_fill_manual(
    
    values = motif_colors,
    
    name   = "ETS motif"
  ) +

  scale_x_discrete(
    
    limits = c(
      "ETS Motif",
      "Associated Gene"
    ),
    
    expand = c(0.08, 0.05)
  ) +

  labs(
    x = NULL,
    
    y = NULL
  ) +

  theme_minimal(base_size = 13) +

  theme(
    
    axis.text.x =
      element_text(
        face = "bold",
        size = 12,
        color = "black"
      ),
    
    axis.text.y =
      element_blank(),
    
    axis.ticks.y =
      element_blank(),
    
    axis.title.y =
      element_blank(),
    
    panel.grid =
      element_blank(),
    
    plot.title =
      element_text(
        face = "bold",
        size = 15
      ),
    
    legend.title =
      element_text(
        face = "bold"
      )
  )

ggsave(
  
  file.path(
    outDir,
    "ETS_Motif_Gene_Sankey.png"
  ),
  
  p_sankey,
  
  width  = 6.69,
  height = 8.5,
  dpi    = 300
)

ggsave(
  
  file.path(
    outDir,
    "ETS_Motif_Gene_Sankey.pdf"
  ),
  
  p_sankey,
  
  width  = 6.69,
  height = 8.5
)

write.csv(
  
  sankey_df,
  
  file.path(
    outDir,
    "ETS_Motif_Gene_Table.csv"
  ),
  
  row.names = FALSE
)

cat("ETS motif Sankey plot saved\n")
