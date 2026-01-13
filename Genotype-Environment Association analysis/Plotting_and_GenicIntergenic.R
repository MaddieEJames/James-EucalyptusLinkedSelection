#### Plotting median results for the BF (in dB) ####

library(ggplot2)
library(dplyr)
library(tibble)
library(grid)

# Files
# The meadian file, calculated from the three replicates
median_file <- "repMedian_summary_betai_reg.out"

# As baypass doesn't retain the chromosome and SNP information,
# this tab-separated file contained the chromosome in the first column and the
# position in the second column, so that the results can be plotted based on this
chrpos_file <- "allPopsMAF0.01_LD-50-5-0.2_chr1-11_CHR-POS.tsv"

# Read data
med <- read.table(median_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
chrpos <- read.table(chrpos_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                     col.names = c("CHROM_true", "POS_true"))

# Attach CHR/POS to median table
med$CHROM_true <- as.character(chrpos$CHROM_true)
med$POS_true   <- chrpos$POS_true

# Contig lengths (bp)
lens <- tribble(
  ~CHROM, ~length_bp,
  "1",  44964737,
  "2",  59526133,
  "3",  83949477,
  "4",  41122621,
  "5",  76232106,
  "6",  57468648,
  "7",  54826121,
  "8",  72515911,
  "9",  39306088,
  "10", 37759200,
  "11", 44831864
)

# Order chromosomes as 1..11 and compute cumulative offsets
chrom_levels <- as.character(1:11)
lens <- lens %>%
  mutate(CHROM = factor(CHROM, levels = chrom_levels)) %>%
  arrange(CHROM) %>%
  mutate(cum_offset = cumsum(dplyr::lag(length_bp, default = 0L)))

# Merge offsets into med and compute global bp coordinate
med <- med %>%
  left_join(lens %>% select(CHROM, length_bp, cum_offset),
            by = c("CHROM_true" = "CHROM")) %>%
  mutate(
    x_bp = cum_offset + POS_true,
    chrOrder = factor(CHROM_true, levels = chrom_levels)
  )

# Boundaries (between chrs 1..10)
chr_boundaries <- lens %>%
  filter(as.numeric(as.character(CHROM)) < 11) %>%
  mutate(boundary = cum_offset + length_bp) %>%
  pull(boundary)

# Midpoints for labels (numbers 1..11)
chr_midpoints <- lens %>%
  mutate(midpoint = cum_offset + length_bp / 2) %>%
  select(CHROM, midpoint)

midpoints_df <- data.frame(
  chromosome = as.numeric(as.character(chr_midpoints$CHROM)),
  midpoint   = chr_midpoints$midpoint
)

# Threshold (10 dB)
thr_db <- 10

# Max position (end of chr 11)
max_pos <- max(med$x_bp, na.rm = TRUE)

# Axis breaks: ticks ONLY at boundaries + start + end
x_breaks <- c(0, chr_boundaries, max_pos)
x_labels <- rep("", length(x_breaks))  # no labels at ticks

# Plot
p <- ggplot(med, aes(x = x_bp, y = BF.dB., color = chrOrder)) +
  geom_point(alpha = 1, size = 2.5) +
  geom_hline(yintercept = thr_db,
             linetype = "solid",
             color = "#E41A1C",
             linewidth = 1.2) +
  
# Removed chromosome number labels
  
  labs(x = "Chromosome", y = expression(BF[iS]~"(dB)")) +
  scale_color_manual(values = rep(c("#4682B4", "#70A8CD"), length.out = 11)) +
  scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(3, 15),
    breaks = c(5, 10, 15),
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0, max_pos),
    ylim = c(3, 15),
    expand = FALSE,
    clip = "on"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    legend.position = "none",
    axis.title.x = element_text(size = 28, margin = margin(t = 30), family = "Helvetica"),
    axis.title.y = element_text(size = 28, family = "Helvetica", margin = margin(r = 15)),
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(size = 22, family = "Helvetica", margin = margin(r = 10)),
    axis.ticks.x = element_line(linewidth = 1.2),
    axis.ticks.y = element_line(linewidth = 1.2),
    axis.line    = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.3, "cm"),
    plot.margin  = margin(0.2, 0.2, 0.5, 0.2, "in"),
    aspect.ratio = NULL
  )

# Save plot
ggsave(
  filename = "BFis_FSTstyle_manhattan.png",
  plot     = p,
  width    = 24,
  height   = 4.5,
  units    = "in",
  dpi      = 300
)

# Summaries

# Number of SNPs greater or equal to 20
sum(med$BF.dB. >= 20)
# 0

# Number of SNPs greater or equal to 15 but less than 20
sum(med$BF.dB. >= 15 & med$BF.dB. < 20)
# 19

# Number of SNPs greater or equal to 10 but less than 15
sum(med$BF.dB. >= 10 & med$BF.dB. < 15)
# 790


#### Plotting BF vs. Effect Size (Magnitude) Scatter Plot ####

library(ggplot2)
library(dplyr)

# Load data
covis.snp.med <- read.table("repMedian_summary_betai_reg.out", 
                            header = TRUE, sep = "\t", 
                            stringsAsFactors = FALSE)

# Calculate absolute effect size
covis.snp.med <- covis.snp.med %>%
  mutate(abs_Beta_is = abs(Beta_is))

# Calculate correlations
pearson_cor <- cor.test(covis.snp.med$abs_Beta_is, covis.snp.med$BF.dB., 
                        method = "pearson")
spearman_cor <- cor.test(covis.snp.med$abs_Beta_is, covis.snp.med$BF.dB., 
                         method = "spearman")

# Print correlation statistics
cat("\n=== Correlation Analysis: BF ~ |Beta| ===\n")
cat("\nPearson's correlation:\n")
cat("  r =", round(pearson_cor$estimate, 4), "\n")
cat("  p-value:", format.pval(pearson_cor$p.value, digits = 3), "\n")
cat("\nSpearman's rank correlation:\n")
cat("  rho =", round(spearman_cor$estimate, 4), "\n")
cat("  p-value:", format.pval(spearman_cor$p.value, digits = 3), "\n\n")

# Create scatter plot
p_scatter <- ggplot(covis.snp.med, aes(x = abs_Beta_is, y = BF.dB.)) +
  geom_point(alpha = 0.4, size = 1.7, color = "black") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = expression("|"*beta[iS]*"| (absolute effect size magnitude)"),
    y = expression(BF[iS]~"(dB)")
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.title.x = element_text(size = 28, margin = margin(t = 15)),
    axis.title.y = element_text(size = 28, margin = margin(r = 15)),
    axis.text.x  = element_text(size = 22, family = "Helvetica"),
    axis.text.y  = element_text(size = 22, family = "Helvetica"),
    axis.line    = element_line(linewidth = 1.2),
    axis.ticks   = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    plot.margin  = margin(0.5, 0.5, 0.5, 0.5, "in")
  )
print(p_scatter)
ggsave(
  filename = "BF_vs_EffectSize_scatter.png",
  plot     = p_scatter,
  width    = 10,
  height   = 8,
  units    = "in",
  dpi      = 300
)


#### Map Eucalyptus outlier SNPs to genes and Arabidopsis orthologues ####

library(GenomicRanges)
library(rtracklayer)
library(dplyr)

# Inputs
# The gff file and that annotation file come from the Eucalyptus grandis genome v2.0 (https://phytozome-next.jgi.doe.gov/info/Egrandis_v2_0)
gff_file <- "Egrandis_297_v2.0.gene.gff3"
ann_file <- "Egrandis_297_v2.0.annotation_info.txt"
out_file <- "outlier_genes_arabi.tsv"

# Select outliers: 10 <= BF < 20
outliers <- med[med$BF.dB. >= 10 & med$BF.dB. < 20, ]

# Load GFF
gff   <- import(gff_file)
genes <- gff[gff$type == "gene"]
mrna  <- gff[gff$type == "mRNA"]

# Map SNP chromosomes to GFF seqlevels
gff_seqlevels <- seqlevels(gff)
gff_chr_num   <- suppressWarnings(as.integer(sub(".*?([0-9]+)$", "\\1", gff_seqlevels)))
chr_map       <- setNames(gff_seqlevels, gff_chr_num)

snp_seqnames <- chr_map[as.integer(outliers$CHROM_true)]
if (any(is.na(snp_seqnames))) stop("SNP chromosomes do not map to GFF seqlevels")

snps_gr <- GRanges(
  seqnames = snp_seqnames,
  ranges   = IRanges(outliers$POS_true, outliers$POS_true),
  SNP_ID   = outliers$SNP
)

# Normalize gene/transcript IDs
normalize_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.v[0-9]+\\.[0-9]+$", "", x)
  x <- sub("\\.p$", "", x)
  x <- sub("\\.t$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}

mcols(genes)$gene_id <- if (!is.null(mcols(genes)$Name)) {
  normalize_id(mcols(genes)$Name)
} else {
  normalize_id(mcols(genes)$ID)
}

# Load Arabidopsis annotations (needed for transcript selection)
ann <- read.delim(ann_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE)
ann$locusName <- trimws(as.character(ann$locusName))

strip_arabi_isoform <- function(x) sub("\\.[0-9]+$", "", x)

ann_gene <- ann %>%
  select(locusName, `Best-hit-arabi-name`) %>%
  filter(!is.na(`Best-hit-arabi-name`) & `Best-hit-arabi-name` != "") %>%
  mutate(arabi_gene_stripped = strip_arabi_isoform(`Best-hit-arabi-name`)) %>%
  distinct(locusName, arabi_gene_stripped) %>%
  group_by(locusName) %>%
  summarize(
    arabi_gene = paste(sort(unique(arabi_gene_stripped)), collapse = ";"),
    .groups    = "drop"
  )

# Get longest annotated transcript per gene (or longest if none annotated)
if (length(mrna) > 0) {
  mcols(mrna)$transcript_id <- normalize_id(mcols(mrna)$ID)
  mcols(mrna)$gene_id       <- normalize_id(mcols(mrna)$Parent)
  mcols(mrna)$tx_length     <- width(mrna)
  
  tx_to_gene <- data.frame(
    transcript_id = mcols(mrna)$transcript_id,
    gene_id       = mcols(mrna)$gene_id,
    tx_length     = mcols(mrna)$tx_length,
    stringsAsFactors = FALSE
  ) %>%
    group_by(transcript_id) %>%
    slice_max(tx_length, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Join annotations early to prioritize annotated transcripts
  tx_with_ann <- tx_to_gene %>%
    left_join(ann_gene, by = c("gene_id" = "locusName"))
  
  # Get longest annotated transcript per gene
  longest_annotated <- tx_with_ann %>%
    filter(!is.na(arabi_gene) & arabi_gene != "") %>%
    group_by(gene_id) %>%
    slice_max(tx_length, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Get longest transcript (any) for genes with no annotated transcripts
  longest_any <- tx_to_gene %>%
    filter(!gene_id %in% longest_annotated$gene_id) %>%
    group_by(gene_id) %>%
    slice_max(tx_length, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Combine: annotated first, then unannotated
  longest_tx_per_gene <- bind_rows(
    longest_annotated,
    longest_any %>% left_join(ann_gene, by = c("gene_id" = "locusName"))
  )
} else {
  longest_tx_per_gene <- tibble(gene_id = character(), transcript_id = character(), tx_length = integer())
}

# Map SNPs to genes
hits_gene    <- findOverlaps(snps_gr, genes, ignore.strand = TRUE)
hit_gene_ids <- unique(mcols(genes)$gene_id[subjectHits(hits_gene)])

snp_gene_map <- data.frame(
  SNP_ID  = mcols(snps_gr)$SNP_ID[queryHits(hits_gene)],
  gene_id = mcols(genes)$gene_id[subjectHits(hits_gene)],
  stringsAsFactors = FALSE
)

cat("Total outlier SNPs:", length(snps_gr), "\n")
cat("Genic SNPs:", length(unique(queryHits(hits_gene))), "\n")
cat("Unique genes hit:", length(hit_gene_ids), "\n")

# Create output tables
snp_annotations <- snp_gene_map %>%
  left_join(longest_tx_per_gene, by = "gene_id")

genes_for_GO <- data.frame(gene_id = hit_gene_ids, stringsAsFactors = FALSE) %>%
  left_join(longest_tx_per_gene, by = "gene_id")

arabi_for_GO <- genes_for_GO %>%
  filter(!is.na(arabi_gene) & arabi_gene != "") %>%
  pull(arabi_gene) %>%
  unique()

# Write outputs
write.table(snp_annotations, "outlier_SNPs_annotations.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(genes_for_GO, out_file, 
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(arabidopsis_gene = arabi_for_GO), 
            "outlier_genes_arabi_GO_input.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


#### Plotting pi chart ####

library(ggplot2)

# Create data frame
data <- data.frame(
  category = c("Intergenic", "Genic"),
  count = c(592, 217)
)

# Calculate percentages
data$percentage <- round(data$count/sum(data$count)*100, 1)

# Create pie chart
p <- ggplot(data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.7) +
  coord_polar("y", start = pi) +
  scale_fill_manual(values = c("Genic" = "#999999", "Intergenic" = "#D3D3D3")) +
  theme_void() +
  theme(
    legend.position = "none"
  )

# Display the plot
print(p)

# Save to PNG
ggsave("snp_distribution_pie.png", plot = p, width = 4, height = 4, dpi = 300)


#### Density Plot: Effect Size for Intergenic vs Genic SNPs (BF >= 10) ####

library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)

# Load and filter data
covis.snp.med <- read.table("repMedian_summary_betai_reg.out", 
                            header = TRUE, sep = "\t", 
                            stringsAsFactors = FALSE)

outlier_snps <- covis.snp.med %>%
  filter(BF.dB. >= 10) %>%
  mutate(abs_Beta_is = abs(Beta_is))

# Attach chromosome/position info
chrpos <- read.table("allPopsMAF0.01_LD-50-5-0.2_chr1-11_CHR-POS.tsv", 
                     header = FALSE, sep = "\t", 
                     stringsAsFactors = FALSE,
                     col.names = c("CHROM_true", "POS_true"))

outlier_indices <- which(covis.snp.med$BF.dB. >= 10)
outlier_snps$CHROM_true <- as.character(chrpos$CHROM_true[outlier_indices])
outlier_snps$POS_true <- chrpos$POS_true[outlier_indices]

# Load GFF and classify SNPs
gff <- import("Egrandis_297_v2.0.gene.gff3")
genes <- gff[gff$type == "gene"]

gff_seqlevels <- seqlevels(gff)
gff_chr_num <- suppressWarnings(as.integer(sub(".*?([0-9]+)$", "\\1", gff_seqlevels)))
chr_map <- setNames(gff_seqlevels, gff_chr_num)

snp_seqnames <- chr_map[as.integer(outlier_snps$CHROM_true)]
valid_snps <- !is.na(snp_seqnames)
outlier_snps_valid <- outlier_snps[valid_snps, ]
snp_seqnames_valid <- snp_seqnames[valid_snps]

snps_gr <- GRanges(
  seqnames = snp_seqnames_valid,
  ranges = IRanges(outlier_snps_valid$POS_true, outlier_snps_valid$POS_true)
)

hits_gene <- findOverlaps(snps_gr, genes, ignore.strand = TRUE)
genic_indices <- unique(queryHits(hits_gene))

outlier_snps_valid$category <- "Intergenic"
outlier_snps_valid$category[genic_indices] <- "Genic"
outlier_snps_valid$category <- factor(outlier_snps_valid$category, 
                                      levels = c("Genic", "Intergenic"))

# Wilcoxon rank-sum test
wilcox_result <- wilcox.test(abs_Beta_is ~ category, 
                             data = outlier_snps_valid,
                             exact = FALSE)
cat("Wilcoxon rank-sum test: W =", wilcox_result$statistic, 
    ", p =", format.pval(wilcox_result$p.value, digits = 3), "\n")

# Create density plot
p_density <- ggplot(outlier_snps_valid, aes(x = abs_Beta_is, fill = category, color = category)) +
  geom_density(data = filter(outlier_snps_valid, category == "Intergenic"),
               alpha = 0.8, linewidth = 1.5) +
  geom_density(data = filter(outlier_snps_valid, category == "Genic"),
               alpha = 0.8, linewidth = 1.5) +
  scale_fill_manual(
    values = c("Genic" = "#999999", "Intergenic" = "#D3D3D3"),
    labels = c("Genic" = "Genic SNPs", "Intergenic" = "Intergenic SNPs")
  ) +
  scale_color_manual(
    values = c("Genic" = "#666666", "Intergenic" = "#999999"),
    labels = c("Genic" = "Genic SNPs", "Intergenic" = "Intergenic SNPs")
  ) +
  labs(
    x = expression("|"*beta[iS]*"| (absolute effect size magnitude)"),
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.title.x = element_text(size = 35, margin = margin(t = 15)),
    axis.title.y = element_text(size = 35, margin = margin(r = 15)),
    axis.text.x = element_text(size = 29, family = "Helvetica"),
    axis.text.y = element_text(size = 29, family = "Helvetica"),
    axis.line = element_line(linewidth = 1.2),
    axis.ticks = element_line(linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.position = c(0.24, 0.87),
    legend.title = element_blank(),
    legend.text = element_text(size = 29, family = "Helvetica"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(1, "cm"),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "in"),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )

print(p_density)

ggsave(
  filename = "EffectSize_Genic_vs_Intergenic_density.png",
  plot = p_density,
  width = 11,
  height = 8,
  units = "in",
  dpi = 300
)

