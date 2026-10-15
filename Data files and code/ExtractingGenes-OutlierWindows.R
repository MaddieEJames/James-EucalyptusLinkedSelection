#!/usr/bin/env Rscript
# Extract genes overlapping genomic outlier windows and map to Arabidopsis orthologs

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
  library(tidyr)
})

# Configuration
MIN_OVERLAP_FRACTION <- 0.25  # 0 = any overlap, 0.5 = 50% of gene, 1.0 = entire gene

# Input files
GFF_FILE <- "Egrandis_297_v2.0.gene.gff3"
ANN_FILE <- "Egrandis_297_v2.0.annotation_info.txt"

# Species pairs to process
SPECIES_PAIRS <- c("dealba_vicina", "moluc_largiflo", "nitida_regnans")

# Helper function: normalize gene/transcript IDs
normalize_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.v[0-9]+\\.[0-9]+$", "", x)
  x <- sub("\\.p$", "", x)
  x <- sub("\\.t$", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x
}

# Helper function: strip Arabidopsis isoform suffixes
strip_arabi_isoform <- function(x) {
  sub("\\.[0-9]+$", "", x)
}

# Helper function: map chromosome names to GFF format
map_chromosome_to_gff <- function(chr_vec, gff_seqlevels) {
  gff_chr_num <- suppressWarnings(as.integer(sub(".*?([0-9]+)$", "\\1", gff_seqlevels)))
  chr_map <- setNames(gff_seqlevels, gff_chr_num)
  chr_numeric <- suppressWarnings(as.integer(sub("Chr", "", chr_vec)))
  mapped <- chr_map[chr_numeric]
  
  if (any(is.na(mapped))) {
    warning("Some chromosomes could not be mapped to GFF seqlevels")
  }
  
  return(mapped)
}

# Load reference data
cat("\nLoading reference data...\n")

# Load GFF annotation
cat("Loading GFF file:", GFF_FILE, "\n")
gff <- import(GFF_FILE)
genes <- gff[gff$type == "gene"]

# Normalize gene IDs
genes$gene_id <- if (!is.null(mcols(genes)$Name)) {
  normalize_id(mcols(genes)$Name)
} else {
  normalize_id(mcols(genes)$ID)
}

cat("Total genes in GFF:", length(genes), "\n")

# Load Arabidopsis annotation
cat("Loading annotation file:", ANN_FILE, "\n")
ann <- read.delim(ANN_FILE, header = TRUE, sep = "\t", quote = "", 
                  check.names = FALSE, stringsAsFactors = FALSE)
ann$locusName <- trimws(as.character(ann$locusName))

# Create gene-to-Arabidopsis mapping
ann_gene <- ann %>%
  select(locusName, `Best-hit-arabi-name`) %>%
  filter(!is.na(`Best-hit-arabi-name`) & `Best-hit-arabi-name` != "") %>%
  mutate(arabi_gene_stripped = strip_arabi_isoform(`Best-hit-arabi-name`)) %>%
  distinct(locusName, arabi_gene_stripped) %>%
  group_by(locusName) %>%
  summarize(
    arabi_gene = paste(sort(unique(arabi_gene_stripped)), collapse = ";"),
    .groups = "drop"
  )

cat("Genes with Arabidopsis orthologs:", nrow(ann_gene), "\n")

# Get chromosome mapping
gff_seqlevels <- seqlevels(gff)

# Function to process outlier windows
process_outlier_windows <- function(data, outlier_column, species_pair) {
  
  cat("\nProcessing", outlier_column, "for", species_pair, "\n")
  
  # Filter for outlier windows (value = 1)
  windows <- data %>%
    filter(!!sym(outlier_column) == 1)
  
  if (nrow(windows) == 0) {
    cat("  No outlier windows found\n")
    return(NULL)
  }
  
  cat("  Number of outlier windows:", nrow(windows), "\n")
  
  # Map chromosome names to GFF format
  windows$gff_seqname <- map_chromosome_to_gff(windows$chromosome, gff_seqlevels)
  
  # Remove windows that couldn't be mapped
  valid_windows <- windows[!is.na(windows$gff_seqname), ]
  
  if (nrow(valid_windows) == 0) {
    cat("  ERROR: No valid windows after chromosome mapping\n")
    return(NULL)
  }
  
  # Ensure window coordinates are valid integers
  valid_windows <- valid_windows %>%
    mutate(
      window_pos_1 = as.integer(window_pos_1),
      window_pos_2 = as.integer(window_pos_2)
    ) %>%
    filter(!is.na(window_pos_1) & !is.na(window_pos_2) & window_pos_1 <= window_pos_2)
  
  if (nrow(valid_windows) == 0) {
    cat("  ERROR: No valid window coordinates\n")
    return(NULL)
  }
  
  # Create GRanges object for windows
  windows_gr <- GRanges(
    seqnames = valid_windows$gff_seqname,
    ranges = IRanges(start = valid_windows$window_pos_1, 
                     end = valid_windows$window_pos_2)
  )
  
  # Find overlaps between windows and genes
  overlaps <- findOverlaps(windows_gr, genes, ignore.strand = TRUE)
  
  if (length(overlaps) == 0) {
    cat("  No genes overlap these windows\n")
    return(NULL)
  }
  
  # Calculate overlap statistics
  window_indices <- queryHits(overlaps)
  gene_indices <- subjectHits(overlaps)
  
  overlap_ranges <- pintersect(windows_gr[window_indices], genes[gene_indices])
  overlap_widths <- width(overlap_ranges)
  gene_widths <- width(genes[gene_indices])
  overlap_fractions <- overlap_widths / gene_widths
  
  # Filter by minimum overlap threshold
  pass_threshold <- overlap_fractions >= MIN_OVERLAP_FRACTION
  
  cat("  Genes passing", MIN_OVERLAP_FRACTION * 100, "% threshold:", 
      length(unique(gene_indices[pass_threshold])), "\n")
  
  if (sum(pass_threshold) == 0) {
    cat("  No genes pass the overlap threshold\n")
    return(NULL)
  }
  
  # Create result data frame
  result <- data.frame(
    window_chromosome = as.character(seqnames(windows_gr[window_indices[pass_threshold]])),
    window_pos_1 = start(windows_gr[window_indices[pass_threshold]]),
    window_pos_2 = end(windows_gr[window_indices[pass_threshold]]),
    gene_id = genes$gene_id[gene_indices[pass_threshold]],
    gene_start = start(genes[gene_indices[pass_threshold]]),
    gene_end = end(genes[gene_indices[pass_threshold]]),
    gene_length = gene_widths[pass_threshold],
    overlap_length = overlap_widths[pass_threshold],
    overlap_fraction = round(overlap_fractions[pass_threshold], 3),
    stringsAsFactors = FALSE
  )
  
  # Add window statistics (FST and dXY residuals)
  windows_for_merge <- valid_windows %>%
    mutate(
      window_pos_1 = as.integer(window_pos_1),
      window_pos_2 = as.integer(window_pos_2)
    ) %>%
    select(gff_seqname, window_pos_1, window_pos_2, 
           any_of(c("fst_residual", "dxy_residual")))
  
  result <- result %>%
    left_join(windows_for_merge,
              by = c("window_chromosome" = "gff_seqname",
                     "window_pos_1" = "window_pos_1",
                     "window_pos_2" = "window_pos_2"),
              relationship = "many-to-many")
  
  # Add Arabidopsis ortholog information
  result <- result %>%
    left_join(ann_gene, by = c("gene_id" = "locusName")) %>%
    replace_na(list(arabi_gene = ""))
  
  # Reorder columns
  stat_cols <- intersect(c("fst_residual", "dxy_residual"), names(result))
  base_cols <- c("window_chromosome", "window_pos_1", "window_pos_2")
  gene_cols <- c("gene_id", "gene_start", "gene_end", "gene_length",
                 "overlap_length", "overlap_fraction", "arabi_gene")
  
  result <- result %>%
    select(all_of(base_cols), any_of(stat_cols), all_of(gene_cols)) %>%
    arrange(window_chromosome, window_pos_1, gene_start)
  
  # Determine output file prefix based on outlier type
  if (outlier_column == "top1pct_residual_overlap") {
    output_prefix <- paste0(species_pair, "_top1pct_overlap_FST_dXY_residuals")
  } else if (outlier_column == "top10pct_overlap") {
    output_prefix <- paste0(species_pair, "_top10pct_overlap")
  }
  
  # Write gene-level output with Arabidopsis
  genes_file <- paste0(output_prefix, "_genes_arabidopsis.tsv")
  write.table(result, genes_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Output written:", genes_file, "\n")
  
  # Write Arabidopsis genes for GO analysis
  # For each Eucalyptus gene, randomly select one Arabidopsis ortholog
  # This prevents counting the same Eucalyptus gene multiple times in GO enrichment
  set.seed(123)  # For reproducibility
  arabi_genes <- result %>%
    filter(arabi_gene != "") %>%
    select(gene_id, arabi_gene) %>%
    distinct() %>%
    mutate(arabi_gene_split = strsplit(arabi_gene, ";")) %>%
    rowwise() %>%
    mutate(selected_arabi = sample(arabi_gene_split, 1)) %>%
    ungroup() %>%
    pull(selected_arabi) %>%
    sort()
  
  if (length(arabi_genes) > 0) {
    go_file <- paste0(output_prefix, "_arabidopsis_GO_input.txt")
    write.table(data.frame(arabidopsis_gene = arabi_genes),
                go_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("  Arabidopsis genes for GO:", go_file, "(", length(arabi_genes), "genes )\n")
  }
  
  return(result)
}

# Process each species pair
cat("\n" , rep("=", 70), "\n", sep = "")
cat("Processing species pairs\n")
cat(rep("=", 70), "\n", sep = "")

for (species_pair in SPECIES_PAIRS) {
  
  input_file <- paste0(species_pair, "_merged_outliers.csv")
  
  if (!file.exists(input_file)) {
    cat("\nWARNING: File not found:", input_file, "- skipping\n")
    next
  }
  
  cat("\n", rep("-", 70), "\n", sep = "")
  cat("Species pair:", species_pair, "\n")
  cat(rep("-", 70), "\n", sep = "")
  
  # Read merged data file
  data <- read.csv(input_file, stringsAsFactors = FALSE)
  cat("Loaded", nrow(data), "windows from", input_file, "\n")
  
  # Process top 1% residual overlap outliers
  if ("top1pct_residual_overlap" %in% names(data)) {
    process_outlier_windows(data, "top1pct_residual_overlap", species_pair)
  } else {
    cat("\nWARNING: Column 'top1pct_residual_overlap' not found\n")
  }
  
  # Process top 10% overlap outliers
  if ("top10pct_overlap" %in% names(data)) {
    process_outlier_windows(data, "top10pct_overlap", species_pair)
  } else {
    cat("\nWARNING: Column 'top10pct_overlap' not found\n")
  }
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("Analysis complete\n")
cat(rep("=", 70), "\n", sep = "")
