# Correlations of genome-wide landscapes among species (20kb windows)

# This script calculates Pearson and Spearman correlations of FST, dXY, pi and
# Tajima's D landscapes across non-overlapping 20kb windows, between the three
# eucalyptus species pairs (input: the merged files from CombiningFiles-allThreePairs.R)

# Set working directory (update this path to your location)
setwd("~/pwd")

# Function: read a species pair file and keep the columns needed for the correlations
read_pair <- function(filename, pair_name, pop1, pop2) {
  
  # Read data
  data <- read.csv(filename, check.names = TRUE)
  
  # Identify window position column (window start)
  window_col <- grep("^window_pos_", names(data), value = TRUE)[1]
  
  # Filter to chromosomes 1-11
  data <- subset(data, chromosome %in% sprintf("Chr%02d", 1:11))
  
  # Remove the extreme dXY/FST outlier window in nitida_regnans (Chr08: 25760001-25780000)
  if (pair_name == "nitida_regnans") {
    outlier_idx <- which(data$chromosome == "Chr08" & data[[window_col]] == 25760001)
    data$fst[outlier_idx] <- NA
    data$dxy[outlier_idx] <- NA
  }
  
  # Keep windows with FST, dXY and pi for both populations
  # (negative FST values are retained; windows with missing FST are dropped)
  pi_cols <- paste0("pi_", c(pop1, pop2))
  data <- data[complete.cases(data[, c("fst", "dxy", pi_cols)]), ]
  
  # Window ID used to match windows among the three species pairs
  out <- data.frame(window = paste(data$chromosome, data[[window_col]], sep = "_"))
  out[[paste0("fst_", pair_name)]] <- data$fst
  out[[paste0("dxy_", pair_name)]] <- data$dxy
  for (pop in c(pop1, pop2)) {
    out[[paste0("pi_", pop)]]   <- data[[paste0("pi_", pop)]]
    out[[paste0("TajD_", pop)]] <- data[[paste0("TajD_", pop)]]
  }
  out
}

# Function: Pearson and Spearman correlation for one comparison
cor_row <- function(x, y, statistic, comparison) {
  pe <- cor.test(x, y, method = "pearson")
  sp <- suppressWarnings(cor.test(x, y, method = "spearman"))
  data.frame(
    Statistic = statistic,
    Comparison = comparison,
    Pearson_r = round(unname(pe$estimate), 2),
    Pearson_p = signif(pe$p.value, 3),
    Spearman_rho = round(unname(sp$estimate), 2),
    Spearman_p = signif(sp$p.value, 3)
  )
}

# Function: correlations between all combinations of the columns in 'cols'
cor_all_combinations <- function(data, cols, statistic, sep = " vs ") {
  rows <- combn(cols, 2, simplify = FALSE)
  do.call(rbind, lapply(rows, function(r) {
    cor_row(data[[r[1]]], data[[r[2]]], statistic, paste(r, collapse = sep))
  }))
}

# Read the three species pairs
dealba_vicina  <- read_pair("dealba_vicina_merged.csv",  "dealba_vicina",  "dealba", "vicina")
moluc_largiflo <- read_pair("moluc_largiflo_merged.csv", "moluc_largiflo", "moluc",  "largiflo")
nitida_regnans <- read_pair("nitida_regnans_merged.csv", "nitida_regnans", "nitida", "regnans")

# Keep windows present in all three species pairs
all_data <- Reduce(function(x, y) merge(x, y, by = "window"),
                   list(dealba_vicina, moluc_largiflo, nitida_regnans))
cat("Number of 20kb windows:", nrow(all_data), "\n")

# Correlations among species pairs (FST and dXY)
pairs <- c("dealba_vicina", "moluc_largiflo", "nitida_regnans")
fst_cors <- cor_all_combinations(all_data, paste0("fst_", pairs), "FST")
dxy_cors <- cor_all_combinations(all_data, paste0("dxy_", pairs), "dXY")

# Correlations among species (pi and Tajima's D)
species <- c("dealba", "vicina", "moluc", "largiflo", "nitida", "regnans")
pi_cors   <- cor_all_combinations(all_data, paste0("pi_", species), "Nucleotide diversity (pi)")
tajd_cors <- cor_all_combinations(all_data, paste0("TajD_", species), "Tajima's D")

# Combine and save
cor_table <- rbind(fst_cors, dxy_cors, pi_cors, tajd_cors)
write.csv(cor_table, "landscape_correlations_20kb.csv", row.names = FALSE)
print(cor_table)