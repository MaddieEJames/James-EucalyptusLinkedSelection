# Genomic islands of differentiation: island features compared with background (20kb windows)

# This script compares each FST island with background windows (all windows not in any island)
# for pi (wet and dry species), dXY, recombination rate, gene density and GC content, using
# two-sided Wilcoxon rank-sum tests with Benjamini-Hochberg correction (within each feature and
# species pair). It then plots the proportion of islands lower than, not different from, or higher
# than the background.
# Input: the *_merged_outliers.csv files (one per species pair)

library(ggplot2)

# Set working directory (update this path to your location)
setwd("~/pwd")

# Species pairs (top to bottom in the figure) with the wet and dry species
pairs <- data.frame(
  pair    = c("moluc_largiflo", "dealba_vicina", "nitida_regnans"),
  pop_wet = c("moluc", "dealba", "nitida"),
  pop_dry = c("largiflo", "vicina", "regnans"),
  label   = c("E. moluccana\nE. largiflorens", "E. dealbata\nE. vicina", "E. nitida\nE. regnans")
)

# Features compared (column name = display name)
features <- c(pi_wet = "π (wet)", pi_dry = "π (dry)", dxy = "dXY",
              r_cM_Mb = "Recombination", gene_count = "Gene density", gc_content = "GC content")

# FST islands (Table S11)
islands <- data.frame(
  pair = c(rep("dealba_vicina", 2), rep("moluc_largiflo", 11), rep("nitida_regnans", 10)),
  chromosome = c("Chr05", "Chr06",
                 "Chr02", "Chr02", "Chr03", "Chr04", "Chr04", "Chr04", "Chr05", "Chr08", "Chr10", "Chr10", "Chr11",
                 "Chr01", "Chr01", "Chr02", "Chr02", "Chr03", "Chr03", "Chr04", "Chr06", "Chr10", "Chr11"),
  start = c(900001, 30260001,
            23320001, 55280001, 40060001, 4760001, 22300001, 39520001, 9880001, 5340001, 9720001, 14080001, 26420001,
            22060001, 23040001, 23240001, 52600001, 39940001, 40600001, 22200001, 43580001, 14800001, 27400001),
  end = c(1580000, 30740000,
          24460000, 55560000, 42500000, 5240000, 22620000, 39880000, 10280000, 5720000, 10000000, 15320000, 26760000,
          22380000, 23500000, 23660000, 53020000, 40220000, 40920000, 22560000, 43960000, 15320000, 28060000)
)

# Test each island against the background, for each species pair and feature
results <- list()

for (i in seq_len(nrow(pairs))) {

  # Read data and create generic pi_wet / pi_dry columns
  data <- read.csv(paste0(pairs$pair[i], "_merged_outliers.csv"), fileEncoding = "UTF-8-BOM")
  data$pi_wet <- data[[paste0("pi_", pairs$pop_wet[i])]]
  data$pi_dry <- data[[paste0("pi_", pairs$pop_dry[i])]]

  # Flag the windows in each island (NA = background window)
  isl <- islands[islands$pair == pairs$pair[i], ]
  data$island <- NA_integer_
  for (k in seq_len(nrow(isl))) {
    in_island <- data$chromosome == isl$chromosome[k] &
      data$window_pos_1 <= isl$end[k] & data$window_pos_2 >= isl$start[k]
    data$island[in_island] <- k
  }

  for (f in names(features)) {

    # Background windows: all windows not in any island
    bg_vals <- data[[f]][is.na(data$island)]
    bg_vals <- bg_vals[is.finite(bg_vals)]

    # Wilcoxon rank-sum test for each island
    res <- do.call(rbind, lapply(seq_len(nrow(isl)), function(k) {
      island_vals <- data[[f]][which(data$island == k)]
      island_vals <- island_vals[is.finite(island_vals)]
      p <- if (length(island_vals) >= 2) suppressWarnings(wilcox.test(island_vals, bg_vals)$p.value) else NA
      data.frame(species_pair = pairs$pair[i], island = k, chromosome = isl$chromosome[k],
                 start = isl$start[k], end = isl$end[k], feature = f,
                 island_median = median(island_vals), background_median = median(bg_vals), p = p)
    }))

    # Benjamini-Hochberg correction across islands (within feature and species pair)
    res$p_BH <- p.adjust(res$p, method = "BH")
    results[[length(results) + 1]] <- res
  }
}

results <- do.call(rbind, results)

# Classify each island: significantly lower, not significant, or significantly higher than background
results$class <- "ns"
results$class[!is.na(results$p_BH) & results$p_BH < 0.05 & results$island_median < results$background_median] <- "lower"
results$class[!is.na(results$p_BH) & results$p_BH < 0.05 & results$island_median > results$background_median] <- "higher"
write.table(results, "island_vs_background_wilcoxon.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Proportion of islands in each class, for each species pair and feature
counts <- as.data.frame(table(species_pair = results$species_pair, feature = results$feature,
                              class = factor(results$class, levels = c("lower", "ns", "higher"))),
                        responseName = "n")
counts$total <- ave(counts$n, counts$species_pair, counts$feature, FUN = sum)
counts$prop  <- counts$n / counts$total
write.table(counts, "island_summary_proportions.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Order species pairs (first pair at the top) and features, and add display labels
counts$species_pair <- factor(counts$species_pair, levels = rev(pairs$pair), labels = rev(pairs$label))
counts$feature      <- factor(counts$feature, levels = names(features), labels = features)
n_labels <- unique(counts[, c("species_pair", "feature", "total")])

# Plot: proportion of islands per species pair and feature
plot <- ggplot(counts, aes(x = prop, y = species_pair, fill = class)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.6, position = position_stack(reverse = TRUE)) +
  geom_text(data = n_labels, aes(x = 1.03, y = species_pair, label = paste0("n=", total)),
            inherit.aes = FALSE, hjust = 0, size = 2.2, colour = "grey40") +
  facet_wrap(~ feature, nrow = 1) +
  scale_fill_manual(values = c(lower = "#78A3AB", ns = "#E8E8E4", higher = "#D8BC7A"),
                    labels = c("Lower than background", "Not significant", "Higher than background"),
                    name = NULL) +
  scale_x_continuous(name = "Proportion of islands", breaks = c(0, 0.5, 1), limits = c(0, 1),
                     expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        axis.title.y = element_blank(), axis.text.y = element_text(size = 8, face = "italic"),
        panel.spacing = unit(22, "pt"))

pdf("island_summary_proportion_bars.pdf", width = 15, height = 3.2)
print(plot)
dev.off()
