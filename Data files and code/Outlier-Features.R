# Association of outlier windows with genome features (20kb windows)

# This script compares outlier windows with background windows for pi (wet and dry
# species), recombination rate, gene density, GC content, FST and dXY, using a
# permutation test, and plots the results as half-violin plots.
# Input: the *_merged_outliers.csv files (one per species pair)

library(ggplot2)
library(patchwork)

# Set working directory (update this path to your location)
setwd("~/pwd")

set.seed(42)
n_perm <- 10000

# Species pairs (top to bottom in the figures), with the wet and dry species and plot colour
pairs <- data.frame(
  pair    = c("moluc_largiflo", "dealba_vicina", "nitida_regnans"),
  pop_wet = c("moluc", "dealba", "nitida"),
  pop_dry = c("largiflo", "vicina", "regnans"),
  colour  = c("#BC83AC", "#9BA2B6", "#BFA37A")
)

# Outlier definitions (columns in the csv files): top 1% residual outliers and top 10% outliers
outlier_defs <- c("top1pct_residual_overlap", "top10pct_overlap")

# Metrics plotted in each figure
figures <- list(
  pi_rec_genedensity = list(
    list(col = "pi_wet",     label = expression(pi ~ "(wet)")),
    list(col = "pi_dry",     label = expression(pi ~ "(dry)")),
    list(col = "r_cM_Mb",    label = "Recombination rate (cM/Mb)"),
    list(col = "gene_count", label = "Gene density (genes / 20 kb)")
  ),
  fst_dxy_gc = list(
    list(col = "fst",        label = expression(italic(F[ST]))),
    list(col = "dxy",        label = expression(italic(d[XY]))),
    list(col = "gc_content", label = "GC content")
  )
)

# Read each species pair file and create generic pi_wet / pi_dry columns
read_pair <- function(i) {
  data <- read.csv(paste0(pairs$pair[i], "_merged_outliers.csv"), fileEncoding = "UTF-8-BOM")
  data$pi_wet <- data[[paste0("pi_", pairs$pop_wet[i])]]
  data$pi_dry <- data[[paste0("pi_", pairs$pop_dry[i])]]
  data
}
pair_data <- lapply(seq_len(nrow(pairs)), read_pair)

# Function: permutation test. The observed statistic is the median of the outlier windows,
# compared against the medians of n_perm random samples of background windows (same number
# of windows as the outliers). Two-sided empirical P-value.
perm_test <- function(out_vals, bg_vals) {
  obs         <- median(out_vals)
  null        <- replicate(n_perm, median(sample(bg_vals, length(out_vals))))
  null_centre <- mean(null)
  p <- (1 + sum(abs(null - null_centre) >= abs(obs - null_centre))) / (n_perm + 1)
  data.frame(n_outlier = length(out_vals), n_background = length(bg_vals),
             obs_median = obs, null_centre = null_centre,
             background_median = median(bg_vals), p_empirical = p)
}

# Function: half-violin polygon (density curve above a baseline) and a tick at the median
half_violin <- function(vals, baseline, height = 0.42) {
  d    <- density(vals)
  dens <- d$y / max(d$y) * height
  med  <- median(vals)
  list(poly = data.frame(x = c(d$x, rev(d$x)), y = c(baseline + dens, rep(baseline, length(d$x)))),
       med  = data.frame(x = med, y0 = baseline, y1 = baseline + approx(d$x, dens, xout = med)$y))
}

# Function: one panel (one metric, one species pair)
make_panel <- function(data, metric, label, outlier_col, colour) {

  # Keep windows with a value for the metric and an outlier status
  data <- data[is.finite(data[[metric]]) & !is.na(data[[outlier_col]]), ]
  out_vals <- data[[metric]][data[[outlier_col]] == 1]
  bg_vals  <- data[[metric]][data[[outlier_col]] == 0]

  # Permutation test
  stats <- perm_test(out_vals, bg_vals)
  p     <- stats$p_empirical
  stars <- if (p < 0.001) " ***" else if (p < 0.01) " **" else if (p < 0.05) " *" else ""
  p_txt <- if (p < 1e-4) "p < 0.0001" else sprintf("p = %.4f", p)
  stats_label <- sprintf("N = %d, %s%s", length(out_vals), p_txt, stars)

  # Violins: outliers (top) and an equal-sized random sample of background windows (bottom)
  v_bg  <- half_violin(sample(bg_vals, length(out_vals)), baseline = 1)
  v_out <- half_violin(out_vals, baseline = 1.65)
  x_ann <- min(v_bg$poly$x, v_out$poly$x) +
    0.55 * (max(v_bg$poly$x, v_out$poly$x) - min(v_bg$poly$x, v_out$poly$x))

  plot <- ggplot() +
    geom_polygon(data = v_bg$poly,  aes(x = x, y = y), fill = "grey95", colour = "black", linewidth = 0.5) +
    geom_polygon(data = v_out$poly, aes(x = x, y = y), fill = colour,   colour = "black", linewidth = 0.5) +
    geom_segment(data = rbind(v_bg$med, v_out$med), aes(x = x, xend = x, y = y0, yend = y1), linewidth = 0.5) +
    annotate("text", x = x_ann, y = 1.5, label = stats_label, hjust = 0, size = 2.6, colour = "grey30") +
    scale_y_continuous(breaks = c(1, 1.65) + 0.21, labels = c("BG", "Outlier")) +
    labs(x = label, y = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.ticks.y = element_blank(), plot.margin = margin(20, 4, 20, 4))

  list(plot = plot, stats = stats)
}

# Make the figures (one PDF per figure and outlier definition) and collect the test results
results <- list()

for (outlier_col in outlier_defs) {
  for (fig_name in names(figures)) {

    rows <- list()
    for (i in seq_len(nrow(pairs))) {
      panels <- list()
      for (m in figures[[fig_name]]) {
        res <- make_panel(pair_data[[i]], m$col, m$label, outlier_col, pairs$colour[i])
        panels[[length(panels) + 1]] <- res$plot
        results[[length(results) + 1]] <- cbind(species_pair = pairs$pair[i], outlier_def = outlier_col,
                                                metric = m$col, res$stats)
      }
      rows[[i]] <- wrap_plots(panels, nrow = 1)
    }

    pdf(paste0(fig_name, "_violinplots_", outlier_col, ".pdf"),
        width = ifelse(fig_name == "pi_rec_genedensity", 13, 11), height = 6)
    print(wrap_plots(rows, ncol = 1))
    dev.off()
  }
}

# Save the permutation test results
results <- do.call(rbind, results)
write.table(results, "outlier_vs_background_permutation.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
