# Genome-wide landscape plots: FST, dXY, π, and Tajima's D

# This script creates stacked landscape plots for population genomic statistics
# across chromosomes for three eucalyptus species pairs

# Load required libraries
library(ggplot2)
library(gtable)
library(grid)

# Set working directory (update this path to your location)
setwd("~/pwd")

# Function: Plot genome-wide landscapes for a species pair
plot_landscape <- function(filename, pair_name) {
  
  # Read data
  data <- read.csv(filename, check.names = TRUE)
  
  # Identify window position column
  window_pos_cols <- grep("^window_pos_", names(data), value = TRUE)
  window_col <- if (length(window_pos_cols) >= 2) window_pos_cols[2] else window_pos_cols[1]
  cat("Using window position column:", window_col, "\n")
  
  # Filter to chromosomes 1-11
  chrom_filter <- sprintf("Chr%02d", 1:11)
  data <- subset(data, chromosome %in% chrom_filter)
  
  # Order chromosomes and sort data
  chrOrder <- sort(unique(data$chromosome))
  data$chrOrder <- factor(data$chromosome, levels = chrOrder)
  data <- data[order(data$chrOrder, data[[window_col]]), ]
  
  # Calculate cumulative positions for x-axis
  chr_lengths <- tapply(data[[window_col]], data$chrOrder, max, na.rm = TRUE)
  offsets <- c(0, cumsum(chr_lengths)[-length(chr_lengths)])
  names(offsets) <- levels(data$chrOrder)
  data$cum_pos <- data[[window_col]] + offsets[as.character(data$chrOrder)]
  
  # Calculate chromosome boundaries and midpoints for plotting
  chr_boundaries <- cumsum(chr_lengths)[-length(chr_lengths)]
  chr_midpoints <- tapply(data$cum_pos, data$chrOrder, function(x) mean(range(x)))
  midpoints_df <- data.frame(chromosome = 1:11, midpoint = chr_midpoints)
  
  # Ensure FST is non-negative
  data$fst <- pmax(data$fst, 0, na.rm = TRUE)
  
  # Define species-specific names and colors
  if (pair_name == "dealba_vicina") {
    pop1 <- "dealba"; pop2 <- "vicina"
    pop1_color <- "#456FAF"; pop2_color <- "#D77C3E"
  } else if (pair_name == "moluc_largiflo") {
    pop1 <- "moluc"; pop2 <- "largiflo"
    pop1_color <- "#456FAF"; pop2_color <- "#D77C3E"
  } else if (pair_name == "nitida_regnans") {
    pop1 <- "nitida"; pop2 <- "regnans"
    pop1_color <- "#456FAF"; pop2_color <- "#D77C3E"
  } else {
    pop1 <- "pop1"; pop2 <- "pop2"
    pop1_color <- "#456FAF"; pop2_color <- "#D77C3E"
  }
  
  # Identify the actual column names for pi and TajD
  pi_col1 <- paste0("pi_", pop1)
  pi_col2 <- paste0("pi_", pop2)
  tajd_col1 <- paste0("TajD_", pop1)
  tajd_col2 <- paste0("TajD_", pop2)
  
  # Smooth data using 1.5 Mb sliding windows
  smooth_window_size <- 1500000
  smoothed_list <- list()
  
  for (chrom in unique(data$chromosome)) {
    chr_data <- data[data$chromosome == chrom, ]
    chr_data <- chr_data[order(chr_data[[window_col]]), ]
    chr_data$smooth_bin <- floor(chr_data[[window_col]] / smooth_window_size) * smooth_window_size
    
    # Create formula dynamically with actual column names
    formula_str <- paste0("cbind(", pi_col1, ", ", pi_col2, ", ", 
                          tajd_col1, ", ", tajd_col2, ", dxy) ~ smooth_bin + chromosome")
    
    smooth_chr <- aggregate(
      as.formula(formula_str),
      data = chr_data, FUN = mean, na.rm = TRUE
    )
    
    smooth_chr$chrOrder <- factor(smooth_chr$chromosome, levels = levels(data$chrOrder))
    smooth_chr$cum_pos <- smooth_chr$smooth_bin + offsets[as.character(smooth_chr$chrOrder)]
    smoothed_list[[chrom]] <- smooth_chr
  }
  
  smoothed_data <- do.call(rbind, smoothed_list)
  
  # Rename smoothed columns to generic names for plotting
  names(smoothed_data)[names(smoothed_data) == pi_col1] <- "pi_pop1_smooth"
  names(smoothed_data)[names(smoothed_data) == pi_col2] <- "pi_pop2_smooth"
  names(smoothed_data)[names(smoothed_data) == tajd_col1] <- "TajD_pop1_smooth"
  names(smoothed_data)[names(smoothed_data) == tajd_col2] <- "TajD_pop2_smooth"
  names(smoothed_data)[names(smoothed_data) == "dxy"] <- "dxy_smooth"
  
  # Create FST landscape plot
  fst_plot <- ggplot(data, aes(cum_pos, fst, color = chrOrder)) +
    geom_point(alpha = 1, size = 2) +
    labs(x = NULL, y = expression(italic(F)[ST])) +
    scale_color_manual(values = rep(c("#3F8E5C", "#7CBF8B"), length.out = 11)) +
    scale_x_continuous(
      breaks = c(0, chr_boundaries, max(data$cum_pos, na.rm = TRUE)),
      labels = NULL, expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max(data$cum_pos, na.rm = TRUE)), ylim = c(0, 1.02), expand = FALSE) +
    theme_classic() +
    theme(
      text = element_text(family = "Helvetica"), 
      legend.position = "none",
      axis.title.y = element_text(size = 32, family = "Helvetica", margin = margin(r = 15)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 26, family = "Helvetica", margin = margin(r = 10)),
      axis.ticks.x = element_line(linewidth = 1.2), 
      axis.ticks.y = element_line(linewidth = 1.2),
      axis.line = element_line(linewidth = 1.2),
      axis.ticks.length = unit(0.3, "cm"),
      plot.margin = margin(0.3, 0.5, 0.1, 0.5, "in"),
      aspect.ratio = 1/12
    )
  
  # Create dXY landscape plot
  dxy_plot <- ggplot(smoothed_data, aes(cum_pos, dxy_smooth, group = chrOrder)) +
    geom_line(linewidth = 0.8, color = "black") +
    labs(x = NULL, y = expression(italic(d)[XY])) +
    scale_x_continuous(
      breaks = c(0, chr_boundaries, max(data$cum_pos, na.rm = TRUE)),
      labels = NULL, expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0.01, 0.045), breaks = c(0.01, 0.02, 0.03, 0.04), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max(data$cum_pos, na.rm = TRUE)), expand = FALSE) +
    theme_classic() +
    theme(
      text = element_text(family = "Helvetica"), 
      legend.position = "none",
      axis.title.y = element_text(size = 32, family = "Helvetica", margin = margin(r = 15)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 26, family = "Helvetica", margin = margin(r = 10)),
      axis.ticks.x = element_line(linewidth = 1.2), 
      axis.ticks.y = element_line(linewidth = 1.2),
      axis.line = element_line(linewidth = 1.2),
      axis.ticks.length = unit(0.3, "cm"),
      plot.margin = margin(0.1, 0.5, 0.1, 0.5, "in"),
      aspect.ratio = 1/12
    )
  
  # Create π landscape plot
  pi_plot <- ggplot() +
    geom_line(data = smoothed_data, aes(x = cum_pos, y = pi_pop2_smooth, color = "pop2", group = chrOrder), linewidth = 0.8) +
    geom_line(data = smoothed_data, aes(x = cum_pos, y = pi_pop1_smooth, color = "pop1", group = chrOrder), linewidth = 0.8) +
    labs(x = NULL, y = expression(pi)) +
    scale_color_manual(values = c("pop1" = pop1_color, "pop2" = pop2_color),
                       labels = c("pop1" = pop1, "pop2" = pop2), name = "Population") +
    scale_x_continuous(
      breaks = c(0, chr_boundaries, max(data$cum_pos, na.rm = TRUE)),
      labels = NULL, expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(0.005, 0.045), breaks = c(0.01, 0.02, 0.03, 0.04), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max(data$cum_pos, na.rm = TRUE)), expand = FALSE) +
    theme_classic() +
    theme(
      text = element_text(family = "Helvetica"), 
      legend.position = "none",
      axis.title.y = element_text(size = 32, family = "Helvetica", margin = margin(r = 15)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 26, family = "Helvetica", margin = margin(r = 10)),
      axis.ticks.x = element_line(linewidth = 1.2), 
      axis.ticks.y = element_line(linewidth = 1.2),
      axis.line = element_line(linewidth = 1.2),
      axis.ticks.length = unit(0.3, "cm"),
      plot.margin = margin(0.1, 0.5, 0.1, 0.5, "in"),
      aspect.ratio = 1/12
    )
  
  # Create boundary data for chromosome labels
  boundary_data <- data.frame(
    x = c(0, chr_boundaries, max(data$cum_pos, na.rm = TRUE)),
    label = ""
  )
  
  # Create Tajima's D landscape plot with increased bottom margin
  tajd_plot <- ggplot() +
    geom_line(data = smoothed_data, aes(x = cum_pos, y = TajD_pop2_smooth, color = "pop2", group = chrOrder), linewidth = 0.8) +
    geom_line(data = smoothed_data, aes(x = cum_pos, y = TajD_pop1_smooth, color = "pop1", group = chrOrder), linewidth = 0.8) +
    labs(x = "Chromosome", y = "Tajima's D") +
    scale_color_manual(values = c("pop1" = pop1_color, "pop2" = pop2_color),
                       labels = c("pop1" = pop1, "pop2" = pop2), name = "Population") +
    scale_x_continuous(
      breaks = c(0, chr_boundaries, max(data$cum_pos, na.rm = TRUE)),
      labels = boundary_data$label, expand = c(0, 0)
    ) +
    scale_y_continuous(limits = c(-1.5, 0.8), breaks = c(-1, -0.5, 0, 0.5), expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max(data$cum_pos, na.rm = TRUE)), ylim = c(-1.1, 0.8), expand = FALSE, clip = "off") +
    annotate("text", x = midpoints_df$midpoint, y = -1.35, label = midpoints_df$chromosome, size = 9, family = "Helvetica") +
    theme_classic() +
    theme(
      text = element_text(family = "Helvetica"), 
      legend.position = "none",
      axis.title.x = element_text(size = 32, margin = margin(t = 30), family = "Helvetica"),
      axis.title.y = element_text(size = 32, family = "Helvetica", margin = margin(r = 15)),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 26, family = "Helvetica", margin = margin(r = 10)),
      axis.ticks.x = element_line(linewidth = 1.2), 
      axis.ticks.y = element_line(linewidth = 1.2),
      axis.line = element_line(linewidth = 1.2),
      axis.ticks.length = unit(0.3, "cm"),
      plot.margin = margin(0.02, 0.5, 0.5, 0.5, "in"),
      aspect.ratio = 1/12
    )
  
  # Convert plots to grobs and align widths
  g1 <- ggplotGrob(fst_plot)
  g2 <- ggplotGrob(dxy_plot)
  g3 <- ggplotGrob(pi_plot)
  g4 <- ggplotGrob(tajd_plot)
  
  maxWidthLandscape <- unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths)
  g1$widths <- maxWidthLandscape
  g2$widths <- maxWidthLandscape
  g3$widths <- maxWidthLandscape
  g4$widths <- maxWidthLandscape
  
  # Stack all plots vertically using rbind from gtable
  landscape_plot <- rbind(g1, g2, g3, g4)
  
  # Save the combined landscape plot
  ggsave(paste0(pair_name, "_landscape.png"),
         plot = landscape_plot, device = "png", dpi = 300, width = 22, height = 16, units = "in")
  
  cat(paste0("Saved landscape plot for ", pair_name, "\n"))
}

# Run for all three species pairs
plot_landscape("dealba_vicina_merged.csv", "dealba_vicina")
plot_landscape("moluc_largiflo_merged.csv", "moluc_largiflo")
plot_landscape("nitida_regnans_merged.csv", "nitida_regnans")
