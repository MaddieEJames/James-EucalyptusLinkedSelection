#### Code to merge pop gen statistics into a single file per species pair ####

# Merge PI, Dxy, Fst, and Tajima's D for three species pairs

library(tidyverse)

# Read PI file
pi <- read_tsv("20000_pi.txt") %>%
  select(pop, chromosome, window_pos_1, window_pos_2, avg_pi)

# Read Dxy file
dxy <- read_tsv("20000_dxy.txt") %>%
  select(pop1, pop2, chromosome, window_pos_1, window_pos_2, avg_dxy)

# Read Fst file
fst <- read_tsv("20000_fst.txt") %>%
  select(pop1, pop2, chromosome, window_pos_1, window_pos_2, avg_wc_fst)

# Read Tajima's D files and add species column
dealba_tajd <- read_tsv("dealba_tajimasD_20000.Tajima.D.txt") %>% mutate(species = "dealba")
largiflo_tajd <- read_tsv("largiflo_tajimasD_20000.Tajima.D.txt") %>% mutate(species = "largiflo")
moluc_tajd <- read_tsv("moluc_tajimasD_20000.Tajima.D.txt") %>% mutate(species = "moluc")
vicina_tajd <- read_tsv("vicina_tajimasD_20000.Tajima.D.txt") %>% mutate(species = "vicina")
nitida_tajd <- read_tsv("nitida_tajimasD_20000.Tajima.D.txt") %>% mutate(species = "nitida")
regnans_tajd <- read_tsv("regnans_tajimasD_20000.Tajima.D.txt") %>% mutate(species = "regnans")

# Combine all Tajima's D data and shift coordinates
tajd <- bind_rows(dealba_tajd, largiflo_tajd, moluc_tajd, vicina_tajd, nitida_tajd, regnans_tajd) %>%
  mutate(BIN_START = BIN_START + 1) %>%
  select(species, CHROM, BIN_START, TajimaD)

# Function to create merged file for a species pair (handles order automatically)
merge_pair <- function(pop1_name, pop2_name) {
  
  # Extract Fst (try both orders)
  fst_pair <- fst %>% 
    filter((pop1 == pop1_name & pop2 == pop2_name) | 
             (pop1 == pop2_name & pop2 == pop1_name))
  
  if(nrow(fst_pair) == 0) {
    warning(paste("No Fst data found for", pop1_name, "vs", pop2_name))
  }
  
  # Extract Dxy (try both orders)
  dxy_pair <- dxy %>% 
    filter((pop1 == pop1_name & pop2 == pop2_name) | 
             (pop1 == pop2_name & pop2 == pop1_name))
  
  if(nrow(dxy_pair) == 0) {
    warning(paste("No Dxy data found for", pop1_name, "vs", pop2_name))
  }
  
  # Join Fst and Dxy (keeping their original pop1/pop2 columns for now)
  merged <- full_join(fst_pair, dxy_pair, 
                      by = c("chromosome", "window_pos_1", "window_pos_2"),
                      suffix = c("_fst", "_dxy"))
  
  # Standardize to ensure pop1 = pop1_name, pop2 = pop2_name
  merged <- merged %>%
    mutate(
      pop1 = pop1_name,
      pop2 = pop2_name
    ) %>%
    select(pop1, pop2, chromosome, window_pos_1, window_pos_2, avg_wc_fst, avg_dxy)
  
  # Add PI for pop1_name
  pi_pop1 <- pi %>% 
    filter(pop == pop1_name) %>%
    select(chromosome, window_pos_1, window_pos_2, avg_pi) %>%
    rename(!!paste0("pi_", pop1_name) := avg_pi)
  
  merged <- left_join(merged, pi_pop1, 
                      by = c("chromosome", "window_pos_1", "window_pos_2"))
  
  # Add PI for pop2_name
  pi_pop2 <- pi %>% 
    filter(pop == pop2_name) %>%
    select(chromosome, window_pos_1, window_pos_2, avg_pi) %>%
    rename(!!paste0("pi_", pop2_name) := avg_pi)
  
  merged <- left_join(merged, pi_pop2, 
                      by = c("chromosome", "window_pos_1", "window_pos_2"))
  
  # Add Tajima's D for pop1_name
  tajd_pop1 <- tajd %>% 
    filter(species == pop1_name) %>%
    select(CHROM, BIN_START, TajimaD) %>%
    rename(!!paste0("TajD_", pop1_name) := TajimaD)
  
  merged <- left_join(merged, tajd_pop1, 
                      by = c("chromosome" = "CHROM", "window_pos_1" = "BIN_START"))
  
  # Add Tajima's D for pop2_name
  tajd_pop2 <- tajd %>% 
    filter(species == pop2_name) %>%
    select(CHROM, BIN_START, TajimaD) %>%
    rename(!!paste0("TajD_", pop2_name) := TajimaD)
  
  merged <- left_join(merged, tajd_pop2, 
                      by = c("chromosome" = "CHROM", "window_pos_1" = "BIN_START"))
  
  # Rename Fst and Dxy for clarity
  merged <- merged %>%
    rename(fst = avg_wc_fst, dxy = avg_dxy)
  
  return(merged)
}

# Create merged files for all three pairs
moluc_largiflo <- merge_pair("moluc", "largiflo")
dealba_vicina <- merge_pair("dealba", "vicina")
nitida_regnans <- merge_pair("nitida", "regnans")

# Write output files
write_csv(moluc_largiflo, "moluc_largiflo_merged.csv")
write_csv(dealba_vicina, "dealba_vicina_merged.csv")
write_csv(nitida_regnans, "nitida_regnans_merged.csv")

message("Done!")
message("moluc vs largiflo: ", nrow(moluc_largiflo), " rows")
message("dealba vs vicina: ", nrow(dealba_vicina), " rows")
message("nitida vs regnans: ", nrow(nitida_regnans), " rows")

