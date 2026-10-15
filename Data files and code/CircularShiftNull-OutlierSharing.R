# Circular-shift permutation test for overlap of outlier windows among species pairs,
# plus the aridity-specific sharing statistic Delta.
#
# Pair 1 = E. moluccana - E. largiflorens   (wet-dry)
# Pair 2 = E. dealbata  - E. vicina         (wet-dry)
# Pair 3 = E. nitida    - E. regnans        (wetter-wet, control for this analysis)
#
# Null: for each comparison, pair A's outlier positions are rotated by a random offset
# within each chromosome (pair B stays fixed) and the overlap is recounted. This keeps
# the clustering of outliers along chromosomes. Expected overlap = mean permuted overlap.

library(readxl)

## ---- 1. Settings (edit these) ---------------------------------------------------
files <- c(
  mol_lar = "moluc_largiflo_merged_outliers.xlsx",  # pair 1
  dea_vic = "dealba_vicina_merged_outliers.xlsx",         # pair 2
  nit_reg = "nitida_regnans_merged_outliers.xlsx"           # pair 3
)
flag_col <- "top1pct_residual_overlap"  # 1 = joint top-1% residual FST and dXY outlier
nperm    <- 10000

# FALSE: rotate over every window in the chromosome (matches the outlier counts in the paper).
# TRUE : for each comparison, keep only windows with data in BOTH pairs before rotating
#        (cleaner null if windows with missing data differ between pairs; outlier counts
#        can drop slightly).
restrict_to_common <- FALSE

set.seed(1)

## ---- 2. Read the three files onto one common window grid -----------------------
read_pair <- function(path, name) {
  d   <- as.data.frame(read_excel(path))
  has <- !is.na(d[[flag_col]])
  out <- data.frame(chromosome   = d$chromosome,
                    window_pos_1 = d$window_pos_1,
                    stringsAsFactors = FALSE)
  out[[paste0(name, "_has")]] <- has
  out[[paste0(name, "_out")]] <- as.integer(has & d[[flag_col]] == 1)
  out
}

pairs <- Map(read_pair, files, names(files))
win <- Reduce(function(x, y) merge(x, y, by = c("chromosome", "window_pos_1"), all = TRUE),
              pairs)

# Windows absent from a file count as "no data", not outliers
for (nm in names(files)) {
  h <- paste0(nm, "_has"); o <- paste0(nm, "_out")
  win[[h]][is.na(win[[h]])] <- FALSE
  win[[o]][is.na(win[[o]])] <- 0L
}
win <- win[order(win$chromosome, win$window_pos_1), ]

cat("Windows on common grid:", nrow(win), "\n")
for (nm in names(files)) {
  cat(sprintf("%s: %d outliers, %d windows with data\n",
              nm, sum(win[[paste0(nm, "_out")]]), sum(win[[paste0(nm, "_has")]])))
}

## ---- 3. Circular-shift permutation ----------------------------------------------
circ_shift_overlap <- function(win, A, B, nperm, restrict) {
  d <- win
  if (restrict) d <- d[d[[paste0(A, "_has")]] & d[[paste0(B, "_has")]], ]
  xa <- d[[paste0(A, "_out")]]
  xb <- d[[paste0(B, "_out")]]
  obs  <- sum(xa * xb)
  perm <- numeric(nperm)
  
  for (idx in split(seq_len(nrow(d)), d$chromosome)) {
    L  <- length(idx)
    ia <- which(xa[idx] == 1)          # positions of A's outliers on this chromosome
    yb <- xb[idx]                      # B's outlier indicator (fixed)
    if (length(ia) == 0 || sum(yb) == 0) next
    shifts <- sample.int(L, nperm, replace = TRUE) - 1L   # independent offset per replicate
    perm <- perm + vapply(shifts,
                          function(k) sum(yb[((ia - 1L + k) %% L) + 1L]),
                          numeric(1))
  }
  
  list(A = A, B = B, a = sum(xa), b = sum(xb), N = nrow(d),
       observed = obs,
       expected = mean(perm),
       excess   = obs - mean(perm),
       ratio    = obs / mean(perm),
       p        = (1 + sum(perm >= obs)) / (nperm + 1),   # one-sided, add-one correction
       perm     = perm)
}

comps <- list(c("mol_lar", "dea_vic"),   # 1-2
              c("mol_lar", "nit_reg"),   # 1-3
              c("dea_vic", "nit_reg"))   # 2-3
res <- lapply(comps, function(p) circ_shift_overlap(win, p[1], p[2], nperm, restrict_to_common))

tab <- data.frame(
  comparison = c("1-2", "1-3", "2-3"),
  a          = sapply(res, `[[`, "a"),
  b          = sapply(res, `[[`, "b"),
  N          = sapply(res, `[[`, "N"),
  observed   = sapply(res, `[[`, "observed"),
  expected   = round(sapply(res, `[[`, "expected"), 3),
  excess     = round(sapply(res, `[[`, "excess"), 3),
  obs_over_exp = round(sapply(res, `[[`, "ratio"), 2),
  P          = round(sapply(res, `[[`, "p"), 4)
)
print(tab)

## ---- 4. Delta -------------------------------------------------------------------
# Delta = E12 - (E13 + E23)/2, where E = observed - expected overlap.
# Measures aridity-specific sharing only if the three pairs are equally prone to share signal.
E     <- tab$excess
delta <- E[1] - (E[2] + E[3]) / 2

# Chance scale of Delta: spread of Delta when every comparison is a circular-shift null
dev        <- lapply(res, function(r) r$perm - mean(r$perm))
delta_null <- dev[[1]] - (dev[[2]] + dev[[3]]) / 2
delta_sd   <- sd(delta_null)
delta_p    <- (1 + sum(abs(delta_null) >= abs(delta))) / (nperm + 1)

cat(sprintf("\nDelta = %.3f (null SD = %.2f; two-sided P = %.3f)\n", delta, delta_sd, delta_p))

## ---- 5. Shared windows, for checking where they sit (islands, recombination) ----
shared <- do.call(rbind, lapply(seq_along(comps), function(i) {
  p <- comps[[i]]
  s <- win[win[[paste0(p[1], "_out")]] == 1 & win[[paste0(p[2], "_out")]] == 1,
           c("chromosome", "window_pos_1")]
  if (nrow(s) > 0) cbind(comparison = tab$comparison[i], s)
}))
print(shared)

write.csv(tab,    "overlap_circular_shift_results.csv", row.names = FALSE)
write.csv(shared, "shared_outlier_windows.csv",         row.names = FALSE)
