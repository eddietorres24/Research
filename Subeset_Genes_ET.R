#Code to subset genes with overlap in H3K27me3 regions

library(GenomicRanges)
library(rtracklayer)
library(dplyr)

# === Load Gene BED ===
# gene_bed_path <- "./bed_files/K27_genes_trimmed_2.bed"
# genes <- import(gene_bed_path, format = "BED")
genes <- read.csv(file = "./csv_files/K27_genes_trimmed.csv")
names(mcols(genes)) <- c("gene_id")

# === Load Peak File (BED or broadPeak format) ===
peak_bed_path <- "peaks.bed"  # or peaks.broadPeak
peaks <- import(peak_bed_path, format = "BED")

# === Function to Calculate Percent Overlap ===
calculate_overlap_fraction <- function(genes, peaks) {
  overlaps <- findOverlaps(genes, peaks)
  
  gene_hits <- genes[queryHits(overlaps)]
  peak_hits <- peaks[subjectHits(overlaps)]
  
  intersect_ranges <- pintersect(gene_hits, peak_hits)
  intersect_df <- data.frame(
    gene_id = mcols(gene_hits)$gene_id,
    width = width(gene_hits),
    intersect_width = width(intersect_ranges)
  )
  
  # Sum overlaps per gene
  overlap_summary <- intersect_df %>%
    group_by(gene_id, width) %>%
    summarise(total_overlap = sum(intersect_width), .groups = "drop") %>%
    mutate(overlap_fraction = total_overlap / width)
  
  return(overlap_summary)
}

# === Run Overlap Calculation ===
overlap_results <- calculate_overlap_fraction(genes, peaks)

# === Bin Genes by Overlap Tier (e.g. 50-59%, 60-69%, ..., 100%) ===
tiers <- seq(0.5, 1.0, by = 0.1)

overlap_tiers <- lapply(tiers, function(t) {
  filtered <- overlap_results %>%
    filter(overlap_fraction >= t & overlap_fraction < t + 0.1)
  return(filtered$gene_id)
})
names(overlap_tiers) <- paste0(seq(50, 100, by = 10), "%")

# === Output Tiered Gene Lists ===
for (tier in names(overlap_tiers)) {
  cat(sprintf("Genes with ≥%s overlap:\n", tier))
  print(overlap_tiers[[tier]])
  cat("\n")
}
