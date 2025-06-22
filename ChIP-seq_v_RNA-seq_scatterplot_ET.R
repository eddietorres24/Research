library(dplyr)
library(rtracklayer)
library(ggplot2)
library(tidyr)

#Read in PRC2 targets
Prc2targets <- read.table("./bed_files/K27_genes_trimmed.bed", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t") 

#Read in All summary files
# WT1 <- read.table("bigwig_summaries/WT_Rep1_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# WT2 <- read.table("bigwig_summaries/WT_Rep2_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# WT3 <- read.table("bigwig_summaries/WT_Rep3_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# 
# cac1_1 <- read.table("bigwig_summaries/cac1_Rep1_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# cac1_2 <- read.table("bigwig_summaries/cac1_Rep2_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# cac1_3 <- read.table("bigwig_summaries/cac1_Rep3_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# 
# cac2_1 <- read.table("bigwig_summaries/cac2_Rep1_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# cac2_2 <- read.table("bigwig_summaries/cac2_Rep2_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# cac2_3 <- read.table("bigwig_summaries/cac2_Rep3_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# 
# cac3_1 <- read.table("bigwig_summaries/cac3_Rep1_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# cac3_2 <- read.table("bigwig_summaries/cac3_Rep2_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# cac3_3 <- read.table("bigwig_summaries/cac3_Rep3_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# 
# cac1_2_1 <- read.table("bigwig_summaries/cac1-2_Rep1_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
# cac1_2_2 <- read.table("bigwig_summaries/cac1-2_Rep2_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))

#assign rownames for later
# rownames(WT1) = WT1$gene

#combine replicates into 1 df (will take column w mean values)
# combined_df <- data.frame(
#   sample1 = WT1$mean,
#   sample2 = WT2$mean,
#   sample3 = WT3$mean,
#   sample4 = cac1_1$mean,
#   sample5 = cac1_2$mean,
#   sample6 = cac1_3$mean,
#   sample7 = cac2_1$mean,
#   sample8 = cac2_2$mean,
#   sample9 = cac2_3$mean,
#   sample10 = cac3_1$mean,
#   sample11 = cac3_2$mean,
#   sample12 = cac3_3$mean,
#   sample13 = cac1_2_1$mean,
#   sample14 = cac1_2_2$mean
# )

#rename rows
# rownames(combined_df) <- rownames(WT1)
# 
# #rename columns
# newnames <- c("WT1", "WT2", "WT3", "cac1_1", "cac1_2", "cac1_3", "cac2_1", "cac2_2", "cac2_3", "cac3_1", "cac3_2", "cac3_3", "cac1_2_1", "cac1_2_2")
# colnames(combined_df) <- newnames
# 
# #Average WT Data and eliminate rows that are below a mean of 4 (< 2 log2FC over background)
# combined_df$WTavg <- rowMeans(combined_df[,1:3])
# combined_df <- subset(combined_df, combined_df$WTavg > 4)
# combined_df <- combined_df[,-15]
# 
# #calculate rowmeans
# WT <- rowMeans(combined_df[, 1:3], na.rm = TRUE)
# cac1 <- rowMeans(combined_df[, 4:6], na.rm = TRUE)
# cac2 <- rowMeans(combined_df[, 7:9], na.rm = TRUE)
# cac3 <- rowMeans(combined_df[, 10:12], na.rm = TRUE)
# cac1_2 <- rowMeans(combined_df[, 13:14], na.rm = TRUE)
# 
# #calculate log2FC
# log2FC_cac1 <- log2(WT / cac1)
# log2FC_cac2 <- log2(WT / cac2)
# log2FC_cac3 <- log2(WT / cac3)
# log2FC_cac1_2 <- log2(WT / cac1_2)
# 
# ##make df for log2FC values
# log2FC_df <- data.frame(
#   log2FC_cac1   = log2FC_cac1,
#   log2FC_cac2   = log2FC_cac2,
#   log2FC_cac3   = log2FC_cac3,
#   log2FC_cac1_2 = log2FC_cac1_2
# )
# 
# ##assign a column for gene names
# log2FC_df$gene <- rownames(log2FC_df)
# ##add anmes to df in addition to NCU column
# log2FC_df$name <- Prc2targets$V10[match(log2FC_df$gene, Prc2targets$gene)]


###Calculate ChIP-seq Signal Averages

# Load your gene/promoter BED file as GRanges
colnames(K27_genes)[1:3] <- c("chr", "start", "end")  # set column names if missing
gene_ranges <- GRanges(
  seqnames = K27_genes$chr,
  ranges = IRanges(
    start = K27_genes$start,
    end = K27_genes$end
  ),
  strand = K27_genes$V6
)
mcols(gene_ranges)$gene_id <- K27_genes$V10

# Step 1: Split gene_ranges by strand
plus_ranges  <- gene_ranges[strand(gene_ranges) == "+"]
minus_ranges <- gene_ranges[strand(gene_ranges) == "-"]

# Step 2: Run dba.count separately
chip_data_plus  <- dba.count(chip_data, peaks = plus_ranges, summits = FALSE, filter = 0, minOverlap = 0)
chip_data_minus <- dba.count(chip_data, peaks = minus_ranges, summits = FALSE, filter = 0, minOverlap = 0)

# Step 3: Extract counts as GRanges
counts_plus  <- dba.peakset(chip_data_plus,  bRetrieve = TRUE)
counts_minus <- dba.peakset(chip_data_minus, bRetrieve = TRUE)

# Step 4: Combine into full GRanges object
full_counts <- c(counts_plus, counts_minus)

# Extract metadata (sample names, counts matrix)
count_matrix <- as.data.frame(mcols(full_counts))

# Extract genomic coordinates
region_info <- data.frame(
  CHR   = as.character(seqnames(full_counts)),
  START = start(full_counts),
  END   = end(full_counts),
  STRAND = as.character(strand(full_counts))
)

# Combine into final count dataframe
norm_counts <- cbind(region_info, count_matrix)

# Ensure coordinate columns match in name

###DEseq2 for ChIP data

# Step 1: Convert full_counts to a count matrix
count_matrix <- as.data.frame(mcols(full_counts))
rownames(count_matrix) <- paste0(seqnames(full_counts), ":", start(full_counts), "-", end(full_counts))

# Step 2: Create sample information (colData)
sample_info <- data.frame(
  sample = colnames(count_matrix),
  strain = gsub("_CS_K27_Rep_\\d+", "", colnames(count_matrix))  # Extract strain names
)
rownames(sample_info) <- sample_info$sample

# Step 3: Make DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),         # DESeq2 expects integer counts
  colData = sample_info,
  design = ~ strain
)

# Step 4: Run DESeq
dds2 <- DESeq(dds)
plotDispEsts(dds2)

# Step 5: Extract results (mutants vs. WT)
res_cac1 <- results(dds2, contrast = c("strain", "cac.1", "WT"))
res_cac2 <- results(dds2, contrast = c("strain", "cac.2", "WT"))
res_cac3 <- results(dds2, contrast = c("strain", "cac.3", "WT"))

#Step 6: Make new df w/ log2FC values

#Extract coordinates from DESeq2 results
coord_split <- do.call(rbind, strsplit(rownames(res_cac1), "[:-]"))
summary_df <- data.frame(
  CHR = coord_split[, 1],
  START = as.numeric(coord_split[, 2]),
  END = as.numeric(coord_split[, 3]),
  log2FC_cac1_vs_WT = res_cac1$log2FoldChange,
  log2FC_cac2_vs_WT = res_cac2$log2FoldChange,
  log2FC_cac3_vs_WT = res_cac3$log2FoldChange
)

#Merge with known gene IDs from norm_counts
summary_df <- merge(
  summary_df,
  norm_counts[, c("CHR", "START", "END", "gene_id")],
  by = c("CHR", "START", "END"),
  all.x = TRUE
)

#Reorder columns
summary_df <- summary_df[, c("gene_id", "CHR", "START", "END",
                             "log2FC_cac1_vs_WT", "log2FC_cac2_vs_WT", "log2FC_cac3_vs_WT")]

###Incorporate Analyzed RNA-seq data

###NEED TO RUN DEseq TO OBTAIN RNA-seq Log2FC VALUES###
##assign RNA-seq Log2FC for each strain after running DEseq (I am not eliminating non-DE genes in this code because i want both DE & non-DE genes present on my plot, because we are correlating DE w/ K27 changes)
cac1_seq <- read.csv("./CAF-1_RNA-seq_Analysis/csv_files/cac1_new_ALL.csv", stringsAsFactors=FALSE, row.names = 1, check.names=FALSE)
cac2_seq <- read.csv("./CAF-1_RNA-seq_Analysis/csv_files/cac2_ALL.csv", stringsAsFactors=FALSE, row.names = 1, check.names=FALSE)
cac3_seq <- read.csv("./CAF-1_RNA-seq_Analysis/csv_files/cac3_ALL.csv", stringsAsFactors=FALSE, row.names = 1, check.names=FALSE)
cac1_cac2_seq <- read.csv("./CAF-1_RNA-seq_Analysis/csv_files/cac1_cac2_ALL.csv", stringsAsFactors=FALSE, row.names = 1, check.names=FALSE)

#append RNA-seq log2FC from DEseq to the end of the ChIP log2FC df
summary_df$cac1RNA <- cac1_seq$log2FoldChange[match(summary_df$gene_id, rownames(cac1_seq))]
summary_df$cac2RNA <- cac2_seq$log2FoldChange[match(summary_df$gene_id, rownames(cac2_seq))]
summary_df$cac3RNA <- cac3_seq$log2FoldChange[match(summary_df$gene_id, rownames(cac3_seq))]
summary_df$cac1_cac2RNA <- cac1_cac2_seq$log2FoldChange[match(summary_df$gene_id, rownames(cac1_cac2_seq))]

#remove any rows w/ NA
summary_df <- na.omit(summary_df)

#statistics
cor_test <- cor.test(summary_df$log2FC_cac3_vs_WT, summary_df$cac3RNA, method = "pearson")
r2 <- round(cor_test$estimate^2, 2)
pval <- signif(cor_test$p.value, 3)

#Scater Plot
ggplot(summary_df, aes(x = log2FC_cac3_vs_WT, y = cac3RNA)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "blue") +
  annotate("text", x = -7.5, y = 12.5, 
           label = paste0("R² = ", r2, ", p = ", pval), hjust = 0) +
  theme_minimal() +
  labs(
    x = "H3K27me3 Signal (cac-3/WT)",
    y = "log2 Fold Change (cac-3/WT)",
    title = "ChIP-seq vs RNA-seq Scatter Plot"
  )


ggsave(filename = "./K27_v_RNA_cac3_Paper.pdf", plot = plot, dpi=600, height=9, width=12, units = "in")



###TO-DO###
#1. color genes that are located within ectopic CAF-1 peaks
#2. 

