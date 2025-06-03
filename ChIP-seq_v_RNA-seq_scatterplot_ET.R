library(dplyr)
library(rtracklayer)
library(ggplot2)
library(tidyr)

#Read in PRC2 targets
Prc2targets <- read.table("./bed_files/K27_genes_trimmed.bed", header=FALSE, stringsAsFactors=FALSE, check.names=FALSE, sep="\t") 

#Read in All summary files
WT1 <- read.table("bigwig_summaries/WT_Rep1_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
WT2 <- read.table("bigwig_summaries/WT_Rep2_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
WT3 <- read.table("bigwig_summaries/WT_Rep3_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))

cac1_1 <- read.table("bigwig_summaries/cac1_Rep1_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
cac1_2 <- read.table("bigwig_summaries/cac1_Rep2_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
cac1_3 <- read.table("bigwig_summaries/cac1_Rep3_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))

cac2_1 <- read.table("bigwig_summaries/cac2_Rep1_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
cac2_2 <- read.table("bigwig_summaries/cac2_Rep2_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
cac2_3 <- read.table("bigwig_summaries/cac2_Rep3_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))

cac3_1 <- read.table("bigwig_summaries/cac3_Rep1_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
cac3_2 <- read.table("bigwig_summaries/cac3_Rep2_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
cac3_3 <- read.table("bigwig_summaries/cac3_Rep3_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))

cac1_2_1 <- read.table("bigwig_summaries/cac1-2_Rep1_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
cac1_2_2 <- read.table("bigwig_summaries/cac1-2_Rep2_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))

#assign rownames for later
rownames(WT1) = WT1$gene

#combine replicates into 1 df (will take column w mean values)
combined_df <- data.frame(
  sample1 = WT1$mean,
  sample2 = WT2$mean,
  sample3 = WT3$mean,
  sample4 = cac1_1$mean,
  sample5 = cac1_2$mean,
  sample6 = cac1_3$mean,
  sample7 = cac2_1$mean,
  sample8 = cac2_2$mean,
  sample9 = cac2_3$mean,
  sample10 = cac3_1$mean,
  sample11 = cac3_2$mean,
  sample12 = cac3_3$mean,
  sample13 = cac1_2_1$mean,
  sample14 = cac1_2_2$mean
)

#rename rows
rownames(combined_df) <- rownames(WT1)

#rename columns
newnames <- c("WT1", "WT2", "WT3", "cac1_1", "cac1_2", "cac1_3", "cac2_1", "cac2_2", "cac2_3", "cac3_1", "cac3_2", "cac3_3", "cac1_2_1", "cac1_2_2")
colnames(combined_df) <- newnames

#Average WT Data and eliminate rows that are below a mean of 4 (< 2 log2FC over background)
combined_df$WTavg <- rowMeans(combined_df[,1:3])
combined_df <- subset(combined_df, combined_df$WTavg > 4)
combined_df <- combined_df[,-15]

#calculate rowmeans
WT <- rowMeans(combined_df[, 1:3], na.rm = TRUE)
cac1 <- rowMeans(combined_df[, 4:6], na.rm = TRUE)
cac2 <- rowMeans(combined_df[, 7:9], na.rm = TRUE)
cac3 <- rowMeans(combined_df[, 10:12], na.rm = TRUE)
cac1_2 <- rowMeans(combined_df[, 13:14], na.rm = TRUE)

#calculate log2FC
log2FC_cac1 <- log2(WT / cac1)
log2FC_cac2 <- log2(WT / cac2)
log2FC_cac3 <- log2(WT / cac3)
log2FC_cac1_2 <- log2(WT / cac1_2)

##make df for log2FC values
log2FC_df <- data.frame(
  log2FC_cac1   = log2FC_cac1,
  log2FC_cac2   = log2FC_cac2,
  log2FC_cac3   = log2FC_cac3,
  log2FC_cac1_2 = log2FC_cac1_2
)

##assign a column for gene names
log2FC_df$gene <- rownames(log2FC_df)
##add anmes to df in addition to NCU column
log2FC_df$name <- Prc2targets$V10[match(log2FC_df$gene, Prc2targets$gene)]

###NEED TO RUN DEseq TO OBTAIN RNA-seq Log2FC VALUES###
##assign RNA-seq Log2FC for each strain after running DEseq (I am not eliminating non-DE genes in this code because i want both DE & non-DE genes present on my plot, because we are correlating DE w/ K27 changes)
cac1_seq <- as.data.frame(cac1_new)
cac2_seq <- as.data.frame(cac2)
cac3_seq <- as.data.frame(cac3)
cac1_cac2_seq <- as.data.frame(cac1_cac2)

#append RNA-seq log2FC from DEseq to the end of the ChIP log2FC df
log2FC_df$cac1RNA <- cac1_seq$log2FoldChange[match(log2FC_df$gene, rownames(cac1_seq))]
log2FC_df$cac2RNA <- cac2_seq$log2FoldChange[match(log2FC_df$gene, rownames(cac2_seq))]
log2FC_df$cac3RNA <- cac3_seq$log2FoldChange[match(log2FC_df$gene, rownames(cac3_seq))]
log2FC_df$cac1_cac2RNA <- cac1_cac2_seq$log2FoldChange[match(log2FC_df$gene, rownames(cac1_cac2_seq))]

#remove any rows w/ NA
log2FC_df <- na.omit(log2FC_df)

#Scater Plot
plot = ggplot(log2FC_df, aes(x = log2FC_cac1_2, y = cac1_cac2RNA)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    x = "H3K27me3 Signal (WT/cac1-2)",
    y = "log2 Fold Change (cac1-2/WT)",
    title = "ChIP-seq vs RNA-seq Scatter Plot"
  ) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "blue")

ggsave(filename = "./K27_v_RNA_cac1-2.pdf", plot = plot, dpi=600, height=9, width=12, units = "in")


#Statistics
model <- lm(log2FC_df[,10] ~ log2FC_df[,4])

# Get R-squared value
summary_model <- summary(model)
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]  # p-value for slope

# Print to console
cat("R-squared:", round(r_squared, 4), "\n")
cat("P-value:", format.pval(p_value, digits = 3), "\n")


###TO-DO###
#1. color genes that are located within ectopic CAF-1 peaks
#2. 

