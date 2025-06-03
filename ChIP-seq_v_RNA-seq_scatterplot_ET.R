library(dplyr)
library(rtracklayer)
library(ggplot2)

dfchip <- read.table("bigwig_summaries/WT_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean"))
dfchip <- dfchip[,-5]
dfchip$log2mean = log2(dfchip$mean)

cac1chip <- read.table("cac1_K27signal_K27genes.tab", header = FALSE, col.names = c("gene", "size", "covered", "sum", "mean0", "mean"))
dfchip$log2meancac1 = log2(cac1chip$mean)

rna_expr <- data.frame(
  gene = rownames(AVERAGE_Prc2targetTPM),
  WT = AVERAGE_Prc2targetTPM[,1],
  cac1 = AVERAGE_Prc2targetTPM[,3],
  cac2 = AVERAGE_Prc2targetTPM[,4],
  cac1cac2 = AVERAGE_Prc2targetTPM[,5],
  cac3 = AVERAGE_Prc2targetTPM[,6],
  set7 = AVERAGE_Prc2targetTPM[,2]
)

merged_df <- merge(dfchip, rna_expr, by = "gene", all.x = TRUE)
merged_df <- merged_df[!is.na(merged_df$log2FoldChange), ]

#deltas
merged_df$deltachipcac1 <- merged_df$log2mean - merged_df$log2meancac1
merged_df$deltaRNAcac1 <- merged_df$cac1 - merged_df$WT

#Scater Plot
plot = ggplot(merged_df, aes(x = deltachipcac1, y = log2FoldChange)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  labs(
    x = "H3K27me3 Signal (WT - cac1)",
    y = "log2 Fold Change (cac1/WT)",
    title = "ChIP-seq vs RNA-seq Scatter Plot"
  ) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "blue")

ggsave(filename = "./K27_v_RNA_cac1.pdf", plot = plot, dpi=600, height=9, width=12, units = "in")


#Statistics
model <- lm(merged_df[,18] ~ merged_df[,15])

# Get R-squared value
summary_model <- summary(model)
r_squared <- summary_model$r.squared
p_value <- summary_model$coefficients[2, 4]  # p-value for slope

# Print to console
cat("R-squared:", round(r_squared, 4), "\n")
cat("P-value:", format.pval(p_value, digits = 3), "\n")
