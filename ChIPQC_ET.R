### Script for ChIPQC

## Load Libraries
library(BiocManager)
library(ChIPQC)

## Set Working Directory
setwd("C:\\Users\\eddie\\Desktop\\Research\\GitHub")

## Load Sample Data
samples <- read.csv('meta/samplesheet_chr12.csv')

## Make ChIPQC Object
chipObj <- ChIPQC(samples) 

## Make QC Report
ChIPQCreport(chipObj, reportName="ChIP_QC_report_cac_H3K27me3", reportFolder="ChIPQCreport")

### DiffBind Script

## Load Libraries
library(DiffBind)
library(tidyverse)

## Run DBA
dbObj <- dba(sampleSheet=samples)

## Plot PCA & Correlation Plot
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj)

## Contrast Samples
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR, minMembers = 2)

## Differential Enrichment Analysis
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
dba.show(dbObj, bContrasts=T)

# Plot
dba.plotPCA(dbObj, contrast=1, method=DBA_DESEQ2, attributes=DBA_FACTOR, label=DBA_ID)

## Visualize Results

# Peaks recognized by edgeR vs. DEseq
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

# MA Plots
dba.plotMA(dbObj, method=DBA_DESEQ2)
dba.plotMA(dbObj, bXY=TRUE)
pvals <- dba.plotBox(dbObj)

## Extract Results
res_deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)

# Write to File
out <- as.data.frame(res_deseq)
write.table(out, file="results/cac1_H3K27me3_DiffBind.txt", sep="\t", quote=F, row.names=F)

# Create bed files for each keeping only significant peaks (p < 0.05)

WT_enrich <- out %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  select(seqnames, start, end)

# Write to file
write.table(WT_enrich, file="WT_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

cac1_enrich <- out %>% 
  filter(FDR < 0.05 & Fold < 0) %>% 
  select(seqnames, start, end)

# Write to file
write.table(cac1_enrich, file="cac1_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

