# antisense_deseq_minimal.R
suppressPackageStartupMessages({
  library(DESeq2)
})

## --- USER PATHS (edit these two) ---
BASE_DIR   <- "C:/Users/eddie/Research/GitHub/Research"
COLDATA_CSV <- "C:/Users/eddie/Research/GitHub/coldata_anti.csv"

## If your copied counts live in BASE_DIR/counts/..., this is fine.
## Otherwise, point COUNTS_DIR to wherever the .antisense.txt files are.
COUNTS_DIR <- file.path(BASE_DIR, "counts")

## Switch between "antisense" and "sense" by changing this:
COUNT_KIND <- "antisense"

#--- helpers ---------------------------------------------------------------

findCountsFile <- function(accession, kind = c("antisense","sense"),
                           counts_dir = COUNTS_DIR) {
  kind <- match.arg(kind)
  # candidate layouts:
  cand <- c(
    file.path(counts_dir, accession, paste0(accession, ".", kind, ".txt")),
    file.path(counts_dir, paste0(accession, ".", kind, ".txt"))
  )
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0) {
    stop(sprintf("Counts file not found for %s (%s): tried\n  - %s\n  - %s",
                 accession, kind, cand[1], cand[2]), call. = FALSE)
  }
  hit[1]
}

read_fc_one <- function(f) {
  df <- read.table(f, header = TRUE, sep = "\t",
                   comment.char = "#", check.names = FALSE, quote = "")
  gene  <- df[[1]]                 # "Geneid" column
  count <- df[[ncol(df)]]          # the only sample column in each file
  data.frame(gene = gene, cnt = as.numeric(count), 
             stringsAsFactors = FALSE)
}

#--- load coldata ----------------------------------------------------------

coldata <- read.csv(COLDATA_CSV, row.names = 1, stringsAsFactors = FALSE)
stopifnot(all(c("condition","type","file_id") %in% colnames(coldata)))
coldata$condition <- factor(coldata$condition)
coldata$type      <- factor(coldata$type)

# Basenames to look up files by
samples <- rownames(coldata)
file_ids <- setNames(coldata$file_id, samples)

#--- build counts matrix ---------------------------------------------------

message("Locating count files...")
files <- vapply(file_ids, findCountsFile, FUN.VALUE = character(1),
                kind = COUNT_KIND, counts_dir = COUNTS_DIR)

message("Reading featureCounts tables...")
lst <- lapply(files, read_fc_one)

# merge by gene ID
merged <- Reduce(function(x,y) merge(x,y, by = "gene", all = TRUE), lst)
rownames(merged) <- merged$gene
merged$gene <- NULL
colnames(merged) <- samples  # label columns by sample_id

# replace NAs with 0 and coerce to integer matrix
merged[is.na(merged)] <- 0
cts <- as.matrix(round(merged, 0))

# ensure order matches coldata
cts <- cts[, rownames(coldata), drop = FALSE]

# quick write-out for inspection
write.csv(data.frame(Geneid = rownames(cts), cts, check.names = FALSE),
          file = file.path(BASE_DIR, paste0("combined_", COUNT_KIND, "_counts.csv")),
          row.names = FALSE)

#--- minimal DESeq2 --------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData   = coldata,
                              design    = ~ condition)

# light prefilter: keep genes with at least 10 total counts
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

# set reference
if ("WT" %in% levels(dds$condition)) {
  dds$condition <- relevel(dds$condition, ref = "WT")
}

dds <- DESeq(dds)

# example results: set7 vs WT (change the contrast as you like)
if ("set7" %in% levels(dds$condition)) {
  res_set7 <- results(dds, contrast = c("condition","set7","WT"))
  write.csv(as.data.frame(res_set7),
            file = file.path(BASE_DIR, "DE_set7_vs_WT_antisense.csv"))
}

# save the DESeq2 object for downstream work
saveRDS(dds, file = file.path(BASE_DIR, paste0("dds_", COUNT_KIND, ".rds")))

##setting p-value cutoff, starting with 0.ant
alpha = 0.ant

#Subsetting results based on p-value
cac1_old_ant <- results(dds, contrast=c("condition", "cac1_old", "WT"))
cac2_ant <- results(dds, contrast=c("condition", "cac2", "WT"))
cac3_ant <- results(dds, contrast=c("condition", "cac3", "WT"))
asf1_ant <- results(dds, contrast=c("condition", "asf1", "WT"))
set7_ant <- results(dds, contrast=c("condition", "set7", "WT"))
cac1_new_ant <- results(dds, contrast=c("condition", "cac1_new", "WT"))
cac1_cac2_ant <- results(dds, contrast=c("condition", "cac1_cac2", "WT"))
cac1_suz12_ant <- results(dds, contrast=c("condition", "cac1_suz12", "WT"))
suz12_ant <- results(dds, contrast=c("condition", "suz12", "WT"))
ash1_ant <- results(dds, contrast=c("condition", "ash1", "WT"))

#Make csvs
write.csv(as.data.frame(cac1_old_ant), file = "cac1_ant.csv")
write.csv(as.data.frame(cac2_ant), file = "cac2_ant.csv")
write.csv(as.data.frame(cac3_ant), file = "cac3_ant.csv")
write.csv(as.data.frame(asf1_ant), file = "asf1_ant.csv")
write.csv(as.data.frame(set7_ant), file = "set7_ant.csv")
write.csv(as.data.frame(cac1_new_ant), file = "cac1_new_ant.csv")
write.csv(as.data.frame(cac1_cac2_ant), file = "cac1_cac2_ant.csv")
write.csv(as.data.frame(cac1_suz12_ant), file = "cac1_suz12_ant.csv")
write.csv(as.data.frame(suz12_ant), file = "suz12_ant.csv")
write.csv(as.data.frame(ash1_ant), file = "ash1_ant.csv")


##############Making IGV Track
#calling list of files
##Note: this will only work if the ONLY csv's in this directory are the ones containing the Data you want on the track (i.e. you need to move any other csv's from this directory, or move the ones you want on the track to a new one and change the path)
list_of_files <- list.files(path = "./Research",
                            recursive = FALSE,
                            pattern = ".csv$")

setwd("C:/Users/eddie/Research/GitHub/Research")

induced_genes <- readr::read_csv(list_of_files, id = "file_name")

##convert to wide format
l2fc_data <- data.frame(pivot_wider(data = induced_genes, id_cols = "...1", names_from = "file_name", values_from = "log2FoldChange"))

#assigning row names and deleting extra column
rownames(l2fc_data) <- l2fc_data[,1]
l2fc_data <- l2fc_data[,-1]

#assigning NA values to 0
l2fc_data[is.na(l2fc_data)] <- 0

#tibble NCU names back to first column
l2fc <- tibble::rownames_to_column(l2fc_data, "VALUE")

#read in file w/ all genes
ngenes = read.csv("../csv_files/neurospora_genes_edit.csv")

#Make IGV Track
##Note: only genes with DE in one or more strains will display data in the track
combine <- transform(merge(ngenes, l2fc, by = "VALUE"))
track = combine %>% select("SequenceID", "FeatureStart", "FeatureEnd", "VALUE", "asf1_ant.csv", "cac1_ant.csv", "cac2_ant.csv", "cac3_ant.csv", "set7_ant.csv", "cac1_new_ant.csv", "cac1_cac2_ant.csv", "cac1_suz12_ant.csv", "suz12_ant.csv", "ash1_ant.csv")

#write track to file
write.table(track, file="newseq.igv", sep="\t", row.names = FALSE, quote = FALSE)
