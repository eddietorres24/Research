#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --time=8:00:00
#SBATCH --output=../RNAseqMap/logs/featureCounts.%j.out
#SBATCH --error=../RNAseqMap/logs/featureCounts.%j.err

#This code will only be used for building a count file across multiple Samples/Strains

# cd $SLURM_SUBMIT_DIR
#
# THREADS=2
#
# #Make & Assign Directories
# BAMDIR="/scratch/evt82290/MappingOutputs/RNAseq/bamFiles"
# B149DIR="/scratch/evt82290/MappingOutputs/Run149/RNA/bamFiles"
# OUTDIR="/scratch/evt82290/MappingOutputs/RNAseq/counts"
#
# # #if output directory doesn't exist, create it
# if [ ! -d $OUTDIR ]
# then
#     mkdir -p $OUTDIR
# fi
# ###
# # /scratch/evt82290/RNAseq/CAF-1_Heatmap/bamFiles/SRR7970630/SRR7970630_Aligned.sortedByCoord.out.bam
# # /scratch/evt82290/RNAseq/CAF-1_Heatmap/bamFiles/SRR7970630/SRR7970630_Aligned.sortedByCoord.out.bam
# #Run featureCounts
# module load Subread/
#
# featureCounts -T $THREADS \
# -p \
# -t CDS \
# -g gene_name \
# -s 0 --primary \
# -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
# -o $OUTDIR/readcounts_All_CAF1paper.txt \
# $BAMDIR/SRR8444037/SRR8444037_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR8444038/SRR8444038_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR8444043/SRR8444043_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR7970629/SRR7970629_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR7970630/SRR7970630_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR7970631/SRR7970631_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR7970598/SRR7970598_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR7970599/SRR7970599_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR7970600/SRR7970600_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR8269830/SRR8269830_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR8269647/SRR8269647_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR8269650/SRR8269650_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR10916318/SRR10916318_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR10916319/SRR10916319_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR10916320/SRR10916320_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR7970603/SRR7970603_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR7970606/SRR7970606_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR7970610/SRR7970610_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR9027727/SRR9027727_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR9027728/SRR9027728_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR9027730/SRR9027730_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR8269825/SRR8269825_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR8269775/SRR8269775_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR8269782/SRR8269782_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR8269810/SRR8269810_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR10916163/SRR10916163_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR10916164/SRR10916164_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR10916165/SRR10916165_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-115_RNA_cac-1__Rep1_S139/149-115_RNA_cac-1__Rep1_S139_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-116_RNA_cac-1__Rep2_S140/149-116_RNA_cac-1__Rep2_S140_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-117_RNA_cac-1__Rep3_S141/149-117_RNA_cac-1__Rep3_S141_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-65_RNA_cac-1_cac-2___S65/149-65_RNA_cac-1_cac-2___S65_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-66_RNA_cac-1_cac-2___S66/149-66_RNA_cac-1_cac-2___S66_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-67_RNA_cac-1_cac-2___S67/149-67_RNA_cac-1_cac-2___S67_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-96_RNA_cac-1-suz-12__Rep1_S96/149-96_RNA_cac-1-suz-12__Rep1_S96_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-97_RNA_cac-1-suz-12__Rep2_S121/149-97_RNA_cac-1-suz-12__Rep2_S121_Aligned.sortedByCoord.out.bam \
# $B149DIR/149-98_RNA_cac-1-suz-12__Rep3_S122/149-98_RNA_cac-1-suz-12__Rep3_S122_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR7690267/SRR7690267_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR7690268/SRR7690268_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR9027658/SRR9027658_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR9027759/SRR9027759_Aligned.sortedByCoord.out.bam \
# $BAMDIR/SRR9027689/SRR9027689_Aligned.sortedByCoord.out.bam

#ChatGPT code

# set -euo pipefail
#
# # --- modules (adjust names if your cluster uses different module keys) ---
# module load SAMtools
# module load Subread
# module load R
#
# cd "$SLURM_SUBMIT_DIR"
#
# # ====================== CONFIG ======================
# THREADS=${SLURM_CPUS_PER_TASK:-4}
#
# # BAM roots to search (recursive)
# BAMDIR="/scratch/evt82290/MappingOutputs/RNAseq/bamFiles"
# B149DIR="/scratch/evt82290/MappingOutputs/Run149/RNA/bamFiles"
#
# # Output
# OUTDIR="/scratch/evt82290/MappingOutputs/RNAseq/counts"
# MERGED="$OUTDIR/readcounts_All_CAF1paper.merged.txt"
#
# # Annotation
# GTF="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf"
# STRAND=0          # 0=unstranded, 1=forward, 2=reverse
#
# # FeatureCounts target (your current choice)
# FTYPE="CDS"       # common alternative for gene-level RNA-seq: "exon"
# GATTR="gene_name" # common alternative: "gene_id"
# # Optional sample renaming map (TSV with columns: file  sample). Leave empty to skip.
# SAMPLE_MAP=""
# # =====================================================
#
# mkdir -p "$OUTDIR" "$OUTDIR/tmp"
#
# # 1) Collect BAMs recursively from both roots
# echo "[INFO] Finding BAMs under: $BAMDIR  and  $B149DIR"
# mapfile -d '' ALL_BAMS < <(find "$BAMDIR" "$B149DIR" -type f -name "*Aligned.sortedByCoord.out.bam" -print0)
# if (( ${#ALL_BAMS[@]} == 0 )); then
#   echo "[ERROR] No BAMs found."
#   exit 1
# fi
# echo "[INFO] Found ${#ALL_BAMS[@]} BAMs."
#
# # 2) Split into PE vs SE via FLAG 0x1
# PE_BAMS=()
# SE_BAMS=()
# for bam in "${ALL_BAMS[@]}"; do
#   if (( $(samtools view -c -f 1 "$bam") > 0 )); then
#     PE_BAMS+=("$bam")
#   else
#     SE_BAMS+=("$bam")
#   fi
# done
# echo "[INFO] Detected ${#PE_BAMS[@]} paired-end BAMs and ${#SE_BAMS[@]} single-end BAMs."
#
# # 3) Run featureCounts separately
# SE_OUT=""
# PE_OUT=""
#
# if (( ${#SE_BAMS[@]} )); then
#   SE_OUT="$OUTDIR/readcounts.SE.txt"
#   echo "[INFO] Running featureCounts (SE) -> $SE_OUT"
#   featureCounts -T "$THREADS" \
#     -t "$FTYPE" -g "$GATTR" -s "$STRAND" --primary \
#     --tmpDir "$OUTDIR/tmp" \
#     -a "$GTF" -o "$SE_OUT" \
#     "${SE_BAMS[@]}"
# fi
#
# if (( ${#PE_BAMS[@]} )); then
#   PE_OUT="$OUTDIR/readcounts.PE.txt"
#   echo "[INFO] Running featureCounts (PE) -> $PE_OUT"
#   featureCounts -T "$THREADS" -p -B -C \
#     -t "$FTYPE" -g "$GATTR" -s "$STRAND" --primary \
#     --tmpDir "$OUTDIR/tmp" \
#     -a "$GTF" -o "$PE_OUT" \
#     "${PE_BAMS[@]}"
# fi
#
# # 4) Merge into one table with a robust R parser that preserves headers
# echo "[INFO] Merging outputs -> $MERGED"
# Rscript - "$SE_OUT" "$PE_OUT" "$MERGED" "$SAMPLE_MAP" <<'RSCRIPT'
# args <- commandArgs(trailingOnly=TRUE)
# se <- args[1]; pe <- args[2]; out <- args[3]; map_path <- if (length(args) >= 4) args[4] else ""
#
# files <- c(se, pe)
# files <- files[file.exists(files) & nzchar(files)]
# if (length(files) == 0) stop("No featureCounts outputs found")
#
# read_fc <- function(f) {
#   x <- read.table(
#     f,
#     header = TRUE,           # force header row
#     sep = "\t",
#     quote = "",
#     comment.char = "#",      # drop comment lines
#     check.names = FALSE,
#     stringsAsFactors = FALSE
#   )
#   # Remove BOM if present
#   colnames(x)[1] <- sub("^\ufeff", "", colnames(x)[1])
#
#   need <- c("Geneid","Chr","Start","End","Strand","Length")
#   if (!all(need %in% colnames(x)[1:6])) {
#     stop(sprintf("Header parse failed for %s. First cols: %s",
#                  f, paste(colnames(x)[1:6], collapse=" | ")))
#   }
#   x
# }
#
# dfs <- lapply(files, read_fc)
#
# # Start with annotation cols from first file
# base <- dfs[[1]][, 1:6]
#
# # Append sample columns from each file
# for (d in dfs) {
#   d[[1]] <- as.character(d[[1]])  # Geneid as character
#   base <- merge(
#     base,
#     d[, c(1, 7:ncol(d)), drop = FALSE],  # Geneid + sample columns
#     by = "Geneid",
#     all = FALSE,
#     sort = FALSE
#   )
# }
#
# # Order columns: annotation first
# ann <- c("Geneid","Chr","Start","End","Strand","Length")
# base <- base[, c(ann, setdiff(colnames(base), ann))]
#
# # Always write colnames
# write.table(base, out, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# RSCRIPT
#
# echo "[OK] Merged counts written: $MERGED"
#
# # 5) (Optional) print header mapping to verify columns quickly
# echo "[INFO] Column header check:"
# awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) printf "%3d  %s\n", i,$i; exit}' "$MERGED"



#!/bin/bash
#SBATCH --job-name=fc_merge_only
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=../RNAseqMap/logs/fc_merge_only.%j.out
#SBATCH --error=../RNAseqMap/logs/fc_merge_only.%j.err

#!/bin/bash
#SBATCH --job-name=fc_merge_only
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=../RNAseqMap/logs/fc_merge_only.%j.out
#SBATCH --error=../RNAseqMap/logs/fc_merge_only.%j.err

set -euo pipefail
module load R

OUTDIR="/scratch/evt82290/MappingOutputs/RNAseq/counts"
SE_OUT="$OUTDIR/readcounts.SE.txt"   # set to your SE featureCounts file, or "" if none
PE_OUT="$OUTDIR/readcounts.PE.txt"   # set to your PE featureCounts file, or "" if none
MERGED="$OUTDIR/readcounts_All_CAF1paper.merged.txt"

Rscript - "$SE_OUT" "$PE_OUT" "$MERGED" <<'RS'
args <- commandArgs(trailingOnly=TRUE)
files <- args[1:2]
files <- files[file.exists(files) & nzchar(files)]
if (length(files) == 0) stop("No input featureCounts files found.")

read_fc <- function(f) {
  x <- read.table(f, header=TRUE, sep="\t", quote="", comment.char="#",
                  check.names=FALSE, stringsAsFactors=FALSE)
  colnames(x)[1] <- sub("^\ufeff", "", colnames(x)[1])  # strip BOM if present
  x
}

dfs <- lapply(files, read_fc)

# Start with the 6 annotation columns from the first file
base <- dfs[[1]][, 1:6]

# Append sample columns from each file
for (d in dfs) {
  base <- merge(base, d[, c(1, 7:ncol(d)), drop=FALSE],
                by="Geneid", all=FALSE, sort=FALSE)
}

# Write merged table with header
write.table(base, args[3], sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
RS
