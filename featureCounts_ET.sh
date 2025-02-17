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

cd $SLURM_SUBMIT_DIR

THREADS=2

#Make & Assign Directories
BAMDIR="/scratch/evt82290/RNAseq/CAF-1_Heatmap/bamFiles"
OUTDIR="/scratch/evt82290/RNAseq/CAF-1_Heatmap/counts"

# #if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
###
# /scratch/evt82290/RNAseq/CAF-1_Heatmap/bamFiles/SRR7970630/SRR7970630_Aligned.sortedByCoord.out.bam
# /scratch/evt82290/RNAseq/CAF-1_Heatmap/bamFiles/SRR7970630/SRR7970630_Aligned.sortedByCoord.out.bam
#Run featureCounts
module load Subread/

featureCounts -T $THREADS \
-t CDS \
-g gene_name \
-s 0 --primary \
-a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
-o $OUTDIR/readcounts_histchap.txt \
$BAMDIR/SRR8444037/SRR8444037_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR8444038/SRR8444038_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR8444043/SRR8444043_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR7970629/SRR7970629_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR7970630/SRR7970630_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR7970631/SRR7970631_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR7970598/SRR7970598_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR7970600/SRR7970600_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR8269830/SRR8269830_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR8269647/SRR8269647_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR8269650/SRR8269650_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR10916318/SRR10916318_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR10916319/SRR10916319_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR10916320/SRR10916320_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR7970603/SRR7970603_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR7970606/SRR7970606_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR7970610/SRR7970610_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR9027727/SRR9027727_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR9027728/SRR9027728_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR9027730/SRR9027730_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR8269825/SRR8269825_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR8269775/SRR8269775_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR8269782/SRR8269782_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR8269810/SRR8269810_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR10916163/SRR10916163_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR10916164/SRR10916164_Aligned.sortedByCoord.out.bam \
$BAMDIR/SRR10916165/SRR10916165_Aligned.sortedByCoord.out.bam
