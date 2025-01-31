#!/bin/bash
#SBATCH --job-name=Call_Peaks_Homer_ET
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../Homer/logs/CallPeak.%j.out
#SBATCH --error=../Homer/logs/CallPeak.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $auto, )

#HOMER

#source config.txt

#Make Output Directory
auto="/scratch/evt82290/Peaks/qa-suz12"
BAMDIR="/scratch/evt82290/MappingOutputs/Run137/bamFiles"
BAMDIR2="/scratch/evt82290/MappingOutputs/Run144/bamFiles"
INPDIR="/scratch/evt82290/MappingOutputs/Run139/bamFiles"
TAGDIR="/scratch/evt82290/Peaks/qa-suz12/homertags"
#PEAKDIR="/scratch/evt82290/Run137/Peaks"

#if output directory doesn't exist, create it
if [ ! -d $auto ]
then
    mkdir -p $auto
fi
###

ml SAMtools
ml BWA
# #
# bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $autotempReps $auto "$bam" -
# samtools index "$bam"
#
# #samtools view -b -q 30 $bam > "$QualityBam"
# #samtools index "$QualityBam"
#
# ############################
# # # #deeptools

ml deepTools

#HOMER
module load Homer/4.11-foss-2022a

# 144-62_ChIP_WT_Input__S62.bam
# 144-69_ChIP_suz12_Input__S69.bam
# 144-6_ChIP_tetO_cac-1_p0_H3K27me3_Rep2_S6.bam
# 144-71_ChIP_qa-suz12_Input__S71.bam
# 144-73_ChIP_qa-suz12_Input__S73.bam
# 144-79_ChIP_qa-suz12_Input__S79.bam
# 144-132_ChIP_qa-suz12_Input__S132.bam

#Make Tag Directories
#Rep 1 qa-suz12
# makeTagDirectory $TAGDIR/WT_0hr_rep1 $BAMDIR/137-66_ChIP_WT_H3K27me3_Rep1_S63.bam
# makeTagDirectory $TAGDIR/qasuz12_0hr_rep1 $BAMDIR/137-67_ChIP_qa-suz12_H3K27me3_Rep1_S64.bam
# makeTagDirectory $TAGDIR/qasuz12_4hr_rep1 $BAMDIR/137-69_ChIP_qa-suz12_H3K27me3_Rep1_S66.bam
# makeTagDirectory $TAGDIR/qasuz12_8hr_rep1 $BAMDIR/137-71_ChIP_qa-suz12_H3K27me3_Rep1_S68.bam
# makeTagDirectory $TAGDIR/qasuz12_12hr_rep1 $BAMDIR/137-73_ChIP_qa-suz12_H3K27me3_Rep1_S70.bam
# makeTagDirectory $TAGDIR/qasuz12_24hr_rep1 $BAMDIR/137-75_ChIP_qa-suz12_H3K27me3_Rep1_S72.bam
# makeTagDirectory $TAGDIR/WT_24hr_rep1 $BAMDIR/137-74_ChIP_WT_H3K27me3_Rep1_S71.bam
# makeTagDirectory $TAGDIR/WT_0hr_input_rep1 $INPDIR/139-29_ChIP_WT_Input__S29.bam
# makeTagDirectory $TAGDIR/qasuz12_0hr_input_rep1 $INPDIR/139-30_ChIP_qa-suz12_Input__S30.bam
# makeTagDirectory $TAGDIR/qasuz12_4hr_input_rep1 $INPDIR/139-32_ChIP_qa-suz12_Input__S32.bam
# makeTagDirectory $TAGDIR/qasuz12_8hr_input_rep1 $INPDIR/139-34_ChIP_qa-suz12_Input__S34.bam
# makeTagDirectory $TAGDIR/qasuz12_12hr_input_rep1 $INPDIR/139-36_ChIP_qa-suz12_Input__S36.bam
# makeTagDirectory $TAGDIR/qasuz12_24hr_input_rep1 $INPDIR/139-38_ChIP_qa-suz12_Input__S38.bam
# makeTagDirectory $TAGDIR/WT_24hr_input_rep1 $INPDIR/139-37_ChIP_WT_Input__S37.bam
#
# #Rep 2 qa-suz12
# makeTagDirectory $TAGDIR/WT_0hr_rep2 $BAMDIR2/144-44_ChIP_WT_0hr_H3K27me3_Rep4_S44.bam
# makeTagDirectory $TAGDIR/qasuz12_0hr_rep2 $BAMDIR2/144-48_ChIP_qa-suz12_0hr_H3K27me3_Rep4_S48.bam
# makeTagDirectory $TAGDIR/qasuz12_4hr_rep2 $BAMDIR2/144-49_ChIP_qa-suz12_4hr_H3K27me3_Rep3_S49.bam
# makeTagDirectory $TAGDIR/qasuz12_8hr_rep2 $BAMDIR2/144-50_ChIP_qa-suz12_8hr_H3K27me3_Rep3_S50.bam
# makeTagDirectory $TAGDIR/qasuz12_12hr_rep2 $BAMDIR2/144-51_ChIP_qa-suz12_12hr_H3K27me3_Rep3_S51.bam
# makeTagDirectory $TAGDIR/qasuz12_24hr_rep2 $BAMDIR2/144-52_ChIP_qa-suz12_24hr_H3K27me3_Rep4_S52.bam
# makeTagDirectory $TAGDIR/WT_24hr_rep2 $BAMDIR2/144-45_ChIP_WT_24hr_H3K27me3_Rep3_S45.bam
# makeTagDirectory $TAGDIR/WT_0hr_input_rep2 $BAMDIR2/144-62_ChIP_WT_Input__S62.bam
# makeTagDirectory $TAGDIR/qasuz12_0hr_input_rep2 $BAMDIR2/144-71_ChIP_qa-suz12_Input__S71.bam
# makeTagDirectory $TAGDIR/qasuz12_4hr_input_rep2 $BAMDIR2/144-73_ChIP_qa-suz12_Input__S73.bam
# makeTagDirectory $TAGDIR/qasuz12_8hr_input_rep2 $BAMDIR2/144-79_ChIP_qa-suz12_Input__S79.bam
# #makeTagDirectory $TAGDIR/qasuz12_12hr_input_rep2 $BAMDIR2/
# makeTagDirectory $TAGDIR/qasuz12_24hr_input_rep2 $BAMDIR2/144-132_ChIP_qa-suz12_Input__S132.bam
#makeTagDirectory $TAGDIR/WT_24hr_input_rep2 $INPDIR/

#call Peaks
#rep 1

findPeaks $TAGDIR/WT_0hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_0hr_input_rep1
findPeaks $TAGDIR/qasuz12_0hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_0hr_input_rep1
findPeaks $TAGDIR/qasuz12_4hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_4hr_input_rep1
findPeaks $TAGDIR/qasuz12_8hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_8hr_input_rep1
findPeaks $TAGDIR/qasuz12_12hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_12hr_input_rep1
findPeaks $TAGDIR/qasuz12_24hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_24hr_input_rep1
findPeaks $TAGDIR/WT_24hr_rep1 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_24hr_input_rep1
#rep 2
findPeaks $TAGDIR/WT_0hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_0hr_input_rep2
findPeaks $TAGDIR/qasuz12_0hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_0hr_input_rep2
findPeaks $TAGDIR/qasuz12_4hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_4hr_input_rep2
findPeaks $TAGDIR/qasuz12_8hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_8hr_input_rep2
findPeaks $TAGDIR/qasuz12_12hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_12hr_input_rep1
findPeaks $TAGDIR/qasuz12_24hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/qasuz12_24hr_input_rep2
findPeaks $TAGDIR/WT_24hr_rep2 -style histone -size 500 -minDist 530 -o auto -i $TAGDIR/WT_24hr_input_rep1

#Find Motifs
# findMotifsGenome.pl qasuz12_4hr_peaks.bed /scratch/evt82290/Foxy_Ncrassa_merged.fasta $automotifs/4hr -size given
# findMotifsGenome.pl qasuz12_12hr_peaks.bed /scratch/evt82290/Foxy_Ncrassa_merged.fasta $automotifs/12hr -size given
# findMotifsGenome.pl qasuz12_8hr_no_telo_macs_peaks.bed /home/evt82290/Research/GCA_000182925.2_NC12_genomic_wTetO_at_his3_CLEAN.fasta $automotifs/nucsites -size given -len 8,9,10,11,12,13,14,15

#bedtools

# module load BEDTools

#Combining all overlapping peaks
# bedtools multiinter -header -i ${auto3}/2024_04_23_WT_peaks.bed \
#                                ${auto3}/2024_04_23_136_abcam_cac-1_peaks.bed \
#                                ${auto3}/2024_04_23_136_abcam_cac-2_peaks.bed \
#                                ${auto3}/2024_04_23_136_abcam_cac-3_peaks.bed \
#                                ${auto3}/2024_04_23_24hr_peaks.bed > ${auto3}/merge_peaks.txt
#
#
# bedtools sort -i ${auto3}/merged_sorted.bed > ${auto3}/merged_sorted_2.bed
# bedtools merge -i ${auto3}/merged_sorted_2.bed > ${auto3}/merged_file.txt

#determining which peaks overlap across peak files
# bedtools intersect -wa -wb \
#     -a ${auto3}/2024_04_23_WT_peaks.bed \
#     -b ${auto3}/2024_04_23_136_abcam_cac-1_peaks.bed ${auto3}/2024_04_23_136_abcam_cac-2_peaks.bed ${auto3}/2024_04_23_136_abcam_cac-3_peaks.bed ${auto3}/2024_04_23_24hr_peaks.bed \
#     -names cac-1 cac-2 cac-3 24hr \
#     -sorted > ${auto3}/intersect_peaks.txt

#Run137 bams
# 137-60_ChIP_WT_HA_Rep1_S57.bam
# 137-61_ChIP_ETX51-cac-1-HA_HA_Rep1_S58.bam
# 137-62_ChIP_ETX52-cac-1-HA_HA_Rep1_S59.bam
# 137-63_ChIP_WT_input__S60.bam
# 137-64_ChIP_ETX51-cac-1-HA_input__S61.bam
# 137-65_ChIP_ETX52-cac-1-HA_input__S62.bam
# 137-66_ChIP_WT_H3K27me3_Rep1_S63.bam
# 137-67_ChIP_qa-suz12_H3K27me3_Rep1_S64.bam
# 137-68_ChIP_WT_H3K27me3_Rep1_S65.bam
# 137-69_ChIP_qa-suz12_H3K27me3_Rep1_S66.bam
# 137-70_ChIP_WT_H3K27me3_Rep1_S67.bam
# 137-71_ChIP_qa-suz12_H3K27me3_Rep1_S68.bam
# 137-72_ChIP_WT_H3K27me3_Rep1_S69.bam
# 137-73_ChIP_qa-suz12_H3K27me3_Rep1_S70.bam
# 137-74_ChIP_WT_H3K27me3_Rep1_S71.bam
# 137-75_ChIP_qa-suz12_H3K27me3_Rep1_S72.bam

#Run139 bam
# 139-29_ChIP_WT_Input__S29.bam
# 139-30_ChIP_qa-suz12_Input__S30.bam
# 139-31_ChIP_WT_Input__S31.bam
# 139-32_ChIP_qa-suz12_Input__S32.bam
# 139-33_ChIP_WT_Input__S33.bam
# 139-34_ChIP_qa-suz12_Input__S34.bam
# 139-35_ChIP_WT_Input__S35.bam
# 139-36_ChIP_qa-suz12_Input__S36.bam
# 139-37_ChIP_WT_Input__S37.bam
# 139-38_ChIP_qa-suz12_Input__S38.bam
# 139-49_ChIP_qa-suz12_H3K27me3_Rep3_S40.bam
# 139-50_ChIP_qa-suz12-cac-1_H3K27me3_Rep2_S41.bam
# 139-51_ChIP_qa-suz12_H3K27me3_Rep1_S42.bam
# 139-52_ChIP_qa-suz12-cac-1_H3K27me3_Rep1_S43.bam
# 139-53_ChIP_qa-suz12_H3K27me3_Rep1_S44.bam
# 139-54_ChIP_qa-suz12-cac-1_H3K27me3_Rep1_S45.bam
# 139-55_ChIP_qa-suz12_H3K27me3_Rep3_S46.bam
# 139-56_ChIP_qa-suz12-cac-1_H3K27me3_Rep2_S47.bam
# 139-57_ChIP_WT_H3K27me3_Rep1_S48.bam
# 139-58_ChIP_qa-suz12_H3K27me3_Rep1_S49.bam
# 139-59_ChIP_qa-suz12-cac-1_H3K27me3_Rep1_S50.bam
# 139-60_ChIP_WT_H3K27me2me3_Rep1_S51.bam
# 139-61_ChIP_qa-suz12_H3K27me2me3_Rep1_S52.bam
# 139-62_ChIP_qa-suz12-cac-1_H3K27me2me3_Rep1_S53.bam
# 139-63_ChIP_qa-suz12_H3K27me2me3_Rep1_S54.bam
# 139-64_ChIP_qa-suz12-cac-1_H3K27me2me3_Rep1_S55.bam
# 139-65_ChIP_qa-suz12_H3K27me2me3_Rep1_S56.bam
# 139-66_ChIP_qa-suz12-cac-1_H3K27me2me3_Rep1_S57.bam
# 139-67_ChIP_qa-suz12_H3K27me2me3_Rep1_S58.bam
# 139-68_ChIP_qa-suz12-cac-1_H3K27me2me3_Rep1_S59.bam
# 139-69_ChIP_WT_H3K27me2me3_Rep1_S60.bam
# 139-70_ChIP_qa-suz12_H3K27me2me3_Rep1_S61.bam

# Run144 bams
# 144-10_ChIP_tetO_p0_H3K27me3_Rep3_S10.bam
# 144-11_ChIP_tetO_p1_H3K27me3_Rep3_S11.bam
# 144-122_ChIP_tetO_p3_input__S122.bam
# 144-124_ChIP_tetO_p1_input__S124.bam
# 144-125_ChIP_tetO_p2_input__S125.bam
# 144-126_ChIP_tetO_p3_input__S126.bam
# 144-127_ChIP_tetO_p0_tet_input__S127.bam
# 144-12_ChIP_tetO_p2_H3K27me3_Rep3_S12.bam
# 144-132_ChIP_qa-suz12_Input__S132.bam
# 144-13_ChIP_tetO_p3_H3K27me3_Rep3_S13.bam
# 144-141_ChIP_qa-suz12_24hr_H3K27me2_Rep1_S141.bam
# 144-14_ChIP_tetO_p0_H3K27me3_Rep4_S14.bam
# 144-15_ChIP_tetO_p1_H3K27me3_Rep4_S15.bam
# 144-16_ChIP_tetO_p2_H3K27me3_Rep4_S16.bam
# 144-17_ChIP_tetO_p3_H3K27me3_Rep4_S17.bam
# 144-18_ChIP_tetO_p0_H3K36me3_Rep1_S18.bam
# 144-19_ChIP_tetO_p1_H3K36me3_Rep1_S19.bam
# 144-1_ChIP_tetO_p0_H3K27me3_Rep2_S1.bam
# 144-20_ChIP_tetO_p2_H3K36me3_Rep1_S20.bam
# 144-21_ChIP_tetO_p3_H3K36me3_Rep1_S21.bam
# 144-22_ChIP_tetO_p0_H3K36me3_Rep2_S22.bam
# 144-23_ChIP_tetO_p1_H3K36me3_Rep2_S23.bam
# 144-24_ChIP_tetO_p2_H3K36me3_Rep2_S24.bam
# 144-25_ChIP_tetO_p3_H3K36me3_Rep2_S25.bam
# 144-26_ChIP_tetO_p0_tet_H3K27me3_Rep2_S26.bam
# 144-27_ChIP_tetO_p0_tet_H3K36me3_Rep2_S27.bam
# 144-28_ChIP_WT_0hr_H3K27me3_Rep3_S28.bam
# 144-29_ChIP_WT_24hr_H3K27me3_Rep2_S29.bam
# 144-2_ChIP_tetO_p1_H3K27me3_Rep2_S2.bam
# 144-30_ChIP_qa-suz12_0hr_H3K27me3_Rep3_S30.bam
# 144-31_ChIP_qa-suz12_4hr_H3K27me3_Rep2_S31.bam
# 144-32_ChIP_qa-suz12_8hr_H3K27me3_Rep2_S32.bam
# 144-33_ChIP_qa-suz12_12hr_H3K27me3_Rep2_S33.bam
# 144-34_ChIP_qa-suz12_24hr_H3K27me3_Rep3_S34.bam
# 144-35_ChIP_qa-suz12_48hr_H3K27me3_Rep2_S35.bam
# 144-37_ChIP_suz12_24hr_H3K27me3_Rep1_S37.bam
# 144-38_ChIP_qa-suz12_cac-1_0hr_H3K27me3_Rep2_S38.bam
# 144-39_ChIP_qa-suz12_cac-1_4hr_H3K27me3_Rep1_S39.bam
# 144-3_ChIP_tetO_p2_H3K27me3_Rep2_S3.bam
# 144-40_ChIP_qa-suz12_cac-1_8hr_H3K27me3_Rep1_S40.bam
# 144-41_ChIP_qa-suz12_cac-1_12hr_H3K27me3_Rep2_S41.bam
# 144-42_ChIP_qa-suz12_cac-1_24hr_H3K27me3_Rep2_S42.bam
# 144-43_ChIP_qa-suz12_cac-1_48hr_H3K27me3_Rep2_S43.bam
# 144-44_ChIP_WT_0hr_H3K27me3_Rep4_S44.bam
# 144-45_ChIP_WT_24hr_H3K27me3_Rep3_S45.bam
# 144-46_ChIP_suz12_0hr_H3K27me3_Rep2_S46.bam
# 144-47_ChIP_suz12_24hr_H3K27me3_Rep2_S47.bam
# 144-48_ChIP_qa-suz12_0hr_H3K27me3_Rep4_S48.bam
# 144-49_ChIP_qa-suz12_4hr_H3K27me3_Rep3_S49.bam
# 144-4_ChIP_tetO_p3_H3K27me3_Rep2_S4.bam
# 144-50_ChIP_qa-suz12_8hr_H3K27me3_Rep3_S50.bam
# 144-51_ChIP_qa-suz12_12hr_H3K27me3_Rep3_S51.bam
# 144-52_ChIP_qa-suz12_24hr_H3K27me3_Rep4_S52.bam
# 144-53_ChIP_WT_0hr_H3K27me2_Rep1_S53.bam
# 144-54_ChIP_WT_24hr_H3K27me2_Rep1_S54.bam
# 144-55_ChIP_suz12_0hr_H3K27me2_Rep1_S55.bam
# 144-56_ChIP_suz12_24hr_H3K27me2_Rep1_S56.bam
# 144-57_ChIP_qa-suz12_0hr_H3K27me2_Rep1_S57.bam
# 144-58_ChIP_qa-suz12_4hr_H3K27me2_Rep1_S58.bam
# 144-59_ChIP_qa-suz12_8hr_H3K27me2_Rep1_S59.bam
# 144-5_ChIP_tetO_p0_tet_H3K27me3_Rep2_S5.bam
# 144-60_ChIP_qa-suz12_12hr_H3K27me2_Rep1_S60.bam
# 144-62_ChIP_WT_Input__S62.bam
# 144-69_ChIP_suz12_Input__S69.bam
# 144-6_ChIP_tetO_cac-1_p0_H3K27me3_Rep2_S6.bam
# 144-71_ChIP_qa-suz12_Input__S71.bam
# 144-73_ChIP_qa-suz12_Input__S73.bam
# 144-79_ChIP_qa-suz12_Input__S79.bam
# 144-7_ChIP_tetO_cac-1_p1_H3K27me3_Rep2_S7.bam
# 144-8_ChIP_tetO_cac-1_p2_H3K27me3_Rep2_S8.bam
# 144-9_ChIP_tetO_cac-1_p3_H3K27me3_Rep2_S9.bam
