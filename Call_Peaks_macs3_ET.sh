#!/bin/bash
#SBATCH --job-name=Call_Peaks_macs3_ET
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../CallPeak.%j.out
#SBATCH --error=../CallPeak.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

#Make Output Directory
OUTDIR="/scratch/evt82290/Run137"
OUTDIR2="/scratch/evt82290/Run136"
OUTDIR3="/scratch/evt82290/Peaks"

#if output directory doesn't exist, create it
# if [ ! -d $OUTDIR ]
# then
#     mkdir -p $OUTDIR
# fi
###

#process reads using trimGalore

#  ml Trim_Galore/0.6.7-GCCcore-11.2.0
#  trim_galore --paired --length 20 --fastqc --gzip -o ${OUTDIR}/TrimmedReads ${FASTQ}/*fastq\.gz
# #
# FILES="${OUTDIR}/TrimmedReads/*R1_001_val_1\.fq\.gz" #Don't forget the *
#
# mkdir "${OUTDIR}/SortedBamFiles"
# mkdir "${OUTDIR}/BigWigs"
# mkdir "${OUTDIR}/Peaks"
# #mkdir "$OUTDIR/HomerTagDirectories"
# #mkdir "$OUTDIR/TdfFiles"
# #
# #Iterate over the files
# for f in $FILES
# do
# #
# # 	#Examples to Get Different parts of the file name
# # 		#See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
# 		#${string//substring/replacement}
# # 		#dir=${f%/*}
#
# 	file=${f##*/}
# 	#remove ending from file name to create shorter names for bam files and other downstream output
# 	name=${file/%_S[1-12]*_L001_R1_001_val_1.fq.gz/}
#
# #
# # 	# File Vars
# # 	#use sed to get the name of the second read matching the input file
# 	read2=$(echo "$f" | sed 's/R1_001_val_1\.fq\.gz/R2_001_val_2\.fq\.gz/g')
# 	#variable for naming bam file
#  	bam="${OUTDIR}/SortedBamFiles/${name}.bam"
# 	#variable name for bigwig output
# 	bigwig="${OUTDIR}/BigWigs/${name}"
# 	#QualityBam="${OUTDIR}/SortedBamFiles/${name}_Q30.bam"
# #

ml SAMtools
ml BWA
# #
# bwa mem -M -v 3 -t $THREADS $GENOME $f $read2 | samtools view -bhSu - | samtools sort -@ $THREADS -T $OUTDIR/SortedBamFiles/tempReps -o "$bam" -
# samtools index "$bam"
#
# #samtools view -b -q 30 $bam > "$QualityBam"
# #samtools index "$QualityBam"
#
# ############################
# # # #deeptools
#
ml deepTools
# #Plot all reads
# bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

#plot mononucleosomes (don't need to do for ChIP)
#bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"

#call Peaks
module load MACS3

#using --nolambda paramenter to call peaks without control
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/6147_136-2_ChIP_cac-1_H3K27me3_abcam_Rep2.bam" -c "${OUTDIR}/SortedBamFiles/6147_136-12_ChIP_cac-1_input.bam" -f BAMPE -n "136-2_ChIP_cac-1_H3K27me3_abcam_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/6147_136-8_ChIP_cac-2_H3K27me3_CS_Rep1_S8_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/6147_136-13_ChIP_cac-2_input.bam" -f BAMPE -n "136-8_ChIP_cac-2_H3K27me3_CS_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
#
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/6147_136-4_ChIP_cac-3_H3K27me3_abcam_Rep2_S4_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/6147_136-14_ChIP_cac-3_input.bam" -f BAMPE -n "136-4_ChIP_cac-3_H3K27me3_abcam_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/6147_136-9_ChIP_cac-3_H3K27me3_CS_Rep1_S9_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/6147_136-14_ChIP_cac-3_input.bam" -f BAMPE -n "136-9_ChIP_cac-3_H3K27me3_CS_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
#
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/6147_136-78_ChIP_WT_H3K27me3_CS_Rep2_S77_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/6147_136-84_ChIP_WT_input_S83_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "136-78_ChIP_WT_H3K27me3_CS_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/6147_136-79_ChIP_cac-1_H3K27me3_CS_Rep2_S78_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/6147_136-85_ChIP_cac-1_input_S84_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "136-79_ChIP_cac-1_H3K27me3_CS_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/6147_136-80_ChIP_cac-2_H3K27me3_CS_Rep2_S79_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/6147_136-89_ChIP_cac-2_input_S88_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "136-80_ChIP_cac-2_H3K27me3_CS_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/6147_136-81_ChIP_cac-3_H3K27me3_CS_Rep2_S80_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/6147_136-14_ChIP_cac-3_input.bam" -f BAMPE -n "136-81_ChIP_cac-3_H3K27me3_CS_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/6147_136-83_ChIP_set-7_H3K27me3_CS_Rep2_S82_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/6147_136-92_ChIP_set-7_input_S91_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "136-83_ChIP_set-7_H3K27me3_CS_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375

#Run137 callpeaks

# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/137-66_ChIP_WT_H3K27me3_Rep1_S63_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR2}/SortedBamFiles/6147_136-11_ChIP_WT_input.bam" -f BAMPE -n "137-66_ChIP_WT_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/137-67_ChIP_qa-suz12_H3K27me3_Rep1_S64_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR2}/SortedBamFiles/6147_136-11_ChIP_WT_input.bam" -f BAMPE -n "137-67_ChIP_qa-suz12_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.05 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/137-69_ChIP_qa-suz12_H3K27me3_Rep1_S66_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR2}/SortedBamFiles/6147_136-11_ChIP_WT_input.bam" -f BAMPE -n "137-69_ChIP_qa-suz12_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.05 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/137-71_ChIP_qa-suz12_H3K27me3_Rep1_S68_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR2}/SortedBamFiles/6147_136-11_ChIP_WT_input.bam" -f BAMPE -n "137-71_ChIP_qa-suz12_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.05 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/137-73_ChIP_qa-suz12_H3K27me3_Rep1_S70_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR2}/SortedBamFiles/6147_136-11_ChIP_WT_input.bam" -f BAMPE -n "137-73_ChIP_qa-suz12_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.05 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/137-75_ChIP_qa-suz12_H3K27me3_Rep1_S72_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR2}/SortedBamFiles/6147_136-11_ChIP_WT_input.bam" -f BAMPE -n "137-75_ChIP_qa-suz12_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.05 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375

#Run129
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/129-38_ChIP_WT_K27me3_AbC_Rep_1_S37_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/129-43_ChIP_WT_input_S42_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "129-38_ChIP_WT_H3K27me3_abcam_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/129-39_ChIP_cac-1_K27me3_AbC_Rep_1_S38_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/129-44_ChIP_cac-1_input_S43_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "129-39_ChIP_cac-1_H3K27me3_abcam_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/129-40_ChIP_cac-2_K27me3_AbC_Rep_1_S39_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/129-45_ChIP_cac-2_input_S44_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "129-40_ChIP_cac-2_H3K27me3_abcam_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/129-41_ChIP_cac-3_K27me3_AbC_Rep_1_S40_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/129-46_ChIP_cac-3_input_S45_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "129-41_ChIP_cac-3_H3K27me3_abcam_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/129-42_ChIP_set-7_K27me3_AbC_Rep_1_S41_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/129-47_ChIP_set-7_input_S46_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "129-42_ChIP_set-7_H3K27me3_abcam_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375

# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/129-38_ChIP_WT_K27me3_AbC_Rep_1_S37_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/129-43_ChIP_WT_input_S42_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "129-38_ChIP_WT_H3K27me3_abcam_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/129-39_ChIP_cac-1_K27me3_AbC_Rep_1_S38_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/129-44_ChIP_cac-1_input_S43_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "129-39_ChIP_cac-1_H3K27me3_abcam_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/129-40_ChIP_cac-2_K27me3_AbC_Rep_1_S39_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/129-45_ChIP_cac-2_input_S44_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "129-40_ChIP_cac-2_H3K27me3_abcam_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/129-41_ChIP_cac-3_K27me3_AbC_Rep_1_S40_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/129-46_ChIP_cac-3_input_S45_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "129-41_ChIP_cac-3_H3K27me3_abcam_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375
# macs3 callpeak -t "${OUTDIR}/SortedBamFiles/129-42_ChIP_set-7_K27me3_AbC_Rep_1_S41_L001_R1_001_val_1.fq.gz.bam" -c "${OUTDIR}/SortedBamFiles/129-47_ChIP_set-7_input_S46_L001_R1_001_val_1.fq.gz.bam" -f BAMPE -n "129-42_ChIP_set-7_H3K27me3_abcam_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375

#bedtools

module load BEDTools

#Combining all overlapping peaks
bedtools multiinter -header -i 2024_04_23_WT_peaks.bed \
                               2024_04_23_136_abcam_cac-1_peaks.bed \
                               2024_04_23_136_abcam_cac-2_peaks.bed \
                               2024_04_23_136_abcam_cac-3_peaks.bed \
                               2024_04_23_24hr_peaks.bed

#determining which peaks overlap across peak files
bedtools intersect -wa -wb \
    -a ${OUTDIR}/2024_04_23_WT_peaks.bed \
    -b 2024_04_23_136_abcam_cac-1_peaks.bed 2024_04_23_136_abcam_cac-2_peaks.bed 2024_04_23_136_abcam_cac-3_peaks.bed 2024_04_23_24hr_peaks.bed \
    -names cac-1 cac-2 cac-3 24hr \
    -sorted


# #Run129
# 124228-2.bam                                                            129-42_ChIP_set-7_K27me3_AbC_Rep_1_S41_L001_R1_001_val_1.fq.gz.bam
# 124228-2.bam.bai                                                        129-42_ChIP_set-7_K27me3_AbC_Rep_1_S41_L001_R1_001_val_1.fq.gz.bam.bai
# 124228-2_S2_L002_R1_001_val_1.fq.gz.bam                                 129-43_ChIP_WT_input_S42_L001_R1_001_val_1.fq.gz.bam
# 124228-2_S2_L002_R1_001_val_1.fq.gz.bam.bai                             129-43_ChIP_WT_input_S42_L001_R1_001_val_1.fq.gz.bam.bai
# 124228-2_S2_L003_R1_001_val_1.fq.gz.bam                                 129-44_ChIP_cac-1_input_S43_L001_R1_001_val_1.fq.gz.bam
# 124228-2_S2_L003_R1_001_val_1.fq.gz.bam.bai                             129-44_ChIP_cac-1_input_S43_L001_R1_001_val_1.fq.gz.bam.bai
# 124228-2_S2_L004_R1_001_val_1.fq.gz.bam                                 129-45_ChIP_cac-2_input_S44_L001_R1_001_val_1.fq.gz.bam
# 124228-2_S2_L004_R1_001_val_1.fq.gz.bam.bai                             129-45_ChIP_cac-2_input_S44_L001_R1_001_val_1.fq.gz.bam.bai
# 129-33_ChIP_WT_K27me3_AM_Rep_2_S32_L001_R1_001_val_1.fq.gz.bam          129-46_ChIP_cac-3_input_S45_L001_R1_001_val_1.fq.gz.bam
# 129-33_ChIP_WT_K27me3_AM_Rep_2_S32_L001_R1_001_val_1.fq.gz.bam.bai      129-46_ChIP_cac-3_input_S45_L001_R1_001_val_1.fq.gz.bam.bai
# 129-34_ChIP_cac-1_K27me3_AM_Rep_2_S33_L001_R1_001_val_1.fq.gz.bam       129-47_ChIP_set-7_input_S46_L001_R1_001_val_1.fq.gz.bam
# 129-34_ChIP_cac-1_K27me3_AM_Rep_2_S33_L001_R1_001_val_1.fq.gz.bam.bai   129-47_ChIP_set-7_input_S46_L001_R1_001_val_1.fq.gz.bam.bai
# 129-35_ChIP_cac-2_K27me3_AM_Rep_2_S34_L001_R1_001_val_1.fq.gz.bam       129-48_ChIP-seq_Dim-5_input_S47_L001_R1_001_val_1.fq.gz.bam
# 129-35_ChIP_cac-2_K27me3_AM_Rep_2_S34_L001_R1_001_val_1.fq.gz.bam.bai   129-48_ChIP-seq_Dim-5_input_S47_L001_R1_001_val_1.fq.gz.bam.bai
# 129-36_ChIP_cac-3_K27me3_AM_Rep_2_S35_L001_R1_001_val_1.fq.gz.bam       129-90_ChIP_WT_H3K36me3_Rep1_S71_L001_R1_001_val_1.fq.gz.bam
# 129-36_ChIP_cac-3_K27me3_AM_Rep_2_S35_L001_R1_001_val_1.fq.gz.bam.bai   129-90_ChIP_WT_H3K36me3_Rep1_S71_L001_R1_001_val_1.fq.gz.bam.bai
# 129-37_ChIP_set-7_K27me3_AM_Rep_2_S36_L001_R1_001_val_1.fq.gz.bam       129-91_ChIP_cac-1_H3K36me3_Rep1_S72_L001_R1_001_val_1.fq.gz.bam
# 129-37_ChIP_set-7_K27me3_AM_Rep_2_S36_L001_R1_001_val_1.fq.gz.bam.bai   129-91_ChIP_cac-1_H3K36me3_Rep1_S72_L001_R1_001_val_1.fq.gz.bam.bai
# 129-38_ChIP_WT_K27me3_AbC_Rep_1_S37_L001_R1_001_val_1.fq.gz.bam         129-92_ChIP_cac-2_H3K36me3_Rep1_S73_L001_R1_001_val_1.fq.gz.bam
# 129-38_ChIP_WT_K27me3_AbC_Rep_1_S37_L001_R1_001_val_1.fq.gz.bam.bai     129-92_ChIP_cac-2_H3K36me3_Rep1_S73_L001_R1_001_val_1.fq.gz.bam.bai
# 129-39_ChIP_cac-1_K27me3_AbC_Rep_1_S38_L001_R1_001_val_1.fq.gz.bam      129-93_ChIP_cac-3_H3K36me3_Rep1_S74_L001_R1_001_val_1.fq.gz.bam
# 129-39_ChIP_cac-1_K27me3_AbC_Rep_1_S38_L001_R1_001_val_1.fq.gz.bam.bai  129-93_ChIP_cac-3_H3K36me3_Rep1_S74_L001_R1_001_val_1.fq.gz.bam.bai
# 129-40_ChIP_cac-2_K27me3_AbC_Rep_1_S39_L001_R1_001_val_1.fq.gz.bam      129-94_ChIP_set-7_H3K36me3_Rep1_S75_L001_R1_001_val_1.fq.gz.bam
# 129-40_ChIP_cac-2_K27me3_AbC_Rep_1_S39_L001_R1_001_val_1.fq.gz.bam.bai  129-94_ChIP_set-7_H3K36me3_Rep1_S75_L001_R1_001_val_1.fq.gz.bam.bai
# 129-41_ChIP_cac-3_K27me3_AbC_Rep_1_S40_L001_R1_001_val_1.fq.gz.bam      MV.bam
# 129-41_ChIP_cac-3_K27me3_AbC_Rep_1_S40_L001_R1_001_val_1.fq.gz.bam.bai  MV.bam.bai

#Run137

# 137-60_ChIP_WT_HA_Rep1_S57_L001_R1_001_val_1.fq.gz.bam                  137-68_ChIP_WT_H3K27me3_Rep1_S65_L001_R1_001_val_1.fq.gz.bam
# 137-60_ChIP_WT_HA_Rep1_S57_L001_R1_001_val_1.fq.gz.bam.bai              137-68_ChIP_WT_H3K27me3_Rep1_S65_L001_R1_001_val_1.fq.gz.bam.bai
# 137-61_ChIP_ETX51-cac-1-HA_HA_Rep1_S58_L001_R1_001_val_1.fq.gz.bam      137-69_ChIP_qa-suz12_H3K27me3_Rep1_S66_L001_R1_001_val_1.fq.gz.bam
# 137-61_ChIP_ETX51-cac-1-HA_HA_Rep1_S58_L001_R1_001_val_1.fq.gz.bam.bai  137-69_ChIP_qa-suz12_H3K27me3_Rep1_S66_L001_R1_001_val_1.fq.gz.bam.bai
# 137-62_ChIP_ETX52-cac-1-HA_HA_Rep1_S59_L001_R1_001_val_1.fq.gz.bam      137-70_ChIP_WT_H3K27me3_Rep1_S67_L001_R1_001_val_1.fq.gz.bam
# 137-62_ChIP_ETX52-cac-1-HA_HA_Rep1_S59_L001_R1_001_val_1.fq.gz.bam.bai  137-70_ChIP_WT_H3K27me3_Rep1_S67_L001_R1_001_val_1.fq.gz.bam.bai
# 137-63_ChIP_WT_input__S60_L001_R1_001_val_1.fq.gz.bam                   137-71_ChIP_qa-suz12_H3K27me3_Rep1_S68_L001_R1_001_val_1.fq.gz.bam
# 137-63_ChIP_WT_input__S60_L001_R1_001_val_1.fq.gz.bam.bai               137-71_ChIP_qa-suz12_H3K27me3_Rep1_S68_L001_R1_001_val_1.fq.gz.bam.bai
# 137-64_ChIP_ETX51-cac-1-HA_input__S61_L001_R1_001_val_1.fq.gz.bam       137-72_ChIP_WT_H3K27me3_Rep1_S69_L001_R1_001_val_1.fq.gz.bam
# 137-64_ChIP_ETX51-cac-1-HA_input__S61_L001_R1_001_val_1.fq.gz.bam.bai   137-72_ChIP_WT_H3K27me3_Rep1_S69_L001_R1_001_val_1.fq.gz.bam.bai
# 137-65_ChIP_ETX52-cac-1-HA_input__S62_L001_R1_001_val_1.fq.gz.bam       137-73_ChIP_qa-suz12_H3K27me3_Rep1_S70_L001_R1_001_val_1.fq.gz.bam
# 137-65_ChIP_ETX52-cac-1-HA_input__S62_L001_R1_001_val_1.fq.gz.bam.bai   137-73_ChIP_qa-suz12_H3K27me3_Rep1_S70_L001_R1_001_val_1.fq.gz.bam.bai
# 137-66_ChIP_WT_H3K27me3_Rep1_S63_L001_R1_001_val_1.fq.gz.bam            137-74_ChIP_WT_H3K27me3_Rep1_S71_L001_R1_001_val_1.fq.gz.bam
# 137-66_ChIP_WT_H3K27me3_Rep1_S63_L001_R1_001_val_1.fq.gz.bam.bai        137-74_ChIP_WT_H3K27me3_Rep1_S71_L001_R1_001_val_1.fq.gz.bam.bai
# 137-67_ChIP_qa-suz12_H3K27me3_Rep1_S64_L001_R1_001_val_1.fq.gz.bam      137-75_ChIP_qa-suz12_H3K27me3_Rep1_S72_L001_R1_001_val_1.fq.gz.bam
# 137-67_ChIP_qa-suz12_H3K27me3_Rep1_S64_L001_R1_001_val_1.fq.gz.bam.bai  137-75_ChIP_qa-suz12_H3K27me3_Rep1_S72_L001_R1_001_val_1.fq.gz.bam.bai


# #Run136
# 6147_136-10_ChIP_cac-1-2_H3K27me3_CS_Rep1.bam                                                6147_136-45_ChIP_set-7_H3K27me3_abcam_6hr_Rep1_S45_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-10_ChIP_cac-1-2_H3K27me3_CS_Rep1.bam.bai                                            6147_136-4_ChIP_cac-3_H3K27me3_abcam_Rep2_S4_L001_R1_001_val_1.fq.gz.bam
# 6147_136-11_ChIP_WT_input.bam                                                                6147_136-4_ChIP_cac-3_H3K27me3_abcam_Rep2_S4_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-11_ChIP_WT_input.bam.bai                                                            6147_136-58_ChIP_WT_H3K27me3_abcam_Rep3_S58_L001_R1_001_val_1.fq.gz.bam
# 6147_136-12_ChIP_cac-1_input.bam                                                             6147_136-58_ChIP_WT_H3K27me3_abcam_Rep3_S58_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-12_ChIP_cac-1_input.bam.bai                                                         6147_136-5_ChIP_cac-1-2_H3K27me3_abcam_Rep2_S5_L001_R1_001_val_1.fq.gz.bam
# 6147_136-13_ChIP_cac-2_input.bam                                                             6147_136-5_ChIP_cac-1-2_H3K27me3_abcam_Rep2_S5_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-13_ChIP_cac-2_input.bam.bai                                                         6147_136-6_ChIP_WT_H3K27me3_CS_Rep1_S6_L001_R1_001_val_1.fq.gz.bam
# 6147_136-14_ChIP_cac-3_input.bam                                                             6147_136-6_ChIP_WT_H3K27me3_CS_Rep1_S6_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-14_ChIP_cac-3_input.bam.bai                                                         6147_136-71_ChIP_cac-1_H3K27me3_abcam_Rep3_S70_L001_R1_001_val_1.fq.gz.bam
# 6147_136-1_ChIP_WT_H3K27me3_abcam_Rep2.bam                                                   6147_136-71_ChIP_cac-1_H3K27me3_abcam_Rep3_S70_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-1_ChIP_WT_H3K27me3_abcam_Rep2.bam.bai                                               6147_136-72_ChIP_cac-2_H3K27me3_abcam_Rep3_S71_L001_R1_001_val_1.fq.gz.bam
# 6147_136-21_ChIP_cac-1-2_input.bam                                                           6147_136-72_ChIP_cac-2_H3K27me3_abcam_Rep3_S71_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-21_ChIP_cac-1-2_input.bam.bai                                                       6147_136-75_ChIP_cac-3_H3K27me3_abcam_Rep3_S74_L001_R1_001_val_1.fq.gz.bam
# 6147_136-22_ChIP_WT_H3K27me3_CS_0hr_Rep1.bam                                                 6147_136-75_ChIP_cac-3_H3K27me3_abcam_Rep3_S74_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-22_ChIP_WT_H3K27me3_CS_0hr_Rep1.bam.bai                                             6147_136-76_ChIP_cac-1-2_H3K27me3_abcam_Rep3_S75_L001_R1_001_val_1.fq.gz.bam
# 6147_136-23_ChIP_qa-suz12_H3K27me3_CS_0hr_Rep1.bam                                           6147_136-76_ChIP_cac-1-2_H3K27me3_abcam_Rep3_S75_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-23_ChIP_qa-suz12_H3K27me3_CS_0hr_Rep1.bam.bai                                       6147_136-77_ChIP_set-7_H3K27me3_abcam_Rep3_S76_L001_R1_001_val_1.fq.gz.bam
# 6147_136-24_ChIP_qa-suz12_cac-1_H3K27me3_CS_0hr_Rep1.bam                                     6147_136-77_ChIP_set-7_H3K27me3_abcam_Rep3_S76_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-24_ChIP_qa-suz12_cac-1_H3K27me3_CS_0hr_Rep1.bam.bai                                 6147_136-78_ChIP_WT_H3K27me3_CS_Rep2_S77_L001_R1_001_val_1.fq.gz.bam
# 6147_136-25_ChIP_WT_H3K27me3_CS_8hr_Rep1.bam                                                 6147_136-78_ChIP_WT_H3K27me3_CS_Rep2_S77_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-25_ChIP_WT_H3K27me3_CS_8hr_Rep1.bam.bai                                             6147_136-79_ChIP_cac-1_H3K27me3_CS_Rep2_S78_L001_R1_001_val_1.fq.gz.bam
# 6147_136-26_ChIP_qa-suz12_H3K27me3_CS_8hr_Rep1.bam                                           6147_136-79_ChIP_cac-1_H3K27me3_CS_Rep2_S78_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-26_ChIP_qa-suz12_H3K27me3_CS_8hr_Rep1.bam.bai                                       6147_136-7_ChIP_cac-1_H3K27me3_CS_Rep1_S7_L001_R1_001_val_1.fq.gz.bam
# 6147_136-27_ChIP_qa-suz12_cac-1_H3K27me3_CS_8hr_Rep1.bam                                     6147_136-7_ChIP_cac-1_H3K27me3_CS_Rep1_S7_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-27_ChIP_qa-suz12_cac-1_H3K27me3_CS_8hr_Rep1.bam.bai                                 6147_136-80_ChIP_cac-2_H3K27me3_CS_Rep2_S79_L001_R1_001_val_1.fq.gz.bam
# 6147_136-28_ChIP_WT_H3K27me3_CS_24hr_Rep1.bam                                                6147_136-80_ChIP_cac-2_H3K27me3_CS_Rep2_S79_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-28_ChIP_WT_H3K27me3_CS_24hr_Rep1.bam.bai                                            6147_136-81_ChIP_cac-3_H3K27me3_CS_Rep2_S80_L001_R1_001_val_1.fq.gz.bam
# 6147_136-29_ChIP_qa-suz12_H3K27me3_CS_24hr_Rep1.bam                                          6147_136-81_ChIP_cac-3_H3K27me3_CS_Rep2_S80_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-29_ChIP_qa-suz12_H3K27me3_CS_24hr_Rep1.bam.bai                                      6147_136-82_ChIP_cac-1-2_H3K27me3_CS_Rep2_S81_L001_R1_001_val_1.fq.gz.bam
# 6147_136-2_ChIP_cac-1_H3K27me3_abcam_Rep2.bam                                                6147_136-82_ChIP_cac-1-2_H3K27me3_CS_Rep2_S81_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-2_ChIP_cac-1_H3K27me3_abcam_Rep2.bam.bai                                            6147_136-83_ChIP_set-7_H3K27me3_CS_Rep2_S82_L001_R1_001_val_1.fq.gz.bam
# 6147_136-39_ChIP_qa-suz12_cac-1_H3K27me3_CS_24hr_Rep1_S39_L001_R1_001_val_1.fq.gz.bam        6147_136-83_ChIP_set-7_H3K27me3_CS_Rep2_S82_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-39_ChIP_qa-suz12_cac-1_H3K27me3_CS_24hr_Rep1_S39_L001_R1_001_val_1.fq.gz.bam.bai    6147_136-84_ChIP_WT_input_S83_L001_R1_001_val_1.fq.gz.bam
# 6147_136-3_ChIP_cac-2_H3K27me3_abcam_Rep2_S3_L001_R1_001_val_1.fq.gz.bam                     6147_136-84_ChIP_WT_input_S83_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-3_ChIP_cac-2_H3K27me3_abcam_Rep2_S3_L001_R1_001_val_1.fq.gz.bam.bai                 6147_136-85_ChIP_cac-1_input_S84_L001_R1_001_val_1.fq.gz.bam
# 6147_136-40_ChIP_WT_H3K27me3_abcam_6hr_Rep1_S40_L001_R1_001_val_1.fq.gz.bam                  6147_136-85_ChIP_cac-1_input_S84_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-40_ChIP_WT_H3K27me3_abcam_6hr_Rep1_S40_L001_R1_001_val_1.fq.gz.bam.bai              6147_136-89_ChIP_cac-2_input_S88_L001_R1_001_val_1.fq.gz.bam
# 6147_136-41_ChIP_qa-suz12_H3K27me3_abcam_6hr_Rep1_S41_L001_R1_001_val_1.fq.gz.bam            6147_136-89_ChIP_cac-2_input_S88_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-41_ChIP_qa-suz12_H3K27me3_abcam_6hr_Rep1_S41_L001_R1_001_val_1.fq.gz.bam.bai        6147_136-8_ChIP_cac-2_H3K27me3_CS_Rep1_S8_L001_R1_001_val_1.fq.gz.bam
# 6147_136-42_ChIP_qa-suz12_H3K27me3_abcam_6hr_Rep1_S42_L001_R1_001_val_1.fq.gz.bam            6147_136-8_ChIP_cac-2_H3K27me3_CS_Rep1_S8_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-42_ChIP_qa-suz12_H3K27me3_abcam_6hr_Rep1_S42_L001_R1_001_val_1.fq.gz.bam.bai        6147_136-91_ChIP_cac-1-2_input_S90_L001_R1_001_val_1.fq.gz.bam
# 6147_136-43_ChIP_qa-suz12_cac-1_H3K27me3_abcam_6hr_Rep1_S43_L001_R1_001_val_1.fq.gz.bam      6147_136-91_ChIP_cac-1-2_input_S90_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-43_ChIP_qa-suz12_cac-1_H3K27me3_abcam_6hr_Rep1_S43_L001_R1_001_val_1.fq.gz.bam.bai  6147_136-92_ChIP_set-7_input_S91_L001_R1_001_val_1.fq.gz.bam
# 6147_136-44_ChIP_qa-suz12_cac-1_H3K27me3_abcam_6hr_Rep1_S44_L001_R1_001_val_1.fq.gz.bam      6147_136-92_ChIP_set-7_input_S91_L001_R1_001_val_1.fq.gz.bam.bai
# 6147_136-44_ChIP_qa-suz12_cac-1_H3K27me3_abcam_6hr_Rep1_S44_L001_R1_001_val_1.fq.gz.bam.bai  6147_136-9_ChIP_cac-3_H3K27me3_CS_Rep1_S9_L001_R1_001_val_1.fq.gz.bam
# 6147_136-45_ChIP_set-7_H3K27me3_abcam_6hr_Rep1_S45_L001_R1_001_val_1.fq.gz.bam               6147_136-9_ChIP_cac-3_H3K27me3_CS_Rep1_S9_L001_R1_001_val_1.fq.gz.bam.bai
