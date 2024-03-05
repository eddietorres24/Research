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
OUTDIR="/scratch/evt82290/Run136"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
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

ml SAMtools/1.16.1-GCC-11.3.0
ml BWA/0.7.17-GCCcore-11.3.0
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
ml deepTools/3.5.2-foss-2022a
# #Plot all reads
# bamCoverage -p $THREADS -bs $BIN --normalizeUsing BPM --smoothLength $SMOOTH -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}Bulk.bw"

#plot mononucleosomes (don't need to do for ChIP)
#bamCoverage -p $THREADS --MNase -bs 1 --normalizeUsing BPM --smoothLength 25 -of bigwig -b "$bam" -o "${bigwig}.bin_${BIN}.smooth_${SMOOTH}_MNase.bw"

#call Peaks
module load MACS3/3.0.0b1-foss-2022a-Python-3.10.4

#using --nolambda paramenter to call peaks without control
macs3 callpeak -t "${OUTDIR}/SortedBamFiles/6147_136-1_ChIP_WT_H3K27me3_abcam_Rep2.bam" -c "${OUTDIR}/SortedBamFiles/6147_136-11_ChIP_WT_input.bam" -f BAMPE -n "136-1_ChIP_WT_H3K27me3_abcam_Rep2_TEST_12" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/Peaks" --min-length 650 --max-gap 375



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
