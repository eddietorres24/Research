#!/bin/bash
#SBATCH --job-name=Call_Peaks_macs3_ET
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../MACS3/logs/CallPeak.%j.out
#SBATCH --error=../MACS3/logs/CallPeak.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($threads, $FASTQ, $OUTDIR, )

source config.txt

#Make Directories
#Output
OUTDIR1="/scratch/evt82290/Peaks/qa-suz12"
OUTDIR2="/scratch/evt82290/Peaks/CAF-1/ATAC"
OUTDIR3="/scratch/evt82290/Peaks/CAF-1/ChIP"
#bams
P126DIR="/scratch/evt82290/MappingOutputs/Run126/bamFiles"
P129DIR="/scratch/evt82290/MappingOutputs/Run129/bamFiles"
P131DIR="/scratch/evt82290/MappingOutputs/Run131/bamFiles"
P133DIR="/scratch/evt82290/MappingOutputs/Run133/bamFiles"
P136DIR="/scratch/evt82290/MappingOutputs/Run136/bamFiles"
P137DIR="/scratch/evt82290/MappingOutputs/Run137/bamFiles"
P138DIR="/scratch/evt82290/MappingOutputs/Run138/bamFiles"
P139DIR="/scratch/evt82290/MappingOutputs/Run139/bamFiles"
P141DIR="/scratch/evt82290/MappingOutputs/Run141/bamFiles"
P144DIR="/scratch/evt82290/MappingOutputs/Run144/bamFiles"
P145DIR="/scratch/evt82290/MappingOutputs/Run145/bamFiles"
P146DIR="/scratch/evt82290/MappingOutputs/Run146/bamFiles"
MISCDIR="/scratch/evt82290/MappingOutputs/iswi_ash1/bamFiles"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR3 ]
then
    mkdir -p $OUTDIR3
fi
###

#load MACS
module load MACS3

#Calling Peaks
###QA-SUZ12###

#masayuki isw peack calls for paper
#macs3 callpeak -t "${MISCDIR}/SRR11806698.bam" -f BAMPE -n "isw_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/isw" --min-length 650 --max-gap 375 --nolambda
#macs3 callpeak -t "${MISCDIR}/SRR11806688.bam" -f BAMPE -n "isw_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/isw" --min-length 650 --max-gap 375 --nolambda

#H3K27me3 24 hr
#Rep 1
# macs3 callpeak -t "${P137DIR}/137-66_ChIP_WT_H3K27me3_Rep1_S63.bam" -c "${P139DIR}/139-29_ChIP_WT_Input__S29.bam" -f BAMPE -n "WT_0hr_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/WT_0hr" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P137DIR}/137-67_ChIP_qa-suz12_H3K27me3_Rep1_S64.bam" -c "${P139DIR}/139-30_ChIP_qa-suz12_Input__S30.bam" -f BAMPE -n "qa-suz12_0hr_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/qa-suz12_0hr" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P137DIR}/137-69_ChIP_qa-suz12_H3K27me3_Rep1_S66.bam" -c "${P139DIR}/139-32_ChIP_qa-suz12_Input__S32.bam" -f BAMPE -n "qa-suz12_4hr_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/qa-suz12_4hr" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P137DIR}/137-71_ChIP_qa-suz12_H3K27me3_Rep1_S68.bam" -c "${P139DIR}/139-34_ChIP_qa-suz12_Input__S34.bam" -f BAMPE -n "qa-suz12_8hr_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/qa-suz12_8hr" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P137DIR}/137-73_ChIP_qa-suz12_H3K27me3_Rep1_S70.bam" -c "${P139DIR}/139-36_ChIP_qa-suz12_Input__S36.bam" -f BAMPE -n "qa-suz12_12hr_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/qa-suz12_12hr" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P137DIR}/137-75_ChIP_qa-suz12_H3K27me3_Rep1_S72.bam" -c "${P139DIR}/139-38_ChIP_qa-suz12_Input__S38.bam" -f BAMPE -n "qa-suz12_24hr_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/qa-suz12_24hr" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P137DIR}/137-74_ChIP_WT_H3K27me3_Rep1_S71.bam" -c "${P139DIR}/139-37_ChIP_WT_Input__S37.bam" -f BAMPE -n "WT_24hr_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/WT_24hr" --min-length 650 --max-gap 375
#
# #Rep 2
# macs3 callpeak -t "${P144DIR}/144-44_ChIP_WT_0hr_H3K27me3_Rep4_S44.bam" -c "${P144DIR}/144-62_ChIP_WT_Input__S62.bam" -f BAMPE -n "WT_0hr_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/WT_0hr" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P144DIR}/144-48_ChIP_qa-suz12_0hr_H3K27me3_Rep4_S48.bam" -c "${P144DIR}/144-71_ChIP_qa-suz12_Input__S71.bam" -f BAMPE -n "qa-suz12_0hr_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/qa-suz12_0hr" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P144DIR}/144-49_ChIP_qa-suz12_4hr_H3K27me3_Rep3_S49.bam" -c "${P144DIR}/144-73_ChIP_qa-suz12_Input__S73.bam" -f BAMPE -n "qa-suz12_4hr_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/qa-suz12_4hr" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P144DIR}/144-50_ChIP_qa-suz12_8hr_H3K27me3_Rep3_S50.bam" -c "${P144DIR}/144-79_ChIP_qa-suz12_Input__S79.bam" -f BAMPE -n "qa-suz12_8hr_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/qa-suz12_8hr" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P144DIR}/144-51_ChIP_qa-suz12_12hr_H3K27me3_Rep3_S51.bam" -c "${P139DIR}/139-36_ChIP_qa-suz12_Input__S36.bam" -f BAMPE -n "qa-suz12_12hr_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/qa-suz12_12hr" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P144DIR}/144-52_ChIP_qa-suz12_24hr_H3K27me3_Rep4_S52.bam" -c "${P144DIR}/144-132_ChIP_qa-suz12_Input__S132.bam" -f BAMPE -n "qa-suz12_24hr_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/qa-suz12_24hr" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P144DIR}/144-45_ChIP_WT_24hr_H3K27me3_Rep3_S45.bam" -c "${P139DIR}/139-38_ChIP_qa-suz12_Input__S38.bam" -f BAMPE -n "WT_24hr_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR}/WT_24hr" --min-length 650 --max-gap 375

#Rep 3
# macs3 callpeak -t "${P146DIR}/146-125_ChIP_qa-suz12_0hr_H3K27me3_Rep1_S144.bam" -c "${P139DIR}/139-30_ChIP_qa-suz12_Input__S30.bam" -f BAMPE -n "qa-suz12_0hr_H3K27me3_Rep3" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/qa-suz12_0hr" --min-length 650 --max-gap 250
# macs3 callpeak -t "${P146DIR}/146-126_ChIP_qa-suz12_4hr_H3K27me3_Rep1_S145.bam" -c "${P139DIR}/139-32_ChIP_qa-suz12_Input__S32.bam" -f BAMPE -n "qa-suz12_4hr_H3K27me3_Rep3" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/qa-suz12_4hr" --min-length 650 --max-gap 250
# macs3 callpeak -t "${P146DIR}/146-127_ChIP_qa-suz12_8hr_H3K27me3_Rep1_S146.bam" -c "${P139DIR}/139-34_ChIP_qa-suz12_Input__S34.bam" -f BAMPE -n "qa-suz12_8hr_H3K27me3_Rep3" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/qa-suz12_8hr" --min-length 650 --max-gap 250
# macs3 callpeak -t "${P146DIR}/146-128_ChIP_qa-suz12_12hr_H3K27me3_Rep1_S147.bam" -c "${P139DIR}/139-36_ChIP_qa-suz12_Input__S36.bam" -f BAMPE -n "qa-suz12_12hr_H3K27me3_Rep3" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/qa-suz12_12hr" --min-length 650 --max-gap 250
# macs3 callpeak -t "${P146DIR}/146-129_ChIP_qa-suz12_24hr_H3K27me3_Rep1_S148.bam" -c "${P139DIR}/139-38_ChIP_qa-suz12_Input__S38.bam" -f BAMPE -n "qa-suz12_24hr_H3K27me3_Rep3" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/qa-suz12_24hr" --min-length 650 --max-gap 250
# macs3 callpeak -t "${P146DIR}/ -c "${P146DIR}/146-18_ChIP_WT_input__S18.bam" -f BAMPE -n "WT_24hr_H3K27me3_Rep3" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/WT_24hr" --min-length 650 --max-gap 250

#H3K27me3 96 hr
#Rep 1
# macs3 callpeak -t "${P146DIR}/146-39_ChIP_qa-suz12_0hr_H3K27me3_Rep1_S39.bam" -c "${P139DIR}/139-30_ChIP_qa-suz12_Input__S30.bam" -f BAMPE -n "qa-suz12_0hr_H3K27me3_Rep4" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/qa-suz12_0hr" --min-length 650 --max-gap 275
# macs3 callpeak -t "${P146DIR}/146-40_ChIP_qa-suz12_24hr_H3K27me3_Rep1_S40.bam" -c "${P139DIR}/139-32_ChIP_qa-suz12_Input__S32.bam" -f BAMPE -n "qa-suz12_24hr_H3K27me3_Rep4" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/qa-suz12_24hr" --min-length 650 --max-gap 275
# macs3 callpeak -t "${P146DIR}/146-41_ChIP_qa-suz12_48hr_H3K27me3_Rep1_S41.bam" -c "${P139DIR}/139-34_ChIP_qa-suz12_Input__S34.bam" -f BAMPE -n "qa-suz12_48hr_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/qa-suz12_48hr" --min-length 650 --max-gap 275
# macs3 callpeak -t "${P146DIR}/146-42_ChIP_qa-suz12_72hr_H3K27me3_Rep1_S42.bam" -c "${P139DIR}/139-36_ChIP_qa-suz12_Input__S36.bam" -f BAMPE -n "qa-suz12_72hr_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/qa-suz12_72hr" --min-length 650 --max-gap 275
# macs3 callpeak -t "${P146DIR}/146-43_ChIP_qa-suz12_96hr_H3K27me3_Rep1_S43.bam" -c "${P139DIR}/139-38_ChIP_qa-suz12_Input__S38.bam" -f BAMPE -n "qa-suz12_96hr_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/qa-suz12_96hr" --min-length 650 --max-gap 275
# macs3 callpeak -t "${P146DIR}/146-44_ChIP_WT_H3K27me3_Rep1_S44.bam" -c "${P146DIR}/146-18_ChIP_WT_input__S18.bam" -f BAMPE -n "WT_96hr_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR1}/WT_96hr" --min-length 650 --max-gap 275





###CAF-1###

#ATAC-seq
# macs3 callpeak -t "${P141DIR}/141-N11_ATAC_WT__Rep1_S100.bam" -n "WT_ATAC_Rep1" --outdir "${OUTDIR2}/WT_ATAC_Rep1" -f BAMPE -g 41037538 -q 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all --nolambda
# macs3 callpeak -t "${P141DIR}/141-N12_ATAC_cac-1__Rep1_S101.bam" -n "cac-1_ATAC_Rep1" --outdir "${OUTDIR2}/cac-1_ATAC_Rep1" -f BAMPE -g 41037538 -q 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all --nolambda
# macs3 callpeak -t "${P141DIR}/141-N13_ATAC_cac-2__Rep1_S102.bam" -n "cac-2_ATAC_Rep1" --outdir "${OUTDIR2}/cac-2_ATAC_Rep1" -f BAMPE -g 41037538 -q 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all --nolambda
# macs3 callpeak -t "${P141DIR}/141-N14_ATAC_cac-3__Rep1_S103.bam" -n "cac-3_ATAC_Rep1" --outdir "${OUTDIR2}/cac-3_ATAC_Rep1" -f BAMPE -g 41037538 -q 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all --nolambda
# macs3 callpeak -t "${P141DIR}/141-N15_ATAC_set-7__Rep1_S104.bam" -n "set-7_ATAC_Rep1" --outdir "${OUTDIR2}/set-7_ATAC_Rep1" -f BAMPE -g 41037538 -q 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all --nolambda
# 146-N1_ATAC_WT__Rep2_S97.bam
# 146-N3_ATAC_cac-1__Rep2_S98.bam
# 146-N4_ATAC_cac-2__Rep2_S99.bam
# 146-N5_ATAC_set-7__Rep2_S100.bam
# 146-N6_ATAC_cac-3__Rep2_S101.bam
# macs3 callpeak -t "${P146DIR}/146-N1_ATAC_WT__Rep2_S97.bam" -n "WT_ATAC_Rep2" --outdir "${OUTDIR2}/WT_ATAC_Rep2" -f BAMPE -g 41037538 -q 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all --nolambda
# macs3 callpeak -t "${P146DIR}/146-N3_ATAC_cac-1__Rep2_S98.bam" -n "cac-1_ATAC_Rep2" --outdir "${OUTDIR2}/cac-1_ATAC_Rep2" -f BAMPE -g 41037538 -q 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all --nolambda
# macs3 callpeak -t "${P146DIR}/146-N4_ATAC_cac-2__Rep2_S99.bam" -n "cac-2_ATAC_Rep2" --outdir "${OUTDIR2}/cac-2_ATAC_Rep2" -f BAMPE -g 41037538 -q 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all --nolambda
# macs3 callpeak -t "${P146DIR}/146-N5_ATAC_set-7__Rep2_S100.bam" -n "set-7_ATAC_Rep2" --outdir "${OUTDIR2}/set-7_ATAC_Rep2" -f BAMPE -g 41037538 -q 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all --nolambda
# macs3 callpeak -t "${P146DIR}/146-N6_ATAC_cac-3__Rep2_S101.bam" -n "cac-3_ATAC_Rep2" --outdir "${OUTDIR2}/cac-3_ATAC_Rep2" -f BAMPE -g 41037538 -q 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all --nolambda

#ChIP-seq

#H3K27me3
#abcam Rep1
# macs3 callpeak -t "${P129DIR}/129-38_ChIP_WT_K27me3_AbC_Rep_1_S37.bam" -c "${P129DIR}/129-43_ChIP_WT_input_S42.bam" -f BAMPE -n "WT_abc_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/WT_abc_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P129DIR}/129-39_ChIP_cac-1_K27me3_AbC_Rep_1_S38.bam" -c "${P129DIR}/129-44_ChIP_cac-1_input_S43.bam" -f BAMPE -n "cac-1_abc_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-1_abc_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P129DIR}/129-40_ChIP_cac-2_K27me3_AbC_Rep_1_S39.bam" -c "${P129DIR}/129-45_ChIP_cac-2_input_S44.bam" -f BAMPE -n "cac-2_abc_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-2_abc_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P129DIR}/129-41_ChIP_cac-3_K27me3_AbC_Rep_1_S40.bam" -c "${P129DIR}/129-46_ChIP_cac-3_input_S45.bam" -f BAMPE -n "cac-3_abc_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-3_abc_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P129DIR}/129-42_ChIP_set-7_K27me3_AbC_Rep_1_S41.bam" -c "${P129DIR}/129-47_ChIP_set-7_input_S46.bam" -f BAMPE -n "set-7_abc_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/set-7_abc_K27" --min-length 650 --max-gap 375
#abcam Rep2
# macs3 callpeak -t "${P136DIR}/6147_136-1_ChIP_WT_H3K27me3_abcam_Rep2_S1.bam" -c "${P136DIR}/6147_136-11_ChIP_WT_input_S11.bam" -f BAMPE -n "WT_abc_H3K27me3_Rep2" -g 41037538 -p 0.01 --outdir "${OUTDIR3}/WT_abc_K27" --min-length 650 --max-gap 360
# macs3 callpeak -t "${P136DIR}/6147_136-2_ChIP_cac-1_H3K27me3_abcam_Rep2_S2.bam" -c "${P136DIR}/6147_136-12_ChIP_cac-1_input_S12.bam" -f BAMPE -n "cac-1_abc_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-1_abc_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P136DIR}/6147_136-3_ChIP_cac-2_H3K27me3_abcam_Rep2_S3.bam" -c "${P136DIR}/6147_136-13_ChIP_cac-2_input_S13.bam" -f BAMPE -n "cac-2_abc_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-2_abc_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P136DIR}/6147_136-4_ChIP_cac-3_H3K27me3_abcam_Rep2_S4.bam" -c "${P136DIR}/6147_136-14_ChIP_cac-3_input_S14.bam" -f BAMPE -n "cac-3_abc_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-3_abc_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P136DIR}/6147_136-77_ChIP_set-7_H3K27me3_abcam_Rep3_S76.bam" -c "${P136DIR}/6147_136-92_ChIP_set-7_input_S91.bam" -f BAMPE -n "set-7_abc_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/set-7_abc_K27" --min-length 650 --max-gap 375
# #CS Rep1
# macs3 callpeak -t "${P136DIR}/6147_136-6_ChIP_WT_H3K27me3_CS_Rep1_S6.bam" -c "${P136DIR}/6147_136-11_ChIP_WT_input_S11.bam" -f BAMPE -n "WT_CS_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/WT_CS_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P136DIR}/6147_136-7_ChIP_cac-1_H3K27me3_CS_Rep1_S7.bam" -c "${P136DIR}/6147_136-12_ChIP_cac-1_input_S12.bam" -f BAMPE -n "cac-1_CS_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-1_CS_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P136DIR}/6147_136-8_ChIP_cac-2_H3K27me3_CS_Rep1_S8.bam" -c "${P136DIR}/6147_136-13_ChIP_cac-2_input_S13.bam" -f BAMPE -n "cac-2_CS_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-2_CS_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P136DIR}/6147_136-9_ChIP_cac-3_H3K27me3_CS_Rep1_S9.bam" -c "${P136DIR}/6147_136-14_ChIP_cac-3_input_S14.bam" -f BAMPE -n "cac-3_CS_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-3_CS_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P126DIR}/126-58_ChIP_set-7_H3K27me3_CS_Rep3_S52.bam" -c "${P136DIR}/6147_136-92_ChIP_set-7_input_S91.bam" -f BAMPE -n "set-7_CS_H3K27me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/set-7_CS_K27" --min-length 650 --max-gap 375
# #CS Rep2
# macs3 callpeak -t "${P136DIR}/6147_136-78_ChIP_WT_H3K27me3_CS_Rep2_S77.bam" -c "${P136DIR}/6147_136-84_ChIP_WT_input_S83.bam" -f BAMPE -n "WT_CS_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/WT_CS_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P136DIR}/6147_136-79_ChIP_cac-1_H3K27me3_CS_Rep2_S78.bam" -c "${P136DIR}/6147_136-85_ChIP_cac-1_input_S84.bam" -f BAMPE -n "cac-1_CS_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-1_CS_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P136DIR}/6147_136-80_ChIP_cac-2_H3K27me3_CS_Rep2_S79.bam" -c "${P136DIR}/6147_136-89_ChIP_cac-2_input_S88.bam" -f BAMPE -n "cac-2_CS_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-2_CS_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P136DIR}/6147_136-81_ChIP_cac-3_H3K27me3_CS_Rep2_S80.bam" -c "${P136DIR}/6147_136-14_ChIP_cac-3_input_S14.bam" -f BAMPE -n "cac-3_CS_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-3_CS_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P136DIR}/6147_136-83_ChIP_set-7_H3K27me3_CS_Rep2_S82.bam" -c "${P138DIR}/138-76_ChIP_set-7_input__6252_S75.bam" -f BAMPE -n "set-7_CS_H3K27me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/set-7_CS_K27" --min-length 650 --max-gap 375
# #CS Rep3
# macs3 callpeak -t "${P138DIR}/138-57_ChIP_WT_H3K27me3_Rep3_6252_S56.bam" -c "${P138DIR}/138-72_ChIP_WT_input__6252_S71.bam" -f BAMPE -n "WT_CS_H3K27me3_Rep3" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/WT_CS_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P138DIR}/138-58_ChIP_cac-1_H3K27me3_Rep3_6252_S57.bam" -c "${P138DIR}/138-73_ChIP_cac-1_input__6252_S72.bam" -f BAMPE -n "cac-1_CS_H3K27me3_Rep3" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-1_CS_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P138DIR}/138-59_ChIP_cac-2_H3K27me3_Rep3_6252_S58.bam" -c "${P138DIR}/138-74_ChIP_cac-2_input__6252_S73.bam" -f BAMPE -n "cac-2_CS_H3K27me3_Rep3" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-2_CS_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P138DIR}/138-60_ChIP_cac-3_H3K27me3_Rep3_6252_S59.bam" -c "${P138DIR}/138-75_ChIP_cac-3_input__6252_S74.bam" -f BAMPE -n "cac-3_CS_H3K27me3_Rep3" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-3_CS_K27" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P138DIR}/138-61_ChIP_set-7_H3K27me3_Rep3_6252_S60.bam" -c "${P138DIR}/138-76_ChIP_set-7_input__6252_S75.bam" -f BAMPE -n "set-7_CS_H3K27me3_Rep3" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/set-7_CS_K27" --min-length 650 --max-gap 375
#
# #H3K36me3
# #Rep1
# macs3 callpeak -t "${P129DIR}/129-90_ChIP_WT_H3K36me3_Rep1_S71.bam" -c "${P138DIR}/138-72_ChIP_WT_input__6252_S71.bam" -f BAMPE -n "WT_H3K36me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/WT_H3K36" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P129DIR}/129-91_ChIP_cac-1_H3K36me3_Rep1_S72.bam" -c "${P138DIR}/138-73_ChIP_cac-1_input__6252_S72.bam" -f BAMPE -n "cac-1_H3K36me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-1_H3K36" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P129DIR}/129-92_ChIP_cac-2_H3K36me3_Rep1_S73.bam" -c "${P138DIR}/138-74_ChIP_cac-2_input__6252_S73.bam" -f BAMPE -n "cac-2_H3K36me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-2_H3K36" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P129DIR}/129-93_ChIP_cac-3_H3K36me3_Rep1_S74.bam" -c "${P138DIR}/138-75_ChIP_cac-3_input__6252_S74.bam" -f BAMPE -n "cac-3_H3K36me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-3_H3K36" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P129DIR}/129-94_ChIP_set-7_H3K36me3_Rep1_S75.bam" -c "${P138DIR}/138-76_ChIP_set-7_input__6252_S75.bam" -f BAMPE -n "set-7_H3K36me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/set-7_H3K36" --min-length 650 --max-gap 375
#Rep2
# macs3 callpeak -t "${P138DIR}/138-62_ChIP_WT_H3K36me3_Rep2_6252_S61.bam" -c "${P138DIR}/138-72_ChIP_WT_input__6252_S71.bam" -f BAMPE -n "WT_H3K36me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/WT_H3K36" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P138DIR}/138-63_ChIP_cac-1_H3K36me3_Rep2_6252_S62.bam" -c "${P138DIR}/138-73_ChIP_cac-1_input__6252_S72.bam" -f BAMPE -n "cac-1_H3K36me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-1_H3K36" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P138DIR}/138-64_ChIP_cac-2_H3K36me3_Rep2_6252_S63.bam" -c "${P138DIR}/138-74_ChIP_cac-2_input__6252_S73.bam" -f BAMPE -n "cac-2_H3K36me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-2_H3K36" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P138DIR}/138-65_ChIP_cac-3_H3K36me3_Rep2_6252_S64.bam" -c "${P138DIR}/138-75_ChIP_cac-3_input__6252_S74.bam" -f BAMPE -n "cac-3_H3K36me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-3_H3K36" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P138DIR}/138-66_ChIP_set-7_H3K36me3_Rep2_6252_S65.bam" -c "${P138DIR}/138-76_ChIP_set-7_input__6252_S75.bam" -f BAMPE -n "set-7_H3K36me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/set-7_H3K36" --min-length 650 --max-gap 375

#H4K20me3
#Rep1
# macs3 callpeak -t "${P133DIR}/133-19_ChIP_WT_H4K20me3_Rep1_S17.bam" -c "${P138DIR}/138-72_ChIP_WT_input__6252_S71.bam" -f BAMPE -n "WT_H4K20me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/WT_H4K20" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P133DIR}/133-20_ChIP_cac-1_H4K20me3_Rep1_S18.bam" -c "${P138DIR}/138-73_ChIP_cac-1_input__6252_S72.bam" -f BAMPE -n "cac-1_H4K20me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-1_H4K20" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P133DIR}/133-21_ChIP_cac-2_H4K20me3_Rep1_S19.bam" -c "${P138DIR}/138-74_ChIP_cac-2_input__6252_S73.bam" -f BAMPE -n "cac-2_H4K20me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-2_H4K20" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P133DIR}/" -c "${P131DIR}/131-40_ChIP_cac-3_input_Rep1_S30.bam" -f BAMPE -n "cac-3_H4K20me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-3_H4K20" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P133DIR}/133-22_ChIP_set-7_H4K20me3_Rep1_S20.bam" -c "${P138DIR}/138-76_ChIP_set-7_input__6252_S75.bam" -f BAMPE -n "set-7_H4K20me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/set-7_H4K20" --min-length 650 --max-gap 375
#Rep2
# macs3 callpeak -t "${P138DIR}/138-67_ChIP_WT_H4K20me3_Rep2_6252_S66.bam" -c "${P138DIR}/138-72_ChIP_WT_input__6252_S71.bam" -f BAMPE -n "WT_H4K20me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/WT_H4K20" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P138DIR}/138-68_ChIP_cac-1_H4K20me3_Rep2_6252_S67.bam" -c "${P138DIR}/138-73_ChIP_cac-1_input__6252_S72.bam" -f BAMPE -n "cac-1_H4K20me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-1_H4K20" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P138DIR}/138-69_ChIP_cac-2_H4K20me3_Rep2_6252_S68.bam" -c "${P138DIR}/138-74_ChIP_cac-2_input__6252_S73.bam" -f BAMPE -n "cac-2_H4K20me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-2_H4K20" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P138DIR}/138-70_ChIP_cac-3_H4K20me3_Rep2_6252_S69.bam" -c "${P138DIR}/138-75_ChIP_cac-3_input__6252_S74.bam" -f BAMPE -n "cac-3_H4K20me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-3_H4K20" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P138DIR}/138-71_ChIP_set-7_H4K20me3_Rep2_6252_S70.bam" -c "${P138DIR}/138-76_ChIP_set-7_input__6252_S75.bam" -f BAMPE -n "set-7_H4K20me3_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/set-7_H4K20" --min-length 650 --max-gap 375

#H3K4me2
# Rep1
# macs3 callpeak -t "${P146DIR}/146-19_ChIP_WT_H3K4me2_Rep1_S19.bam" -c "${P146DIR}/146-34_ChIP_WT_input__S34.bam" -f BAMPE -n "WT_H3K4me2_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/WT_H3K4" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P146DIR}/146-20_ChIP_cac-1_H3K4me2_Rep1_S20.bam" -c "${P146DIR}/146-35_ChIP_cac-1_input__S35.bam" -f BAMPE -n "cac-1_H3K4me2_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-1_H3K4" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P146DIR}/146-21_ChIP_cac-2_H3K4me2_Rep1_S21.bam" -c "${P146DIR}/146-36_ChIP_cac-2_input__S36.bam" -f BAMPE -n "cac-2_H3K4me2_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-2_H3K4" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P146DIR}/146-22_ChIP_cac-3_H3K4me2_Rep1_S22.bam" -c "${P146DIR}/146-37_ChIP_cac-3_input__S37.bam" -f BAMPE -n "cac-3_H3K4me2_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-3_H3K4" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P146DIR}/146-23_ChIP_set-7_H3K4me2_Rep1_S23.bam" -c "${P146DIR}/146-38_ChIP_set-7_input__S38.bam" -f BAMPE -n "set-7_H3K4me2_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/set-7_H3K4" --min-length 650 --max-gap 375
# #Rep2
# macs3 callpeak -t "${P146DIR}/146-29_ChIP_WT_H3K4me2_Rep2_S29.bam" -c "${P146DIR}/146-34_ChIP_WT_input__S34.bam" -f BAMPE -n "WT_H3K4me2_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/WT_H3K4" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P146DIR}/146-30_ChIP_cac-1_H3K4me2_Rep2_S30.bam" -c "${P146DIR}/146-35_ChIP_cac-1_input__S35.bam" -f BAMPE -n "cac-1_H3K4me2_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-1_H3K4" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P146DIR}/146-31_ChIP_cac-2_H3K4me2_Rep2_S31.bam" -c "${P146DIR}/146-36_ChIP_cac-2_input__S36.bam" -f BAMPE -n "cac-2_H3K4me2_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-2_H3K4" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P146DIR}/146-32_ChIP_cac-3_H3K4me2_Rep2_S32.bam" -c "${P146DIR}/146-37_ChIP_cac-3_input__S37.bam" -f BAMPE -n "cac-3_H3K4me2_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-3_H3K4" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P146DIR}/146-33_ChIP_set-7_H3K4me2_Rep2_S33.bam" -c "${P146DIR}/146-38_ChIP_set-7_input__S38.bam" -f BAMPE -n "set-7_H3K4me2_Rep2" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/set-7_H3K4" --min-length 650 --max-gap 375

#H3K9me3
#Rep1
# Rep2
# macs3 callpeak -t "${P133DIR}/133-23_ChIP_WT_H3K9me3_Rep2_S21.bam" -c "${P146DIR}/146-34_ChIP_WT_input__S34.bam" -f BAMPE -n "WT_H3K9me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/WT_H3K9" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P133DIR}/133-24_ChIP_cac-1_H3K9me3_Rep2_S22.bam" -c "${P146DIR}/146-35_ChIP_cac-1_input__S35.bam" -f BAMPE -n "cac-1_H3K9me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-1_H3K9" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P133DIR}/133-25_ChIP_cac-2_H3K9me3_Rep2_S23.bam" -c "${P146DIR}/146-36_ChIP_cac-2_input__S36.bam" -f BAMPE -n "cac-2_H3K9me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/cac-2_H3K9" --min-length 650 --max-gap 375
# macs3 callpeak -t "${P133DIR}/133-26_ChIP_set-7_H3K9me3_Rep2_S24.bam" -c "${P146DIR}/146-38_ChIP_set-7_input__S38.bam" -f BAMPE -n "set-7_H3K9me3_Rep1" --broad -g 41037538 --broad-cutoff 0.01 --outdir "${OUTDIR3}/set-7_H3K9" --min-length 650 --max-gap 375

#Convert Broadpeaks to bed format
#Rep 1
# cut -f 1-6 $OUTDIR1/WT_0hr_H3K27me3_Rep1_peaks.broadPeak > $OUTDIR1/WT_0hr_H3K27me3_Rep1_peaks.bed
# cut -f 1-6 $OUTDIR1/qa-suz12_0hr_H3K27me3_Rep1_peaks.broadPeak > $OUTDIR1/qa-suz12_0hr_H3K27me3_Rep1_peaks.bed
# cut -f 1-6 $OUTDIR1/qa-suz12_4hr_H3K27me3_Rep1_peaks.broadPeak > $OUTDIR1/qa-suz12_4hr_H3K27me3_Rep1_peaks.bed
# cut -f 1-6 $OUTDIR1/qa-suz12_8hr_H3K27me3_Rep1_peaks.broadPeak > $OUTDIR1/qa-suz12_8hr_H3K27me3_Rep1_peaks.bed
# cut -f 1-6 $OUTDIR1/qa-suz12_12hr_H3K27me3_Rep1_peaks.broadPeak > $OUTDIR1/qa-suz12_12hr_H3K27me3_Rep1_peaks.bed
# cut -f 1-6 $OUTDIR1/qa-suz12_24hr_H3K27me3_Rep1_peaks.broadPeak > $OUTDIR1/qa-suz12_24hr_H3K27me3_Rep1_peaks.bed
# cut -f 1-6 $OUTDIR1/WT_24hr_H3K27me3_Rep1_peaks.broadPeak > $OUTDIR1/WT_24hr_H3K27me3_Rep1_peaks.bed
#Rep 3
# cut -f 1-6 WT_H4K20me3_Rep2_peaks.broadPeak > WT_H4K20me3_peaks.bed
# cut -f 1-6 cac-1_H4K20me3_Rep2_peaks.broadPeak > cac-1_H4K20me3_peaks.bed
# cut -f 1-6 cac-2_H4K20me3_Rep2_peaks.broadPeak > cac-2_H4K20me3_peaks.bed
# cut -f 1-6 WT_H3K36me3_Rep2_peaks.broadPeak > WT_H3K36me3_peaks.bed
# cut -f 1-6 cac-1_H3K36me3_Rep2_peaks.broadPeak > cac-1_H3K36me3_peaks.bed
# cut -f 1-6 cac-2_H3K36me3_Rep2_peaks.broadPeak > cac-2_H3K36me3_peaks.bed

#bedtools
# module load BEDTools


# bedtools intersect -a -b -wa #overlap b/w a and b, keep sequence in file a
# bedtools intersect -a -b -v #no overlap b/w a and b, keep sequence in file a
# bedtools intersect -a qa-suz12_no_WT.bed -b WT_macs_24hr_rep2.bed -v > qa-suz12_no_WT_0_or_24.bed
# bedtools intersect -a internal_K27_cac3_peaks.bed -b WT_macs_0hr_rep2.bed -v > internal_K27_ectopic.bed
# bedtools intersect -a neurospora.bed -b WT_K27_narrow_rep2.bed -wa -f 0.9 > K27_mraked_genes.bed
# bedtools intersect -a CAF-1_ATAC_Peaks_merge_2.bed -b WT_ATAC_peaks.bed -wa > CAF1_ATAC_WT.bed
# bedtools intersect -a subtelomeric_K27_no_cac-3.bed -b WT_macs_0hr_rep2.bed -wa > subtelomeric_K27_normal.bed

# bedtools intersect -a CAF-1_K27_sorted.bed -b WT_CS_H3K27me3_Rep1_peaks.bed -v > CAF-1_Ectopic_K27.bed
# bedtools intersect -a CAF1_ATAC_WT.bed -b cac-2_ATAC_peaks.bed -v > WT_ATAC_NoCAF2.bed
# bedtools intersect -a CAF1_ATAC_WT.bed -b cac-3_ATAC_peaks.bed -v > WT_ATAC_NoCAF3.bed
# bedtools intersect -a CAF1_ATAC_WT.bed -b set-7_ATAC_peaks.bed -v > WT_ATAC_NoCAF4.bed

#Combining all overlapping peaks & merging
#CAF-1
# bedtools multiinter -header -i ${OUTDIR3}/2024_04_23_WT_peaks.bed \
#                                ${OUTDIR3}/2024_04_23_136_abcam_cac-1_peaks.bed \
#                                ${OUTDIR3}/2024_04_23_136_abcam_cac-2_peaks.bed > ${OUTDIR3}/merge_peaks.txt


#   cat cac-1_CS_H3K27me3_Rep1_peaks.bed \
#       cac-1_CS_H3K27me3_Rep1_peaks.bed \
#       cac-3_CS_H3K27me3_Rep1_peaks.bed > CAF-1_K27_Ectopic.bed
#
# #
# sort -k1,1 -k2,2n CAF-1_K27_Ectopic.bed > CAF-1_K27_Ectopic_sort.bed
# bedtools sort -i CAF-1_K27_Ectopic_sort.bed > CAF-1_K27_Ectopic_bed_sort.bed
#
# bedtools merge -i CAF-1_K27_Ectopic_bed_sort.bed -d 350 > CAF-1_K27_Ectopic_sorted.bed


#qa-suz12
#Rep1
# bedtools multiinter-bams ${OUTDIR1}/WT_0hr_H3K27me3_Rep1_peaks_sorted.bed \
#                                ${OUTDIR1}/qa-suz12_4hr_H3K27me3_Rep1_peaks_sorted.bed \
#                                ${OUTDIR1}/qa-suz12_8hr_H3K27me3_Rep1_peaks_sorted.bed \
#                                ${OUTDIR1}/qa-suz12_12hr_H3K27me3_Rep1_peaks_sorted.bed \
#                                ${OUTDIR1}/qa-suz12_24hr_H3K27me3_Rep1_peaks_sorted.bed \
#                                ${OUTDIR1}/WT_24hr_H3K27me3_Rep1_peaks_sorted.bed > ${OUTDIR1}/qa-suz12_rep1_overlap_peaks.bed

# cat WT_macs_0hr_rep2.bed \
#     qa-suz12_24hr_H3K27me3_Rep3_peaks.bed \
#     qa-suz12_12hr_H3K27me3_Rep3_peaks.bed > qa-suz12_WT_K27.bed
#
# sort -k1,1 -k2,2n qa-suz12_WT_K27.bed > qa-suz12_WT_K27_sorted.bed
# bedtools sort -i qa-suz12_WT_K27_sorted.bed > qa-suz12_WT_K27_bed_sorted.bed
#
# bedtools merge -i qa-suz12_WT_K27_bed_sorted.bed > qa-suz12_K27_merge_peaks.bed

# bedtools sort -i ${OUTDIR1}/merged_sorted.bed > ${OUTDIR1}/merged_sorted_2.bed
# bedtools merge -i ${OUTDIR1}/merged_sorted_2.bed > ${OUTDIR1}/merged_file.txt

########FILES#############
#Run129 bams
# 129-33_ChIP_WT_K27me3_AM_Rep_2_S32.bam
# 129-34_ChIP_cac-1_K27me3_AM_Rep_2_S33.bam
# 129-35_ChIP_cac-2_K27me3_AM_Rep_2_S34.bam
# 129-36_ChIP_cac-3_K27me3_AM_Rep_2_S35.bam
# 129-37_ChIP_set-7_K27me3_AM_Rep_2_S36.bam
# 129-38_ChIP_WT_K27me3_AbC_Rep_1_S37.bam
# 129-39_ChIP_cac-1_K27me3_AbC_Rep_1_S38.bam
# 129-40_ChIP_cac-2_K27me3_AbC_Rep_1_S39.bam
# 129-41_ChIP_cac-3_K27me3_AbC_Rep_1_S40.bam
# 129-42_ChIP_set-7_K27me3_AbC_Rep_1_S41.bam
# 129-43_ChIP_WT_input_S42.bam
# 129-44_ChIP_cac-1_input_S43.bam
# 129-45_ChIP_cac-2_input_S44.bam
# 129-46_ChIP_cac-3_input_S45.bam
# 129-47_ChIP_set-7_input_S46.bam
# 129-48_ChIP-seq_Dim-5_input_S47.bam
# 129-90_ChIP_WT_H3K36me3_Rep1_S71.bam
# 129-91_ChIP_cac-1_H3K36me3_Rep1_S72.bam
# 129-92_ChIP_cac-2_H3K36me3_Rep1_S73.bam
# 129-93_ChIP_cac-3_H3K36me3_Rep1_S74.bam
# 129-94_ChIP_set-7_H3K36me3_Rep1_S75.bam

#Run131 bams
# 131-22_ChIP_WT_H3K27me2me3_Rep1_S20.bam
# 131-23_ChIP_cac-1_H3K27me2me3_Rep1_S21.bam
# 131-24_ChIP_cac-2_H3K27me2me3_Rep1_S22.bam
# 131-25_ChIP_cac-3_H3K27me2me3_Rep1_S23.bam
# 131-26_ChIP_set-7_H3K27me2me3_Rep1_S24.bam
# 131-27_ChIP_WT_H3K36me2_Rep1_S25.bam
# 131-28_ChIP_cac-1_H3K36me2_Rep1_S26.bam
# 131-37_ChIP_WT_input_Rep1_S27.bam
# 131-38_ChIP_cac-1_input_Rep1_S28.bam
# 131-39_ChIP_cac-2_input_Rep1_S29.bam
# 131-40_ChIP_cac-3_input_Rep1_S30.bam
# 131-41_ChIP_set-7_input_Rep1_S31.bam
# 131-82_ChIP_WT_H3K27me3_Rep1_S70.bam
# 131-83_ChIP_cac-1_H3K27me3_Rep1_S71.bam
# 131-84_ChIP_cac-2_H3K27me3_Rep1_S72.bam
# 131-85_ChIP_cac-3_H3K27me3_Rep1_S73.bam
# 131-86_ChIP_set-7_H3K27me3_Rep1_S74.bam
# 131-87_ChIP_WT_H3K27me2me3_Rep3_S75.bam
# 131-88_ChIP_cac-1_H3K27me2me3_Rep3_S76.bam
# 131-89_ChIP_cac-2_H3K27me2me3_Rep3_S77.bam
# 131-90_ChIP_cac-3_H3K27me2me3_Rep3_S78.bam
# 131-91_ChIP_set-7_H3K27me2me3_Rep3_S79.bam

#Run133 bams
# 133-10_ChIP_cac-1_H3K27me2me3_17hr_Rep1_S10.bam
# 133-11_ChIP_cac-2_H3K27me2me3_17hr_Rep1_S11.bam
# 133-14_ChIP_set-7_H3K27me2me3_17hr_Rep1_S12.bam
# 133-15_ChIP_WT_H3K27me2me3_24hr_Rep1_S13.bam
# 133-16_ChIP_cac-1_H3K27me2me3_24hr_Rep1_S14.bam
# 133-17_ChIP_cac-2_H3K27me2me3_24hr_Rep1_S15.bam
# 133-18_ChIP_set-7_H3K27me2me3_24hr_Rep1_S16.bam
# 133-19_ChIP_WT_H4K20me3_Rep1_S17.bam
# 133-1_ChIP_WT_H3K27me2me3_4hr_Rep1_S1.bam
# 133-20_ChIP_cac-1_H4K20me3_Rep1_S18.bam
# 133-21_ChIP_cac-2_H4K20me3_Rep1_S19.bam
# 133-22_ChIP_set-7_H4K20me3_Rep1_S20.bam
# 133-23_ChIP_WT_H3K9me3_Rep2_S21.bam
# 133-24_ChIP_cac-1_H3K9me3_Rep2_S22.bam
# 133-25_ChIP_cac-2_H3K9me3_Rep2_S23.bam
# 133-26_ChIP_set-7_H3K9me3_Rep2_S24.bam
# 133-27_ChIP_WT_H3K4me2_Rep1_S25.bam
# 133-28_ChIP_cac-1_H3K4me2_Rep2_S26.bam
# 133-29_ChIP_cac-2_H3K4me2_Rep1_S27.bam
# 133-2_ChIP_cac-1_H3K27me2me3_4hr_Rep1_S2.bam
# 133-30_ChIP_set-7_H3K4me2_Rep1_S28.bam
# 133-31_ChIP_WT_H3K27me3_4hr_Rep1_S29.bam
# 133-32_ChIP_cac-1_H3K27me3_4hr_Rep1_S30.bam
# 133-33_ChIP_cac-2_H3K27me3_4hr_Rep1_S31.bam
# 133-34_ChIP_set-7_H3K27me3_4hr_Rep1_S32.bam
# 133-35_ChIP_WT_H3K27me3_8hr_Rep1_S33.bam
# 133-36_ChIP_cac-1_H3K27me3_8hr_Rep1_S34.bam
# 133-37_ChIP_cac-2_H3K27me3_8hr_Rep1_S35.bam
# 133-38_ChIP_set-7_H3K27me3_8hr_Rep1_S36.bam
# 133-3_ChIP_cac-2_H3K27me2me3_4hr_Rep1_S3.bam
# 133-40_ChIP_cac-1_H3K27me3_17hr_Rep1_S37.bam
# 133-41_ChIP_cac-2_H3K27me3_17hr_Rep1_S38.bam
# 133-42_ChIP_set-7_H3K27me3_17hr_Rep1_S39.bam
# 133-43_ChIP_WT_H3K27me3_24hr_Rep1_S40.bam
# 133-44_ChIP_cac-1_H3K27me3_24hr_Rep1_S41.bam
# 133-45_ChIP_cac-2_H3K27me3_24hr_Rep1_S42.bam
# 133-46_ChIP_set-7_H3K27me3_24hr_Rep1_S43.bam
# 133-47_ChIP_WT_Acetylated-Lysine_Rep1_S44.bam
# 133-48_ChIP_cac-1_Acetylated-Lysine_Rep1_S45.bam
# 133-49_ChIP_cac-2_Acetylated-Lysine_Rep1_S46.bam
# 133-4_ChIP_set-7_H3K27me2me3_4hr_Rep1_S4.bam
# 133-5_ChIP_WT_H3K27me2me3_8hr_Rep1_S5.bam
# 133-6_ChIP_cac-1_H3K27me2me3_8hr_Rep1_S6.bam
# 133-7_ChIP_cac-2_H3K27me2me3_8hr_Rep1_S7.bam
# 133-8_ChIP_set-7_H3K27me2me3_8hr_Rep1_S8.bam
# 133-95_ChIP_cac-3_Acetylated-Lysine_Rep1_S92.bam
# 133-96_ChIP_set-7_Acetylated-Lysine_Rep1_S93.bam
# 133-9_ChIP_WT_H3K27me2me3_17hr_Rep1_S9.bam

#Run136 bams
# 6147_136-10_ChIP_cac-1-2_H3K27me3_CS_Rep1_S10.bam
# 6147_136-11_ChIP_WT_input_S11.bam
# 6147_136-12_ChIP_cac-1_input_S12.bam
# 6147_136-13_ChIP_cac-2_input_S13.bam
# 6147_136-14_ChIP_cac-3_input_S14.bam
# 6147_136-1_ChIP_WT_H3K27me3_abcam_Rep2_S1.bam
# 6147_136-21_ChIP_cac-1-2_input_S21.bam
# 6147_136-22_ChIP_WT_H3K27me3_CS_0hr_Rep1_S22.bam
# 6147_136-23_ChIP_qa-suz12_H3K27me3_CS_0hr_Rep1_S23.bam
# 6147_136-24_ChIP_qa-suz12_cac-1_H3K27me3_CS_0hr_Rep1_S24.bam
# 6147_136-25_ChIP_WT_H3K27me3_CS_8hr_Rep1_S25.bam
# 6147_136-26_ChIP_qa-suz12_H3K27me3_CS_8hr_Rep1_S26.bam
# 6147_136-27_ChIP_qa-suz12_cac-1_H3K27me3_CS_8hr_Rep1_S27.bam
# 6147_136-28_ChIP_WT_H3K27me3_CS_24hr_Rep1_S28.bam
# 6147_136-29_ChIP_qa-suz12_H3K27me3_CS_24hr_Rep1_S29.bam
# 6147_136-2_ChIP_cac-1_H3K27me3_abcam_Rep2_S2.bam
# 6147_136-39_ChIP_qa-suz12_cac-1_H3K27me3_CS_24hr_Rep1_S39.bam
# 6147_136-3_ChIP_cac-2_H3K27me3_abcam_Rep2_S3.bam
# 6147_136-40_ChIP_WT_H3K27me3_abcam_6hr_Rep1_S40.bam
# 6147_136-41_ChIP_qa-suz12_H3K27me3_abcam_6hr_Rep1_S41.bam
# 6147_136-42_ChIP_qa-suz12_H3K27me3_abcam_6hr_Rep1_S42.bam
# 6147_136-43_ChIP_qa-suz12_cac-1_H3K27me3_abcam_6hr_Rep1_S43.bam
# 6147_136-44_ChIP_qa-suz12_cac-1_H3K27me3_abcam_6hr_Rep1_S44.bam
# 6147_136-45_ChIP_set-7_H3K27me3_abcam_6hr_Rep1_S45.bam
# 6147_136-4_ChIP_cac-3_H3K27me3_abcam_Rep2_S4.bam
# 6147_136-58_ChIP_WT_H3K27me3_abcam_Rep3_S58.bam
# 6147_136-5_ChIP_cac-1-2_H3K27me3_abcam_Rep2_S5.bam
# 6147_136-6_ChIP_WT_H3K27me3_CS_Rep1_S6.bam
# 6147_136-71_ChIP_cac-1_H3K27me3_abcam_Rep3_S70.bam
# 6147_136-72_ChIP_cac-2_H3K27me3_abcam_Rep3_S71.bam
# 6147_136-75_ChIP_cac-3_H3K27me3_abcam_Rep3_S74.bam
# 6147_136-76_ChIP_cac-1-2_H3K27me3_abcam_Rep3_S75.bam
# 6147_136-77_ChIP_set-7_H3K27me3_abcam_Rep3_S76.bam
# 6147_136-78_ChIP_WT_H3K27me3_CS_Rep2_S77.bam
# 6147_136-79_ChIP_cac-1_H3K27me3_CS_Rep2_S78.bam
# 6147_136-7_ChIP_cac-1_H3K27me3_CS_Rep1_S7.bam
# 6147_136-80_ChIP_cac-2_H3K27me3_CS_Rep2_S79.bam
# 6147_136-81_ChIP_cac-3_H3K27me3_CS_Rep2_S80.bam
# 6147_136-82_ChIP_cac-1-2_H3K27me3_CS_Rep2_S81.bam
# 6147_136-83_ChIP_set-7_H3K27me3_CS_Rep2_S82.bam
# 6147_136-84_ChIP_WT_input_S83.bam
# 6147_136-85_ChIP_cac-1_input_S84.bam
# 6147_136-89_ChIP_cac-2_input_S88.bam
# 6147_136-8_ChIP_cac-2_H3K27me3_CS_Rep1_S8.bam
# 6147_136-91_ChIP_cac-1-2_input_S90.bam
# 6147_136-92_ChIP_set-7_input_S91.bam
# 6147_136-9_ChIP_cac-3_H3K27me3_CS_Rep1_S9.bam

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

#Run138 bams
# 138-57_ChIP_WT_H3K27me3_Rep3_6252_S56.bam
# 138-58_ChIP_cac-1_H3K27me3_Rep3_6252_S57.bam
# 138-59_ChIP_cac-2_H3K27me3_Rep3_6252_S58.bam
# 138-60_ChIP_cac-3_H3K27me3_Rep3_6252_S59.bam
# 138-61_ChIP_set-7_H3K27me3_Rep3_6252_S60.bam
# 138-62_ChIP_WT_H3K36me3_Rep2_6252_S61.bam
# 138-63_ChIP_cac-1_H3K36me3_Rep2_6252_S62.bam
# 138-64_ChIP_cac-2_H3K36me3_Rep2_6252_S63.bam
# 138-65_ChIP_cac-3_H3K36me3_Rep2_6252_S64.bam
# 138-66_ChIP_set-7_H3K36me3_Rep2_6252_S65.bam
# 138-67_ChIP_WT_H4K20me3_Rep2_6252_S66.bam
# 138-68_ChIP_cac-1_H4K20me3_Rep2_6252_S67.bam
# 138-69_ChIP_cac-2_H4K20me3_Rep2_6252_S68.bam
# 138-70_ChIP_cac-3_H4K20me3_Rep2_6252_S69.bam
# 138-71_ChIP_set-7_H4K20me3_Rep2_6252_S70.bam
# 138-72_ChIP_WT_input__6252_S71.bam
# 138-73_ChIP_cac-1_input__6252_S72.bam
# 138-74_ChIP_cac-2_input__6252_S73.bam
# 138-75_ChIP_cac-3_input__6252_S74.bam
# 138-76_ChIP_set-7_input__6252_S75.bam
# 138-77_ChIP_WT_0hr_H3K27me3_Rep1_6252_S76.bam
# 138-78_ChIP_tetO-cac-1_0hr_H3K27me3_Rep1_6252_S77.bam
# 138-79_ChIP_tetO-cac-2_0hr_H3K27me3_Rep1_6252_S78.bam
# 138-80_ChIP_tetO-cac-3_0hr_H3K27me3_Rep1_6252_S79.bam
# 138-81_ChIP_tetO_0hr_H3K27me3_Rep1_6252_S80.bam
# 138-82_ChIP_WT_6hr_H3K27me3_Rep1_6252_S81.bam
# 138-83_ChIP_tetO-cac-1_6hr_H3K27me3_Rep1_6252_S82.bam
# 138-84_ChIP_tetO-cac-2_6hr_H3K27me3_Rep1_6252_S83.bam
# 138-85_ChIP_tetO-cac-3_6hr_H3K27me3_Rep1_6252_S84.bam
# 138-86_ChIP_tetO_6hr_H3K27me3_Rep1_6252_S85.bam
# 138-87_ChIP_WT_12hr_H3K27me3_Rep1_6252_S86.bam
# 138-88_ChIP_tetO-cac-1_12hr_H3K27me3_Rep1_6252_S87.bam
# 138-89_ChIP_tetO-cac-2_12hr_H3K27me3_Rep1_6252_S88.bam
# 138-90_ChIP_tetO-cac-3_12hr_H3K27me3_Rep1_6252_S89.bam
# 138-91_ChIP_tetO_12hr_H3K27me3_Rep1_6252_S90.bam
# 138-92_ChIP_WT_24hr_H3K27me3_Rep1_6252_S91.bam
# 138-93_ChIP_tetO-cac-1_24hr_H3K27me3_Rep1_6252_S92.bam
# 138-94_ChIP_tetO-cac-2_24hr_H3K27me3_Rep1_6252_S93.bam
# 138-95_ChIP_tetO-cac-3_24hr_H3K27me3_Rep1_6252_S94.bam
# 138-96_ChIP_tetO_24hr_H3K27me3_Rep1_6252_S95.bam

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

#Run141 bams
# 141-58_ChIP_WT_H3K23me1_Rep1_S52.bam
# 141-59_ChIP_cac-1_H3K23me1_Rep1_S53.bam
# 141-60_ChIP_cac-2_H3K23me1_Rep1_S54.bam
# 141-61_ChIP_cac-3_H3K23me1_Rep1_S55.bam
# 141-62_ChIP_set-7_H3K23me1_Rep1_S56.bam
# 141-63_ChIP_WT_teto_P0_H3K27me3_Rep1_S57.bam
# 141-64_ChIP_WT_teto_P1_H3K27me3_Rep1_S58.bam
# 141-65_ChIP_WT_teto_P2_H3K27me3_Rep1_S59.bam
# 141-66_ChIP_WT_teto_P3_H3K27me3_Rep1_S60.bam
# 141-67_ChIP_cac-1_teto_P0_H3K27me3_Rep1_S61.bam
# 141-68_ChIP_cac-1_teto_P1_H3K27me3_Rep1_S62.bam
# 141-69_ChIP_cac-1_teto_P2_H3K27me3_Rep1_S63.bam
# 141-70_ChIP_cac-1_teto_P3_H3K27me3_Rep1_S64.bam
# 141-88_ChIP_WT_H3K27me2_Rep1_S82.bam
# 141-89_ChIP_cac-1_H3K27me2_Rep1_S83.bam
# 141-90_ChIP_cac-2_H3K27me2_Rep1_S84.bam
# 141-91_ChIP_cac-3_H3K27me2_Rep1_S85.bam
# 141-92_ChIP_set-7_H3K27me2_Rep1_S86.bam
# 141-93_ChIP_WT_H4K12ac_Rep1_S87.bam
# 141-94_ChIP_cac-1_H4K12ac_Rep1_S88.bam
# 141-95_ChIP_cac-2_H4K12ac_Rep1_S89.bam
# 141-96_ChIP_set-7_H4K12ac_Rep1_S90.bam
# 141-N11_ATAC_WT__Rep1_S100.bam
# 141-N12_ATAC_cac-1__Rep1_S101.bam
# 141-N13_ATAC_cac-2__Rep1_S102.bam
# 141-N14_ATAC_cac-3__Rep1_S103.bam
# 141-N15_ATAC_set-7__Rep1_S104.bam


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

#Run145 bams
# 145-10_ChIP_csr1_LGVI_2_nuc_H3K27me3_Rep1_S10.bam
# 145-11_ChIP_WT_H3K27me3_nuc_S11.bam
# 145-15_ChIP_tetO_0hr_H3K27me3_Rep1_S12.bam
# 145-16_ChIP_tetO_6hr_H3K27me3_Rep1_S13.bam
# 145-17_ChIP_tetO_12hr_H3K27me3_Rep1_S14.bam
# 145-18_ChIP_tetO_18hr_H3K27me3_Rep1_S15.bam
# 145-19_ChIP_tetO_24hr_H3K27me3_Rep1_S16.bam
# 145-1_ChIP_csr1_LGIII_2_nuc_H3K27me3_Rep1_S1.bam
# 145-20_ChIP_tetO_30hr_H3K27me3_Rep1_S17.bam
# 145-21_ChIP_tetO_36hr_H3K27me3_Rep1_S18.bam
# 145-22_ChIP_tetO_cac-1_0hr_H3K27me3_Rep1_S19.bam
# 145-23_ChIP_tetO_cac-1_6hr_H3K27me3_Rep1_S20.bam
# 145-24_ChIP_tetO_cac-1_12hr_H3K27me3_Rep1_S21.bam
# 145-25_ChIP_tetO_cac-1_18hr_H3K27me3_Rep1_S22.bam
# 145-26_ChIP_tetO_cac-1_24hr_H3K27me3_Rep1_S23.bam
# 145-27_ChIP_tetO_cac-1_30hr_H3K27me3_Rep1_S24.bam
# 145-28_ChIP_tetO_cac-1_36hr_H3K27me3_Rep1_S25.bam
# 145-2_ChIP_csr1_LGIII_2_nuc_H3K27me3_Rep1_S2.bam
# 145-3_ChIP_csr1_LGV_1_nuc_H3K27me3_Rep1_S3.bam
# 145-4_ChIP_csr1_LGVI_1_nuc_H3K27me3_Rep1_S4.bam
# 145-5_ChIP_csr1_LGIII_1_nuc_H3K27me3_Rep1_S5.bam
# 145-60_ChIP_WT_0hr_input__S57.bam
# 145-6_ChIP_csr1_LGIII_2_nuc_H3K27me3_Rep1_S6.bam
# 145-7_ChIP_csr1_LGV_1_nuc_H3K27me3_Rep1_S7.bam
# 145-8_ChIP_csr1_LGVI_1_nuc_H3K27me3_Rep1_S8.bam
# 145-9_ChIP_csr1_LGVI_1_nuc_H3K27me3_Rep1_S9.bam

#Run146 bams
# 146-106_ChIP_epr-1GFP_GFP-trap_Rep1_S125.bam
# 146-10_ChIP_qa-suz12_72hr_H3K27me2_Rep1_S10.bam
# 146-117_ChIP_qa-suz12_96hr_H3K27me2_Rep1_S136.bam
# 146-118_ChIP_epr-1_GFP_H3K27me3_Rep1_S137.bam
# 146-11_ChIP_qa-suz12_96hr_H3K27me2_Rep1_S11.bam
# 146-123_ChIP_WT_H3K27me2_Rep1_S142.bam
# 146-124_ChIP_suz12_0hr_H3K27me3_Rep1_S143.bam
# 146-125_ChIP_qa-suz12_0hr_H3K27me3_Rep1_S144.bam
# 146-126_ChIP_qa-suz12_4hr_H3K27me3_Rep1_S145.bam
# 146-127_ChIP_qa-suz12_8hr_H3K27me3_Rep1_S146.bam
# 146-128_ChIP_qa-suz12_12hr_H3K27me3_Rep1_S147.bam
# 146-129_ChIP_qa-suz12_24hr_H3K27me3_Rep1_S148.bam
# 146-12_ChIP_WT_H3K27me2_Rep1_S12.bam
# 146-130_ChIP_suz12_0hr_H3K36me3_Rep1_S149.bam
# 146-131_ChIP_qa-suz12_0hr_H3K36me3_Rep1_S150.bam
# 146-132_ChIP_qa-suz12_4hr_H3K36me3_Rep1_S151.bam
# 146-133_ChIP_qa-suz12_8hr_H3K36me3_Rep1_S152.bam
# 146-134_ChIP_qa-suz12_12hr_H3K36me3_Rep1_S153.bam
# 146-135_ChIP_qa-suz12_24hr_H3K36me3_Rep1_S154.bam
# 146-136_ChIP_cac-2-epr-1-GFP_GFP-trap_Rep1_S155.bam
# 146-137_ChIP_cac-2-epr-1-GFP_GFP-trap_Rep1_S156.bam
# 146-13_ChIP_qa-suz12_0hr_input__S13.bam
# 146-14_ChIP_qa-suz12_24hr_input__S14.bam
# 146-15_ChIP_qa-suz12_48hr_input__S15.bam
# 146-16_ChIP_qa-suz12_72hr_input__S16.bam
# 146-17_ChIP_qa-suz12_96hr_input__S17.bam
# 146-18_ChIP_WT_input__S18.bam
# 146-19_ChIP_WT_H3K4me2_Rep1_S19.bam
# 146-1_ChIP_qa-suz12_0hr_H3K27me3_Rep4_S1.bam
# 146-20_ChIP_cac-1_H3K4me2_Rep1_S20.bam
# 146-21_ChIP_cac-2_H3K4me2_Rep1_S21.bam
# 146-22_ChIP_cac-3_H3K4me2_Rep1_S22.bam
# 146-23_ChIP_set-7_H3K4me2_Rep1_S23.bam
# 146-24_ChIP_WT_H3K27me2_Rep2_S24.bam
# 146-25_ChIP_cac-1_H3K27me2_Rep2_S25.bam
# 146-26_ChIP_cac-2_H3K27me2_Rep2_S26.bam
# 146-27_ChIP_cac-3_H3K27me2_Rep2_S27.bam
# 146-28_ChIP_set-7_H3K27me2_Rep2_S28.bam
# 146-29_ChIP_WT_H3K4me2_Rep2_S29.bam
# 146-2_ChIP_qa-suz12_24hr_H3K27me3_Rep4_S2.bam
# 146-30_ChIP_cac-1_H3K4me2_Rep2_S30.bam
# 146-31_ChIP_cac-2_H3K4me2_Rep2_S31.bam
# 146-32_ChIP_cac-3_H3K4me2_Rep2_S32.bam
# 146-33_ChIP_set-7_H3K4me2_Rep2_S33.bam
# 146-34_ChIP_WT_input__S34.bam
# 146-35_ChIP_cac-1_input__S35.bam
# 146-36_ChIP_cac-2_input__S36.bam
# 146-37_ChIP_cac-3_input__S37.bam
# 146-38_ChIP_set-7_input__S38.bam
# 146-39_ChIP_qa-suz12_0hr_H3K27me3_Rep1_S39.bam
# 146-3_ChIP_qa-suz12_48hr_H3K27me3_Rep1_S3.bam
# 146-40_ChIP_qa-suz12_24hr_H3K27me3_Rep1_S40.bam
# 146-41_ChIP_qa-suz12_48hr_H3K27me3_Rep1_S41.bam
# 146-42_ChIP_qa-suz12_72hr_H3K27me3_Rep1_S42.bam
# 146-43_ChIP_qa-suz12_96hr_H3K27me3_Rep1_S43.bam
# 146-44_ChIP_WT_H3K27me3_Rep1_S44.bam
# 146-45_ChIP_qa-suz12_0hr_H3K27me2_Rep1_S45.bam
# 146-46_ChIP_qa-suz12_24hr_H3K27me2_Rep1_S46.bam
# 146-47_ChIP_qa-suz12_48hr_H3K27me2_Rep1_S47.bam
# 146-48_ChIP_qa-suz12_72hr_H3K27me2_Rep1_S48.bam
# 146-4_ChIP_qa-suz12_72hr_H3K27me3_Rep1_S4.bam
# 146-5_ChIP_qa-suz12_96hr_H3K27me3_Rep1_S5.bam
# 146-6_ChIP_WT_H3K27me3_Rep1_S6.bam
# 146-7_ChIP_qa-suz12_0hr_H3K27me2_Rep1_S7.bam
# 146-8_ChIP_qa-suz12_24hr_H3K27me2_Rep1_S8.bam
# 146-9_ChIP_qa-suz12_48hr_H3K27me2_Rep1_S9.bam
# 146-N1_ATAC_WT__Rep2_S97.bam
# 146-N3_ATAC_cac-1__Rep2_S98.bam
# 146-N4_ATAC_cac-2__Rep2_S99.bam
# 146-N5_ATAC_set-7__Rep2_S100.bam
# 146-N6_ATAC_cac-3__Rep2_S101.bam
