#!/bin/bash
#SBATCH --job-name=cac_GATK
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --time=08:00:00
#SBATCH --output=../GATK/logs/%x.out
#SBATCH --error=../GATK/logs/%x.err

#change directroy & load GATK
cd $SLURM_SUBMIT_DIR
ml GATK
ml SAMtools

#set working & output directories
OUTDIR="/scratch/evt82290/Run136/GATK"
BAMDIR="/scratch/evt82290/MappingOutputs/Run136/bamFiles"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#Create dict & index files for genome (can comment out once done)
gatk CreateSequenceDictionary -R /home/evt82290/Research/GCA_000182925.2_NC12_genomic_wTetO_at_his3_CLEAN.fasta

samtools faidx /home/evt82290/Research/GCA_000182925.2_NC12_genomic_wTetO_at_his3_CLEAN.fasta

#HaplotypeCaller (-T chooses what tool you want to use)
gatk HaplotypeCaller -R /home/evt82290/Research/GCA_000182925.2_NC12_genomic_wTetO_at_his3_CLEAN.fasta -I ${BAMDIR}/6147_136-11_ChIP_WT_input.bam -O ${OUTDIR}/WT.vcf.gz\

#6147_136-11_ChIP_WT_input.b
#6147_136-12_ChIP_cac-1_input.bam
#6147_136-13_ChIP_cac-2_input.bam
#6147_136-14_ChIP_cac-3_input.bam
#147_136-92_ChIP_set-7_input_S91_L001_R1_001_val_1.fq.gz.bam


#2022_Run124_ET  2023_Run131_ET  2024_Run137_ET  2024_Run141_ET
#2022_Run126_ET  2023_Run133_ET  2024_Run138_ET  2024_Run144_ET
#2022_Run129_ET  2024_Run136_ET  2024_Run139_ET  2024_Run145_ET
