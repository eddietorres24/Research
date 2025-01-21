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

#set working & output directories
OUTDIR="/scratch/evt82290/Run136/GATK"
BAMDIR="/scratch/evt82290/Run136/SortedBamFiles"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#HaplotypeCaller (-T chooses what tool you want to use)
gatk HaplotypeCaller -R /home/evt82290/Research/Foxy_Ncrassa_Genome/Foxy_Ncrassa_merged.fasta -I ${BAMDIR}/6147_136-11_ChIP_WT_input.bam -O ${OUTDIR}/WT.vcf.gz\

#6147_136-11_ChIP_WT_input.bam
#6147_136-12_ChIP_cac-1_input.bam
#6147_136-13_ChIP_cac-2_input.bam
#6147_136-14_ChIP_cac-3_input.bam
#147_136-92_ChIP_set-7_input_S91_L001_R1_001_val_1.fq.gz.bam
