#!/bin/bash
#SBATCH --job-name=BCFtools
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --time=08:00:00
#SBATCH --output=../BCFtools/logs/%x.out
#SBATCH --error=../BCFtools/logs/%x.err

#change directory & load BCFtools
cd $SLURM_SUBMIT_DIR
ml BCFtools

#set working & output directories
OUTDIR="/scratch/evt82290/BCFtools/cac_analysis"
BAMDIR="/scratch/evt82290/MappingOutputs/Run136/bamFiles"

#if output directory doesn't exist, create it
if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

#Run bcftools
bcftools view --call -s -O z -o WT_v_cac
