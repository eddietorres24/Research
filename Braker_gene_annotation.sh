#!/bin/bash
#SBATCH --job-name=ET_Braker.%j.job
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50gb
#SBATCH --time=48:00:00
#SBATCH --output=../Braker.%j.out
#SBATCH --error=../Braker.%j.err

cd $SLURM_SUBMIT_DIR

#read in variables from the config file ($GENOME, fastqpath)

source config.txt

#load modules
module load BRAKER/3.0.8-foss-2022a
module load GeneMark-ET/4.72-GCCcore-12.3.0
module load AUGUSTUS
module load SAMtools

#Run BRAKER
braker.pl --genome ${GENOME} --bam /scratch/evt82290/RNAseq/CAF-1_Heatmap/bamFiles/SRR7970598/SRR7970598_Aligned.sortedByCoord.out.bam --softmasking --gff3 --species Neurospora_crassa_cac3 --cores 8 --workingdir /scratch/evt82290/BRAKER
