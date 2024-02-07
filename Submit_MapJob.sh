#!/bin/bash

#check for required command line argument

if [[ -z "$1" ]]
then
		echo "
		################
		ERROR: You must include the path to your accession file in the command line call.
		eg. sh submit.sh Path/To/Your/AccessionFile.txt
		##################
		"

		exit
fi

#iterates through list of accessions and passes to mapping script

fastqPath="/scratch/evt82290/histone_chaperone_plot_FINAL" #fastq directory generate by https://github.com/UGALewisLab/downloadSRA.git
outdir="/scratch/evt82290/RNAseqIGV"

mkdir ${outdir}
mkdir ${outdir}/logs

#make output file folders
mkdir "${outdir}/TrimmedFastQs"
mkdir "${outdir}/bamFiles"
mkdir "${outdir}/counts"
mkdir "${outdir}/bigWig"

while read -r line

	do
	sleep 10
	echo "$line mapping job submitted"
	sbatch --export=ALL,accession="${line}",fastqPath="${fastqPath}",outdir="${outdir}" MapRNAseq.sh & done <"$1"
