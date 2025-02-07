#!/bin/bash
#SBATCH --job-name=RNAseq_Map
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --time=8:00:00
#SBATCH --output=../RNAseqMap/logs/MapRNAseq.%j.out
#SBATCH --error=../RNAseqMap/logs/MapRNAseq.%j.err

cd $SLURM_SUBMIT_DIR

THREADS=2

##ADD a source file with path to FastqFiles
#variables imported from submission script
#accession=SRR10916163
fastqPath="/scratch/evt82290/SRA/FastqFiles"
outdir="/scratch/evt82290/RNAseq/cac_aberrant_transcripts"

# #if output directory doesn't exist, create it
if [ ! -d $outdir ]
then
    mkdir -p $outdir
fi
# ###


###################
#start
####################################################

#input file variables
  read1=${fastqPath}/${accession}/${accession}_1.fastq.gz
  read2=${fastqPath}/${accession}/${accession}_2.fastq.gz
  unpaired=${fastqPath}/${accession}/${accession}.fastq.gz

###################################

#make output file folders
trimmed="${outdir}/TrimmedFastQs/${accession}"
mkdir $trimmed

bamdir="${outdir}/bamFiles/${accession}"
mkdir "${bamdir}"

countsdir="${outdir}/counts/${accession}"
mkdir "${countsdir}"

bwDir="${outdir}/bigWig/${accession}"
mkdir "${bwDir}"

#pipeaccession summary: trim reads, map with STAR, get Counts

#notes
#generated a STAR genome index with the following call:
#STAR --runMode genomeGenerate --runThreadN 1 --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR --genomeFastaFiles /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_00182925.2plusHphplusBarplusTetO.fna --sjdbGTFfile /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_WithExtras_GFFtoGTFconversion.gtf
#need to rerun with normal genome assembly. The his-3 duplicated region will create multi-mappers


# #make output file folders
# trimmed= mkdir "${outdir}/TrimmedFastQs/${accession}"
#
# bamdir= mkdir "${outdir}/bamFiles/${accession}"
#
# countsdir= mkdir "${outdir}/counts/${accession}"
#
# bwDir= mkdir "${outdir}/bigWig/${accession}"

#make variables for output file names

bam="${bamdir}/${accession}_"
fowcounts="${countsdir}/${accession}_counts_fow.txt"
fowbw="${bwDir}/${accession}_fow.bw"
revcounts="${countsdir}/${accession}_counts_rev.txt"
revbw="${bwDir}/${accession}_rev.bw"

############# Read Trimming ##############
#remove adaptors, trim low quality reads (default = phred 20), length > 25

##fastq files from the ebi link are in folders that either have one file with a SRR##.fastq.gz or a SRR##_1.fastq.gz ending, or have two files with a SRR##_1.fastq.gz ending or a SRR##_2.fastq.gz ending
#or have three files with a SRR##_1.fastq.gz, SRR##_2.fastq.gz and SRR##.fastq.gz ending. In this case, the third file corresponds to unpaired reads that the depositers mapped.

#if read1 file does not exist, do single-end trimming using the only file in the folder i.e. SRR##.fastq.gz filename format
##This entire section can be simplified for JGI data

if [ ! -f $read1 ]; then
#trim reads
  echo "${line} running as unpaired file only"

  module load Trim_Galore/0.6.7-GCCcore-11.2.0

  trim_galore --illumina --fastqc --length 25 --basename ${accession} --gzip -o $trimmed $unpaired

  wait

#map with STAR
  module load STAR/2.7.10b-GCC-11.3.0

  STAR --runMode alignReads \
  --runThreadN $THREADS \
  --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
  --outFileNamePrefix ${bam} \
  --readFilesIn $trimmed/${accession}_trimmed.fq.gz \
  --readFilesCommand zcat \
  --alignIntronMax 10000 \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMsortingBinsN 100 \
  --outSAMunmapped Within \
  --outSAMattributes Standard \
  --limitBAMsortRAM 19990000000

#MAPPING READS BASED ON STRAND
#Load samtools
module load SAMtools/1.16.1-GCC-11.3.0

# Forward strand.
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
samtools view -b -f 128 -F 16 ${bam}Aligned.sortedByCoord.out.bam > ${bam}fow1.bam
samtools index ${bam}fow1.bam

samtools view -b -f 80 ${bam}Aligned.sortedByCoord.out.bam > ${bam}fow2.bam
samtools index ${bam}fow2.bam

# Combine alignments that originate on the forward strand.
samtools merge -f ${bam}fwd.bam ${bam}fow1.bam ${bam}fow2.bam
samtools index ${bam}fwd.bam

# Reverse strand
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
samtools view -b -f 144 ${bam}Aligned.sortedByCoord.out.bam > ${bam}rev1.bam
samtools index ${bam}rev1.bam

samtools view -b -f 64 -F 16 ${bam}Aligned.sortedByCoord.out.bam > ${bam}rev2.bam
samtools index ${bam}rev2.bam

# Combine alignments that originate on the reverse strand.
#
samtools merge -f ${bam}rev.bam rev1.bam ${bam}rev2.bam
samtools index ${bam}rev.bam

  #create index
  # module load SAMtools/1.16.1-GCC-11.3.0
  # samtools view -b -f 0x40 ${bam}Aligned.sortedByCoord.out.bam > ${bam}forward.bam
  # samtools view -b -f 0x80 ${bam}Aligned.sortedByCoord.out.bam > ${bam}reverse.bam
  # samtools index "${bam}forward.bam"
  # samtools index "${bam}reverse.bam"

  ##quantify with featureCounts
  module load Subread/2.0.6-GCC-11.3.0

  featureCounts -T $THREADS \
  -p \
  -t CDS \
  -g gene_name \
  -s 0 --primary \
  -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
  -o ${fowcounts} \
  ${bam}fwd.bam

  featureCounts -T $THREADS \
  -p \
  -t CDS \
  -g gene_name \
  -s 0 --primary \
  -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
  -o ${revcounts} \
  ${bam}rev.bam



  ##Plot reads to visualize tracks if needed
       module load deepTools/3.5.2-foss-2022a
       #Plot all reads
       bamCoverage -p $THREADS -bs 50 --normalizeUsing BPM -of bigwig -b "${bam}fwd.bam" -o "${fowbw}"
       bamCoverage -p $THREADS -bs 50 --normalizeUsing BPM -of bigwig -b "${bam}rev.bam" -o "${revbw}"


#elseif read2 exists, do paired-end Trimming and PE mapping
elif [ -f $read2 ]; then
  echo "${line} running as PE file only"


  ##################
  #Trimming
  #################
  	  # module load Trim_Galore/0.6.7-GCCcore-11.2.0
      #
  	  # trim_galore --illumina --fastqc --paired --length 25 --basename ${accession} --gzip -o $trimmed $read1 $read2
  	  # wait

#map with STAR
  # module load STAR/2.7.10b-GCC-11.3.0
  #
  # STAR --runMode alignReads \
  # --runThreadN $THREADS \
  # --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
  # --outFileNamePrefix ${bam} \
  # --readFilesIn $trimmed/${accession}_val_1.fq.gz $trimmed/${accession}_val_2.fq.gz \
  # --readFilesCommand zcat \
  # --alignIntronMax 10000 \
  # --outSAMtype BAM SortedByCoordinate \
  # --outBAMsortingBinsN 100 \
  # --outSAMunmapped Within \
  # --outSAMattributes Standard \
  # --limitBAMsortRAM 19990000000

  #MAPPING READS BASED ON STRAND
  #Load samtools
  module load SAMtools/1.16.1-GCC-11.3.0

  # Forward strand.
  # 1. alignments of the second in pair if they map to the forward strand
  # 2. alignments of the first in pair if they map to the reverse  strand
  samtools view -b -f 128 -F 16 ${bam}Aligned.sortedByCoord.out.bam > ${bam}fow1.bam
  samtools index ${bam}fow1.bam

  samtools view -b -f 80 ${bam}Aligned.sortedByCoord.out.bam > ${bam}fow2.bam
  samtools index ${bam}fow2.bam

  # Combine alignments that originate on the forward strand.
  samtools merge -f ${bam}fwd.bam ${bam}fow1.bam ${bam}fow2.bam
  samtools index ${bam}fwd.bam

  # Reverse strand
  # 1. alignments of the second in pair if they map to the reverse strand
  # 2. alignments of the first in pair if they map to the forward strand
  samtools view -b -f 144 ${bam}Aligned.sortedByCoord.out.bam > ${bam}rev1.bam
  samtools index ${bam}rev1.bam

  samtools view -b -f 64 -F 16 ${bam}Aligned.sortedByCoord.out.bam > ${bam}rev2.bam
  samtools index ${bam}rev2.bam

  # Combine alignments that originate on the reverse strand.
  #
  samtools merge -f ${bam}rev.bam rev1.bam ${bam}rev2.bam
  samtools index ${bam}rev.bam

  #create index
  # module load SAMtools/1.16.1-GCC-11.3.0
  # samtools view -b -f 0x40 ${bam}Aligned.sortedByCoord.out.bam > ${bam}forward.bam
  # samtools view -b -f 0x80 ${bam}Aligned.sortedByCoord.out.bam > ${bam}reverse.bam
  # samtools index "${bam}forward.bam"
  # samtools index "${bam}reverse.bam"

  ##quantify with featureCounts
  module load Subread/2.0.6-GCC-11.3.0

  featureCounts -T $THREADS \
  -p \
  -t CDS \
  -g gene_name \
  -s 0 --primary \
  -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
  -o ${fowcounts} \
  ${bam}fwd.bam

  featureCounts -T $THREADS \
  -p \
  -t CDS \
  -g gene_name \
  -s 0 --primary \
  -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
  -o ${revcounts} \
  ${bam}rev.bam



  ##Plot reads to visualize tracks if needed
       module load deepTools/3.5.2-foss-2022a
       #Plot all reads
       bamCoverage -p $THREADS -bs 50 --normalizeUsing BPM -of bigwig -b "${bam}fwd.bam" -o "${fowbw}"
       bamCoverage -p $THREADS -bs 50 --normalizeUsing BPM -of bigwig -b "${bam}rev.bam" -o "${revbw}"

#in rare cases there will only be a SRR##_1.fastq.gz format. Use this if nothing else exists.
else

    echo "${accesion} running as Read1 file only"

       trim_galore --illumina --fastqc --length 25 --basename ${accession} --gzip -o $trimmed $read1

       #map with STAR
         STAR --runMode alignReads \
         --runThreadN $THREADS \
         --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
         --outFileNamePrefix ${accession} \
         --readFilesIn ${accession}_trimmed.fq.gz  \
         --readFilesCommand zcat \
         --alignIntronMax 10000 \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --outBAMsortingBinsN 100 \
         --outSAMunmapped Within \
         --outSAMattributes Standard \
         --limitBAMsortRAM 19990000000


         #create index
         module load SAMtools/1.16.1-GCC-11.3.0
         samtools index "${bam}Aligned.sortedByCoord.out.bam"

         ##quantify with featureCounts
         module load Subread/2.0.6-GCC-11.3.0

         featureCounts -T $THREADS \
         -t CDS \
         -g gene_name \
         -s 0 --primary \
         -p \
         -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
         -o $counts \
         ${bam}Aligned.sortedByCoord.out.bam


         ##Plot reads to visualize tracks if needed
         	    module load deepTools/3.5.2-foss-2022a
         	    #Plot all reads
         	    bamCoverage -p $THREADS -bs 50 --normalizeUsing BPM -of bigwig -b "${bam}Aligned.sortedByCoord.out.bam" -o "${bw}"


fi


fi
