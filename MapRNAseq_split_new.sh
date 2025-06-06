#!/bin/bash
#SBATCH --job-name=RNAseq_Map
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
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
Run149path="/scratch/evt82290/FastqFiles/2025_Run149_ET/RNA"
outdir="/scratch/evt82290/RNAseq/cac_aberrant_transcripts"

# #if output directory doesn't exist, create it
if [ ! -d $outdir ]
then
    mkdir -p $outdir
fi

###################
#start
####################################################

#input file variables
  read1=${fastqPath}/${accession}_1.fastq.gz
  read2=${fastqPath}/${accession}_2.fastq.gz
  unpaired=${fastqPath}/${accession}.fastq.gz

# Check if at least one of the expected files exists in the primary path
  if [[ ! -f "$read1" && ! -f "$read2" && ! -f "$unpaired" ]]; then
    # Fallback to alternate path
    read1="${Run149path}/${accession}*R1_001.fastq.gz"
    read2="${Run149path}/${accession}*R2_001.fastq.gz"
    unpaired="${Run149path}/${accession}.fastq.gz"
  fi

###################################

#make output file folders
trimmed="${outdir}/TrimmedFastQs/${accession}"
mkdir $trimmed

bamdir="${outdir}/bamFiles/${accession}"
mkdir "${bamdir}"

countsdir="${outdir}/counts/${accession}"
mkdir "${countsdir}"

countsdir2="${outdir}/counts2/${accession}"
mkdir "${countsdir2}"

bwDir="${outdir}/bigWig/${accession}"
mkdir "${bwDir}"

bedDir="${outdir}/beds/${accession}"
mkdir "${bedDir}"

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

sencounts="${countsdir2}/${accession}_counts_sen.txt"
senbw="${bwDir}/${accession}_sen.bw"
antcounts="${countsdir2}/${accession}_counts_ant.txt"
antbw="${bwDir}/${accession}_ant.bw"

############# Read Trimming ##############
#remove adaptors, trim low quality reads (default = phred 20), length > 25

##fastq files from the ebi link are in folders that either have one file with a SRR##.fastq.gz or a SRR##_1.fastq.gz ending, or have two files with a SRR##_1.fastq.gz ending or a SRR##_2.fastq.gz ending
#or have three files with a SRR##_1.fastq.gz, SRR##_2.fastq.gz and SRR##.fastq.gz ending. In this case, the third file corresponds to unpaired reads that the depositers mapped.

#if read1 file does not exist, do single-end trimming using the only file in the folder i.e. SRR##.fastq.gz filename format
##This entire section can be simplified for JGI data

if [ ! -f $read1 ]; then
#trim reads
  echo "${line} running as unpaired file only"

  # module load Trim_Galore/0.6.7-GCCcore-11.2.0
  #
  # trim_galore --illumina --fastqc --length 25 --basename ${accession} --gzip -o $trimmed $unpaired
  #
  # wait

#map with STAR
  module load STAR/2.7.10b-GCC-11.3.0

#OG STAR mapping script
  # STAR --runMode alignReads \
  # --runThreadN $THREADS \
  # --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
  # --outFileNamePrefix ${bam} \
  # --readFilesIn $trimmed/${accession}_trimmed.fq.gz \
  # --readFilesCommand zcat \
  # --alignIntronMax 10000 \
  # --outSAMtype BAM SortedByCoordinate \
  # --outBAMsortingBinsN 100 \
  # --outSAMunmapped Within \
  # --outSAMattributes Standard \
  # --limitBAMsortRAM 19990000000

#Altered STAR script for antisense mapping w/ help from ChatGPT
  # STAR --runMode alignReads \
  #   --runThreadN $THREADS \
  #   --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
  #   --readFilesIn $trimmed/${accession}_trimmed.fq.gz \
  #   --readFilesCommand zcat \
  #   --outFileNamePrefix ${bam} \
  #   --outSAMtype BAM SortedByCoordinate \
  #   --twopassMode Basic \
  #   --outSAMstrandField intronMotif \
  #   --alignIntronMax 10000 \
  #   --outBAMsortingBinsN 100 \
  #   --outSAMunmapped Within \
  #   --outSAMattributes Standard \
  #   --limitBAMsortRAM 19990000000 \
  #   --quantMode TranscriptomeSAM GeneCounts

    #ChatGPT helping me make bigwigs and count files with antisnese reads only

    #load modules
    ml BEDTools
    module load SAMtools/1.16.1-GCC-11.3.0
    module load StringTie
    module load deepTools/3.5.2-foss-2022a

          # Set input and output paths
          BAM=${bam}Aligned.sortedByCoord.out.bam
          ANNOT=/home/evt82290/Research/annotation.bed    # BED12 gene annotation file
          GENOME="/home/ad45368/chrom_sizes.txt"   # Chromosome sizes

          # 1. Convert BAM to BED (splice-aware)
          # bedtools bamtobed -split -i $BAM > $bedDir/${accession}_all_reads.bed

          # 2. Extract antisense reads (opposite strand of gene annotation)
          # bedtools intersect -S -wa -a $bedDir/${accession}_all_reads.bed -b $ANNOT >  $bedDir/${accession}_antisense_reads.bed
          #
          # # 3. Convert antisense BED to BAM
          # bedtools bedtobam -i $bedDir/${accession}_antisense_reads.bed -g $GENOME > $bamdir/${accession}_antisense.bam

# 2. Extract reads antisense to genes on the + strand
bedtools intersect -s -wa -abam $BAM -b <(awk '$6 == "+"' $ANNOT) \
    | samtools view -f 16 -b - > $bamdir/${accession}_antisense_from_plus.bam

# 3. Extract reads antisense to genes on the - strand
bedtools intersect -s -wa -abam $BAM -b <(awk '$6 == "-"' $ANNOT) \
    | samtools view -f 0 -b - > $bamdir/${accession}_antisense_from_minus.bam

# 4. Merge and sort
samtools merge -f $bamdir/${accession}_antisense.bam \
    $bamdir/${accession}_antisense_from_plus.bam \
    $bamdir/${accession}_antisense_from_minus.bam

samtools sort -o $bamdir/${accession}_antisense.sorted.bam $bamdir/${accession}_antisense.bam
samtools index $bamdir/${accession}_antisense.sorted.bam


          # 5. Generate strand-specific bigWigs for visualization
          # bamCoverage -b $bamdir/${accession}_antisense.sorted.bam \
          #   -o $bwDir/${accession}_antisense_forward.bw \
          #   --filterRNAstrand forward \
          #   --normalizeUsing CPM \
          #   --binSize 10 \
          #   --numberOfProcessors 8
          #
          # bamCoverage -b $bamdir/${accession}_antisense.sorted.bam \
          #   -o $bwDir/${accession}_antisense_reverse.bw \
          #   --filterRNAstrand reverse \
          #   --normalizeUsing CPM \
          #   --binSize 10 \
          #   --numberOfProcessors 8

            # Get only antisense reads moving in the forward direction (FLAG 16 = reverse strand read)
  bamCoverage -b $bamdir/${accession}_antisense.sorted.bam \
    -o $bwDir/${accession}_antisense_forward.bw \
    --samFlagExclude 16 \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 8

  # Get only antisense reads moving in the reverse direction (FLAG 16 = reverse strand read)
  bamCoverage -b $bamdir/${accession}_antisense.sorted.bam \
    -o $bwDir/${accession}_antisense_reverse.bw \
    --samFlagInclude 16 \
    --normalizeUsing CPM \
    --binSize 10 \
    --numberOfProcessors 8


#quantify w/ featureCounts
  # featureCounts -T $THREADS \
  # -p \
  # -t CDS \
  # -g gene_name \
  # -s 0 --primary \
  # -a /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf \
  # -o ${fowcounts} \
  # ${bam}fwd.bam

#elseif read2 exists, do paired-end Trimming and PE mapping
elif [ -f $read2 ]; then
  echo "${line} running as PE file only"

  ##################
  #Trimming
  #################
  	  # module load Trim_Galore/0.6.7-GCCcore-11.2.0
      #
  	  # trim_galore --illumina --fastqc --paired --length 25 --basename ${accession} --gzip -o $trimmed $read1 $read2
      #
      # wait

#Map with STAR
  module load STAR/2.7.10b-GCC-11.3.0

#OG STAR mapping script
  # STAR --runMode alignReads \
  # --runThreadN $THREADS \
  # --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
  # --outFileNamePrefix ${bam} \
  # --readFilesIn $trimmed/${accession}_val_1.fq.gz $trimmed/${accession}_val_2.fq.gz \
  # --readFilesCommand zcat \
  # --outFilterType BySJout \
  # --alignIntronMax 10000 \
  # --outSAMtype BAM SortedByCoordinate \
  # --outBAMsortingBinsN 100 \
  # --outSAMunmapped Within \
  # --outSAMattributes Standard \
  # --limitBAMsortRAM 19990000000

  #Altered STAR script for antisense mapping w/ help from ChatGPT
    STAR --runMode alignReads \
      --runThreadN $THREADS \
      --genomeDir /home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR \
      --readFilesIn $trimmed/${accession}_val_1.fq.gz $trimmed/${accession}_val_2.fq.gz \
      --readFilesCommand zcat \
      --outFileNamePrefix ${bam} \
      --outSAMtype BAM SortedByCoordinate \
      --twopassMode Basic \
      --outFilterType BySJout \
      --outSAMstrandField intronMotif \
      --alignIntronMax 10000 \
      --outBAMsortingBinsN 100 \
      --outSAMunmapped Within \
      --outSAMattributes Standard \
      --limitBAMsortRAM 19990000000 \
      --quantMode TranscriptomeSAM GeneCounts

#ChatGPT helping me make bigwigs and count files with antisnese reads only

#load modules
ml BEDTools
module load SAMtools/1.16.1-GCC-11.3.0
module load StringTie
module load deepTools/3.5.2-foss-2022a

      # Set input and output paths
      BAM=${bam}Aligned.sortedByCoord.out.bam
      ANNOT=/home/evt82290/Research/annotation.bed    # BED12 gene annotation file
      GENOME="/home/ad45368/chrom_sizes.txt"   # Chromosome sizes

      # 1. Convert BAM to BED (splice-aware)
      bedtools bamtobed -split -i $BAM > $bedDir/${accession}_all_reads.bed

      # 2. Extract antisense reads (opposite strand of gene annotation)
      bedtools intersect -S -wa -a  $bedDir/${accession}_all_reads.bed -b $ANNOT >  $bedDir/${accession}_antisense_reads.bed

      # 3. Convert antisense BED to BAM
      bedtools bedtobam -i  $bedDir/${accession}_antisense_reads.bed -g $GENOME > $bamdir/${accession}_antisense.bam

      # 4. Sort and index the antisense BAM
      samtools sort -o $bamdir/${accession}_antisense.sorted.bam $bamdir/${accession}_antisense.bam
      samtools index $bamdir/${accession}_antisense.sorted.bam

      # 5. Generate strand-specific bigWigs for visualization
      bamCoverage -b $bamdir/${accession}_antisense.sorted.bam \
        -o $bwDir/${accession}_antisense_forward.bw \
        --filterRNAstrand forward \
        --normalizeUsing CPM \
        --binSize 10 \
        --numberOfProcessors 8

      bamCoverage -b $bamdir/${accession}_antisense.sorted.bam \
        -o $bwDir/${accession}_antisense_reverse.bw \
        --filterRNAstrand reverse \
        --normalizeUsing CPM \
        --binSize 10 \
        --numberOfProcessors 8

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
