#!/usr/bin/env bash
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

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

############################
# Required env from submitter
############################
: "${accession:?Missing accession (exported by submit script)}"
: "${fastqPath:?Missing fastqPath (exported by submit script)}"
: "${outdir:?Missing outdir (exported by submit script)}"

# Auto-detect available threads from SLURM, default to 4
THREADS="${SLURM_CPUS_PER_TASK:-4}"

############################
# Reference / tools (edit if needed)
############################
STAR_INDEX="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR"
GTF="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf"
CHROMSIZES="/home/ad45368/chrom_sizes.txt"

# Library strandedness:
#   "reverse" = NEB Ultra II Directional (dUTP) (most common)
#   "forward" for forward-stranded kits
LIB_STRAND="reverse"

# Modules
module load Trim_Galore
module load STAR
module load SAMtools
module load Subread              # featureCounts
module load BEDTools
module load deepToolsF
# UCSC tools if available (else leave as plain names if on $PATH)
BWTOBEDGRAPH=${BWTOBEDGRAPH:-bigWigToBedGraph}
BEDGRAPHTOBW=${BEDGRAPHTOBW:-bedGraphToBigWig}

############################
# Output dirs (your structure)
############################
TRIMDIR="${outdir}/TrimmedFastQs/${accession}"
BAMDIR="${outdir}/bamFiles/${accession}"
COUNTDIR="${outdir}/counts/${accession}"
BWDIR="${outdir}/bigWig/${accession}"
BEDGRAPHDIR="${outdir}/bedGraph/${accession}"
BEDDIR="${outdir}/beds/${accession}"
TMPDIR="${outdir}/tmp/${accession}"

mkdir -p "$TRIMDIR" "$BAMDIR" "$COUNTDIR" "$BWDIR" "$BEDGRAPHDIR" "$BEDDIR" "$TMPDIR"

############################
# FASTQ discovery
# - Supports multiple dirs via fastqPath as:
#     /dir1:/dir2:/dir3
#   or "/dir1 /dir2 /dir3"
#   or "/dir1,/dir2,/dir3"
############################
normalize_path_list() {
  # Convert commas and colons to spaces; squeeze repeats
  echo "$1" | tr ',:' ' ' | xargs -n999 echo
}

find_fastqs() {
  local paths=($(normalize_path_list "$fastqPath"))
  local d
  shopt -s nullglob
  for d in "${paths[@]}"; do
    # SRA paired: acc_1/acc_2
    if [[ -f "${d}/${accession}_1.fastq.gz" && -f "${d}/${accession}_2.fastq.gz" ]]; then
      R1="${d}/${accession}_1.fastq.gz"; R2="${d}/${accession}_2.fastq.gz"; MODE="PE"; echo "[INFO] Found PE (SRA) in $d"; return 0
    fi
    # Illumina paired: *R1_001/*R2_001
    mapfile -t candR1 < <(ls -1 "${d}/${accession}"*R1*_001.fastq.gz 2>/dev/null || true)
    mapfile -t candR2 < <(ls -1 "${d}/${accession}"*R2*_001.fastq.gz 2>/dev/null || true)
    if [[ ${#candR1[@]} -gt 0 && ${#candR2[@]} -gt 0 ]]; then
      R1="${candR1[0]}"; R2="${candR2[0]}"; MODE="PE"; echo "[INFO] Found PE (Illumina) in $d"; return 0
    fi
    # Single-end: acc.fastq.gz
    if [[ -f "${d}/${accession}.fastq.gz" ]]; then
      RU="${d}/${accession}.fastq.gz"; MODE="SE"; echo "[INFO] Found SE in $d"; return 0
    fi
    # SE using R1-only
    mapfile -t candR1only < <(ls -1 "${d}/${accession}"*R1*_001.fastq.gz 2>/dev/null || true)
    if [[ ${#candR1only[@]} -gt 0 ]]; then
      R1="${candR1only[0]}"; MODE="SE"; echo "[INFO] Found SE (R1-only) in $d"; return 0
    fi
  done
  shopt -u nullglob
  return 1
}

R1=""; R2=""; RU=""; MODE=""
if ! find_fastqs; then
  echo "[ERROR] Could not locate FASTQs for ${accession} under: ${fastqPath}" >&2
  exit 1
fi

############################
# Trimming
############################
if [[ "$MODE" == "PE" ]]; then
  echo "[INFO] Trimming PE"
  trim_galore --illumina --paired --length 25 --basename "${accession}" --gzip \
              -o "$TRIMDIR" "$R1" "$R2"
  R1T="${TRIMDIR}/${accession}_val_1.fq.gz"
  R2T="${TRIMDIR}/${accession}_val_2.fq.gz"
else
  echo "[INFO] Trimming SE"
  SE_IN="${RU:-$R1}"
  trim_galore --illumina --length 25 --basename "${accession}" --gzip \
              -o "$TRIMDIR" "$SE_IN"
  RUT="${TRIMDIR}/${accession}_trimmed.fq.gz"
fi

############################
# STAR alignment
############################
OUTPFX="${BAMDIR}/${accession}_"
STAR_COMMON=(
  --runThreadN "$THREADS"
  --genomeDir "$STAR_INDEX"
  --outFileNamePrefix "$OUTPFX"
  --readFilesCommand zcat
  --outSAMtype BAM SortedByCoordinate
  --twopassMode Basic
  --outFilterType BySJout
  --alignIntronMax 10000
  --outSAMunmapped Within
  --outSAMattributes Standard
  --outSAMstrandField intronMotif
  --quantMode GeneCounts
)
if [[ "$MODE" == "PE" ]]; then
  STAR "${STAR_COMMON[@]}" --readFilesIn "$R1T" "$R2T"
else
  STAR "${STAR_COMMON[@]}" --readFilesIn "$RUT"
fi
BAM="${OUTPFX}Aligned.sortedByCoord.out.bam"
samtools index "$BAM"

############################
# Stranded bigWigs (CPM + log2(CPM+1))
############################
if [[ "$LIB_STRAND" == "reverse" ]]; then
  FSTRAND="reverse"; RSTRAND="forward"
else
  FSTRAND="forward"; RSTRAND="reverse"
fi

BIN=10
# Linear CPM
bamCoverage -b "$BAM" -o "${BWDIR}/${accession}.${FSTRAND}.cpm.bw" \
  --normalizeUsing CPM --binSize $BIN --filterRNAstrand "$FSTRAND" -p "$THREADS" --skipNonCoveredRegions
bamCoverage -b "$BAM" -o "${BWDIR}/${accession}.${RSTRAND}.cpm.bw" \
  --normalizeUsing CPM --binSize $BIN --filterRNAstrand "$RSTRAND" -p "$THREADS" --skipNonCoveredRegions

# Log2(CPM+1) via bedGraph roundtrip
log_transform_bw () {
  local inbw="$1" outbw="$2" base tmpbg tmpbg2
  base="$(basename "$inbw" .bw)"
  tmpbg="${TMPDIR}/${base}.bedGraph"
  tmpbg2="${TMPDIR}/${base}.log2p1.bedGraph"
  $BWTOBEDGRAPH "$inbw" "$tmpbg"
  awk 'BEGIN{OFS="\t"}
       NF>=4{v=$4+0; if(v<0)v=0; print $1,$2,$3,(log(v+1)/log(2)); next}
       {print}' "$tmpbg" > "$tmpbg2"
  sort -k1,1 -k2,2n "$tmpbg2" > "${tmpbg2}.sorted"
  $BEDGRAPHTOBW "${tmpbg2}.sorted" "$CHROMSIZES" "$outbw"
  # optional keep bedGraphs:
  cp "$tmpbg2" "${BEDGRAPHDIR}/${base}.log2p1.bedGraph" || true
  rm -f "$tmpbg" "$tmpbg2" "${tmpbg2}.sorted"
}
log_transform_bw "${BWDIR}/${accession}.${FSTRAND}.cpm.bw" "${BWDIR}/${accession}.${FSTRAND}.log2p1.bw"
log_transform_bw "${BWDIR}/${accession}.${RSTRAND}.cpm.bw" "${BWDIR}/${accession}.${RSTRAND}.log2p1.bw"

############################
# Gene sense vs antisense counts
############################
FLIPGTF="${BAMDIR}/$(basename "$GTF" .gtf).flipped.gtf"
if [[ ! -s "$FLIPGTF" ]]; then
  awk 'BEGIN{OFS="\t"} $3=="exon"{if($7=="+")$7="-"; else if($7=="-")$7="+"} {print}' \
    "$GTF" > "$FLIPGTF"
fi

FC_STRAND_ARG="-s 2"; [[ "$LIB_STRAND" == "forward" ]] && FC_STRAND_ARG="-s 1"

featureCounts -T "$THREADS" $FC_STRAND_ARG -t exon -g gene_id \
  -a "$GTF"     -o "${COUNTDIR}/${accession}.sense.txt"     "$BAM"
featureCounts -T "$THREADS" $FC_STRAND_ARG -t exon -g gene_id \
  -a "$FLIPGTF" -o "${COUNTDIR}/${accession}.antisense.txt" "$BAM"

echo "[DONE] ${accession}"
