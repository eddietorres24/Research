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

THREADS="${SLURM_CPUS_PER_TASK:-4}"

############################
# Reference / tools
############################
STAR_INDEX="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/STAR"
GTF="/home/zlewis/Genomes/Neurospora/Nc12_RefSeq/GCA_000182925.2_NC12_genomic_GFFtoGTFconversion.gtf"
CHROMSIZES="/home/ad45368/chrom_sizes.txt"

# Library strandedness: "reverse" (NEB dUTP) or "forward"
LIB_STRAND="reverse"

# Modules
module load Trim_Galore
module load STAR
module load SAMtools/1.16.1-GCC-11.3.0
module load Subread
module load BEDTools
module load deepTools/3.5.2-foss-2022a
# UCSC tools (use module if available or ensure on PATH)
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
# [STRICT FASTQ DISCOVERY + DEBUG]
############################
echo "[DEBUG] accession=$accession"
echo "[DEBUG] fastqPath=$fastqPath"

normalize_path_list() {
  # Accept "dir1 dir2", "dir1:dir2", or "dir1,dir2"
  echo "$1" | tr ',:' ' ' | xargs -n999 echo
}

basename_matches_accession () {
  # Accept only files whose basename STARTS with the accession (before first underscore)
  local f="$1" base root
  base="$(basename "$f")"
  root="${base%%_*}"
  [[ "$root" == "$accession" ]]
}

find_fastqs() {
  local paths=($(normalize_path_list "$fastqPath"))
  local d r1 r2 ru
  shopt -s nullglob
  for d in "${paths[@]}"; do
    echo "[DEBUG] scanning dir: $d"
    ls -1 "${d}/${accession}"* 2>/dev/null | head -n 5 | sed 's/^/[DEBUG] cand: /' || true

    # 1) SRA paired: accession_1.fastq.gz / accession_2.fastq.gz (exact-root)
    r1="${d}/${accession}_1.fastq.gz"
    r2="${d}/${accession}_2.fastq.gz"
    if [[ -f "$r1" && -f "$r2" ]] && basename_matches_accession "$r1" && basename_matches_accession "$r2"; then
      R1="$r1"; R2="$r2"; MODE="PE"
      echo "[INFO] Found PE (SRA)"; echo "[INFO] R1=$R1"; echo "[INFO] R2=$R2"
      shopt -u nullglob; return 0
    fi

    # 2) Illumina paired: accession*R1_001.fastq.gz / accession*R2_001.fastq.gz (exact-root)
    mapfile -t candR1 < <(ls -1 "${d}/${accession}"*R1*_001.fastq.gz 2>/dev/null || true)
    mapfile -t candR2 < <(ls -1 "${d}/${accession}"*R2*_001.fastq.gz 2>/dev/null || true)
    candR1=($(for f in "${candR1[@]:-}"; do basename_matches_accession "$f" && echo "$f"; done))
    candR2=($(for f in "${candR2[@]:-}"; do basename_matches_accession "$f" && echo "$f"; done))
    if [[ ${#candR1[@]} -gt 0 && ${#candR2[@]} -gt 0 ]]; then
      R1="${candR1[0]}"; R2="${candR2[0]}"; MODE="PE"
      echo "[INFO] Found PE (Illumina)"; echo "[INFO] R1=$R1"; echo "[INFO] R2=$R2"
      shopt -u nullglob; return 0
    fi

    # 3) Single-end exact: accession.fastq.gz
    ru="${d}/${accession}.fastq.gz"
    if [[ -f "$ru" ]] && basename_matches_accession "$ru"; then
      RU="$ru"; MODE="SE"
      echo "[INFO] Found SE"; echo "[INFO] RU=$RU"
      shopt -u nullglob; return 0
    fi

    # 4) SE via R1-only Illumina (exact-root)
    mapfile -t candR1only < <(ls -1 "${d}/${accession}"*R1*_001.fastq.gz 2>/dev/null || true)
    candR1only=($(for f in "${candR1only[@]:-}"; do basename_matches_accession "$f" && echo "$f"; done))
    if [[ ${#candR1only[@]} -gt 0 ]]; then
      R1="${candR1only[0]}"; MODE="SE"
      echo "[INFO] Found SE (R1-only)"; echo "[INFO] R1=$R1"
      shopt -u nullglob; return 0
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

# Safety: paired-end must not point to the same file
if [[ "$MODE" == "PE" ]]; then
  if [[ -z "${R1:-}" || -z "${R2:-}" ]]; then
    echo "[ERROR] Paired-end mode but R1 or R2 is empty." >&2; exit 1
  fi
  if [[ "$R1" == "$R2" ]]; then
    echo "[ERROR] R1 and R2 resolved to the same file: $R1" >&2; exit 1
  fi
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
  # Keep the transformed bedGraph for optional inspection
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
