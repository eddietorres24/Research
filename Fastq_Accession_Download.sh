#!/bin/bash
#SBATCH --job-name=fqdump_gzip
#SBATCH --partition=batch
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/evt82290/SRA/FastqFiles/logs/fqdump.%j.out
#SBATCH --error=/scratch/evt82290/SRA/FastqFiles/logs/fqdump.%j.err

# Reads SRR accessions (one per line) from JGI_RNA.txt in the submit directory,
# dumps paired FASTQs to OUTDIR, then gzips them.

cd "$SLURM_SUBMIT_DIR"

ACC_FILE="$SLURM_SUBMIT_DIR/JGI_RNA.txt"          # accession list file (one SRR per line)
OUTDIR="/scratch/evt82290/SRA/FastqFiles"         # destination for FASTQs
CACHEDIR="$OUTDIR/sra_cache"                      # local SRA cache for this job

mkdir -p "$OUTDIR" "$CACHEDIR" "$OUTDIR/logs"

# Load SRA Toolkit if your cluster uses modules (non-fatal if not present)
module load sratoolkit >/dev/null 2>&1 || true

echo "Starting: $(date)"
echo "Using accession list: $ACC_FILE"
echo "Output directory: $OUTDIR"
echo "Cache directory: $CACHEDIR"

# Iterate accessions
# - Skips blank lines and lines starting with '#'
while IFS=$'\r' read -r acc || [ -n "$acc" ]; do
  [[ -z "$acc" || "$acc" =~ ^# ]] && continue

  echo "==> Processing $acc"

  # 1) Prefetch to cache (writes .sra/.sralite or a directory under $CACHEDIR)
  if ! prefetch -O "$CACHEDIR" "$acc"; then
    echo "[WARN] prefetch failed for $acc — skipping"
    continue
  fi

  # Resolve the cached object path (covers .sralite or directory)
  obj=( "$CACHEDIR"/"$acc"* )
  if [[ ! -e "${obj[0]}" ]]; then
    echo "[WARN] No cached object found for $acc in
