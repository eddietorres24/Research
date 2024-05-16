#!/bin/bash
#SBATCH --job-name=ET_ChIP_Heatmap_edge
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=6:00:00
#SBATCH --output=../HeatMapEdge.%j.out
#SBATCH --error=../HeatMapEdge.%j.err

OUTDIR2="/scratch/evt82290/Run136/BigWigs"

module load deepTools/3.5.2-foss-2022a

computeMatrix reference-point --referencePoint TSS -p 12 -b 3000 -a 3000 \
 -S  $OUTDIR2/6147_136-1_ChIP_WT_H3K27me3_abcam_Rep2.bin_25.smooth_75Bulk.bw \
     $OUTDIR2/6147_136-2_ChIP_cac-1_H3K27me3_abcam_Rep2.bin_25.smooth_75Bulk.bw \
     $OUTDIR2/6147_136-3_ChIP_cac-2_H3K27me3_abcam_Rep2_S3_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
     $OUTDIR2/6147_136-4_ChIP_cac-3_H3K27me3_abcam_Rep2_S4_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
     $OUTDIR2/6147_136-83_ChIP_set-7_H3K27me3_CS_Rep2_S82_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
     -R Figure2_K27regions_Scaledcenter_FileToCheckOrderFINAL_ZL.txt \
     -o /Research/edge_matrix.matrix \
     --sortRegions keep \
     --missingDataAsZero \
     -bs 10






#https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html
plotHeatmap -m /Research/edge_matrix.matrix -out 2024_H3K27me3_CAF_edges.png \
  --sortRegions keep \
  --refPointLabel "5' edge" \
  --samplesLabel "WT K27me3" "cac-1 K27me3" "cac-2 K27me3" "cac-3 K27me3"  "set-7 K27me3" \
  --zMax 20 20 20 20 20 \
  --colorMap 'Blues'
