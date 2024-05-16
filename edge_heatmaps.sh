#!/bin/bash
#SBATCH --job-name=ET_ChIP_Heatmap
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evt82290@uga.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=50gb
#SBATCH --time=6:00:00
#SBATCH --output=../HeatMap.%j.out
#SBATCH --error=../HeatMap.%j.err

computeMatrix reference-point -p 12 -R /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/Peaks/K9_Peaks.txt -S \
      /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/BigWigs/119-49_CUT_RUN_WT_H3K9__S37_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
      /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/BigWigs/133-75_CUTnRUN_WT_H3K27m2m3__S72_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
      /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/BigWigs/119-41_CUT_RUN_WT_K27me3_.bin_25.smooth_75Bulk.bw \
      /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/BigWigs/119-50_CUT_RUN_hda-1_H3K9__S38_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
      /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/BigWigs/6147_136-18_CUTnRUN_HDA-1_H3K27me2m3.bin_25.smooth_75Bulk.bw \
      /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/BigWigs/119-42_CUT_RUN_hda-1_K27me3__S30_L001_R1_001_val_1.fq.gz.bin_25.smooth_75Bulk.bw \
      --referencePoint TSS \
      -b 3000 -a 3000 \
      -o /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/MatrixFiles/Figure2_K9regions_FivePrime.matrix \
      --outFileNameMatrix Figure2_K9regions_FivePrime.matrix.txt \
      --sortRegions keep \
   --missingDataAsZero -bs 10
 # --sortUsing sum --sortUsingSamples 5 2
      #--sortRegions descend --sortUsing sum \
       #  --sortUsingSamples 3


#https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html
plotHeatmap -m /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/MatrixFiles/test.matrix -o /scratch/zlewis/HDA1_MS/Figure2/Figure2_Output/HeatmapFiles/Figure2_K9peaks_FivePrime.png --sortRegions keep --refPointLabel "5' edge" --samplesLabel "WT K9" "WT K27me2/3" "WT K27me3" "hda-1 K9" "hda-1 K27me2/3"  "hda-1 K27me3" --zMax 10 10 12 10 10 12 --colorMap 'Blues' 'Greens' 'Greens' 'Blues' 'Greens' 'Greens'
