#Code to Run csaw ChIP-seq analysis
## Eddie Torres
### 7/3/24

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("csaw")

#Loading in BAM files
library(csaw)
param <- readParam(minq=20)
data <- windowCounts(bam.files, ext=110, width=10, param=param)
