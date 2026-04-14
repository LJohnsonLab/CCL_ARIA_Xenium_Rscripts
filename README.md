**Date created:** 3/11/26    
**Last Updated:** 4/14/26    
**Author:** Chloe Lucido (parts of code adapted from Jose Arbones-Mainar)    
**Project:** mouse low dose ARIA (Akhil's project)  
**Purpose:** Analysis of low dose ARIA Xenium Spatial Transcriptomics data  
**- Part 1:** Merging files, attaching metadata and pre-processing   
**- Part 2:** RCTD and SPLIT  
**- Part 3:** manually annotating all clusters   
**- Part 4:** subclustering, cleaning and manually annotating microglia   
**- Part 5:** fixed slide and sample ID misalignment in Kai's microglia object, extracted final seurat objs metadata for squidpy and spatial data analyses, and code for publishable figures  

**Important Notes:**    
- Kai subclustered the microglia, cleaned and annotated the subclusters so the object used to make figures and for subsequent analysis uses the microglia object she generated with updated slide and sample ID alignment 

**R session and packages Info**    
R version 4.4.3 (2025-02-28)     

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] cowplot_1.2.0         pheatmap_1.0.13       scCustomize_3.2.4     dittoSeq_1.18.0       SCP_0.5.6             speckle_1.6.0         Polychrome_1.5.4     
 [8] RColorBrewer_1.1-3    qs2_0.1.7             SPLIT_0.1.3           sceasy_0.0.7          reticulate_1.45.0     patchwork_1.3.2       readxl_1.4.5         
[15] lubridate_1.9.5       forcats_1.0.1         stringr_1.6.0         purrr_1.2.1           tidyr_1.3.2           tibble_3.3.1          tidyverse_2.0.0      
[22] dplyr_1.2.0           readr_2.1.6           glmGamPoi_1.18.0      presto_1.0.0          data.table_1.18.2.1   Rcpp_1.1.1            hdf5r_1.3.12         
[29] arrow_23.0.0.1        ggplot2_4.0.2         future_1.69.0         SeuratDisk_0.0.0.9021 spacexr_2.2.1         Seurat_5.4.0          SeuratObject_5.3.0   
[36] sp_2.2-1 
