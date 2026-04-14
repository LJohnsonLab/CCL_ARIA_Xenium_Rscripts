########################################################################
# Name: CCL_ST_v1_ARIA_03.R
# Project: mouse low dose ARIA (Akhil's project)
# Purpose: Analysis of low dose ARIA Xenium Spatial Transcriptomics data
#         - Part 1: Merging files, attaching metadata and pre-processing 
#         - Part 2: RCTD and SPLIT
#         - Part 3: manually annotating all clusters 
#         - Part 4: subclustering, cleaning and manually annotating microglia 
#         - Part 5: fixed slide and sample ID misalignment in Kai's microglia object, extracted final seurat objs metadata for squidpy and spatial data analyses, and code for publishable figures  
# Date created: 2/20/26
# Last updated: 4/14/26
# Author: Chloe Lucido
########################################################################

# Load Libraries ----
BiocManager::install("spacexr")
library(spacexr)
library(Seurat)
library(SeuratDisk)
library(future)
library(ggplot2)
library(arrow)
library(hdf5r)
library(presto)
library(glmGamPoi)
library(readr)
library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)
library(patchwork)
library(sceasy)
library(reticulate)
library(SPLIT)
library(qs2) 
library(RColorBrewer)
library(Polychrome)
library(purrr)
options(future.globals.maxSize = 100 *1024^3)

### 10. load xenium.aria obj with RCTD and SPLIT SHIFT info ----
xenium.aria <- LoadSeuratRds("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260209_aria_RCTDSPLIT.rds")


# CHECKING FOR BATCH EFFECTS ----
# used "Grin2a" as a marker to look for batch effects as this marker/neurons shouldnt change with treatment 
FeaturePlot(xenium.aria, split.by = "Slide", features = "Grin2a", reduction = "umap_SPLIT") 




### 11. finding clusters ----
# CHANGE DEFAULT ASSAY TO SCT_SPLIT

DefaultAssay(xenium.aria) <- "SCT_SPLIT"

# determining how many pcs to use
ElbowPlot(xenium.aria, ndims = 30, reduction = "pca_SPLIT")

xenium.aria <- FindNeighbors(
  xenium.aria,
  reduction = "pca_SPLIT",
  dims = 1:30
)

xenium.aria <- FindClusters(xenium.aria, resolution = 0.3, graph.name = "SCT_snn")


# find markers 
aria.markers <- FindAllMarkers(xenium.aria, group.by = "seurat_clusters")

aria.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_aria

write_csv(top5_aria, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260210_top5all_v1.csv")

table(xenium.aria$seurat_clusters)

#   0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    21 
#70324 46161 39421 34335 23500 22851 21175 15026 13050 11299  9878  7483  5927  5074  4155  4129  4091  3604  1994  1814  1013   624 



### 12. MANUAL ANNOTATION ----
FeaturePlot(xenium.aria, features = "Bgn", reduction = "umap_SPLIT") # VLMC marker

FeaturePlot(xenium.aria, features = "Slc47a1", reduction = "umap_SPLIT") # VLMC marker

FeaturePlot(xenium.aria, features = "Acta2", reduction = "umap_SPLIT") # VSMC marker

FeaturePlot(xenium.aria, features = "Tmem119", reduction = "umap_SPLIT") # microglia marker

FeaturePlot(xenium.aria, features = "Mrc1", reduction = "umap_SPLIT") # macrophage marker

FeaturePlot(xenium.aria, features = "Flt1", reduction = "umap_SPLIT") # endo marker 

FeaturePlot(xenium.aria, features = "Kcnj8", reduction = "umap_SPLIT") # pericyte marker 

FeaturePlot(xenium.aria, features = "Dnah11", reduction = "umap_SPLIT") # epend marker 

FeaturePlot(xenium.aria, features = "Car12", reduction = "umap_SPLIT") # CP marker 

# subclustering cluster 8 (v1) ----

Idents(xenium.aria) <- "seurat_clusters"

subclus_v1 <- FindSubCluster(xenium.aria, cluster = "8", resolution = 0.1, graph.name = "SCT_snn")

# find markers 
subclus_v1.markers <- FindAllMarkers(subclus_v1, group.by = "sub.cluster")

subclus_v1.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_subclus_v1

write_csv(top5_subclus_v1, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260210_top5subclus_v1.csv")

table(subclus_v1$sub.cluster)

#.  0     1    10    11    12    13    14    15    16    17    18    19     2    20    21     3     4     5     6     7   8_0   8_1   8_2 
#70324 46161  9878  7483  5927  5074  4155  4129  4091  3604  1994  1814 39421  1013   624 34335 23500 22851 21175 15026  6179  2536  1719 
# 8_3   8_4   8_5    9 
#1584   658   374 11299 




## subset immune clusters ----

Idents(subclus_v1) <- "sub.cluster"

immune_v1 <- subset(subclus_v1, idents = c("4", "8_3", "8_5", "9", "19", "21"))

# standard seurat pipeline
immune_v1 <- SCTransform(immune_v1, assay = "Xenium_SPLIT")

immune_v1 <- RunPCA(immune_v1, npcs = 30, features = rownames(immune_v1))

immune_v1 <- RunUMAP(immune_v1, dims = 1:30)


table(immune_v1$sub.cluster)

# 19    21     4   8_3   8_5     9 
#1814   624 23500  1584   374 11299 

# SANITY CHECK
dim(immune_v1)
# 474 39195


immune_v1 <- FindNeighbors(immune_v1, reduction = "pca", dims = 1:30)

immune_v1 <- FindClusters(immune_v1, resolution = 0.1, graph.name = "SCT_snn")

immune_v1.markers <- FindAllMarkers(immune_v1, group.by = "seurat_clusters")

immune_v1.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_immune_v1

write_csv(top5_immune_v1, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260210_top5immune_v1.csv") # version 1

# Immune obj: remove other cell-type gene markers ----
# lots of gene marker contamination, so removing these other gene markers 
genes_to_remove <- c("Ermn", "Cldn11", "Mog", # oligos 
                     "Ccdc153", "Dnah11", "Tmem212", # epend
                     "Pdgfra", "Opcml", "Tnr", # OPC
                     "Mgp", "Bgn", "Slc47a1", # VLMC
                     "Acta2", "Tagln", # VSMC
                     "Vtn", "Kcnj8", "Atp13a5", # pericytes
                     "Flt1", "Emcn", "Cldn5", # endo
                     "Mki67", "Car12", "Ttr", # CP
                     "Col1a2", "Lum", # fibroblast
                     "Dcx", "Pax6", # NPC
                     "Klrb1c", # NK cells
                     "Tpsab1", "Tpsb2", "Fcer1a", # mast cell
                     "Plac8", "Msr1", # monocyte
                     "Ly6g", "Camp", # neutrophil
                     "Snap25", "Grin2a", "Grin2b", # pan-neuronal
                     "Cux2", "Otof", "Stard8", "Lypd1", "Lrg1", "Adamts2", "Macc1", # L2/3 IT
                     "Fam84b", "Rapgef3", # L5 PT
                     "Rorb", "Tcap", "Hsd11b1", "Rspo1", "Whrn", # L5 IT
                     "Tunar", "Osr1", "Ppapdc1a", "Oprk1", # L6 IT
                     "Pou3f1", # L5 ET
                     "Tshz2", "Tox2", # L5/6 NP
                     "Syt6", "Trh", # L6 CT
                     "Lamp5", "Ndnf", # Lamp5
                     "Pvalb", "Nts", "Tac1", # pvalb
                     "Sncg", # sncg
                     "Sst", "Calb1", # sst
                     "Vip", "Calb2", "Pthlh", "Crh", # VIP
                     "Gad1", "Serpinf1", "Cnr1", "Lhx6", "Slc32a1", # inhibitory 
                     "Slc17a7", "Fezf2", "Htr2c", "Dpp4", "Slc17a6", # excitatory 
                     "Aldh1l1", "Aldoc", "Gfap", "Slc7a10", "Aqp4" # astrocytes
                     
                     
                     #"Cd3d", "Gzmb", "Pdcd1", # T cell
                     #"Cd209a", "Clec10a", # DC 
                     #"Cd247", # PBMC
                     #"Gja1", "Cpe", "Ubc", "Tubb2b", "Atp1b2", "Ptgds", "Ckb", "Psap", "Scg3", "Htra1", "Wif1" # LARA
)

genes_to_keep <- setdiff(rownames(immune_v1), genes_to_remove)


immune_v2 <- subset(immune_v1, features = genes_to_keep) # version 2

# SANITY CHECK
dim(immune_v2)
# 379 39195

# repeat pipeline with other marker genes subset out 

immune_v2 <- SCTransform(immune_v2, assay = "Xenium_SPLIT")

immune_v2 <- RunPCA(immune_v2, npcs = 30, features = rownames(immune_v2))

immune_v2 <- RunUMAP(immune_v2, dims = 1:30)

immune_v2 <- FindNeighbors(immune_v2, reduction = "pca", dims = 1:30)

immune_v2 <- FindClusters(immune_v2, resolution = 0.1, graph.name = "SCT_snn")

immune_v2.markers <- FindAllMarkers(immune_v2, group.by = "seurat_clusters")

immune_v2.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_immune_v2

write_csv(top5_immune_v2, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260210_top5immune_v2.csv") # version 2

table(immune_v2$seurat_clusters)

#    0     1     2     3     4     5     6 
# 11649 11558  8504  3062  2686  1189   547 


# subclustering clusters 3 and 4
Idents(immune_v2) <- "seurat_clusters"

immune_v3 <- FindSubCluster(immune_v2, cluster = "4", resolution = 0.1, graph.name = "SCT_snn")

dim(immune_v3)
#395 39195

# find markers 
immune_v3.markers <- FindAllMarkers(immune_v3, group.by = "sub.cluster")

immune_v3.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_immune_v3

write_csv(top5_immune_v3, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260211_top5immune_v3.csv")

Idents(immune_v3) <- "sub.cluster"

immune_v4 <- FindSubCluster(immune_v3, cluster = "4_0", resolution = 0.1, graph.name = "SCT_snn")

immune_v4.markers <- FindAllMarkers(immune_v4, group.by = "sub.cluster")

immune_v4.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_immune_v4

write_csv(top5_immune_v4, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260211_top5immune_v4.csv")

dim(immune_v4)
# 395 39195

# rename idents in immune obj

Idents(immune_v4) <- "sub.cluster"

# remove cluster 5 and 4_0_0 (reactive astros: cluster 19_0 and 8_3_0 in whole obj) 

immune_v5 <- subset(immune_v4, idents = c("4_0_0", "5"), invert = T)

# check idents before renaming them 

immune_v5 <- RenameIdents(immune_v5, 
                          "0" = "Microglia", # homeostatic micro
                          "1" = "Microglia", # LDAM
                          "2" = "Microglia", # DAM
                          "3" = "Microglia", # spp1 DAM
                          "4_0_1" = "BAM", # BAMs near astros 
                          "4_1" = "BAM", # BAM
                          "4_2" = "BAM", # AP micro
                          "6" = "T Cell"
) 

immune_v5$annotated_clusters <- Idents(immune_v5)

# original clusters are stored in metadata as sub.cluster

# RE-RUN SEURAT PIPELINE AFTER REMOVING THESE CLUSTERS 
immune_v5 <- SCTransform(immune_v5, assay = "Xenium_SPLIT")

immune_v5 <- RunPCA(immune_v5, npcs = 30, features = rownames(immune_v5))

immune_v5 <- RunUMAP(immune_v5, dims = 1:30)

immune_v5 <- FindNeighbors(immune_v5, reduction = "pca", dims = 1:30)



## BACK TO OVERALL UMAP CLUSTERING ----
# subcluster 19 ----
# after subclustering immune populations, we found reactive astros within the sub clusters from the og obj 8_3 and 19

subclus_v2 <- FindSubCluster(subclus_v1, cluster = "19", resolution = 0.1, graph.name = "SCT_snn")

dim(subclus_v2)
# 474 358399

subclus_v2.markers <- FindAllMarkers(subclus_v2, group.by = "sub.cluster")

subclus_v2.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_subclus_v2

write_csv(top5_subclus_v2, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260211_top5subclus_v2.csv")

# subcluster 8_3 ----
Idents(subclus_v2) <- "sub.cluster"

subclus_v3 <- FindSubCluster(subclus_v2, cluster = "8_3", resolution = 0.1, graph.name = "SCT_snn")

dim(subclus_v3)
# 474 358399

subclus_v3.markers <- FindAllMarkers(subclus_v3, group.by = "sub.cluster")

subclus_v3.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 5) |>
  ungroup() -> top5_subclus_v3

write_csv(top5_subclus_v3, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260211_top5subclus_v3.csv")


### 13. RENAMING CLUSTERS ----

Idents(subclus_v3) <- "sub.cluster"

levels(Idents(subclus_v3))
#[1] "5"     "6"     "1"     "12"    "0"     "14"    "8_2"   "21"    "8_0"   "4"     "9"     "7"     "10"    "2"     "19_1"  "3"     "16"   
#[18] "8_1"   "20"    "8_3_0" "8_3_1" "13"    "8_4"   "11"    "8_5"   "19_0"  "19_2"  "17"    "15"    "18"   

subclus_v3 <- RenameIdents(subclus_v3, 
                           "0" = "Oligodendrocyte", 
                           "1" = "Astrocyte", 
                           "2" = "Glutamatergic Neuron 1", 
                           "3" = "Glutamatergic Neuron 2", 
                           "4" = "Microglia", 
                           "5" = "Endothelial", 
                           "6" = "GABAergic Neuron 1", 
                           "7" = "GABAergic Neuron 2", 
                           "8_0" = "Fibroblast", 
                           "8_1" = "VSMC", 
                           "8_2" = "Endothelial", 
                           "8_3_0" = "Astrocyte", # reactive astro
                           "8_3_1" = "BAM", 
                           "8_4" = "VLMC", 
                           "8_5" = "BAM", 
                           "9" = "Microglia", 
                           "10" = "OPC", 
                           "11" = "Glutamatergic Neuron 3", 
                           "12" = "Pericyte", 
                           "13" = "GABAergic Neuron 3", 
                           "14" = "Glutamatergic Neuron 4", 
                           "15" = "CP", 
                           "16" = "GABAergic Neuron 4", 
                           "17" = "Ependymal", 
                           "18" = "Glutamatergic Neuron 5", 
                           "19_0" = "Astrocyte", # reactive astro
                           "19_1" = "T Cell", 
                           "19_2" = "Astrocyte", # reactive astro
                           "20" = "Endothelial", 
                           "21" = "BAM"
)


DimPlot(subclus_v3, reduction = "umap_SPLIT", cols = "polychrome")

# storing idents 
subclus_v3$annotatedclusters <- Idents(subclus_v3)

table(subclus_v3$annotatedclusters)

#       Oligodendrocyte              Astrocyte Glutamatergic Neuron 1 Glutamatergic Neuron 2              Microglia            Endothelial 
#70324                  48489                  39421                  34335                  34799                  25583 
#GABAergic Neuron 1     GABAergic Neuron 2             Fibroblast                   VSMC                    BAM                   VLMC 
#21175                  15026                   6179                   2536                   1477                    658 
#OPC Glutamatergic Neuron 3               Pericyte     GABAergic Neuron 3 Glutamatergic Neuron 4                     CP 
#9878                   7483                   5927                   5074                   4155                   4129 
#GABAergic Neuron 4              Ependymal Glutamatergic Neuron 5                 T Cell 
#4091                   3604                   1994                    591 




# Saving objects ----

### SAVE THE OG OBJ WITH ALL OF THE SPLIT METADATA BUT CLEAN THE IMMUNE OBJ OF THE METADATA COLS THAT WE DONT NEED 

# final subclus obj
SaveSeuratRds(subclus_v3, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260213_Part1_subclusv3.rds")

# final immune obj
SaveSeuratRds(immune_v5, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260213_Part1_immunev5.rds")
