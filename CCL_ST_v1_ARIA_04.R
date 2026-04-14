########################################################################
# Name: CCL_ST_v1_ARIA_04.R
# Project: low dose ARIA (Akhil's project)
# Purpose: Analysis of low dose ARIA Xenium Spatial Transcriptomics data
#         - Part 1: Merging files, attaching metadata and pre-processing 
#         - Part 2: RCTD and SPLIT
#         - Part 3: manually annotating all clusters 
#         - Part 4: subclustering, cleaning and manually annotating microglia 
#         - Part 5: code for publishable figures  
# Input Files: 
# Output Files: 
# Date created : 2/20/26
# Last updated: 3/23/26
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
library(stringr)
options(future.globals.maxSize = 100 *1024^3)


### OLD IMMUNE AND MICROGLIA MANUAL ANNOTATION ----

### 01. Load objs from part 1 ----
ann_aria <- LoadSeuratRds("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260213_Part1_subclusv3.rds")
clean_immune <- LoadSeuratRds("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260213_Part1_immunev5.rds")


dim(ann_aria)
# 474 358399

dim(clean_immune)
# 395 36932



### 02. subset and re-cluster microglia ----

Idents(clean_immune) <- "seurat_clusters"

sub_immune <- subset(clean_immune, idents = c("0", "1", "2", "3"))

sub_immune <- SCTransform(sub_immune, assay = "Xenium_SPLIT")

sub_immune <- RunPCA(sub_immune, npcs = 30, features = rownames(sub_immune))

sub_immune <- RunUMAP(sub_immune, dims = 1:30)

sub_immune <- FindNeighbors(sub_immune, reduction = "pca", dims = 1:30)

sub_immune <- FindClusters(sub_immune, resolution = 0.2, graph.name = "SCT_snn")

sub_immune.markers <- FindAllMarkers(sub_immune, group.by = "seurat_clusters", only.pos = T, logfc.threshold = 0.25)

sub_immune.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  ungroup() -> top10_subimmune

write_csv(top10_subimmune, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260216_top10_subimmune.csv") # version 1


#### USING SUB_IMMUNE OBJ AS MICROGLIA ONLY OBJ

sub_immune$seurat_clusters <- Idents(sub_immune)

# cluster marker file is "20260216_top5_subimmune"
sub_immune <- RenameIdents(sub_immune, 
                           "0" = "AAM", # astro associated microglia
                           "1" = "DAM 1", 
                           "2" = "HM", 
                           "3" = "ARM", 
                           "4" = "LAM", # lipid-associated micros
                           "5" = "spp1 DAM", 
                           "6" = "APM", 
                           "7" = "IRM", 
                           "8" = "IFN"
)

sub_immune$annotated_subclusters <- Idents(sub_immune)


### 03. Mapping microglia subclusters to immune obj ----
# after annotating microglia subclusters we want to map them to immune obj, which already has all of the microglia cells 

# making sure all cell names match
all(colnames(sub_immune) %in% colnames(clean_immune)) # should output TRUE if they all match 

# extracting microglia subcluster annotations
Idents(sub_immune) <- sub_immune$annotated_subclusters

microglia.subclusters <- as.character(Idents(sub_immune))

common_cells <- intersect(colnames(sub_immune), colnames(clean_immune))

clean_immune$microglia_subclusters <- NA

# match common cells btwn microglia and immune objs
clean_immune$microglia_subclusters[common_cells] <-
  as.character(Idents(sub_immune))[match(
    common_cells,
    colnames(sub_immune)
  )]



DefaultAssay(clean_immune) <- "SCT"

# creating metadata column with all immune clusters and microglial subclusters
clean_immune$annotated_subclusters <- as.character(clean_immune$annotated_clusters)

micro_idx <- !is.na(clean_immune$microglia_subclusters)

clean_immune$annotated_subclusters[micro_idx] <- 
  as.character(clean_immune$microglia_subclusters[micro_idx])

# convert back to factor 
clean_immune$annotated_subclusters <- factor(clean_immune$annotated_subclusters)

# dim plot to ensure it worked 
DimPlot(clean_immune, group.by = "annotated_subclusters")

DimPlot(clean_immune, group.by = "annotated_clusters")


# save as qs
qs_save(sub_immune, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260217_Part2_microglia.qs2") # with updated annotated subcluster name

# save as qs 
qs_save(clean_immune, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260217_Part2_immune.qs2")



# getting rid of all unnecessary metadata cols from full obj, UPDATE COLS AS NEEDED ----
# RE-ALIGN METADATA BEFORE REMOVING COLUMNS!!

ann_aria@meta.data <- ann_aria@meta.data[, !colnames(ann_aria@meta.data) %in% c("Notes...10", 
                                                                                         "Notes...12", 
                                                                                         "min_score",
                                                                                         "conv_all",
                                                                                         "conv_doublet",
                                                                                         "max_doublet_weight",                          
                                                                                         "n_candidates",                              
                                                                                         "rctd_weights_entropy",                        
                                                                                         "weight_first_type",                          
                                                                                         "weight_second_type", 
                                                                                         "annot_min_singlet_score",                     
                                                                                         "annot_max_weight",                           
                                                                                         "annot_max_doublet_weight",                    
                                                                                         "w1_larger_w2", 
                                                                                         "total_neighbors_N",                          
                                                                                         "annotated_neighbors_N",                       
                                                                                         "total_singlets_neighbors_N",                 
                                                                                         "second_type_neighbors_N",                     
                                                                                         "second_type_neighbors_no_reject_N",          
                                                                                         "second_type_singlets_neighbors_N",            
                                                                                         "second_type_class_neighbors_N",             
                                                                                         "first_type_neighbors_N",                       
                                                                                         "first_type_singlets_neighbors_N",            
                                                                                         "first_type_class_neighbors_N",                
                                                                                         "same_second_type_neighbors_N",               
                                                                                         "neighborhood_weights_first_type",             
                                                                                         "neighborhood_weights_second_type",           
                                                                                         "sum_w1_second_type_in_neighborhood",          
                                                                                         "sum_w2_second_type_in_neighborhood",         
                                                                                         "max_weight_of_spilling_type_in_neighborhood", 
                                                                                         "sum_nCount_neighborhood",                   
                                                                                         "sum_nCount_neighborhood_spilling_type",       
                                                                                         "first_type_neighborhood",                    
                                                                                         "first_type_neighborhood_agreement",           
                                                                                         "first_type_class_neighborhood",              
                                                                                         "first_type_class_neighborhood_agreement",     
                                                                                         "second_type_neighborhood",                   
                                                                                         "second_type_neighborhood_agreement",          
                                                                                         "second_type_class_neighborhood",             
                                                                                         "second_type_class_neighborhood_agreement",    
                                                                                         "first_type_neighborhood_certainty",          
                                                                                         "first_type_class_neighborhood_certainty",     
                                                                                         "second_type_neighborhood_certainty",         
                                                                                         "second_type_class_neighborhood_certainty"
)]

# FIX TO RE-ALIGN CELLS AFTER TAKING OUT METADATA COLUMNS

meta_fix <- ann_aria@meta.data[, c("Slide", "sample_ID")]

# align explicitly
meta_fix <- meta_fix[colnames(ann_aria), ]

ann_aria$Slide <- meta_fix$Slide
ann_aria$sample_ID <- meta_fix$sample_ID


#qs_save(ann_aria, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260217_Part2_fullobj.qs2") OLD OBJ W/ MISALIGNED SLIDE AND SAMPLE_ID METADATA COLS

# save new obj with aligned slide and sample ID columns
qs_save(ann_aria, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260320_fullobj.qs2")
















#### 02/19/26: IMMUNE CLUSTER ANNOTATION SESSION W/ KAI ----

### 01. Load objs ----

full.obj <- qs_read("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260217_Part2_fullobj.qs2")


# SANITY CHECK
dim(full.obj) # Output: 474 358399


### 02. Microglia Subcluster ----

microglia <- subset(full.obj, idents = c("Microglia"))

microglia <- SCTransform(microglia, assay = "Xenium_SPLIT") # Xenium_SPLIT assay contains SPLIT purified counts 

microglia <- RunPCA(microglia, npcs = 30, features = rownames(microglia))

microglia <- RunUMAP(microglia, dims = 1:30)

microglia <- FindNeighbors(microglia, reduction = "pca", dims = 1:30)

microglia <- FindClusters(microglia, resolution = 0.2, graph.name = "SCT_snn")

DimPlot(microglia, group.by = "seurat_clusters", label = T)

microglia.markers <- FindAllMarkers(microglia, group.by = "seurat_clusters", only.pos = T, logfc.threshold = 0.25)

microglia.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  ungroup() -> top10_microglia

write_csv(top10_microglia, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260219_top10_microglia_res03.csv") # version 1

dim(microglia) # output: 474 34799

# dirty clusters 0, 4 & 6

FeaturePlot(microglia, features = c("Tmem119", "P2ry12"))

DotPlot(microglia, features = c("Tmem119", "P2ry12")) # for sure removing 0 and 4, unsure abt 6


microglia_v1 <- subset(microglia, idents = c("0", "4", "6"), invert = T)

microglia_v1 <- SCTransform(microglia_v1, assay = "Xenium_SPLIT") # Xenium_SPLIT assay contains SPLIT purified counts 

microglia_v1 <- RunPCA(microglia_v1, npcs = 30, features = rownames(microglia_v1))

microglia_v1 <- RunUMAP(microglia_v1, dims = 1:30)

microglia_v1 <- FindNeighbors(microglia_v1, reduction = "pca", dims = 1:30)

microglia_v1<- FindClusters(microglia_v1, resolution = 0.2, graph.name = "SCT_snn")

DimPlot(microglia_v1, group.by = "seurat_clusters", label = T)

microglia_v1.markers <- FindAllMarkers(microglia_v1, group.by = "seurat_clusters", only.pos = T, logfc.threshold = 0.25)

microglia_v1.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  ungroup() -> top10_microglia_v1

write_csv(top10_microglia_v1, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260219_top10_microglia_res02_clus6removed.csv") # cluster 6 removed this time

> dim(microglia_v1)
[1]   474 23003





### 04. REMOVING GENES ----


genes_to_remove <- c("Ermn", "Cldn11", "Mog", "Mag", "Pllp", # oligos/oligo-related 
                     "Ccdc153", "Dnah11", "Tmem212", # epend
                     "Pdgfra", "Tnr", # OPC
                     "Mgp", "Bgn", "Slc47a1", # VLMC
                     "Acta2", "Tagln", # VSMC
                     "Vtn", "Kcnj8", "Atp13a5", # pericytes
                     "Flt1", "Emcn", "Cldn5", # endo
                     "Car12", "Ttr", # CP
                     "Col1a2", "Lum", # fibroblast
                     "Dcx", "Pax6", # NPC
                     "Hbb-bs", #RBC-related
                     "Snap25", "Grin2a", "Grin2b", # pan-neuronal
                     "Cux2", "Otof", "Stard8", "Lypd1", "Lrg1", "Adamts2", "Macc1", # L2/3 IT
                     "Rapgef3", # L5 PT
                     "Rorb", "Tcap", "Hsd11b1", "Rspo1", "Whrn", # L5 IT
                     "Tunar", "Osr1", "Oprk1", # L6 IT
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
                     "Aldh1l1", "Aldoc", "Gfap", "Slc7a10", "Aqp4", # astrocytes
                     "Cd247", # PBMC
                     "Gja1", "Cpe", "Ubc", "Tubb2b", "Atp1b2", "Ptgds", "Ckb", "Psap", "Scg3", "Htra1", "Wif1",  # LARA
                     "Serpina3n", "Tjp1", # IMMUNOLOGY/COMPLEMENT/BBB
                     "Mt1", "Mt2", "Mt3", # Lance's list
                     "Opcml" # neuron-related 
                     
                     #"Plac8", "Msr1", # monocyte
                     #"Klrb1c", # NK cells
                     #"Tpsab1", "Tpsb2", "Fcer1a", # mast cell
                     #"Ly6g", "Camp", # neutrophil
                     #"Cd3d", "Gzmb", "Pdcd1", # T cell
                     #"Cd209a", "Clec10a" # DC
                     #"Mki67" # proliferation marker
                      
)

genes_to_keep <- setdiff(rownames(microglia_v1), genes_to_remove)


microglia_v2 <- subset(microglia_v1, features = genes_to_keep) # version 2


# SANITY CHECK
dim(microglia_v2) # output: 377 23003

microglia_v2 <- SCTransform(microglia_v2, assay = "Xenium_SPLIT") # Xenium_SPLIT assay contains SPLIT purified counts 

microglia_v2 <- RunPCA(microglia_v2, npcs = 25, features = rownames(microglia_v2))

microglia_v2<- RunUMAP(microglia_v2, dims = 1:25)

microglia_v2<- FindNeighbors(microglia_v2, reduction = "pca", dims = 1:25)

microglia_v2<- FindClusters(microglia_v2, resolution = 0.2, graph.name = "SCT_snn")

DimPlot(microglia_v2, group.by = "seurat_clusters", label = T)

microglia_v2.markers <- FindAllMarkers(microglia_v2, group.by = "seurat_clusters", only.pos = T, logfc.threshold = 0.25)

microglia_v2.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  ungroup() -> top10_microglia_v2

write_csv(top10_microglia_v2, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260223_top10_microglia_v2_euclidean.csv")



# save microglia obj ----
qs_save(microglia_v2, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260223_microglia_25pcs_res02.qs2") # final microglia obj




# 2/23/26 Immune subcluster from full obj
### 01. Load objs ----

full.obj <- qs_read("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260217_Part2_fullobj.qs2")


# SANITY CHECK
dim(full.obj) # Output: 474 358399


### 02. failed immune ----

immune <- subset(full.obj, idents = c("BAM", "T Cell", "Microglia"))

dim(immune) # output: 474 36867

immune <- SCTransform(immune, assay = "Xenium_SPLIT") # Xenium_SPLIT assay contains SPLIT purified counts 

immune <- RunPCA(immune, npcs = 30, features = rownames(immune))

immune<- RunUMAP(immune, dims = 1:30)

immune<- FindNeighbors(immune, reduction = "pca", dims = 1:30)

immune<- FindClusters(immune, resolution = 0.3, graph.name = "SCT_snn")

DimPlot(immune, group.by = "seurat_clusters", label = T)

immune.markers <- FindAllMarkers(immune, group.by = "seurat_clusters", only.pos = T, logfc.threshold = 0.25)

immune.markers |>
  group_by(cluster) |>
  dplyr::filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC)) |>
  slice_head(n = 10) |>
  ungroup() -> top10_immune


write_csv(top10_immune, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/cluster_markers/20260223_top10_immune.csv") # 20 pcs





genes_to_remove <- c("Ermn", "Cldn11", "Mog", "Mag", "Pllp", # oligos/oligo-related 
                     "Ccdc153", "Dnah11", "Tmem212", # epend
                     "Pdgfra", "Tnr", # OPC
                     "Mgp", "Bgn", "Slc47a1", # VLMC
                     "Acta2", "Tagln", # VSMC
                     "Vtn", "Kcnj8", "Atp13a5", # pericytes
                     "Flt1", "Emcn", "Cldn5", # endo
                     "Car12", "Ttr", # CP
                     "Col1a2", "Lum", # fibroblast
                     "Dcx", "Pax6", # NPC
                     "Hbb-bs", #RBC-related
                     "Snap25", "Grin2a", "Grin2b", # pan-neuronal
                     "Cux2", "Otof", "Stard8", "Lypd1", "Lrg1", "Adamts2", "Macc1", # L2/3 IT
                     "Rapgef3", # L5 PT
                     "Rorb", "Tcap", "Hsd11b1", "Rspo1", "Whrn", # L5 IT
                     "Tunar", "Osr1", "Oprk1", # L6 IT
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
                     "Aldh1l1", "Aldoc", "Gfap", "Slc7a10", "Aqp4", # astrocytes
                     "Cd247", # PBMC
                     "Gja1", "Cpe", "Ubc", "Tubb2b", "Atp1b2", "Ptgds", "Ckb", "Psap", "Scg3", "Htra1", "Wif1",  # LARA
                     "Serpina3n", "Tjp1", # IMMUNOLOGY/COMPLEMENT/BBB
                     "Mt1", "Mt2", "Mt3", # Lance's list
                     "Opcml" # neuron-related 
                     
                     #"Plac8", "Msr1", # monocyte
                     #"Klrb1c", # NK cells
                     #"Tpsab1", "Tpsb2", "Fcer1a", # mast cell
                     #"Ly6g", "Camp", # neutrophil
                     #"Cd3d", "Gzmb", "Pdcd1", # T cell
                     #"Cd209a", "Clec10a" # DC
                     #"Mki67" # proliferation marker
                     
)

genes_to_keep <- setdiff(rownames(BAM_Tcells), genes_to_remove)


BAM_Tcells <- subset(BAM_Tcells, features = genes_to_keep) # version 2


BAM_Tcells <- UpdateSeuratObject(BAM_Tcells)

DefaultAssay(BAM_Tcells) <- "Xenium_SPLIT"
DefaultAssay(microglia_v2) <- "Xenium_SPLIT"

BAM_Tcells <- JoinLayers(BAM_Tcells)

immune <- merge(microglia_v2, BAM_Tcells)





