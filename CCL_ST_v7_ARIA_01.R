########################################################################
# Name: CCL_ST_v7_ARIA_01.R
# Project: low dose ARIA (Akhil's project)
# Purpose: Analysis of low dose ARIA Xenium Spatial Transcriptomics data
#         - Part 1: Merging files, attaching metadata and pre-processing 
#         - Part 2: RCTD and SPLIT
#         - Part 3: manually annotating all clusters 
#         - Part 4: subclustering, cleaning and manually annotating microglia 
#         - Part 5: code for publishable figures 
# Input Files: Xenium slide output folders 
# Output Files: 20260211_Part1_subclusv3.rds == final subclustered and fully annotated seurat obj
#             - 20260211_Part1_immunev4.rds == final subclustered (unannotated) immune seurat obj
# UPDATE FOR SOURCE CODE:
#             - change all "SaveSeuratRds" to "qsave" (saves files as qs rather than rds)
# Date created: 2/11/26
# Last updated: 2/20/27
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


# Vignette: https://satijalab.org/seurat/articles/seurat5_spatial_vignette_2 


### 01. Loading in Xenium data ----
# path to xenium output folder 

#slide 1
path1 <- "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/Run2_20250606__202944__20250606_AgingXMetabolism_2/Slide1_output-XETG00118__0021991__Region_1__20250606__202953"
# slide 2
path2 <- "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/Run2_20250606__202944__20250606_AgingXMetabolism_2/Slide2_output-XETG00118__0021998__Region_1__20250606__202953"

# load xenium data
xenium.obj1 <- LoadXenium(path1, fov = "fov", segmentations = "cell") 

xenium.obj2 <- LoadXenium(path2, fov = "fov", segmentations = "cell")

# adding slide number to cell names 

xenium.obj1 <- RenameCells(xenium.obj1, add.cell.id = "slide1")

xenium.obj2 <- RenameCells(xenium.obj2, add.cell.id = "slide2")

# Notes:
  # fov = field of view (spatial info loaded into slots of seurat obj), all data is loaded into fov, later you can crop the fov to the region of interest

xenium.obj1$Slide <- "Slide1"
xenium.obj2$Slide <- "Slide2"


### 02. Merging Xenium Files and saving merged file ----

merged.obj <- merge(xenium.obj1, xenium.obj2)

merged.obj <- JoinLayers(merged.obj, add.prefix = FALSE)

# SANITY CHECK
dim(merged.obj)

# saving merged object 

SaveSeuratRds(merged.obj, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Files/20251124_mergedobj_noQC.rds")

### 03. Attaching metadata (Xenium Explorer annotations and csv file) ----

# load merged object

merged.obj <- LoadSeuratRds("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Files/20251124_mergedobj_noQC.rds")

## 03a. adding sample ID XE annotations ----

files <- list(
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_ExplorerAnnotations/slide1_KK4_465_cells_stats.csv",
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_ExplorerAnnotations/slide1_KK4_492_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_ExplorerAnnotations/slide1_KK4_504_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_ExplorerAnnotations/slide2_KK4_464_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_ExplorerAnnotations/slide2_KK4_496_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_ExplorerAnnotations/slide2_KK4_502_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_ExplorerAnnotations/slide2_M1_36wks_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_ExplorerAnnotations/slide2_F1_16wks_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_ExplorerAnnotations/slide2_M2_16wks_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_ExplorerAnnotations/slide1_M2_36wks_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_ExplorerAnnotations/slide1_M1_16wks_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_ExplorerAnnotations/slide1_F1_36wks_cells_stats.csv")


# reading and tagging each csv cell stats file with corresponding sample name 

annotation_list <- map(files, function(x) { 
  
  # reading each csv file, skipping first 2 metadata lines
  
  read.csv(x, skip = 2, blank.lines.skip = TRUE)
  
})


# defining which samples are in each xenium object for proper suffix placement ("_1" or "_2")

prefix_map <- c(
  "KK4_465" = "slide1_", 
  "KK4_492" = "slide1_", 
  "KK4_504" = "slide1_", 
  "KK4_464" = "slide2_", 
  "KK4_496" = "slide2_", 
  "KK4_502" = "slide2_", 
  "M1_36wks" = "slide2_", 
  "F1_16wks" = "slide2_", 
  "M2_16wks" = "slide2_", 
  "M2_36wks" = "slide1_", 
  "M1_16wks" = "slide1_", 
  "F1_36wks" = "slide1_"
) 

# changing names of each element in annotation_list to its corresponding sample id 

names(annotation_list) <- names(prefix_map)


# iterating through annotation_list and adding the proper suffix to the cell ids

annotation_list <- imap(annotation_list, function(df, name) {
  
  df <- annotation_list[[name]] # extracting each dataframe from annotation_list
  
  # appending suffix to cell ids
  
  df <- df %>% 
    mutate(Cell.ID = paste0(prefix_map[[name]], Cell.ID), 
           
           # adding sample id column
           Sample.ID = name)  
  
  return(df)
})



# combine all csv files in list to one dataframe 

annotation_df <- bind_rows(annotation_list)

# set rownames to cell ID so that they match seurat obj cell names 

rownames(annotation_df) <- annotation_df$Cell.ID

# reordering rows to match seurat object

annotation_df <- annotation_df[Cells(merged.obj), , drop = F]

# adding cell ID and sample ID dataframe (annotation_df) to merged obj

merged.obj <- AddMetaData(merged.obj, annotation_df$Sample.ID, col.name = "sample_ID")



## 03b. adding regional XE annotations ----

files <- list(
  # KK4_464 regional annotations
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_464_regions/Slide0021998_KK4_464_Cortex_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_464_regions/Slide0021998_KK4_464_Hippocampus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_464_regions/Slide0021998_KK4_464_Hypothalamus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_464_regions/Slide0021998_KK4_464_Thalamus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_464_regions/Slide0021998_KK4_464_WM_cells_stats.csv",
  # KK4_465 regional annotations
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_465_regions/Slide0021991_KK4_465_Cortex_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_465_regions/Slide0021991_KK4_465_Hippocampus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_465_regions/Slide0021991_KK4_465_Hypothalamus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_465_regions/Slide0021991_KK4_465_Thalamus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_465_regions/Slide0021991_KK4_465_WM_cells_stats.csv", 
  # KK4_492 regional annotations
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_492_regions/Slide0021991_KK4_492_Cortex_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_492_regions/Slide0021991_KK4_492_Hippocampus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_492_regions/Slide0021991_KK4_492_Hypothalamus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_492_regions/Slide0021991_KK4_492_Thalamus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_492_regions/Slide0021991_KK4_492_WM_cells_stats.csv", 
  # KK4_496 regional annotations
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_496_regions/Slide0021998_KK4_496_WM_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_496_regions/Slide0021998_KK4_496_Cortex_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_496_regions/Slide0021998_KK4_496_Hippocampus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_496_regions/Slide0021998_KK4_496_Hypothalamus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_496_regions/Slide0021998_KK4_496_Thalamus_cells_stats.csv", 
  # KK4_502 regional annotations
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_502_regions/Slide0021998_KK4_502_Cortex_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_502_regions/Slide0021998_KK4_502_Hippocampus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_502_regions/Slide0021998_KK4_502_Hypothalamus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_502_regions/Slide0021998_KK4_502_Thalamus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_502_regions/Slide0021998_KK4_502_WM_cells_stats.csv", 
  # KK4_504 regional annotations
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_504_regions/Slide0021991_KK4_504_Cortex_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_504_regions/Slide0021991_KK4_504_Hippocampus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_504_regions/Slide0021991_KK4_504_Hypothalamus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_504_regions/Slide0021991_KK4_504_Thalamus_cells_stats.csv", 
  "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/ARIA_ExplorerAnnotations/KK4_504_regions/Slide0021991_KK4_504_WM_cells_stats.csv"
)


annotation_list <- map(files, function(path) {
  
  # read file
  df <- read.csv(path, skip = 2, blank.lines.skip = TRUE)
  
  # extract filename
  fname <- basename(path)
  
  # Extract sample ID (KK4_###)
  sample_id <- str_extract(fname, "KK4_\\d+")
  
  # Extract region (word after KK4_###_)
  region <- str_extract(fname, "KK4_\\d+_([A-Za-z]+)") |>
    str_remove("KK4_\\d+_")
  
  # determine slide prefix
  slide_prefix <- case_when(
    sample_id %in% c("KK4_465","KK4_492","KK4_504") ~ "slide1_",
    sample_id %in% c("KK4_464","KK4_496","KK4_502") ~ "slide2_"
  )
  
  df %>%
    mutate(
      Cell.ID  = paste0(slide_prefix, Cell.ID),
      Sample_ID = sample_id,
      Region    = region
    )
})

annotation_df <- bind_rows(annotation_list)

annotation_df_unique <- annotation_df %>%
  distinct(Cell.ID, .keep_all = TRUE)

rownames(annotation_df_unique) <- annotation_df_unique$Cell.ID

annotation_df_unique <- annotation_df_unique[Cells(merged.obj), , drop = FALSE]


merged.obj <- AddMetaData(merged.obj, metadata = annotation_df_unique["Region"])



## 03c. adding metadata from excel file ----

ARIA_SC_Prep_Info <- data.frame(
  read_excel("/Users/cclu223/Desktop/Single_Cell_Analysis/ARIA_dataset/Raw_Data/ARIA SC Prep Info.xlsx"))

# extract metadata from seurat

meta <- merged.obj@meta.data

# changing mouse.id column name to sample.id 

colnames(ARIA_SC_Prep_Info)[colnames(ARIA_SC_Prep_Info) == "Mouse.ID"] <- "sample_ID"

# convert rownames to columns to preserve cell IDs

meta <- meta %>%
  rownames_to_column(var = "cell_id")

# joining metadata from excel sheet to metadata in merged seurat obj

meta <- meta %>%
  left_join(ARIA_SC_Prep_Info, by = "sample_ID"
            #, keep = T
            )

# convert column back to rownames to preserve cell IDs 

meta <- column_to_rownames(meta, var = "cell_id")

# add metadata back to merged obj

merged.obj@meta.data <- meta

# SANITY CHECK
dim(merged.obj)

### 04. Preprocessing ----

merged.obj <- UpdateSeuratObject(merged.obj)

#removing cells with less than 5 counts

Idents(merged.obj) <- "sample_ID"

xenium.aria <- subset(merged.obj, subset = sample_ID %in% c("KK4_465", "KK4_492", "KK4_504", "KK4_502", "KK4_464", "KK4_496"))

xenium.aria <- subset(xenium.aria, subset = nCount_Xenium > 5)

QC_vln <- VlnPlot(xenium.aria, features = c("nFeature_Xenium", "nCount_Xenium"), 
                  raster = F,  
                  ncol = 2, pt.size = 0)

ggsave("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/20251201_mergedQC_vln.png", plot = QC_vln)



### 05. Standard Seurat pipeline ----

xenium.aria <- SCTransform(xenium.aria, assay = "Xenium")

xenium.aria <- RunPCA(xenium.aria, npcs = 30, features = rownames(xenium.aria))

xenium.aria <- RunUMAP(xenium.aria, dims = 1:30)

xenium.aria <- FindNeighbors(xenium.aria, reduction = "pca", dims = 1:30)

xenium.aria <- FindClusters(xenium.aria, resolution = 0.3)

DimPlot(xenium.aria, label = T)

# CHECK FOR BATCH EFFECTS BY FEATURE PLOTTING A CELL MARKER THAT SHOULDNT CHANGE WITH TREATMENT OR AGE OR SOMETHING 

saveRDS(xenium.aria, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Files/20260126_mergedobj_QC.rds")



