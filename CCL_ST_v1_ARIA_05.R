########################################################################
# Name: CCL_ST_v1_ARIA_05.R
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
# Last updated: 2/20/26
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
library(speckle) # for propeller function
options(future.globals.maxSize = 100 *1024^3)

# libraries specific to figure-making
library(SCP) # for pretty UMAPs
library(dittoSeq) # for stacked bar plots
library(scCustomize)

# Load objs ----

ann_aria <- qs_read("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260217_Part2_fullobj.qs2")

microglia <- qs_read("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/aria_adu_KAI/mg_kai.qs2")


# Extracting metadata from full obj for SpatialData obj ----

meta <- ann_aria@meta.data

rownames(meta) <- colnames(ann_aria)

write_csv(meta, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Squidpy/fullobj_metadata.csv")

# Extracting metadata from microglia obj for SpatialData obj ----

mg_meta <- microglia@meta.data

rownames(mg_meta) <- colnames(microglia)

write_csv(mg_meta, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Squidpy/microglia_metadata.csv")

# 01. Formatting levels and selecting colors ----

## 01a. Re-ordering levels ----

# treatment group order 
ann_aria$Treatment.Group <- factor(
  ann_aria$Treatment.Group, 
  levels = c("IgG", "Adu")
)

microglia$Treatment.Group <- factor(
  microglia$Treatment.Group, 
  levels = c("IgG", "Adu")
)

# cluster order
ann_aria$annotatedclusters <- factor(
  ann_aria$annotatedclusters,
  levels = c(
    # vascular
    "CP", 
    "Endothelial",
    "Ependymal", 
    "Fibroblast",
    "Pericyte", 
    "VLMC", 
    "VSMC", 
    # immune
    "BAM", 
    "T Cell",
    "Microglia", 
    # other glia
    "Astrocyte", 
    "Oligodendrocyte", 
    "OPC", 
    # neurons
    "GABAergic Neuron 1", 
    "GABAergic Neuron 2", 
    "GABAergic Neuron 3", 
    "GABAergic Neuron 4", 
    "Glutamatergic Neuron 1", 
    "Glutamatergic Neuron 2", 
    "Glutamatergic Neuron 3",
    "Glutamatergic Neuron 4", 
    "Glutamatergic Neuron 5"
  )
)

# updating active idents to reflect the new order 
Idents(ann_aria) <- "annotatedclusters"

## 01b. Define cluster palette ----
cluster_palette <- c(
  "CP" = "#A6CEE3", 
  "Endothelial" = "#5FA0CA", 
  "Ependymal" = "#257CB2", 
  "Fibroblast" = "#72B29C", 
  "Pericyte" = "#A5D981", 
  "VLMC" = "#63B84F", 
  "VSMC" = "#4FA435", 
  "BAM" = "#B9B458", 
  "T Cell" = "#1F9FB5",
  "Microglia" = "#FDB259", 
  "Astrocyte" = "#FE8524", 
  "Oligodendrocyte" = "#FB9374", 
  "OPC" = "#F47575", 
  "GABAergic Neuron 1" = "#E73233", 
  "GABAergic Neuron 2" = "#DA4C59", 
  "GABAergic Neuron 3" = "#CD9CBB", 
  "GABAergic Neuron 4" = "#A585BF", 
  "Glutamatergic Neuron 1" = "#73489F", 
  "Glutamatergic Neuron 2" = "#A99099", 
  "Glutamatergic Neuron 3" = "#F7F599", 
  "Glutamatergic Neuron 4" = "#D9AF63", 
  "Glutamatergic Neuron 5" = "#B15928"
) 


# 02. Creating figures ----

## 02a. UMAPs (using SCP) ----
p1_UMAP <- CellDimPlot(
  srt = ann_aria, 
  group.by = "annotatedclusters",
  reduction = "umap_SPLIT", 
  theme_use = "theme_blank",
  pt.size = 0.6
) + scale_color_manual(values = cluster_palette) +
  guides(color = guide_legend(title = NULL, # remove legend title
                                ncol = 1, # keep legend in 1 column
                                override.aes = list(size = 3))) # legend dot size

p2_UMAP <- CellDimPlot(
  srt = microglia, 
  group.by = "celltype",
  reduction = "umap", 
  theme_use = "theme_blank",
  palette = "Paired",
  pt.size = 0.6
) + guides(color = guide_legend(title = NULL, # remove legend title
                              ncol = 1, # keep legend in 1 column
                              override.aes = list(size = 3))) # legend dot size

### extract colors from microglia umap ----
# get colour data from the ggplot build
scale_obj <- p2_UMAP$scales$get_scales("colour")

# set the names for each color to the corresponding celltype 
cluster_palette<- setNames(
  scale_obj$palette(length(levels(microglia$celltype))),
  levels(microglia$celltype)
  )

# format output from cluster palette into a dataframe
mg_cluster_palette <- c(
  "Hm1" = "#A6CEE3", 
  "Hm2" = "#1F78B4", 
  "Hm3" = "#B2DF8A", 
  "DAM" = "#33A02C", 
  "Gpnmb+" = "#FDBF6F", 
  "Spp1+" = "#FF7F00", 
  "ARM" = "#FB9A99", 
  "IRM" = "#E31A1C", 
  "Il1b+" = "#CAB2D6", 
  "Proliferating" = "#6A3D9A"
) 

## 02b. Propeller cell-type proportion analysis ----
# adapted from Jose Arbones-Mainar

### for fullobj ----
all_prop <- propeller(
  x = ann_aria,
  clusters = ann_aria$annotatedclusters,
  sample = ann_aria$sample_ID,
  group = ann_aria$Treatment.Group,
  trend = FALSE,
  robust = TRUE,
  transform = "logit")

### for microglia obj -----
mg_prop <- propeller(
  x = microglia,
  clusters = microglia$celltype, # CHECK THIS!!
  sample = microglia$sample_ID,
  group = microglia$Treatment.Group,
  trend = FALSE,
  robust = TRUE,
  transform = "logit")


### microglia by region ----

# updating seurat object to prevent error in map() loop
microglia <- UpdateSeuratObject(microglia)

# Wrap propeller so that regions where the model fails return NULL instead of an error
safe_propeller <- possibly(propeller, otherwise = NULL)

# Get all annotated regions, dropping NA entries
regions <- unique(microglia$Region) |> na.omit() |> as.character()

# Run propeller separately for each brain region
prop_by_region <- map(regions, \(reg) {
  # Subset to cells from this region only
  sub <- subset(microglia, Region == reg)
  
  # Skip if either treatment group has < 2 samples (propeller requires replication)
  n_per_group <- table(
    distinct(sub@meta.data, sample_ID, Treatment.Group)$Treatment.Group
  )
  if (any(n_per_group < 2)) return(NULL)
  
  # Test for proportion differences between Adu and IgG within this region
  safe_propeller(
    x         = sub,
    clusters  = Idents(sub),
    sample    = sub$sample_ID,
    group     = sub$Treatment.Group,
    trend     = FALSE,
    robust    = TRUE,       # down-weights outlier samples
    transform = "logit"     # variance-stabilising transform for bounded proportions
  )
}) |>
  # Name list elements by region, then drop NULLs (skipped/failed regions)
  set_names(regions) |>
  compact()

# Combine per-region results into a single data frame with a Region column
prop_region_df <- imap(prop_by_region, \(tbl, reg) {
  as_tibble(tbl) |> mutate(Region = reg)
}) |>
  list_rbind()


### heatmap for microglia subcluster proportion stats by region ----
prop_heatmap <- prop_region_df |>
  mutate(
    # Convert PropRatio (Adu/IgG) to log2 scale; set Inf/NaN (zero in one group) to NA
    log2ratio = case_when(
      is.nan(log2(PropRatio)) | is.infinite(log2(PropRatio)) ~ NA_real_,
      TRUE ~ log2(PropRatio)
    ),
    # Overlay significance stars on coloured tiles, and × on grey (zero-proportion) tiles
    stars = case_when(
      is.na(log2ratio)   ~ "\u00d7",   # cross for 0-proportion tiles
      P.Value < 0.001    ~ "***",
      P.Value < 0.01     ~ "**",
      P.Value < 0.05     ~ "*",
      TRUE               ~ ""
    )
  )

# Order cell types by mean log2ratio across regions (Adu-enriched at top)
celltype_order <- prop_heatmap |>
  group_by(BaselineProp.clusters) |>
  summarise(mean_r = mean(log2ratio, na.rm = TRUE)) |>
  arrange(desc(mean_r)) |>
  pull(BaselineProp.clusters)

prop_heatmap <- prop_heatmap |>
  mutate(
    CellType  = factor(BaselineProp.clusters, levels = rev(celltype_order)),
    Region    = factor(Region),
    # Use a darker grey for × marks so they remain visible on the grey background
    cross_col = ifelse(is.na(log2ratio), "grey30", "grey15")
  )

# Symmetric colour scale limit (rounded up to nearest integer)
clim <- ceiling(max(abs(prop_heatmap$log2ratio), na.rm = TRUE))

mg_regions_prop_heatmap <- ggplot(prop_heatmap, aes(x = Region, y = CellType, fill = log2ratio)) +
  # Tile fill encodes direction and magnitude of the proportion shift
  geom_tile(colour = "white", linewidth = 0.5) +
  coord_equal() +
  # Significance stars (or ×) overlaid on each tile
  geom_text(aes(label = stars, colour = cross_col), size = 10, vjust = 0.75) +
  # Pass colour strings directly without mapping to a scale
  scale_colour_identity() +
  # Diverging palette: red = Adu-enriched, blue = IgG-enriched, grey = undefined
  scale_fill_distiller(
    palette   = "RdBu",
    direction = -1,
    limits    = c(-clim, clim),
    na.value  = "grey88",
    name      = "log2(Adu/IgG)"
  ) +
  labs(
    title    = "Cell type proportions by brain region",
    subtitle = "Fill: log2(PropRatio Adu/IgG);  * p<0.05  ** p<0.01",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 18),
    axis.text.y     = element_text(size = 18),
    panel.grid      = element_blank(),
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 18),
    plot.subtitle   = element_text(size = 16, colour = "grey40")
  )



## 02c. Stacked bar plots (using dittoseq.) ----

b1_barplot <- dittoBarPlot(ann_aria, 
                           var = "annotatedclusters", 
                           group.by = "sample_ID",
                           split.by = "Treatment.Group",
                           retain.factor.levels = TRUE,  
                           color.panel = cluster_palette, 
                           legend.title = NULL, 
                           main = NULL
                           ) + facet_wrap(~Treatment.Group, scales = "free_x", drop = TRUE)




b2_barplot <- dittoBarPlot(ann_aria, 
             var = "annotatedclusters", 
             group.by = "Treatment.Group", 
             color.panel = cluster_palette,
             retain.factor.levels = TRUE, 
             legend.title = NULL, 
             main = NULL)


b3_barplot <- dittoBarPlot(ann_aria, 
                           var = "annotatedclusters", 
                           group.by = "Treatment.Group",
                           split.by = "Region",
                           retain.factor.levels = TRUE,  
                           color.panel = cluster_palette, 
                           legend.title = NULL, 
                           main = NULL
) 


b4_barplot <- dittoBarPlot(ann_aria, 
                           var = "annotatedclusters", 
                           group.by = "Region",
                           retain.factor.levels = TRUE,  
                           color.panel = cluster_palette, 
                           legend.title = NULL, 
                           main = NULL
) 

b5_barplot <- dittoBarPlot(microglia, 
                           var = "celltype", 
                           group.by = "Treatment.Group",
                           split.by = "Region",
                           retain.factor.levels = TRUE,
                           legend.title = NULL, 
                           main = NULL, 
                           color.panel = mg_cluster_palette
) 

b6_barplot <- dittoBarPlot(microglia, 
                           var = "celltype", 
                           group.by = "Treatment.Group",
                           retain.factor.levels = TRUE,  
                           legend.title = NULL, 
                           main = NULL, 
                           color.panel = mg_cluster_palette
) 

b7_barplot <- dittoBarPlot(microglia, 
                           var = "celltype", 
                           group.by = "Region",
                           retain.factor.levels = TRUE,  
                           legend.title = NULL, 
                           main = NULL, 
                           color.panel = mg_cluster_palette
)

b8_barplot <- dittoBarPlot(microglia, 
                           var = "celltype", 
                           group.by = "sample_ID",
                           split.by = "Treatment.Group",
                           retain.factor.levels = TRUE,  
                           legend.title = NULL, 
                           main = NULL, 
                           color.panel = mg_cluster_palette
) + facet_wrap(~Treatment.Group, scales = "free_x", drop = TRUE)


## 02d. ImageDimPlots ----
imagedim_fov1 <- ImageDimPlot(ann_aria, 
                              group.by = "annotatedclusters", 
                              cols = cluster_palette, 
                              dark.background = F, 
                              fov = "fov") + labs(fill = NULL)



imagedim_fov2 <- ImageDimPlot(ann_aria, 
                              group.by = "annotatedclusters", 
                              cols = cluster_palette, 
                              dark.background = F, 
                              fov = "fov.2") + labs(fill = NULL)



imagedim_split_fov1 <- ImageDimPlot(ann_aria, 
                                    split.by = "annotatedclusters", 
                                    cols = cluster_palette, 
                                    dark.background = F, 
                                    fov = "fov") + labs(fill = NULL)



imagedim_split_fov2 <- ImageDimPlot(ann_aria, 
                                    split.by = "annotatedclusters", 
                                    cols = cluster_palette, 
                                    dark.background = F, 
                                    fov = "fov.2") + labs(fill = NULL)

# microglia obj
imagedim_mg_fov1 <- ImageDimPlot(microglia, 
                              group.by = "celltype",
                              cols = mg_cluster_palette,
                              dark.background = F, 
                              fov = "fov") + labs(fill = NULL)



imagedim_mg_fov2 <- ImageDimPlot(microglia, 
                              group.by = "celltype", 
                              cols = mg_cluster_palette,
                              dark.background = F, 
                              fov = "fov.2") + labs(fill = NULL)



imagedim_mg_split_fov1 <- ImageDimPlot(microglia, 
                                    split.by = "celltype",                              
                                    cols = mg_cluster_palette,
                                    dark.background = F, 
                                    fov = "fov") + labs(fill = NULL)



imagedim_mg_split_fov2 <- ImageDimPlot(microglia, 
                                    split.by = "celltype",
                                    cols = mg_cluster_palette,
                                    dark.background = F, 
                                    fov = "fov.2") + labs(fill = NULL)


# 03. Saving plots ----

# umaps
ggsave(filename = "20260306_UMAP_all.svg", plot = p1_UMAP, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = svglite::svglite, dpi = 600)
ggsave(filename = "20260311_mg_UMAP.pdf", plot = p2_UMAP, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)

# heatmap 
ggsave(filename = "20260308_heatmap_mg_region_props.pdf", plot = mg_regions_prop_heatmap, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600, width = 8, height = 7)

# stacked bar plots
ggsave(filename = "20260306_barplot_bysamples.pdf", plot = b1_barplot, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)
ggsave(filename = "20260306_barplot_bytreatment.pdf", plot = b2_barplot, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)
ggsave(filename = "20260306_barplot_byregion_and_treatment.pdf", plot = b3_barplot, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600, width = 8, height = 11)
ggsave(filename = "20260306_barplot_byregion.pdf", plot = b4_barplot, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)

ggsave(filename = "20260311_mg_barplot_bysamples.pdf", plot = b8_barplot, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)
ggsave(filename = "20260311_mg_barplot_bytreatment.pdf", plot = b6_barplot, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)
ggsave(filename = "20260311_mg_barplot_byregion_and_treatment.pdf", plot = b5_barplot, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)
ggsave(filename = "20260311_mg_barplot_byregion.pdf", plot = b7_barplot, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)


# ImageDimPlots
ggsave(filename = "20260306_imagedimplot_clusters_fov1.pdf", plot = imagedim_fov1, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)
ggsave(filename = "20260306_imagedimplot_clusters_fov2.pdf", plot = imagedim_fov2, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)
ggsave(filename = "20260306_imagedimplot_splitclusters_fov1.pdf", plot = imagedim_split_fov1, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600, width = 11, height = 8)
ggsave(filename = "20260306_imagedimplot_splitclusters_fov2.pdf", plot = imagedim_split_fov2, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600, width = 11, height = 8)

ggsave(filename = "20260311_imagedimplot_mg_clusters_fov1.pdf", plot = imagedim_mg_fov1, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)
ggsave(filename = "20260311_imagedimplot_mg_clusters_fov2.pdf", plot = imagedim_mg_fov2, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)
ggsave(filename = "20260311_imagedimplot_mg_splitclusters_fov1.pdf", plot = imagedim_mg_split_fov1, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600, width = 11, height = 8)
ggsave(filename = "20260311_imagedimplot_mg_splitclusters_fov2.pdf", plot = imagedim_mg_split_fov2, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600, width = 11, height = 8)

