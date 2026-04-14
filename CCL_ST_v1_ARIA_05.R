########################################################################
# Name: CCL_ST_v1_ARIA_05.R
# Project: mouse low dose ARIA (Akhil's project)
# Purpose: Analysis of low dose ARIA Xenium Spatial Transcriptomics data
#         - Part 1: Merging files, attaching metadata and pre-processing 
#         - Part 2: RCTD and SPLIT
#         - Part 3: manually annotating all clusters 
#         - Part 4: subclustering, cleaning and manually annotating microglia 
#         - Part 5: fixed slide and sample ID misalignment in Kai's microglia object, extracted final seurat objs metadata for squidpy and spatial data analyses, and code for publishable figures  
# Date created : 2/20/26
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
library(stringr)
library(speckle) # for propeller function
options(future.globals.maxSize = 100 *1024^3)

# libraries specific to figure-making
library(SCP) # for pretty UMAPs
library(dittoSeq) # for stacked bar plots
library(scCustomize)
library(pheatmap)
library(cowplot) # for imagedimplot with inset and cropped region


# Load objs ----

ann_aria <- qs_read("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260320_fullobj.qs2") # NEW SLIDE AND SAMPLE_ID ALIGNED OBJ

microglia <- qs_read("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/aria_adu_KAI/mg_kai.qs2")

# correcting Slide and sample_ID metadata column misalignment using obj with correct metadata (also used this object on 20260320_fullobj.qs1)
old_aria <- LoadSeuratRds("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260213_Part1_subclusv3.rds")

all(
  microglia$Treatment.Group ==
    old_aria$Treatment.Group[colnames(microglia)]
)


# FIX TO RE-ALIGN SLIDE AND SAMPLE_ID METADATA COLUMNS IN MICROGLIA OBJ
meta_clean <- old_aria@meta.data[, c("Slide", "sample_ID", "Treatment.Group")]

# make sure cells exist
all(colnames(microglia) %in% rownames(meta_clean))

# align properly
meta_clean <- meta_clean[colnames(microglia), ]

# overwrite corrupted metadata
microglia$Slide <- meta_clean$Slide
microglia$sample_ID <- meta_clean$sample_ID
microglia$Treatment.Group <- meta_clean$Treatment.Group



qs_save(microglia, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/aria_adu_KAI/20260323_updated_mg_kai.qs2") # UPDATED OBJ

# read in new microglia obj
microglia <- qs_read("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/aria_adu_KAI/20260323_updated_mg_kai.qs2")

# Extracting metadata from full obj for SpatialData obj ----

meta <- ann_aria@meta.data

rownames(meta) <- colnames(ann_aria)

write_csv(meta, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Squidpy/20260323_updated_fullobj_metadata.csv")

# Extracting metadata from microglia obj for SpatialData obj ----

mg_meta <- microglia@meta.data

rownames(mg_meta) <- colnames(microglia)

write_csv(mg_meta, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Squidpy/20260323_updated_microglia_metadata.csv")

# Extracting metadata from Banksy obj (obj from Jose)

banksy <- qs_read("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/rds_files/20260402_fullobj_correctedChloe_sct_SPLIT_microglia_cleaned.qs2_banksy_agf=TRUE_sct.qs2")

banksy_meta <- banksy@meta.data

rownames(banksy_meta) <- colnames(banksy)

write_csv(banksy_meta, "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Squidpy/20260413_banksy_metadata.csv")

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

ann_aria_clean <- subset(ann_aria, !is.na(annotatedclusters))

b1_barplot <- dittoBarPlot(ann_aria_clean, 
                           var = "annotatedclusters", 
                           group.by = "sample_ID",
                           split.by = "Treatment.Group",
                           retain.factor.levels = TRUE,  
                           color.panel = cluster_palette, 
                           legend.title = NULL, 
                           main = NULL
                           ) + facet_wrap(~Treatment.Group, scales = "free_x", drop = TRUE)




b2_barplot <- dittoBarPlot(ann_aria_clean, 
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


### rep image dim plot with inset and zoomed in portion ----






# find spatial coordinates to subset one brain from fov (PLOT DATA)
data_full <- ggplot_build(
  ImageDimPlot(ann_aria, fov = "fov", dark.background = FALSE)
)$data[[1]]

# extract full x and y range coords 
full_xrange <- diff(range(data_full$x)) # total width of main plot in data coords 
full_yrange <- diff(range(data_full$y)) # total height of main plot in data coords 

# --- Preview crop region before committing to object ---


# Filter to proposed crop region
xlim_preview <- c(5150, 12000)
ylim_preview <- c(0, 8000)

# Get the spatial coordinates from the FOV
coords <- FetchData(ann_aria, 
                    vars = c("annotatedclusters", "x", "y"),
                    layer = "fov")

coords_cropped <- coords[
  coords$x >= xlim_preview[1] & coords$x <= xlim_preview[2] &
    coords$y >= ylim_preview[1] & coords$y <= ylim_preview[2], 
]


# Quick preview plot before adding new cropped fov to obj
ggplot(coords_cropped, aes(x = x, y = y, color = annotatedclusters)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = cluster_palette) +
  coord_fixed() +   # preserves true aspect ratio
  labs(title = paste0("Preview: x [", xlim_preview[1], "-", xlim_preview[2], 
                      "], y [", ylim_preview[1], "-", ylim_preview[2], "]")) +
  guides(color = guide_legend(override.aes = list(size = 3)))



# NOTE:
# ImageDimPlot swaps x and y coordinates relative to GetTissueCoordinates
# so when I define crop regions I will swap the x and y limits in the GetTissueCoordinates space 
# EXAMPLE: x_plot = y_tissue, y_plot = x_tissue

# crop to 1 brain (KK4_492 in this case)
KK4_492_crop <- Crop(ann_aria[["fov"]], 
                   x = ylim_preview,    # FLIPPED (per note above)
                   y = xlim_preview, # FLIPPED (per note above)
                   coords = "plot")


# add cropped fov to obj
ann_aria[["KK4_492_fov"]] <- KK4_492_crop


# plot cropped brain only 
p_main <- ImageDimPlot(ann_aria, 
                       group.by = "annotatedclusters", 
                       cols = cluster_palette,
                       dark.background = F, 
                       fov = "KK4_492_fov", 
                       size = 2)



# pick x and y limits for cropped region (based on ImageDimPlot)

xlim_zoom <- c(1500, 3500)
ylim_zoom <- c(8000, 10000)

# generate x and y ranges from cropped/zoomed region
zoom_xrange <- diff(xlim_zoom)
zoom_yrange <- diff(ylim_zoom)



# generate zoomed in plot with fixed aspect ratio
p_zoom <- p_main +
  coord_cartesian(
    xlim = xlim_zoom,
    ylim = ylim_zoom
  ) +
  theme_void() +   # removes axes for cleaner inset
  theme(legend.position = "none")

# annotation (box) on main plot 
p_main_box <- p_main +
  annotate(
    "rect",
    xmin = xlim_zoom[1], xmax = xlim_zoom[2],
    ymin = ylim_zoom[1], ymax = ylim_zoom[2],
    color = "black", fill = NA, linewidth = 0.8
  ) + theme(legend.position = "none")


## 02e. Neighborhood Enrichment Heatmaps ----

# read in neighborhood enrichment file on FULL BRAINS
nhood_enrich <- read.csv("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Neighborhood_Enrichment/nhood_enrichment_Adu-IgG.csv", row.names = 1)


# convert to data matrix for pheatmap compatibility 
nhood_enrich_mat <- data.matrix(nhood_enrich)

# Define symmetric breaks around 0
max_val <- max(abs(nhood_enrich_mat), na.rm = TRUE)

breaks <- seq(-max_val, max_val, length.out = 101)


full_brain_heatmap <- pheatmap(nhood_enrich_mat, color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100), 
         breaks = breaks,
         display_numbers = TRUE, 
         cluster_rows = FALSE,      # turn off row dendrogram
         cluster_cols = FALSE, # turn off column dendrogram
         # Tweak labels
         angle_col = 45,           # rotate column labels 45° (diagonal)
         fontsize_row = 10,        # adjust row font size
         fontsize_col = 10,        # adjust column font size
         labels_row = rownames(nhood_enrich_mat),  # ensure row names are correct
         labels_col = colnames(nhood_enrich_mat)) 


# Plaque niche 
# read in neighborhood enrichment file on FULL BRAINS
nhood_enrich_plaque <- read.csv("/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Neighborhood_Enrichment/plaque_niche/plaque_niche_nhoodenrich_Adu-IgG.csv", row.names = 1)


# convert to data matrix for pheatmap compatibility 
nhood_enrich_plaque_mat <- data.matrix(nhood_enrich_plaque)

# Define symmetric breaks around 0
plaque_max_val <- max(abs(nhood_enrich_plaque_mat), na.rm = TRUE)

plaque_breaks <- seq(-plaque_max_val, plaque_max_val, length.out = 101)


plaque_heatmap <- pheatmap(nhood_enrich_plaque_mat, color = colorRampPalette(c("#3B4CC0", "white", "#B40426"))(100), 
                               breaks = plaque_breaks,
                               display_numbers = TRUE, 
                               cluster_rows = FALSE,      # turn off row dendrogram
                               cluster_cols = FALSE, # turn off column dendrogram
                               # Tweak labels
                               angle_col = 45,           # rotate column labels 45° (diagonal)
                               fontsize_row = 10,        # adjust row font size
                               fontsize_col = 10,        # adjust column font size
                               labels_row = rownames(nhood_enrich_plaque_mat),  # ensure row names are correct
                               labels_col = colnames(nhood_enrich_plaque_mat)) 





# 03. Saving plots ----

# umaps
ggsave(filename = "20260306_UMAP_all.svg", plot = p1_UMAP, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = svglite::svglite, dpi = 600)
ggsave(filename = "20260311_mg_UMAP.pdf", plot = p2_UMAP, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)

# heatmap 
ggsave(filename = "20260308_heatmap_mg_region_props.pdf", plot = mg_regions_prop_heatmap, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600, width = 8, height = 7)

# stacked bar plots
ggsave(filename = "20260413_barplot_bysamples.pdf", plot = b1_barplot, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600, width = 8, height = 11)
ggsave(filename = "20260413_barplot_bytreatment.pdf", plot = b2_barplot, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600)
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
ggsave(filename = "20260413_imagedimplot_KK4_492_box.pdf", plot = p_main_box, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600, width = 11, height = 8)
ggsave(filename = "20260413_imagedimplot_KK4_492_zoomedregion.pdf", plot = p_zoom, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600, width = 11, height = 8)

# heatmaps
ggsave(filename = "nhoodenrich_heatmap.pdf", plot = full_brain_heatmap, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600, width = 11, height = 8)
ggsave(filename = "plaque_niche_nhoodenrich_heatmap.pdf", plot = plaque_heatmap, path = "/Users/cclu223/Desktop/Xenium_ST_Analysis/Aging_Metabolism_Runs/ARIA_Analysis/Figures", device = "pdf", dpi = 600, width = 11, height = 8)
