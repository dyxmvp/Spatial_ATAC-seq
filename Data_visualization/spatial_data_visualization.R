ibrary(ArchR)
library(Seurat)
library(grid)
library(ggplot2)
library(patchwork)
library(dplyr)

source('getGeneScore_ArchR.R')
source('SpatialPlot_new.R')

## Prepare meta data
meta.data <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

## Create Seurat object 
proj_in_tissue <- addImputeWeights(proj_in_tissue)
gene_score <- getGeneScore_ArchR(ArchRProj = proj_in_tissue, name = markerGenes, imputeWeights = getImputeWeights(proj_in_tissue))

data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = gene_score, assay = assay, meta.data = meta.data)

image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object


## Plot results
p <- SpatialPlot_new(spatial.obj, features = "nFrags",  pt.size.factor = 3, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p

n_clusters <- length(unique(projCUTA$Clusters))
cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
names(cols) <- paste0('C', seq_len(n_clusters))
cols

p <- SpatialPlot(spatial.obj, label = FALSE, label.size = 3, group.by = 'Clusters', pt.size.factor = 3, cols = cols, image.alpha = 0, stroke = 0)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p

feature <- 'Rarg'
p <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 3, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p