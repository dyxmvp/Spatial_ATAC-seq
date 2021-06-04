library(ArchR)
library(Seurat)
library(grid)

source('SpatialDimPlot_new.R')

## Prepare MOCA data
MOCA_dir <- "./ref_data/MOCA/"

meta.data.RNA <- read.csv(file = paste0(MOCA_dir, 'cell_annotate.csv'), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gene.ANN.RNA <- read.csv(file = paste0(MOCA_dir, 'gene_annotate.csv'), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
gene.ANN.RNA <- gene.ANN.RNA[, 'gene_short_name', drop = FALSE]

cds <- readRDS(paste0(MOCA_dir, 'gene_count_cleaned_sampled_100k.RDS'))

MOCA <- CreateSeuratObject(counts = cds, project = 'MOCA')
meta.data.RNA <- meta.data.RNA[colnames(MOCA), ]
meta.data.RNA <- meta.data.RNA[, c('Main_cell_type', 'development_stage')]

MOCA <- AddMetaData(object = MOCA, metadata = meta.data.RNA)
MOCA_E11 <- subset(MOCA, development_stage == 11.5)
MOCA_E11.raw.data <- as.matrix(GetAssayData(MOCA_E11, slot = 'counts'))
MOCA_E11.raw.data <- as.data.frame(MOCA_E11.raw.data)
MOCA_E11.raw.data <- merge(gene.ANN.RNA, MOCA_E11.raw.data, by=0, all=TRUE)
which(is.na(MOCA_E11.raw.data$gene_short_name))

tt <- table(MOCA_E11.raw.data$gene_short_name)
name_rep <- names(which(tt > 1))
row_del_fun <- function(x){
  rows <- which(MOCA_E11.raw.data$gene_short_name == x)
  return(rows[2:length(rows)] )
}
row_del <- unlist(lapply(name_rep, row_del_fun))
MOCA_E11.raw.data <- MOCA_E11.raw.data[-row_del, ]

row.names(MOCA_E11.raw.data) <- MOCA_E11.raw.data$gene_short_name
MOCA_E11.raw.data <- MOCA_E11.raw.data[, -c(1:2), drop=FALSE]
MOCA_E11 <- CreateSeuratObject(counts = MOCA_E11.raw.data, project = 'MOCA_E11', meta.data = MOCA_E11@meta.data)


## Integration with ArchR oject
proj_in_tissue <- addGeneIntegrationMatrix(
  ArchRProj = proj_in_tissue, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = MOCA_E11,
  addToArrow = TRUE,
  groupRNA = "Main_cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)


## Plot results
meta.data.integration <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))[, c('predictedCell', 'predictedGroup', 'predictedScore')]
new_row_names <- row.names(meta.data.integration)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data.integration) <- new_row_names

spatial.obj <- AddMetaData(object = spatial.obj, metadata = meta.data.integration)

Idents(spatial.obj) <- 'predictedGroup'

ids.highlight <- names(table(spatial.obj$predictedGroup))[1]
ids.highlight

p <- SpatialDimPlot_new(spatial.obj, cells.highlight = CellsByIdentities(object = spatial.obj, idents = ids.highlight), 
                        facet.highlight = TRUE, pt.size.factor = 2.5, alpha = c(1,0), stroke = 0)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p