rm(list=ls())
library(ArchR)
library(Seurat)
library(grid)

threads = 8
addArchRThreads(threads = threads)

addArchRGenome("mm10")

inputFiles <- './fragments.tsv.gz'
sampleNames <- 'ME11'

## Create ArchRProject
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 0,
  filterFrags = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  TileMatParams = list(tileSize = 5000)
)
ArrowFiles

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleNames,
  copyArrows = TRUE
)
proj


## Select pixels in tissue
meta.data <- as.data.frame(getCellColData(ArchRProj = proj))
meta.data['cellID_archr'] <- row.names(meta.data)
data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"
image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), 
                       filter.matrix = filter.matrix)
meta.data.spatial <- meta.data[row.names(image@coordinates), ]
proj_in_tissue <- proj[meta.data.spatial$cellID_archr, ]
proj_in_tissue


## Data normalization and dimensionality reduction 
proj_in_tissue <- addIterativeLSI(
  ArchRProj = proj_in_tissue,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

proj_in_tissue <- addClusters(
  input = proj_in_tissue,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5,
  force = TRUE
)

proj_in_tissue <- addUMAP(
  ArchRProj = proj_in_tissue, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

plotEmbedding(ArchRProj = proj_in_tissue, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", size = 1.5)

proj_in_tissue <- addImputeWeights(proj_in_tissue)


## Identify the marker genes for each cluster 
markersGS <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  testMethod = "wilcoxon"
)

markerList_pos <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")

markerGenes <- list()
for (i in seq_len(length(markerList_pos))) {
  markerGenes <- c(markerGenes, markerList_pos[[i]]$name)
}
markerGenes <- unlist(markerGenes)


## Call peaks
proj_in_tissue <- addGroupCoverages(ArchRProj = proj_in_tissue, groupBy = "Clusters")

pathToMacs2 <- findMacs2()

proj_in_tissue <- addReproduciblePeakSet(
  ArchRProj = proj_in_tissue, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

proj_in_tissue <- addPeakMatrix(proj_in_tissue)

getAvailableMatrices(proj_in_tissue)
getPeakSet(proj_in_tissue)

if("Motif" %ni% names(proj_in_tissue@peakAnnotation)){
  proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif", force = TRUE)
}

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_in_tissue,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.1")


## ChromVAR Deviatons Enrichment
proj_in_tissue <- addBgdPeaks(proj_in_tissue, force = TRUE)

proj_in_tissue <- addDeviationsMatrix(
  ArchRProj = proj_in_tissue, 
  peakAnnotation = "Motif",
  force = TRUE
)

markersMotifs <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "MotifMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = 'z'
)


## bulk sample (ENCODE) projection
projBulk <- projectBulkATAC(proj_in_tissue, seATAC = fragment_counts, n = 250)

encode_cell_types <- function(x) {
  x <- tools::file_path_sans_ext(x)
  meta_data_sub[x,]$Biosample.term.name
}
projBulk[[1]][,3] <- sapply(projBulk[[1]]$Type, encode_cell_types)

plotProj <- rbind(projBulk[[2]], projBulk[[1]])
pal <- paletteDiscrete(unique(as.vector(plotProj[,3])))
pal["spatial_ATAC"] <- "lightgrey"

p <- ggPoint(plotProj[,1], plotProj[,2], as.vector(plotProj[,3]), rastr = TRUE, pal = pal)
p


## Pseudotime analysis 
source('SpatialPlot_traj.R')

trajectory <- c("Radial glia", "Postmitotic premature neurons", "Excitatory neurons")

projCUTA <- addTrajectory(
  ArchRProj = proj_in_tissue, 
  name = "Neuron_U", 
  groupBy = "predictedGroup_Un",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE
)

meta.data.integration <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))[, c('Neuron_U'), drop=FALSE]
new_row_names <- row.names(meta.data.integration)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data.integration) <- new_row_names
all(row.names(meta.data.integration) == colnames(spatial.obj))

spatial.obj <- AddMetaData(object = spatial.obj, metadata = meta.data.integration)

p <- SpatialPlot_traj(spatial.obj, features = "Neuron_U",  pt.size.factor = 4, image.alpha = 0, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p
