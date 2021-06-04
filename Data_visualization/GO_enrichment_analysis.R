rm(list=ls())

library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)

## Read markers genes
data_dir <- './data/ME11/markers_list/ME11_C1_markers.txt'
markerList <- read.table(data_dir, header = TRUE, stringsAsFactors = FALSE)
markers <- markerList$name

## Run GO enrichment analysis
ego <- enrichGO(gene         = markers,
                OrgDb         = org.Mm.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

## Plot results
p <- dotplot(ego, showCategory=30)
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
p

png(filename = './data/ME11/ME11_C1_markers_GO.png', width = 2400, height = 3600, res = 300)
p
dev.off()

