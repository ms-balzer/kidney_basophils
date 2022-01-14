library(SCENIC)
library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
library(patchwork)
library(RcisTarget)
library(KernSmooth)
library(AUCell)
library(RColorBrewer)
library (pheatmap)
#====================== Load Seurat Object & Extract the Expression Matrix ===================
#Load PT Subset of UUO Samples
UUO.PT = readRDS ("UUO.PT.rds")
#Extract Expression Matrix
exprMat<-GetAssayData(UUO.PT, slot = "counts")
exprMat<-as.matrix(exprMat)
cellInfo <- data.frame(seuratCluster=Idents(UUO.PT))

#Run SCENIC
org <- "mgi"
scenicOptions <- initializeScenic(org="mgi", dbDir="YOURDIRECTORY", nCores=16)


scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]


runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)

#Run the remaining steps
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123

# Build and score the Regulons
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

# Binarize activity of Regulons
runSCENIC_4_aucell_binarize(scenicOptions)

#Plot AUC
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)

#Average Regulon Activity by cluster 
cellInfo <- readRDS("int/cellInfo.Rds")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
write.csv(topRegulators, 'output/Step4_regulonActivity_byEachCellType_ScaledVersion.csv')

#Visulization of Top Regulons
UUO.PT.TopRegulons = read.csv("Step4_regulonActivity_byEachCellType_ScaledVersion.csv")
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)
tiff('UUO.PT.TopRegulons.tiff', units="in", width=7, height=30, res=300, compression = 'lzw')
pheatmap(UUO.PT.TopRegulons, cluster_rows = T,
+          cluster_cols = F,
+          scale="row",
+          color=my_palette,
+          treeheight_row=30, treeheight_col=30, border_color="white",
+          show_rownames=T,
+          annotation_names_row=F,
+          cellwidth=15, cellheight=10)
> dev.off ()


#Visualization as UMAP
#----------- UMAP embeddings from Seurat Object---------------
dr_coords <- UUO.PT@reductions$umap@cell.embeddings
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]



# Bcl3, Cebp, Stat3, Klf6, Stat5a, Ddit3 Visulaization
regulons <- loadInt(scenicOptions, "regulons")
GOI <- regulons$Bcl3#[1:length(regulons$Bcl3)]
Cairo::CairoPDF("output/Step4_UMAP_Bcl3_targets.pdf", height=20, width=20)
par(mfrow=c(4,4))
AUCell::AUCell_plotTSNE(dr_coords,
                        exprMat = exprMat_filtered,
                        cellsAUC=selectRegulons(regulonAUC, GOI), 
                        plots = c("histogram", "binaryAUC", "AUC", "expression"),
                        exprCols = c("#440154FF", "#238A8DFF", "#FDE725FF"))
dev.off()

GOI <- regulons$Cebpb#[1:length(regulons$Cebpb)]
Cairo::CairoPDF("output/Step4_UMAP_Cebpb_targets.pdf", height=20, width=20)
par(mfrow=c(4,4))
AUCell::AUCell_plotTSNE(dr_coords,
                        exprMat = exprMat_filtered,
                        cellsAUC=selectRegulons(regulonAUC, GOI), 
                        plots = c("histogram", "binaryAUC", "AUC", "expression"),
                        exprCols = c("#440154FF", "#238A8DFF", "#FDE725FF"))
dev.off()

GOI <- regulons$Stat3#[1:length(regulons$Stat3)]
Cairo::CairoPDF("output/Step4_UMAP_Stat3_targets.pdf", height=20, width=20)
par(mfrow=c(4,4))
AUCell::AUCell_plotTSNE(dr_coords,
                        exprMat = exprMat_filtered,
                        cellsAUC=selectRegulons(regulonAUC, GOI), 
                        plots = c("histogram", "binaryAUC", "AUC", "expression"),
                        exprCols = c("#440154FF", "#238A8DFF", "#FDE725FF"))
dev.off()

GOI <- regulons$Klf6#[1:length(regulons$Klf6)]
Cairo::CairoPDF("output/Step4_seuratUMAP_Klf6_targets.pdf", height=20, width=20)
par(mfrow=c(4,4))
AUCell::AUCell_plotTSNE(dr_coords,
                        exprMat = exprMat_filtered,
                        cellsAUC=selectRegulons(regulonAUC, GOI), 
                        plots = c("histogram", "binaryAUC", "AUC", "expression"),
                        exprCols = c("#440154FF", "#238A8DFF", "#FDE725FF"))
dev.off()

GOI <- regulons$Stat5a#[1:length(regulons$Stat5a)]
Cairo::CairoPDF("output/Step4_seuratUMAP_Stat5a_targets.pdf", height=20, width=20)
par(mfrow=c(4,4))
AUCell::AUCell_plotTSNE(dr_coords,
                        exprMat = exprMat_filtered,
                        cellsAUC=selectRegulons(regulonAUC, GOI), 
                        plots = c("histogram", "binaryAUC", "AUC", "expression"),
                        exprCols = c("#440154FF", "#238A8DFF", "#FDE725FF"))
dev.off()

GOI <- regulons$Ddit3#[1:length(regulons$Ddit3)]
Cairo::CairoPDF("output/Step4_seuratUMAP_Ddit3_targets.pdf", height=20, width=20)
par(mfrow=c(4,4))
AUCell::AUCell_plotTSNE(dr_coords,
                        exprMat = exprMat_filtered,
                        cellsAUC=selectRegulons(regulonAUC, GOI), 
                        plots = c("histogram", "binaryAUC", "AUC", "expression"),
                        exprCols = c("#440154FF", "#238A8DFF", "#FDE725FF"))
dev.off()

