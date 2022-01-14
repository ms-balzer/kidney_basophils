library (Seurat)
library (dplyr)
library (ggplot2)
library (harmony)
library (patchwork)

#================================Load Matrix===============================#
Control1 <- Read10X(data.dir = "/outs/filtered_feature_bc_matrix")
Control2 <- Read10X(data.dir = "/outs/filtered_feature_bc_matrix")
Control3 <- Read10X(data.dir = "/outs/filtered_feature_bc_matrix")
Control4 <- Read10X(data.dir = "/outs/filtered_feature_bc_matrix")
Control5 <- Read10X(data.dir = "/outs/filtered_feature_bc_matrix")
Control6 <- Read10X(data.dir = "/outs/filtered_feature_bc_matrix")
UUO1 <- Read10X(data.dir = "/outs/filtered_feature_bc_matrix")
UUO2 <- Read10X(data.dir = "/outs/filtered_feature_bc_matrix")

#================================Create Seurat Object======================================================#
Control1.S <- CreateSeuratObject(counts = Control1, project = "Control1", min.cells = 3, min.features = 300)
Control2.S <- CreateSeuratObject(counts = Control2, project = "Control2", min.cells = 3, min.features = 300)
Control3.S <- CreateSeuratObject(counts = Control3, project = "Control3", min.cells = 3, min.features = 300)
Control4.S <- CreateSeuratObject(counts = Control4, project = "Control4", min.cells = 3, min.features = 300)
Control5.S <- CreateSeuratObject(counts = Control5, project = "Control5", min.cells = 3, min.features = 300)
Control6.S <- CreateSeuratObject(counts = Control6, project = "Control6", min.cells = 3, min.features = 300)
UUO1.S <- CreateSeuratObject(counts = UUO1, project = "UUO1", min.cells = 3, min.features = 300)
UUO2.S <- CreateSeuratObject(counts = UUO2, project = "UUO2", min.cells = 3, min.features = 300)

#============================Merge Datssets======================================================#
UUO.Control.Project = merge (Control1.S, y = c (Control2.S, Control3.S, Control4.S, Control5.S, Control6.S, UUO1.S, UUO2.S))

#============================Add percent.mt to dataset======================================================#
UUO.Control.Project[["percent.mt"]] <- PercentageFeatureSet(UUO.Control.Project, pattern = "^Mt-")

#============================Pre Processing======================================================#
UUO.Control.Project <- subset(UUO.Control.Project, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 50)

UUO.Control.Project <- FindVariableFeatures(object = UUO.Control.Project,selection.method = "vst")
UUO.Control.Project <- NormalizeData(UUO.Control.Project)
UUO.Control.Project <- ScaleData(UUO.Control.Project, vars.to.regress = c ("nCount_RNA", "percent.mt"))
UUO.Control.Project <- RunPCA(UUO.Control.Project, assay = 'RNA', npcs = 30, features = VariableFeatures(object = UUO.Control.Project), 
          verbose = TRUE, ndims.print = 1:5, nfeatures.print = 10)
		  
#============= Run Harmony =============
UUO.Control.Project <- RunHarmony(UUO.Control.Project, group.by.vars = "orig.ident")
UUO.Control.Project <- RunUMAP(UUO.Control.Project, reduction = "harmony", dims = 1:30)
UUO.Control.Project <- FindNeighbors(UUO.Control.Project, reduction = "harmony", dims = 1:30) %>% FindClusters()				   
DimPlot(UUO.Control.Project, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#================================== DEG List==============================#
UUO.Control.Project.Markers <- FindAllMarkers(UUO.Control.Project, only.pos = TRUE, logfc.threshold = 0.25)

#========================== Sub-Clustering after remove PT========================#
UUO.Control.Project.Subset <- subset(UUO.Control.Project, idents = "PT", invert = TRUE)

UUO.Control.Project.Subset <- FindVariableFeatures(object = UUO.Control.Project.Subset,selection.method = "vst")
UUO.Control.Project.Subset <- NormalizeData(UUO.Control.Project.Subset)
UUO.Control.Project.Subset <- ScaleData(UUO.Control.Project.Subset, vars.to.regress = c ("nCount_RNA", "percent.mt"))
UUO.Control.Project.Subset <- RunPCA(UUO.Control.Project.Subset, assay = 'RNA', npcs = 30, features = VariableFeatures(object = UUO.Control.Project.Subset), 
          verbose = TRUE, ndims.print = 1:5, nfeatures.print = 10)
		  
#===================== Run Harmony & Clustering ===============================#
UUO.Control.Project.Subset <- RunHarmony(UUO.Control.Project.Subset, group.by.vars = "orig.ident")
UUO.Control.Project.Subset <- RunUMAP(UUO.Control.Project.Subset, reduction = "harmony", dims = 1:30)
UUO.Control.Project.Subset <- FindNeighbors(UUO.Control.Project.Subset, reduction = "harmony", dims = 1:30) %>% FindClusters()				   
DimPlot(UUO.Control.Project.Subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#================================== DEG List=====================================#
UUO.Control.Project.Subset.Markers <- FindAllMarkers(UUO.Control.Project.Subset, only.pos = TRUE, logfc.threshold = 0.25, test.use = "MAST")


#================================== Human Single Cell Data Processing and Clustering=====================================#
Control <- Read10X(data.dir = "/outs/filtered_feature_bc_matrix")
CKD <- Read10X(data.dir = "/outs/filtered_feature_bc_matrix")
Control.S <- CreateSeuratObject(counts = Control, project = "Control", min.cells = 3, min.features = 300)
CKD.S <- CreateSeuratObject(counts = CKD, project = "CKD", min.cells = 3, min.features = 300)
Human.Control.CKD = merge (Control.S, y = c (CKD.S))
Human.Control.CKD[["percent.mt"]] <- PercentageFeatureSet(Human.Control.CKD, pattern = "^MT-")
Human.Control.CKD <- FindVariableFeatures(object = Human.Control.CKD,selection.method = "vst")
Human.Control.CKD <- NormalizeData(Human.Control.CKD)
Human.Control.CKD <- ScaleData(Human.Control.CKD, vars.to.regress = c ("nCount_RNA", "percent.mt"))
Human.Control.CKD <- RunPCA(Human.Control.CKD, assay = 'RNA', npcs = 30, features = VariableFeatures(object = Human.Control.CKD), 
          verbose = TRUE, ndims.print = 1:5, nfeatures.print = 10)
Human.Control.CKD <- RunHarmony(Human.Control.CKD, group.by.vars = "orig.ident")
Human.Control.CKD <- RunUMAP(Human.Control.CKD, reduction = "harmony", dims = 1:30)
Human.Control.CKD <- FindNeighbors(Human.Control.CKD, reduction = "harmony", dims = 1:30) %>% FindClusters()				   
DimPlot(Human.Control.CKD, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
Human.Control.CKD.Markers <- FindAllMarkers(Human.Control.CKD, only.pos = TRUE, logfc.threshold = 0.25, test.use = "MAST")

#================================== Sub-Clusterin PT Cells in UUO Samples=====================================#
UUO.Project.PT <- subset(UUO.Control.Project, idents = "PT")
UUO.PT = subset(UUO.Project.PT, subset = Groups == "UUO")
UUO.PT <- FindVariableFeatures(object = UUO.PT,selection.method = "vst")
UUO.PT <- NormalizeData(UUO.PT)
UUO.PT <- ScaleData(UUO.PT, vars.to.regress = c ("nCount_RNA", "percent.mt"))
UUO.PT <- RunPCA(UUO.PT, assay = 'RNA', npcs = 30, features = VariableFeatures(object = UUO.PT), 
          verbose = TRUE, ndims.print = 1:5, nfeatures.print = 10)
		  
#===================== Run Harmony & Clustering on PT Subset ===============================#
UUO.PT <- RunHarmony(UUO.PT, group.by.vars = "orig.ident")
UUO.PT <- RunUMAP(UUO.PT, reduction = "harmony", dims = 1:30)
UUO.PT <- FindNeighbors(UUO.PT, reduction = "harmony", dims = 1:30) %>% FindClusters()				   
DimPlot(UUO.PT, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#================================== DEG List of PT Subset=====================================#
UUO.PT.Markers <- FindAllMarkers(UUO.PT, only.pos = TRUE, logfc.threshold = 0.25, test.use = "MAST")
