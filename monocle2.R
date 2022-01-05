library(dplyr)
library(Seurat)
library(data.table)
library(ggplot2)
library(monocle)
set.seed(123)


#=========================================================================
#=========================================================================
#========= LOADING DATA, CREATING CDS & ADDING METADATA ==================
#=========================================================================
#=========================================================================

###LOAD CELLS
seurat_object_subset1 <- readRDS(file="/FILEPATH/seurat_object_subset__1.rds")

###create CDS object
## Building the necessary parts for a basic cds
# part one, gene annotations
gene_annotation <- as.data.frame(seurat_object_subset1@assays[["RNA"]]@counts@Dimnames[[1]], 
                                 row.names = seurat_object_subset1@assays[["RNA"]]@counts@Dimnames[[1]])
gene_annotation <- AnnotatedDataFrame(gene_annotation)
colnames(gene_annotation) <- "gene_short_name"

# part two, cell information
cell_metadata <- as.data.frame(seurat_object_subset1@assays[["RNA"]]@counts@Dimnames[[2]], 
                               row.names = seurat_object_subset1@assays[["RNA"]]@counts@Dimnames[[2]])
cell_metadata <- AnnotatedDataFrame(cell_metadata)
colnames(cell_metadata) <- "barcode"

# part three, counts sparse matrix (from fresh sample)
expression_matrix <- seurat_object_subset1@assays[["RNA"]]@counts

## Construct the basic cds object
cds <- newCellDataSet(as(expression_matrix, "sparseMatrix"),
                      phenoData = cell_metadata,
                      featureData = gene_annotation,
                      expressionFamily=negbinomial.size())



#=================================================
#=================================================
#================= PREPROCESSING =================
#=================================================
#=================================================
#Estimate size factors and dispersions
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

#Filtering low-quality cells
cds <- detectGenes(cds, min_expr = 0.1)



#==================================================================
#==================================================================
#================= Classifying and Counting Cells =================
#==================================================================
#==================================================================
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
cds <- reduceDimension(cds, max_components = 2, num_dim = 15,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 3)



#=========================================================================
#=========================================================================
#================= Constructing Single Cell Trajectories =================
#=========================================================================
#=========================================================================
#========= Trajectory step 1: choose genes that define a cell's progress
#we use the fully unsupervised approach and do not do feature selection here



#========= Trajectory step 2: reduce data dimensionality
cds <- reduceDimension(cds, 
                       max_components = 2,
                       reduction_method = 'DDRTree', 
                       verbose = T)



#========= Trajectory step 3: order cells along the trajectory
cds <- orderCells(cds, reverse = T, root_state=7)



#========= plot_cell_trajectory
#Cluster
plot_cell_trajectory(cds, 
                     color_by = "Cluster",
                     show_branch_points = F,
                     theta=180)

#State
plot_cell_trajectory(cds, 
                     color_by = "State",
                     show_branch_points = F,
                     theta=180)

#Pseudotime
library(viridis)
plot_cell_trajectory(cds, 
                     color_by = "Pseudotime",
                     show_branch_points = F,
                     theta=180) + scale_color_viridis(option="B")



#====================================================================
#====================================================================
#================= Differential Expression Analysis =================
#====================================================================
#====================================================================

#================= Step 1: Basic Differential Analysis =================
marker_genes <- row.names(fData(cds))



#================= Step 2: Finding Genes that Distinguish Cell Type or State =================
#DEG across State
diff_test_res_State <- differentialGeneTest(cds,
                                            fullModelFormulaStr = "~State")
diff_test_res_State[,c("gene_short_name", "pval", "qval")]
sig_gene_names_State <- row.names(subset(diff_test_res_State, qval < 0.1))
sig_gene_names_State0.01 <- row.names(subset(diff_test_res_State, qval < 0.01))



#===== plot_pseudotime_heatmap
plot_pseudotime_heatmap(cds[sig_gene_names_State0.01,],
                        num_clusters = 2,
                        cores = 8,
                        show_rownames = F)



#=========================================================================
#=========================================================================
#================= Analyzing Branches in Single-Cell Trajectories ========
#=========================================================================
#=========================================================================
BEAM_res <- BEAM(cds, branch_point = 2, cores = 8)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]



#===== plot_genes_branched_heatmap
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                 qval < 1e-4)),],
                            branch_point = 2,
                            num_clusters = 2,
                            cores = 8,
                            use_gene_short_name = T,
                            show_rownames = T)


