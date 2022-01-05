import os
import scvelo as scv
import matplotlib.pyplot as plt
import tqdm as tqdm
import ipywidgets as ipywidgets
import numpy as np
import pandas as pd
import scanpy as sc
import pyreadr
scv.set_figure_params('scvelo')



#================== Read in data ==================
adata = scv.read("FILEPATH/YOUR.loom", cache=True)
result = pyreadr.read_r("FILEPATH/barcodes_seurat_object_subset.rds")
df = result[None]
barcodes = df[None].tolist()
adata_sub = adata[adata.obs.index.isin(barcodes)]
seurat_h5 = sc.read_h5ad("/home/balzermi/data/velocyto/UUOnew_subset_UUO1only/seurat_object_subset.h5ad")
adata_sub.obs['cluster'] = seurat_h5.obs['cluster']

scv.pl.proportions(adata_sub)



#================== Basic preprocessing ==================
scv.pp.filter_and_normalize(adata_sub)
scv.pp.moments(adata_sub)



#================== compute UMAP, diffmap ==================
scv.tl.umap(adata_sub)
scv.tl.louvain(adata_sub)
sc.tl.diffmap(adata_sub)



#================== Velocity Tools ==================
scv.tl.velocity(adata_sub, mode='stochastic')
scv.tl.velocity_graph(adata_sub)

scv.pl.velocity_graph(adata_sub, basis='umap')



#================== Project the velocities ==================
palette = ["#FF7F00", "#984EA3", "#377EB8", "#4DAF4A", "#E41A1C"]
scv.pl.velocity_embedding_stream(adata_sub, 
                                 basis='X_diffmap', 
                                 color='cluster', 
                                 palette=palette, 
                                 legend_loc='none', 
                                 arrow_size=1.5, 
                                 arrow_style='-|>', 
                                 linewidth=1.5, 
                                 title='',
                                 size=175)


