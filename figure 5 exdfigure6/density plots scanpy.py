import scanpy as sc
import pandas as pd
import seaborn as sns


adata = sc.read_h5ad("myeloid.h5ad")
adata

sc.settings.set_figure_params(scanpy = True, dpi=80, dpi_save = 300, frameon=False, figsize=(10, 10), facecolor='white')


sc.tl.embedding_density(adata, groupby='condition', basis='ref.umap')


sc.pl.embedding_density(adata, basis = 'ref.umap', key = "ref.umap_density_condition")
