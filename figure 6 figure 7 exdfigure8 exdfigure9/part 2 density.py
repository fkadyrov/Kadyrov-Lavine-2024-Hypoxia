import scanpy as sc
import pandas as pd
import seaborn as sns



adata = sc.read_h5ad("myeloid.h5ad")
adata

sc.settings.set_figure_params(scanpy= True, dpi=300, dpi_save= 300, frameon=False, figsize=(3, 3.5), facecolor='white')

sc.tl.embedding_density(adata, groupby='condition', basis = 'umap')


sc.pl.embedding_density(adata, fg_dotsize=5, bg_dotsize=5, basis = 'umap', key='umap_density_condition')



sc.tl.embedding_density(adata, groupby='zsgr', basis = 'umap')


sc.pl.embedding_density(adata, fg_dotsize=5, bg_dotsize=5, basis = 'umap', key='umap_density_zsgr')