import palantir
import scanpy as sc
import numpy as np
import os
import tkinter
import seaborn as sns
sns.set_style("whitegrid")
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

# Plotting 
import matplotlib
import matplotlib.pyplot as plt

# Inline plotting
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')

# Reset random seed
import random
random.seed(10)

# Load data
import pandas as pd
sct_normalized_matrix = './myeloid_SCT_normalized.txt'
counts = pd.read_csv(sct_normalized_matrix, sep=',', index_col=0).transpose()
counts = counts.sort_index()
counts

meta_data = pd.read_csv('./myeloid_meta.csv', index_col=0)
meta_data = meta_data.sort_index()
meta_data

#pca projections
pca_projections, _ = palantir.utils.run_pca(sc.AnnData(counts), use_hvg=False)
pca_projections.shape
pca_projections

# Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(pca_projections)
dm_res

plt.scatter(np.arange(10), dm_res['EigenValues'])
plt.savefig('./eigenvalues')

pd.Series(np.abs(np.diff(dm_res['EigenValues']))).sort_values(ascending=False).index

ms_data = palantir.utils.determine_multiscale_space(dm_res, n_eigs = 3)
ms_data



#visualization
import harmony
fdl = harmony.plot.force_directed_layout(dm_res['kernel'], list(counts.index))
fig, ax = palantir.plot.plot_tsne(fdl)
plt.savefig('./fdl.png')


#magic imputation
imp_df = palantir.utils.run_magic_imputation(counts, dm_res)

palantir.plot.plot_cell_clusters(fdl, meta_data["condition"])
plt.savefig('./fdl_condition.png')

palantir.plot.plot_cell_clusters(fdl, meta_data["predicted.celltype"])
plt.savefig('./fdl_predicted.celltype.png')


meta_data_mono = meta_data[meta_data["predicted.celltype"] == "Mono"]
meta_data_mono

#running palantir
start_cell = "CCACTACGTACATCCA-1_1"
palantir.plot.highlight_cells_on_tsne(fdl, [start_cell])
plt.savefig('./fdl_start_cell.png')

pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=500, use_early_cell_as_start=True)


#identification of terminal states
pr_res.branch_probs.columns
lst = []
for i in list(pr_res.branch_probs.columns):
    print(meta_data["predicted.celltype"][i])
    lst.append(meta_data["predicted.celltype"][i])
    
pr_res.branch_probs.columns = lst
pr_res.branch_probs = pr_res.branch_probs.loc[:, lst]    
lst


#visualization of palantir results
palantir.plot.plot_palantir_results(pr_res, fdl)
plt.savefig('./fdl_palantir_summary.png')


#save meta data for r
df = pd.DataFrame()
df = pr_res.branch_probs.copy()
df["pseudotime"] = np.array(pr_res.pseudotime)
df["entropy"] = np.array(pr_res.entropy)
df["ClusterName"] = list(meta_data["predicted.celltype"])
df["Group"] = list(meta_data["condition"])
df["ClusterName"] = df["ClusterName"].astype("category")

pca_projections.to_csv('./pca_projections.csv', index=True)
ms_data.to_csv('./ms_data.csv', index=True)
fdl.to_csv('./fdl.csv', index=True)
df.to_csv('./palantir_meta_data.csv', index=True)
