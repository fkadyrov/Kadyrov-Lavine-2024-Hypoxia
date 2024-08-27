import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import seaborn as sns
import scipy.stats
import anndata
import os
import igraph as ig
import louvain

adata = sc.read_h5ad("./merged_postQC.h5ad")
adata

def bh(pvalues):
    '''
    Computes the Benjamini-Hochberg FDR correction.

    Input:
        * pvals - vector of p-values to correct
    '''
    n = int(pvalues.shape[0])
    new_pvalues = np.empty(n)
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in range(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues


####### main
meta = adata.obs
sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)

scorenames = ['scrublet_score','scrublet_cluster_score','bh_pval']
os.makedirs('scrublet-scores')

###

for condition in meta['condition'].unique():
    print(condition)
    #import data
    adata_condition = adata[adata.obs['condition'] == condition]
    
    #set up and run Scrublet
    scrub = scr.Scrublet(adata_condition.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    adata_condition.obs['scrublet_score'] = doublet_scores
    #overcluster prep. run turbo basic scanpy pipeline
    sc.pp.filter_genes(adata_condition, min_cells=3)
    sc.pp.normalize_per_cell(adata_condition, counts_per_cell_after=1e4)
    sc.pp.log1p(adata_condition)
    sc.pp.highly_variable_genes(adata_condition, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_condition = adata_condition[:, adata_condition.var['highly_variable']]
    sc.pp.scale(adata_condition, max_value=10)
    sc.tl.pca(adata_condition, svd_solver='arpack')
    sc.pp.neighbors(adata_condition)
    #eoverclustering proper - do basic clustering first, then cluster each cluster
    sc.tl.louvain(adata_condition)
    for clus in np.unique(adata_condition.obs['louvain']):
        sc.tl.louvain(adata_condition, restrict_to=('louvain',[clus]))
        adata_condition.obs['louvain'] = adata_condition.obs['louvain_R']
    #compute the cluster scores - the median of Scrublet scores per overclustered cluster
    for clus in np.unique(adata_condition.obs['louvain']):
        adata_condition.obs.loc[adata_condition.obs['louvain']==clus, 'scrublet_cluster_score'] = \
            np.median(adata_condition.obs.loc[adata_condition.obs['louvain']==clus, 'scrublet_score'])
    #now compute doublet p-values. figure out the median and mad (from above-median values) for the distribution
    med = np.median(adata_condition.obs['scrublet_cluster_score'])
    mask = adata_condition.obs['scrublet_cluster_score']>med
    mad = np.median(adata_condition.obs['scrublet_cluster_score'][mask]-med)
    #let's do a one-sided test. the Bertie write-up does not address this but it makes sense
    pvals = 1-scipy.stats.norm.cdf(adata_condition.obs['scrublet_cluster_score'], loc=med, scale=1.4826*mad)
    adata_condition.obs['bh_pval'] = bh(pvals)
    #create results data frame for single sample and copy stuff over from the adata object
    scrublet_condition = pd.DataFrame(0, index=adata_condition.obs_names, columns=scorenames)
    for meta in scorenames:
        scrublet_condition[meta] = adata_condition.obs[meta]
    #write out complete sample scores
    scrublet_condition.to_csv('scrublet-scores/'+condition+'.csv')
    
    adata.obs['condition'].unique()
    
s1 = pd.read_csv("./scrublet-scores/WT1.csv", index_col=0)
s2 = pd.read_csv("./scrublet-scores/WT2.csv", index_col=0)
s3 = pd.read_csv("./scrublet-scores/KO1.csv", index_col=0)
s4 = pd.read_csv("./scrublet-scores/KO2.csv", index_col=0)

df_list = [s1,s2,s3,s4]
scrub_all = pd.concat(df_list)

scrub_all.to_csv('scrublet-scores/all.csv')