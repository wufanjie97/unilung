#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3
dir = 'lung_cancer_sp'
os.chdir(dir)

ref =  sc.read_h5ad('cell2loc/ref_luca_sub_latest.h5ad')
ref_lusc = ref[ref.obs.disease == "squamous cell lung carcinoma", :]

lusc_raw = sc.read_h5ad('lusc_latest_all.h5ad')
adata_vis = lusc_raw.copy()


adata_vis.var['SYMBOL'] = adata_vis.var_names
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]



# Read data
adata_ref = ref_lusc.copy()
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()

# check the float value, transfer to int
if np.any(adata_ref.X.data % 1 != 0):
    adata_ref.X = np.rint(adata_ref.X).astype(int)

# Estimation of reference cell type signatures (NB regression)
# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='dataset',
                        # cell type, covariate used for constructing signatures
                        labels_key='cell_type_major',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['platform']
                       )
# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()
mod.train(max_epochs=200, use_gpu=False)

# the ELBO loss history plot, have a decreasing trend and level off by the end of training. If it is still decreasing, increase .max_epochs
mod.plot_history(20)
plt.show()
# export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
)

# Save model and adata_ref anndata object
mod.save('cell2loc/ref_sig_lusc', overwrite=True)
adata_ref.write('cell2loc/map_lusc/sc.h5ad')






# load adata_ref and model
adata_ref = sc.read_h5ad('cell2loc/map_lusc/sc.h5ad')
mod = cell2location.models.RegressionModel.load('cell2loc/ref_sig_lusc', adata_ref)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]


# Cell2location: spatial mapping
# find shared genes and subset both anndata and reference signatures
# the var_name of adata_vis is symbol, match with ensemble in adata_ref
gene_mapping = pd.read_csv('total_gene_list_43878.txt', sep='\t', header=None, names=['GeneSymbol', 'LocusGroup' ,'GeneID'])
gene_symbol_to_id = dict(zip(gene_mapping['GeneSymbol'], gene_mapping['GeneID']))
adata_vis.var['Symbol'] = adata_vis.var.index
adata_vis.var['GeneID'] = adata_vis.var['Symbol'].map(gene_symbol_to_id)
adata_vis.var['GeneID'] = adata_vis.var['GeneID'].fillna(adata_vis.var['Symbol'])
adata_vis.var_names = adata_vis.var['GeneID']

intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="patient")
# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    N_cells_per_location=30,
    detection_alpha=20
)
mod.view_anndata_setup()
# Train
mod.train(max_epochs=20000,
          batch_size=None,
          train_size=1,
          use_gpu=False,
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training'])
plt.show()

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis2 = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
)

# Save model and adata_vis
mod.save("cell2loc/map_lusc", overwrite=True)
adata_vis2.write("cell2loc/map_lusc/sp.h5ad")



# Visualising cell abundance in spatial coordinates
# load adata_vis and model
adata_vis = sc.read_h5ad("cell2loc/map_lusc/sp.h5ad")
mod = cell2location.models.Cell2location.load("cell2loc/map_lusc", adata_vis)

adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
)
# Compute expected expression per cell type
expected_dict = mod.module.model.compute_expected_per_cell_type(
    mod.samples["post_sample_q05"], mod.adata_manager
)

# Add to anndata layers
for i, n in enumerate(mod.factor_names_):
    adata_vis.layers[n] = expected_dict['mu'][i]

adata_vis.write("cell2loc/map_lusc/sp_gene_expr.h5ad")


from cell2location import run_colocation
res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='CoLocatedGroupsSklearnNMF',
    train_args={
      'n_fact': np.arange(3, 10), # IMPORTANT: use a wider range of the number of factors (5-30)
      'sample_name_col': 'sample', # columns in adata_vis.obs that identifies sample
      'n_restarts': 3 # number of training restarts
    },
    export_args={'path': 'cell2loc/nmf_lusc/'}
)


