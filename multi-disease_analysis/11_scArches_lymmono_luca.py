#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2024/9/2 0:01

# @Author : Fanjie
"""
import scanpy as sc
import os
import pandas as pd
import torch
import scarches as sca
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp


sc.set_figure_params(dpi_save=300, frameon=False, figsize=(8, 10))
torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)
sc.logging.print_header()
dir = 'multi-disease_analysis'
os.chdir(dir)


lymmono = sc.read_h5ad('scarches/lymmono.h5ad')
ref_luca_sub =  sc.read_h5ad('scarches/ref_luca_sub_latest.h5ad')
lymmono_ref_luca = ref_luca_sub[~ref_luca_sub.obs_names.isin(lymmono.obs.cell_ID), :]


## train reference model for ref_luca
condition_key = 'dataset'
cell_type_key = 'cell_type'
vae_epochs = 500
surgery_epochs = 500

lymmono_ref_luca.layers["counts"] = lymmono_ref_luca.X.copy()
sc.pp.highly_variable_genes(
        lymmono_ref_luca,
        layer="counts",
        flavor="seurat_v3",
        n_top_genes=3000,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_bins=20,
        subset=True
    )

sca.models.SCVI.setup_anndata(lymmono_ref_luca, batch_key=condition_key)
vae = sca.models.SCVI(lymmono_ref_luca)
vae.train(max_epochs=vae_epochs)

reference_latent = sc.AnnData(vae.get_latent_representation())
reference_latent.obs["cell_type"] = lymmono_ref_luca.obs[cell_type_key].tolist()
reference_latent.obs["dataset"] = lymmono_ref_luca.obs[condition_key].tolist()

sc.pp.neighbors(reference_latent, n_neighbors=8)
sc.tl.leiden(reference_latent)
sc.tl.umap(reference_latent)
vae.save('scarches/ref_model/lymmono_ref_luca', overwrite=True)


## Perform surgery on reference model and train on query dataset without cell type labels
adata_to_map = lymmono.copy()

# the gene id of model is ensemble id, so transfer the gene id of target data
gene_mapping = pd.read_csv('scarches/total_gene_list_43878.txt', sep='\t', header=None, names=['GeneSymbol', 'LocusGroup' ,'GeneID'])
gene_symbol_to_id = dict(zip(gene_mapping['GeneSymbol'], gene_mapping['GeneID']))
adata_to_map.var['GeneSymbol'] = adata_to_map.var.index
adata_to_map.var['GeneID'] = adata_to_map.var_names.map(gene_symbol_to_id)

# the gene of target data and query data need to be the same (element and sort)
# contact the adata for same gene and the new adata for missing gene
genes_to_extract = lymmono_ref_luca.var_names
common_genes = adata_to_map.var_names.intersection(genes_to_extract)
adata_to_map_filter = adata_to_map[:, adata_to_map.var['GeneID'].isin(common_genes.tolist())]
adata_to_map_filter.var_names = adata_to_map_filter.var['GeneID']
common_adata = sc.AnnData(
    X=adata_to_map_filter.X,
    obs=adata_to_map_filter.obs.iloc[:, :1],
    var=adata_to_map_filter.var_names.tolist(),
)
common_adata.layers["counts"] = common_adata.X
common_adata.var_names = adata_to_map_filter.var_names

missing_genes = lymmono_ref_luca.var_names.difference(common_genes)
missing_adata = sc.AnnData(
    X=sp.sparse.csr_matrix(np.zeros(shape=(adata_to_map.n_obs, len(missing_genes))), dtype="float32"),
    obs=adata_to_map.obs.iloc[:, :1],
    var=missing_genes.tolist(),
)
missing_adata.layers["counts"] = missing_adata.X
missing_adata.var_names = missing_genes
adata_to_map_augmented = sc.concat([common_adata, missing_adata], join='outer', axis=1, index_unique=None, merge="unique")

# gene sort by reference dataset (lymmono_ref_luca)
adata_to_map_augmented = adata_to_map_augmented[:, lymmono_ref_luca.var.index.tolist()].copy()
(adata_to_map_augmented.var.index == lymmono_ref_luca.var.index).all()
adata_to_map_augmented.obs['dataset'] = adata_to_map.obs['donor_status']

# train model for target data
model = sca.models.SCVI.load_query_data(
    adata_to_map_augmented,
    'scarches/ref_model/lymmono_ref_luca',
    freeze_dropout = True,
)
model.train(max_epochs=surgery_epochs)
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = lymmono.obs['mono_anno'].tolist()

sc.pp.neighbors(query_latent)
sc.tl.leiden(query_latent)
sc.tl.umap(query_latent)
model.save('scarches/query_model/lymmono_ref_luca', overwrite=True)


## Get latent representation of reference + query dataset and compute UMAP
lymmono_ref_luca.obs['ref_query'] = 'refdata'
adata_to_map_augmented.obs['ref_query'] = 'querydata'
adata_to_map_augmented.var.columns=['GeneID']
adata_full = lymmono_ref_luca.concatenate(adata_to_map_augmented)
full_latent = sc.AnnData(model.get_latent_representation(adata=adata_full))
full_latent.obs['cell_type'] = adata_full.obs[cell_type_key].tolist()
full_latent.obs['ref_query'] = adata_full.obs['ref_query'].tolist()
full_latent.obs['cell_type_major'] = adata_full.obs['cell_type_major'].tolist()
full_latent.obs['ann_fine'] = adata_full.obs['ann_fine'].tolist()

sc.pp.neighbors(full_latent)
sc.tl.leiden(full_latent)
sc.tl.umap(full_latent)
sc.pl.umap(full_latent,
           color=['cell_type','ref_query','cell_type_major','ann_fine'],
           frameon=False,
           wspace=0.6,
           ncols=1
           )
plt.savefig('scarches/lymmono_luca.pdf', dpi=300, format='pdf', bbox_inches='tight')
full_latent.write_h5ad('scarches/lymmono_luca_full.h5ad')