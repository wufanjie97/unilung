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


tbcell = sc.read_h5ad('scarches/tbcell.h5ad')
ref_panb_sub =  sc.read_h5ad('scarches/ref_panb_sub_latest.h5ad')
tbcell_ref_panb = ref_panb_sub[~ref_panb_sub.obs_names.isin(tbcell.obs.cell_ID), :]


## train reference model for ref_panb
condition_key = 'cancer'
cell_type_key = 'celltype'
vae_epochs = 500
surgery_epochs = 500

tbcell_ref_panb.layers["counts"] = tbcell_ref_panb.X.copy()
sc.pp.highly_variable_genes(
        tbcell_ref_panb,
        layer="counts",
        flavor="seurat_v3",
        n_top_genes=3000,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_bins=20,
        subset=True
    )

sca.models.SCVI.setup_anndata(tbcell_ref_panb, batch_key=condition_key)
vae = sca.models.SCVI(tbcell_ref_panb)
vae.train(max_epochs=vae_epochs)

reference_latent = sc.AnnData(vae.get_latent_representation())
reference_latent.obs["cell_type"] = tbcell_ref_panb.obs[cell_type_key].tolist()
reference_latent.obs["cancer"] = tbcell_ref_panb.obs[condition_key].tolist()

sc.pp.neighbors(reference_latent, n_neighbors=8)
sc.tl.leiden(reference_latent)
sc.tl.umap(reference_latent)
vae.save('scarches/ref_model/tbcell_ref_panb', overwrite=True)


## Perform surgery on reference model and train on query dataset without cell type labels
adata_to_map = tbcell.copy()

# the gene of target data and query data need to be the same (element and sort)
# contact the adata for same gene and the new adata for missing gene
adata_to_map.var['GeneSymbol'] = adata_to_map.var_names
genes_to_extract = tbcell_ref_panb.var_names
common_genes = adata_to_map.var_names.intersection(genes_to_extract)
adata_to_map_filter = adata_to_map[:, adata_to_map.var['GeneSymbol'].isin(common_genes.tolist())]
common_adata = sc.AnnData(
    X=adata_to_map_filter.X,
    obs=adata_to_map_filter.obs.iloc[:, :1],
    var=adata_to_map_filter.var_names.tolist(),
)
common_adata.layers["counts"] = common_adata.X
common_adata.var_names = adata_to_map_filter.var_names

missing_genes = tbcell_ref_panb.var_names.difference(common_genes)
missing_adata = sc.AnnData(
    X=sp.sparse.csr_matrix(np.zeros(shape=(adata_to_map.n_obs, len(missing_genes))), dtype="float32"),
    obs=adata_to_map.obs.iloc[:, :1],
    var=missing_genes.tolist(),
)
missing_adata.layers["counts"] = missing_adata.X
missing_adata.var_names = missing_genes
adata_to_map_augmented = sc.concat([common_adata, missing_adata], join='outer', axis=1, index_unique=None, merge="unique")

# gene sort by reference dataset (tbcell_ref_panb)
adata_to_map_augmented = adata_to_map_augmented[:, tbcell_ref_panb.var.index.tolist()].copy()
(adata_to_map_augmented.var.index == tbcell_ref_panb.var.index).all()
adata_to_map_augmented.obs['cancer'] = adata_to_map.obs['donor_status']

# train model for target data
model = sca.models.SCVI.load_query_data(
    adata_to_map_augmented,
    'scarches/ref_model/tbcell_ref_panb',
    freeze_dropout = True,
)
model.train(max_epochs=surgery_epochs)
query_latent = sc.AnnData(model.get_latent_representation())
query_latent.obs['cell_type'] = tbcell.obs['bcell_anno'].tolist()

sc.pp.neighbors(query_latent)
sc.tl.leiden(query_latent)
sc.tl.umap(query_latent)
model.save('scarches/query_model/tbcell_ref_panb', overwrite=True)


## Get latent representation of reference + query dataset and compute UMAP
tbcell_ref_panb.obs['ref_query'] = 'refdata'
adata_to_map_augmented.obs['ref_query'] = 'querydata'
adata_to_map_augmented.var.columns=['GeneSymbol']
adata_full = tbcell_ref_panb.concatenate(adata_to_map_augmented)
full_latent = sc.AnnData(model.get_latent_representation(adata=adata_full))
full_latent.obs['cell_type'] = adata_full.obs[cell_type_key].tolist()
full_latent.obs['ref_query'] = adata_full.obs['ref_query'].tolist()
full_latent.obs['cancer'] = adata_full.obs['cancer'].tolist()
full_latent.obs['Phase'] = adata_full.obs['Phase'].tolist()


sc.pp.neighbors(full_latent)
sc.tl.leiden(full_latent)
sc.tl.umap(full_latent)
sc.pl.umap(full_latent,
           color=['cell_type','ref_query','cancer','Phase'],
           frameon=False,
           wspace=0.6,
           ncols=1
           )
plt.savefig('scarches/tbcell_panb.pdf', dpi=300, format='pdf', bbox_inches='tight')
full_latent.write_h5ad('scarches/tbcell_panb_full.h5ad')