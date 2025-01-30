#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scanpy as sc
import numpy as np
import pandas as pd
import pickle
import os

import matplotlib.pyplot as plt
import scanpy as sc
import os

def run_raw(in_h5ad, out_h5ad, out_umap, out_unintegrate, top_genes):
    """
    :param in_h5ad: input adata
    :param out_h5ad: output adata
    :param out_umap: output umap.pdf
    :param out_unintegrate: output embed.pdf
    :param top_genes: number of HVGs
    :return:
    """
    print('Start unintegrate')
    sc.set_figure_params(dpi_save=300, frameon=False, figsize=(8, 10))
    adata = sc.read_h5ad(in_h5ad)
    adata.var_names_make_unique()
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=top_genes,
        layer="counts",
        batch_key='batch',
        subset=True
    )
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.embedding(
        adata, basis="umap", color=['batch', 'map_ann_level1',  'map_ann_level2'],
        ncols=1, title=['unintegrate_batch', 'unintegrate_ann_level1', 'unintegrate_ann_level2'],
        palette=None, cmap=None, frameon=False
    )
    plt.savefig(out_unintegrate, dpi=300, format='pdf', bbox_inches='tight')
    sc.pl.umap(
        adata, color=['batch', 'map_ann_level1', 'map_ann_level2'],
        ncols=1, title=['unintegrate_batch', 'unintegrate_ann_level1', 'unintegrate_ann_level2']
    )
    plt.savefig(out_umap, dpi=300, format='pdf', bbox_inches='tight')
    adata.obsm['X_umapraw'] = adata.obsm['X_umap']
    print("Save output")
    adata.write_h5ad(out_h5ad)
    print("Done unintegerated")
    return adata

dir = 'inte_bench'
outdir = 'inte_bench/inte_data'
figdir = 'inte_bench/fig'
in_h5ad = os.path.join(dir, 'core_raw.h5ad')
out_h5ad = os.path.join(outdir, 'unintegrate.h5ad')
out_umap = os.path.join(figdir, 'uninteg_umap.pdf')
out_integrate = os.path.join(figdir, 'uninteg_embed.pdf')

run_raw(in_h5ad, out_h5ad, out_umap, out_integrate, 2000)
