#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2024/1/16 15:21

# @Author : Fanjie
"""

import gc
import scanpy as sc
import os
import numpy as np
import pandas as pd
import scvi
import matplotlib.pyplot as plt

sc.set_figure_params(dpi_save=300, frameon=False, figsize=(8, 10))
dir = 'multi-disease_analysis'


all = sc.read_h5ad('multi-disease_analysis/alldisease_raw.h5ad')
all2 = all.copy()
sc.pp.highly_variable_genes(
        all,
        layer="counts",
        flavor="seurat_v3",
        n_top_genes=3000,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_bins=20,
        subset=True
    )
sc.pp.normalize_total(all, target_sum=1e4)
sc.pp.log1p(all)
all.raw = all
scvi.model.SCANVI.setup_anndata(all, layer="counts", labels_key='final_celltype', unlabeled_category="Unknown", batch_key='donor_ID')
vae = scvi.model.SCANVI(all)
vae.train(max_epochs=20, n_samples_per_label=100)
all.obsm["X_scANVI"] = vae.get_latent_representation()

sc.pp.highly_variable_genes(
        all2,
        layer="counts",
        flavor="seurat_v3",
        n_top_genes=3000,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_bins=20,
        subset=False
    )
sc.pp.normalize_total(all2, target_sum=1e4)
sc.pp.log1p(all2)
all2.raw = all2
all2.obsm["X_scANVI"] = all.obsm['X_scANVI']
sc.pp.neighbors(all2, use_rep="X_scANVI", key_added='scANVI', n_neighbors=30)
sc.tl.umap(all2, neighbors_key='scANVI', min_dist=0.3)
sc.tl.leiden(all2, neighbors_key='scANVI', resolution=0.1, key_added='cluster')
sc.tl.rank_genes_groups(all2, 'cluster', method='t-test')


marker_genes = {
    'Endothelial cell': ['PECAM1', 'CLDN5', 'VWF'],
    'Epithelial cell': ['EPCAM', 'KRT19', 'SNTN', 'SLPI'],
    'Fibroblast': ['DCN', 'LUM', 'COL1A2', 'COL1A1'],
    'Lymphocyte': ['PTPRC', 'CXCR4', 'NKG7'],
    'Myeloid cell': ['LYZ', 'SRGN', 'MARCO'],
}
sc.pl.dotplot(all2, marker_genes, groupby='first_anno')

## add cell annotation
first_anno = {
    "0": "Myeloid cell", "1": "Lymphocyte", "2": "Epithelial cell", "3": "Epithelial cell",
    "4": "Fibroblast", "5": "Endothelial cell", "6": "Epithelial cell", "7": "Fibroblast",
    "8": "Myeloid cell", "9": "Myeloid cell", "10": "Lymphocyte", "11": "Myeloid cell",
    "12": "Epithelial cell", "13": "Myeloid cell", "14": "Fibroblast"
}
all2.obs["first_anno"] = all2.obs.cluster.map(first_anno)

sc.pl.umap(
    all2, neighbors_key='scANVI', color=['first_anno'],
    ncols=1, title=['first_anno']
)
plt.savefig('fir_anno.pdf', dpi=300, format='pdf', bbox_inches='tight')



## subset mye cells
all = sc.read_h5ad('multi-disease_analysis/alldisease_raw.h5ad')
all.obs['first_anno'] = all2.obs.first_anno
mye = all[all.obs['first_anno'] == 'Myeloid cell']
mye.layers["newcounts"] = mye.X.copy()
mye.write_h5ad('mye_raw.h5ad')
mye2 = mye.copy()
sc.pp.highly_variable_genes(
        mye,
        layer="newcounts",
        flavor="seurat_v3",
        n_top_genes=3000,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_bins=20,
        subset=True
    )
sc.pp.normalize_total(mye, target_sum=1e4)
sc.pp.log1p(mye)
mye.raw = mye
scvi.model.SCANVI.setup_anndata(mye, layer="newcounts", labels_key='final_celltype', unlabeled_category="Unknown", batch_key='donor_ID')
vae = scvi.model.SCANVI(mye)
vae.train(max_epochs=20, n_samples_per_label=100)
mye.obsm["X_scANVI"] = vae.get_latent_representation()


sc.pp.highly_variable_genes(
        mye2,
        layer="newcounts",
        flavor="seurat_v3",
        n_top_genes=3000,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_bins=20,
        subset=False
    )
sc.pp.normalize_total(mye2, target_sum=1e4)
sc.pp.log1p(mye2)
mye2.raw = mye2
mye2.obsm["X_scANVI"] = vae.get_latent_representation()
sc.pp.neighbors(mye2, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(mye2, min_dist=0.3)
sc.tl.leiden(mye2, resolution=0.2, key_added='second_cluster')
sc.tl.rank_genes_groups(mye2, 'second_cluster', method='t-test')


second_anno = {
    "0": "Macrophage", "1": "Monocyte", "2": "Macrophage", "3": "Macrophage",
    "4": "Macrophage", "5": "Mast cell", "6": "Macrophage"
}
mye2.obs["second_anno"] = mye2.obs.second_cluster.map(second_anno)
marker_genes = {
    'Macrophage': ['CD163', 'MRC1', 'MARCO'],
    'Mast cell': ['GATA2', 'MS4A2', 'TPSAB1', 'TPSB2'],
    'Monocyte': ['FCN1', 'VCAN', 'SRGN'],
}
sc.pl.dotplot(mye2, marker_genes, groupby='second_anno')
plt.savefig(dir+'/sec_mye_dotplot.pdf', dpi=300, format='pdf', bbox_inches='tight')

sc.pl.umap(
    mye2, neighbors_key='scANVI', color=['donor_status', 'second_cluster', 'second_anno'],
    ncols=1, title=['donor_status', 'second_cluster', 'second_anno']
)
plt.savefig(os.path.join(dir, 'sec_anno_mye.pdf'), dpi=300, format='pdf', bbox_inches='tight')
mye2.write_h5ad('mye.h5ad')


## subset lym cells
all = sc.read_h5ad('multi-disease_analysis/alldisease_raw.h5ad')
lym = all[all.obs['first_anno'] == 'Lymphocyte']
lym.layers["counts"] = lym.X.copy()
lym.write_h5ad('lym_raw.h5ad')
lym2 = lym.copy()
sc.pp.highly_variable_genes(
        lym,
        layer="newcounts",
        flavor="seurat_v3",
        n_top_genes=3000,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_bins=20,
        subset=True
    )
sc.pp.normalize_total(lym, target_sum=1e4)
sc.pp.log1p(lym)
lym.raw = lym
scvi.model.SCANVI.setup_anndata(lym, layer="newcounts", labels_key='final_celltype', unlabeled_category="Unknown", batch_key='donor_ID')
vae = scvi.model.SCANVI(lym)
vae.train(max_epochs=20, n_samples_per_label=100)

sc.pp.highly_variable_genes(
        lym2,
        layer="newcounts",
        flavor="seurat_v3",
        n_top_genes=3000,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_bins=20,
        subset=False
    )
sc.pp.normalize_total(lym2, target_sum=1e4)
sc.pp.log1p(lym2)
lym2.raw = lym2
lym2.obsm["X_scANVI"] = vae.get_latent_representation()
sc.pp.neighbors(lym2, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(lym2, min_dist=0.3)
sc.tl.leiden(lym2, resolution=0.2, key_added='second_cluster')
sc.tl.rank_genes_groups(lym2, 'second_cluster', method='t-test')


second_anno = {
    "0": "B cell", "1": "NK cell", "2": "T cell", "5": "T cell", "7": "B cell",
    "8": "NK cell", "11": "B cell", "13": "B cell"
}
lym2.obs["second_anno"] = lym2.obs.second_cluster.map(second_anno)
marker_genes = {
    'B cell': ['CD74', 'CD79B', 'CD80'],
    'NK cell': ['KLRD1', 'NKG7', 'GNLY'],
    'T cell': ['CD3E', 'CD3G', 'CD3D']
}
sc.pl.dotplot(lym2, marker_genes, groupby='second_anno')
plt.savefig(dir+'/sec_lym_dotplot.pdf', dpi=300, format='pdf', bbox_inches='tight')

sc.pl.umap(
    lym2, neighbors_key='scANVI', color=['donor_status', 'second_cluster', 'second_anno'],
    ncols=1, title=['donor_status', 'second_cluster', 'second_anno']
)
plt.savefig(os.path.join(dir, 'sec_anno_lym.pdf'), dpi=300, format='pdf', bbox_inches='tight')
lym2.write_h5ad('lym.h5ad')




def build_cellsub(fathercell, celltype, cellsub):
    dir = 'multi-disease_analysis'
    ad = sc.read_h5ad(dir + fathercell + '_raw.h5ad')
    ad2 = sc.read_h5ad(dir + fathercell + '.h5ad')
    ad.obs['second_anno'] = ad2.obs.second_anno
    sub = ad[ad.obs['second_anno'] == celltype]
    sub.layers["newcounts"] = sub.X.copy()
    sub.write_h5ad(dir + cellsub + '_raw.h5ad')
    sc.pp.highly_variable_genes(
        sub,
        layer="newcounts",
        flavor="seurat_v3",
        n_top_genes=3000,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_bins=20,
        subset=True
    )
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    sub.raw = sub
    scvi.model.SCANVI.setup_anndata(sub, layer="newcounts",
                                    labels_key='final_celltype',
                                    unlabeled_category="Unknown",
                                    batch_key='donor_ID')
    vae = scvi.model.SCANVI(sub)
    vae.train(max_epochs=20, n_samples_per_label=100)
    sub.obsm["X_scANVI"] = vae.get_latent_representation()
    sub.write_h5ad(dir + cellsub + '_anvi.h5ad')
    sub.obs.to_csv(os.path.join(dir, cellsub + '_meta.csv'), sep=",", index=False)


build_cellsub('mye', 'Monocyte', 'mono')
build_cellsub('lym', 'B cell', 'bcell')






