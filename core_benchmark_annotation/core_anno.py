#!/usr/bin/env python
# -*- coding: utf-8 -*-
import gc
import matplotlib.pyplot as plt
from scvi.model.utils import mde
import pickle
import scanpy as sc
import scvi
import os
import numpy as np
import pandas as pd

sc.set_figure_params(dpi_save=300, frameon=False, figsize=(8, 10))
sc.logging.print_header()
dir = 'core_anno'

adata = sc.read_h5ad('raw_r2.h5ad')
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=3000,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    n_bins=20,
    batch_key='batch',
)
adata.raw = adata
adata = adata[:, adata.var['highly_variable']==True]

# scANVI
adata = adata.copy()
scvi.model.SCANVI.setup_anndata(adata, layer="counts", labels_key='map_celltype', unlabeled_category="Unknown", batch_key='batch')
vae = scvi.model.SCANVI(adata)
vae.train(max_epochs=20, n_samples_per_label=100)

adata.obsm["X_scANVI"] = vae.get_latent_representation()
adata.obsm["X_mde_scANVI"] = mde(adata.obsm["X_scANVI"])
sc.pl.embedding(
    adata, basis="X_mde_scANVI", color=['batch', 'map_ann_level1',  'map_ann_level2'],
    ncols=1, title=['scANVI_batch', 'scANVI_ann_level1', 'scANVI_ann_level2'],
    palette=None, cmap=None, frameon=False
)

## first annotation
sc.pp.neighbors(adata, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=0.1, key_added='first_cluster')
sc.tl.rank_genes_groups(adata, 'first_cluster', method='t-test')

## add cell annotation
first_anno = {
    "0": "Myeloid cell", "1": "Lymphocyte", "2": "Epithelial cell", "3": "Epithelial cell",
    "4": "Fibroblast", "5": "Epithelial cell", "6": "Endothelial cell", "7": "Lymphocyte",
    "8": "Epithelial cell", "9": "Lymphocyte", "10": "Myeloid cell", "11": "Lymphocyte",
    "12": "Epithelial cell", "13": "Fibroblast", "14": "Epithelial cell", "15": "Epithelial cell"
}
adata_tmp.obs["unilung_ann_level1"] = adata_tmp.obs.first_cluster.map(first_anno)

sc.pl.umap(
    adata, color=['first_cluster', 'unilung_ann_level1', 'map_ann_level1'],
    ncols=1, title=['first_cluster', 'unilung_ann_level1', 'map_ann_level1']
)

marker = {
    'Endothelial cell': ['PECAM1', 'CLDN5', 'VWF' ,'CAV1'],
    'Epithelial cell': ['EPCAM', 'KRT19', 'CAPS', 'KRT7'],
    'Fibroblast': ['DCN', 'LUM', 'COL1A2', 'COL1A1'],
    'Lymphocyte': ['PTPRC', 'TRAC', 'CD3D', 'CXCR4'],
    'Myeloid cell': ['LYZ', 'CD14', 'MARCO', 'SRGN']
}

sc.pl.dotplot(adata, marker, groupby='first_annotation')




## second annotation
## endo
endo = adata[adata.obs.unilung_ann_level1 == "Endothelial cell", :]
sc.pp.filter_genes(endo, min_cells=3)
sc.pp.filter_cells(endo, min_genes=200)
sc.pp.neighbors(endo, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(endo, min_dist=0.3)
sc.tl.leiden(endo, resolution=0.2, key_added='second_cluster')
sc.tl.rank_genes_groups(endo, 'second_cluster', method='t-test')

## add cell annotation
second_anno = {
    "0": "Vein endothelial cell", "1": "Lymphatic endothelial cell", "2": "Vein endothelial cell",
    "3": "Capillary endothelial cell", "4": "Capillary endothelial cell", "5": "Capillary endothelial cell",
    "7": "Vein endothelial cell",
}

endo.obs["unilung_ann_level2"] = endo.obs.second_cluster.map(second_anno)
sc.pl.umap(
    endo, color=['second_cluster', 'unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['second_cluster', 'unilung_ann_level2', 'map_ann_level2']
)

marker = {
    'Artery endothelial cell': ['DKK2'],
    'Capillary endothelial cell': ['CA4', 'RAMP3'],
    'Lymphatic endothelial cell': ['PROX1', 'PDPN', 'RELN'],
    'Vein endothelial cell': ['ACKR1', 'EMCN', 'BACE2']
}
sc.pl.dotplot(endo, marker, groupby='unilung_ann_level2')




## epi
epi = adata[adata.obs.unilung_ann_level1 == "Epithelial cell", :]
sc.pp.filter_genes(epi, min_cells=3)
sc.pp.filter_cells(epi, min_genes=200)
sc.pp.neighbors(epi, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(epi, min_dist=0.3)
sc.tl.leiden(epi, resolution=0.2, key_added='second_cluster')
sc.tl.rank_genes_groups(epi, 'second_cluster', method='t-test')

## add cell annotation
second_anno = {
    "0": "Ciliated cell", "1": "Alveolar cell", "2": "Secretory cell", "3": "Basal cell",
    "4": "Alveolar cell", "5": "Basal cell", "6": "Alveolar cell", "7": "Alveolar cell",
    "8": "Alveolar cell", "9": "Submucosal gland cell", "10": "Basal cell", "11": "Ciliated cell"
}

epi.obs["unilung_ann_level2"] = epi.obs.second_cluster.map(second_anno)
sc.pl.umap(
    epi, color=['second_cluster', 'unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['second_cluster', 'unilung_ann_level2', 'map_ann_level2']
)

marker = {
    'Alveolar cell': ['SFTPC', 'CAV1', 'SFTPB'],
    'Basal cell': ['KRT5', 'KRT15', 'S100A2', 'S100A2'],
    'Ciliated cell': ['FOXJ1', 'PIFO', 'TPPP3'],
    'Secretory cell': ['SCGB1A1', 'SCGB3A1', 'BPIFB1', 'SERPINB3'],
    'Squamous cell': ['SCEL', 'SPRR1B'],
    'Submucosal gland cell': ['LTF', 'SLPI', 'LYZ', 'MUC5B']
}

sc.pl.dotplot(epi, marker, groupby='unilung_ann_level2')




## fibro
fib = adata[adata.obs.unilung_ann_level1 == "Fibroblast", :]
sc.pp.filter_genes(fib, min_cells=3)
sc.pp.filter_cells(fib, min_genes=200)
sc.pp.neighbors(fib, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(fib, min_dist=0.3)
sc.tl.leiden(fib, resolution=0.2, key_added='second_cluster')
sc.tl.rank_genes_groups(fib, 'second_cluster', method='t-test')

## add cell annotation
second_anno = {
    "0": "Alveolar adventitial fibroblast", "1": "Alveolar type I fibroblast cell",
    "2": "Myofibroblast", "3": "Myofibroblast", "4": "Myofibroblast", "5": "Alveolar adventitial fibroblast", "6": "Alveolar adventitial fibroblast", "7": "Myofibroblast"
}

fib.obs["unilung_ann_level2"] = fib.obs.second_cluster.map(second_anno)
sc.pl.umap(
    fib2, color=['second_cluster', 'unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['second_cluster', 'unilung_ann_level2', 'map_ann_level2']
)

marker = {
    'Alveolar adventitial fibroblast': ['SFRP2', 'SERPINF1', 'DPT', 'APOD'],
    'Alveolar type I fibroblast cell': ['PLIN2', 'LUM', 'GPC3'],
    'Myofibroblast': ['ASPN', 'ACTA2'],
    'Pathological fibroblast': ['CTHRC1', 'COL5A2', 'COL1A1']
}

sc.pl.dotplot(fib2, marker, groupby='unilung_ann_level2')





## lym
lym = adata[adata.obs.unilung_ann_level1 == "Lymphocyte", :]
sc.pp.filter_genes(lym, min_cells=3)
sc.pp.filter_cells(lym, min_genes=200)
sc.pp.neighbors(lym, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(lym, min_dist=0.3)
sc.tl.leiden(lym, resolution=0.2, key_added='second_cluster')
sc.tl.rank_genes_groups(lym, 'second_cluster', method='t-test')

## add cell annotation
second_anno = {
    "0": "T cell", "1": "T cell", "2": "T cell", "3": "B cell",
    "4": "NK cell", "5": "B cell", "7": "T cell"
}

lym.obs["unilung_ann_level2"] = lym.obs.second_cluster.map(second_anno)
sc.pl.umap(
    lym, color=['second_cluster', 'unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['second_cluster', 'unilung_ann_level2', 'map_ann_level2']
)

marker = {
    'B cell': ['CD79A', 'CD79B', 'MS4A1', 'CD74'],
    'NK cell': ['KLRD1', 'NKG7', 'FCGR3A'],
    'T cell': ['CD3D', 'CD3E', 'CD3G', 'IL7R']
}

sc.pl.dotplot(lym, marker, groupby='unilung_ann_level2')





## mye
mye = adata[adata.obs.unilung_ann_level1 == "Myeloid cell", :]
sc.pp.filter_genes(mye, min_cells=3)
sc.pp.filter_cells(mye, min_genes=200)
sc.pp.neighbors(mye, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(mye, min_dist=0.3)
sc.tl.leiden(mye, resolution=0.2, key_added='second_cluster')
sc.tl.rank_genes_groups(mye, 'second_cluster', method='t-test')

## add cell annotation
second_anno = {
    "0": "Macrophage", "1": "Monocyte", "2": "Monocyte", "3": "Macrophage",
    "4": "Neutrophilic granulocyte", "5": "Mast cell"
}

mye.obs["unilung_ann_level2"] = mye.obs.second_cluster.map(second_anno)
sc.pl.umap(
    mye, color=['second_cluster', 'unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['second_cluster', 'unilung_ann_level2', 'map_ann_level2']
)

marker = {
    'Dendritic cell': ['MS4A6A', 'CLEC10A', 'GPR183', 'CST3'],
    'Macrophage': ['SPP1', 'MRC1', 'MARCO', 'CTSB'],
    'Mast cell': ['GATA2', 'MS4A2', 'TPSAB1', 'TPSB2'],
    'Monocyte': ['CD14', 'LYZ', 'VCAN'],
    'Neutrophilic granulocyte': ['CSF3R', 'FCGR3B']
}
sc.pl.dotplot(mye, marker, groupby='unilung_ann_level2')







## third annotation
## Capillary endothelial cell
cap = adata[adata.obs.unilung_ann_level2 == "Capillary endothelial cell", :]
sc.pp.filter_genes(cap, min_cells=3)
sc.pp.filter_cells(cap, min_genes=200)
sc.pp.neighbors(cap, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(cap, min_dist=0.3)
sc.tl.leiden(cap, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(cap, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "Alveolar capillary cell", "1": "Alveolar capillary cell", "2": "General capillary cell",
    "3": "Alveolar capillary cell", "4": "Alveolar capillary cell", "5": "Alveolar capillary cell"
}
cap.obs["unilung_ann_level3"] = cap.obs.third_cluster.map(third_anno)
sc.pl.umap(
    cap, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'Alveolar capillary cell': ['EDNRB', 'TBX2', 'HPGD'],
    'General capillary cell': ['APLNR', 'IL7R', 'GPIHBP1']
}
sc.pl.dotplot(cap, marker, groupby='unilung_ann_level3')




## Alveolar cell
ale = adata[adata.obs.unilung_ann_level2 == "Alveolar cell", :]
sc.pp.filter_genes(alv, min_cells=3)
sc.pp.filter_cells(alv, min_genes=200)
sc.pp.neighbors(ale, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(ale, min_dist=0.3)
sc.tl.leiden(ale, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(ale, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "Type II alveolar cell", "1": "Type I alveolar cell", "2": "Type II alveolar cell",
    "3": "Type II alveolar cell", "4": "Type II alveolar cell", "5": "Type II alveolar cell",
    "6": "Type II alveolar cell", "7": "Type I alveolar cell", "8": "Type II alveolar cell",
    "9": "Type II alveolar cell", "10": "Type II alveolar cell"
}

ale.obs["unilung_ann_level3"] = ale.obs.third_cluster.map(third_anno)
sc.pl.umap(
    ale, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'Type I alveolar cell': ['AGER', 'CAV1', 'RTKN2'],
    'Type II alveolar cell': ['SFTPA1', 'SFTPC', 'SFTPB', 'ABCA3']
}
sc.pl.dotplot(ale, marker, groupby='unilung_ann_level3')



## Submucosal gland cell
smg = adata[adata.obs.unilung_ann_level2 == "Submucosal gland cell", :]
sc.pp.filter_genes(smg, min_cells=3)
sc.pp.filter_cells(smg, min_genes=200)
sc.pp.neighbors(smg, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(smg, min_dist=0.3)
sc.tl.leiden(smg, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(smg, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "SMG mucous cell", "1": "SMG serous cell", "2": "SMG serous cell",
    "3": "SMG duct cell", "4": "SMG serous cell", "5": "SMG mucous cell"
}
smg.obs["unilung_ann_level3"] = smg.obs.third_cluster.map(third_anno)
sc.pl.umap(
    smg, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'SMG duct cell': ['MGST1', 'KRT19', 'KLF5'],
    'SMG mucous cell': ['MUC5B', 'BPIFB2', 'MSMB', 'BPIFB1'],
    'SMG serous cell': ['LTF', 'LYZ', 'SAA2']
}
sc.pl.dotplot(smg, marker, groupby='unilung_ann_level3')




## Basal cell
bas = adata[adata.obs.unilung_ann_level2 == "Basal cell", :]
sc.pp.filter_genes(bas, min_cells=3)
sc.pp.filter_cells(bas, min_genes=200)
sc.pp.neighbors(bas, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(bas, min_dist=0.3)
sc.tl.leiden(bas, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(bas, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "Basal resting cell", "1": "Suprabasal cell", "2": "Suprabasal cell",
    "3": "Suprabasal cell", "4": "Basal resting cell", "5": "Suprabasal cell",
    "6": "Suprabasal cell"
}
bas.obs["unilung_ann_level3"] = bas.obs.third_cluster.map(third_anno)
sc.pl.umap(
    bas, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'Basal resting cell': ['KRT15', 'KRT17', 'TP63'],
    'Suprabasal cell': ['KRT19', 'KRT13', 'NOTCH3', 'SERPINB4']
}
sc.pl.dotplot(bas, marker, groupby='unilung_ann_level3')




## Dendritic cell
den = adata[adata.obs.unilung_ann_level2 == "Dendritic cell", :]
sc.pp.filter_genes(den, min_cells=3)
sc.pp.filter_cells(den, min_genes=200)
sc.pp.neighbors(den, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(den, min_dist=0.3)
sc.tl.leiden(den, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(den, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "Conventional dendritic cell", "1": "Conventional dendritic cell",
    "2": "Conventional dendritic cell", "3": "Plasmacytoid dendritic cell",
    "4": "Conventional dendritic cell", "5": "Plasmacytoid dendritic cell"
}

den.obs["unilung_ann_level3"] = den.obs.third_cluster.map(third_anno)
sc.pl.umap(
    den, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'Conventional dendritic cell': ['THBD', 'CD1C', 'CD1E', 'FCGR2B', 'CLEC10A'],
    'Plasmacytoid dendritic cell': ['CSF3R', 'IL3RA', 'TCF4', 'IRF7', 'IRF7'],
}
sc.pl.dotplot(den, marker, groupby='unilung_ann_level3')




## Monocyte
mon = adata[adata.obs.unilung_ann_level2 == "Monocyte", :]
sc.pp.filter_genes(mono, min_cells=3)
sc.pp.filter_cells(mono, min_genes=200)
sc.pp.neighbors(mon, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(mon, min_dist=0.3)
sc.tl.leiden(mon, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(mon, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "Classical monocyte", "1": "Classical monocyte",
    "2": "Non-classical monocyte", "3": "Non-classical monocyte",
    "4": "Non-classical monocyte", "5": "Promonocyte"
}
mon.obs["unilung_ann_level3"] = mon.obs.third_cluster.map(third_anno)
sc.pl.umap(
    mon, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'Classical monocyte': ['CD14', 'S100A9', 'S100A12', 'FCN1'],
    'Non-classical monocyte': ['FCGR3A', 'CSF1R', 'LILRB1', 'LILRB2'],
    'Promonocyte': ['MPO', 'VCAN', 'S100A8']
}
sc.pl.dotplot(mon, marker, groupby='unilung_ann_level3')




## Macrophage
mac = adata[adata.obs.unilung_ann_level2 == "Macrophage", :]
sc.pp.filter_genes(mac, min_cells=3)
sc.pp.filter_cells(mac, min_genes=200)
sc.pp.neighbors(mac, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(mac, min_dist=0.3)
sc.tl.leiden(mac, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(mac, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "Interstitial macrophage", "1": "Alveolar macrophage",
    "2": "Alveolar macrophage", "3": "Alveolar macrophage",
    "4": "Alveolar macrophage"
}
mac.obs["unilung_ann_level3"] = mac.obs.third_cluster.map(third_anno).astype("category")
sc.pl.umap(
    mac, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'Alveolar macrophage': ['FABP4', 'TREM2', 'MARCO', 'ALOX5', 'SPI1', 'ABCG1'],
    'Interstitial macrophage': ['MERTK', 'C1QC', 'C1QB', 'IL1B', 'C3AR1', 'NPL']
}
sc.pl.dotplot(mac, marker, groupby='unilung_ann_level3')




## b cell
bcell = adata[adata.obs.unilung_ann_level2 == "B cell", :]
sc.pp.filter_genes(bcell, min_cells=3)
sc.pp.filter_cells(bcell, min_genes=200)
sc.pp.neighbors(bcell, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(bcell, min_dist=0.3)
sc.tl.leiden(bcell, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(bcell, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "Mature B cell", "1": "Plasma cell",
    "2": "Mature B cell", "3": "Mature B cell",
    "4": "Plasma cell", "5": "Mature B cell",
    "6": "Plasma cell", "7": "Naive B cell",
    "8": "Plasma cell"
}

bcell.obs["unilung_ann_level3"] = bcell.obs.third_cluster.map(third_anno)
sc.pl.umap(
    bcell, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'Mature B cell': ['CD19', 'CD5', 'MS4A2', 'PAX5', 'CD79A'],
    'Naive B cell': ['MS4A1', 'CD22', 'CD37', 'CD19', 'TCL1A'],
    'Plasma cell': ['JCHAIN', 'MZB1', 'IGHG1', 'IGKC', 'XBP1']
}
sc.pl.dotplot(bcell, marker, groupby='unilung_ann_level3')




## t cell
tcell = adata[adata.obs.unilung_ann_level2 == "T cell", :]
sc.pp.filter_genes(tcell, min_cells=3)
sc.pp.filter_cells(tcell, min_genes=200)
sc.pp.neighbors(tcell, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(tcell, min_dist=0.3)
sc.tl.leiden(tcell, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(tcell, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "CD4 T cell", "1": "NKT cell", "2": "Cytotoxic T cell",
    "3": "NKT cell", "4": "CD4 T cell", "5": "CD8 T cell",
    "6": "CD8 T cell", "7": "CD4 T cell", "8": "CD8 T cell"
}

tcell.obs["unilung_ann_level3"] = tcell.obs.third_cluster.map(third_anno)
sc.pl.umap(
    tcell, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'CD4 T cell': ['CD4', 'CCR7', 'CD3E', 'IL7R', 'LEF1'],
    'CD8 T cell': ['CD8A', 'CD8B', 'CD3E'],
    'Cytotoxic T cell': ['PRF1', 'CD3D', 'CD8A'],
    'NKT cell': ['KLRF1', 'NKG7', 'GNLY', 'GZMK', 'GZMA', 'CCL5'],
    'Proliferating T cell': ['TK1', 'CENPW']
}
sc.pl.dotplot(tcell, marker, groupby='unilung_ann_level3')








## fourth annotation
## cd4t cell
cd4 = adata[adata.obs.unilung_ann_level3 == "CD4 T cell", :]
sc.pp.filter_genes(cd4, min_cells=3)
sc.pp.filter_cells(cd4, min_genes=200)
sc.pp.neighbors(cd4, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(cd4, min_dist=0.3)
sc.tl.leiden(cd4, resolution=0.5, key_added='fourth_cluster')
sc.tl.rank_genes_groups(cd4, 'fourth_cluster', method='t-test')

## add cell annotation
fourth_anno = {
    "0": "Naive CD4 T cell", "1": "Memory CD4 T cell", "2": "Naive CD4 T cell",
    "3": "Memory CD4 T cell", "4": "Naive CD4 T cell", "5": "Naive CD4 T cell",
    "6": "Naive CD4 T cell"
}

cd4.obs["unilung_ann_level4"] = cd4.obs.fourth_cluster.map(fourth_anno)
sc.pl.umap(
    cd4, color=['fourth_cluster', 'unilung_ann_level4', 'map_ann_level4'],
    ncols=1, title=['fourth_cluster', 'unilung_ann_level4', 'map_ann_level4']
)

marker = {
    'Memory CD4 T cell': ['IL2RB', 'CD44', 'TRAC', 'S100A11'],
    'Naive CD4 T cell': ['CCR7', 'CD27', 'CD28', 'TMSB4X']
}

sc.pl.dotplot(cd4, marker, groupby='unilung_ann_level4')




## cd8t cell
cd8 = adata[adata.obs.unilung_ann_level3 == "CD8 T cell", :]
sc.pp.filter_genes(cd8, min_cells=3)
sc.pp.filter_cells(cd8, min_genes=200)
sc.pp.neighbors(cd8, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(cd8, min_dist=0.3)
sc.tl.leiden(cd8, resolution=0.5, key_added='fourth_cluster')
sc.tl.rank_genes_groups(cd8, 'fourth_cluster', method='t-test')

## add cell annotation
fourth_anno = {
    "0": "Memory CD8 T cell", "1": "Memory CD8 T cell", "2": "Naive CD8 T cell",
    "3": "Effector CD8 T cell", "4": "Naive CD8 T cell", "5": "Effector CD8 T cell",
    "6": "Memory CD8 T cell", "7": "Naive CD8 T cell", "8": "Cycling CD8 T cell"
}
cd8.obs["unilung_ann_level4"] = cd8.obs.fourth_cluster.map(fourth_anno)
sc.pl.umap(
    cd8, color=['fourth_cluster', 'unilung_ann_level4', 'map_ann_level4'],
    ncols=1, title=['fourth_cluster', 'unilung_ann_level4', 'map_ann_level4']
)

marker = {
    'Cycling CD8 T cell': ['PCNA', 'MKI67'],
    'Effector CD8 T cell': ['CD8A', 'CD8B', 'PRF1', 'GZMB', 'GZMH'],
    'Memory CD8 T cell': ['CD27', 'CD44', 'CXCR3', 'CD28'],
    'Naive CD8 T cell': ['CCR7', 'CXCR4', 'LEF1', 'TMSB4X'],
}
sc.pl.dotplot(cd8, marker, groupby='unilung_ann_level4')


