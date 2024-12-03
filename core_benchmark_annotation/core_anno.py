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
    n_top_genes=2000,
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
sc.tl.leiden(adata, resolution=0.2, key_added='first_cluster')
sc.tl.rank_genes_groups(adata, 'first_cluster', method='t-test')

## add cell annotation
first_anno = {
    "0": "Lymphocyte", "1": "Myeloid cell", "2": "Myeloid cell", "3": "Epithelial cell",
    "4": "Fibroblast", "5": "Epithelial cell", "6": "Epithelial cell", "7": "Endothelial cell",
    "8": "Lymphocyte", "9": "Epithelial cell", "10": "Myeloid cell", "11": "Neuroendocrine cell",
    "12": "Myeloid cell", "13": "Neuroendocrine cell", "14": "Submucosal gland cell",
    "15": "Lymphocyte", "16": "Lymphocyte", "17": "Myeloid cell", "18": "Smooth muscle cell"
}
adata.obs["first_annotation"] = adata.obs.first_cluster.map(first_anno)
adata.obs["unilung_ann_level1"] = adata.obs.first_cluster.map(first_anno)
sc.pl.umap(
    adata, color=['first_cluster', 'unilung_ann_level1', 'map_ann_level1'],
    ncols=1, title=['first_cluster', 'unilung_ann_level1', 'map_ann_level1']
)

marker = {
    'Endothelial cell': ['PECAM1', 'CLDN5', 'VWF'],
    'Epithelial cell': ['EPCAM', 'KRT19', 'CAPS', 'KRT7'],
    'Fibroblast': ['DCN', 'LUM', 'COL1A2', 'COL1A1'],
    'Lymphocyte': ['PTPRC', 'CCL5', 'CD3D'],
    'Myeloid cell': ['LYZ', 'CD14', 'MARCO'],
    'Neuroendocrine cell': ['CHGA', 'TUBA1A', 'UCHL1'],
    'Smooth muscle cell': ['ACTA2', 'CALD1', 'TPM2'],
    'Submucosal gland cell': ['LTF', 'SLPI', 'LYZ']
}
sc.pl.dotplot(adata, marker, groupby='first_annotation')




## second annotation
## endo
endo = adata[adata.obs.unilung_ann_level1 == "Endothelial cell", :]
sc.pp.neighbors(endo, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(endo, min_dist=0.3)
sc.tl.leiden(endo, resolution=0.3, key_added='second_cluster')
sc.tl.rank_genes_groups(endo, 'second_cluster', method='t-test')

## add cell annotation
second_anno = {
    "0": "Capillary endothelial cell", "1": "Vein endothelial cell", "2": "Capillary endothelial cell",
    "3": "Capillary endothelial cell", "4": "Vein endothelial cell", "5": "Lymphatic endothelial cell",
    "6": "Capillary endothelial cell", "7": "Vein endothelial cell",
}
endo.obs["unilung_ann_level2"] = endo.obs.second_cluster.map(second_anno)
sc.pl.umap(
    endo, color=['second_cluster', 'unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['second_cluster', 'unilung_ann_level2', 'map_ann_level2']
)

marker = {
    'Capillary endothelial cell': ['CA4'],
    'Vein endothelial cell': ['ACKR1'],
    'Lymphatic endothelial cell': ['PROX1', 'PDPN'],
    'Artery endothelial cell': ['GJA5', 'BMX']
}
sc.pl.dotplot(endo, marker, groupby='unilung_ann_level2')




## epi
epi = adata[adata.obs.unilung_ann_level1 == "Epithelial cell", :]
sc.pp.neighbors(epi, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(epi, min_dist=0.3)
sc.tl.leiden(epi, resolution=0.3, key_added='second_cluster')
sc.tl.rank_genes_groups(epi, 'second_cluster', method='t-test')

## add cell annotation
second_anno = {
    "0": "Ciliated cell", "1": "Alveolar cell", "2": "Basal cell", "3": "Alveolar cell",
    "4": "Secretory cell", "5": "Alveolar cell", "6": "Alveolar cell", "7": "Alveolar cell",
    "8": "Alveolar cell", "9": "Alveolar cell", "10": "Ciliated cell", "11": "Alveolar cell",
    "12": "Squamous cell", "13": "Ciliated cell", "14": "Alveolar cell",
    "15": "Basal cell", "16": "Deuterosomal cell"
}
epi.obs["unilung_ann_level2"] = epi.obs.second_cluster.map(second_anno)
sc.pl.umap(
    epi, color=['second_cluster', 'unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['second_cluster', 'unilung_ann_level2', 'map_ann_level2']
)

marker = {
    'Secretory cell': ['SCGB3A1', 'BPIFB1'],
    'Alveolar cell': ['AQP4', 'CLDN18'],
    'Ciliated cell': ['FOXJ1', 'PIFO', 'TPPP3'],
    'Basal cell': ['KRT5', 'KRT15', 'S100A2', 'TP63'],
    'Squamous cell': ['SCEL', 'SPRR1B'],
    'Deuterosomal cell': ['HES6'],
}
sc.pl.dotplot(epi, marker, groupby='unilung_ann_level2')




## fibro
fib = adata[adata.obs.unilung_ann_level1 == "Fibroblast", :]
sc.pp.neighbors(fib, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(fib, min_dist=0.3)
sc.tl.leiden(fib, resolution=0.3, key_added='second_cluster')
sc.tl.rank_genes_groups(fib, 'second_cluster', method='t-test')

## add cell annotation
second_anno = {
    "0": "Myofibroblast", "1": "Airway fibroblast", "2": "Basal cell", "3": "Alveolar fibroblast",
    "4": "Adventitial fibroblast", "5": "Adventitial fibroblast", "6": "Secretory cell", "7": "Alveolar fibroblast"
}
fib.obs["unilung_ann_level2"] = fib.obs.second_cluster.map(second_anno)
sc.pl.umap(
    fib2, color=['second_cluster', 'unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['second_cluster', 'unilung_ann_level2', 'map_ann_level2']
)

marker = {
    'Myofibroblast': ['MYLK', 'COL6A3'],
    'Alveolar fibroblast': ['CFD', 'FGFR4'],
    'Adventitial fibroblast': ['SFRP2', 'SERPINF1'],
    'Airway fibroblast': ['S100A4', 'FGF7', 'TNC'],
}
sc.pl.dotplot(fib2, marker, groupby='unilung_ann_level2')





## lym
lym = adata[adata.obs.unilung_ann_level1 == "Lymphocyte", :]
sc.pp.neighbors(lym, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(lym, min_dist=0.3)
sc.tl.leiden(lym, resolution=0.3, key_added='second_cluster')
sc.tl.rank_genes_groups(lym, 'second_cluster', method='t-test')

## add cell annotation
second_anno = {
    "0": "T cell", "1": "T cell", "2": "NK cell", "3": "T cell",
    "4": "B cell", "5": "T cell", "6": "B cell", "7": "T cell",
    "8": "B cell", "9": "NK cell", "10": "Erythrocyte", "11": "B cell"
}
lym.obs["unilung_ann_level2"] = lym.obs.second_cluster.map(second_anno)
sc.pl.umap(
    lym, color=['second_cluster', 'unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['second_cluster', 'unilung_ann_level2', 'map_ann_level2']
)

marker = {
    'B cell': ['CD79A', 'CD79B', 'MS4A1'],
    'T cell': ['CD3D', 'CD3E', 'CD3G'],
    'NK cell': ['KLRD1', 'NKG7', 'FCGR3A'],
}
sc.pl.dotplot(lym, marker, groupby='unilung_ann_level2')





## mye
mye = adata[adata.obs.unilung_ann_level1 == "Myeloid cell", :]
sc.pp.neighbors(mye, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(mye, min_dist=0.3)
sc.tl.leiden(mye, resolution=0.3, key_added='second_cluster')
sc.tl.rank_genes_groups(mye, 'second_cluster', method='t-test')

## add cell annotation
second_anno = {
    "0": "Macrophage", "1": "Macrophage", "2": "Monocyte", "3": "Macrophage", "4": "Macrophage",
    "5": "Mast cell", "6": "Neutrophilic granulocyte", "7": "Macrophage", "8": "Megakaryocyte"
}

mye.obs["unilung_ann_level2"] = mye.obs.second_cluster.map(second_anno)
sc.pl.umap(
    mye, color=['second_cluster', 'unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['second_cluster', 'unilung_ann_level2', 'map_ann_level2']
)

marker = {
    'Mast cell': ['GATA2', 'MS4A2', 'TPSAB1', 'TPSB2'],
    'Neutrophilic granulocyte': ['CSF3R', 'CXCR2', 'FCGR3B'],
    'Dendritic cell': ['MS4A6A', 'CLEC10A', 'GPR183'],
    'Monocyte': ['CD14', 'LYZ', 'VCAN'],
    'Macrophage': ['CD163', 'CD68', 'MRC1', 'MARCO'],
    'Megakaryocyte': ['PPBP', 'PF4', 'ITGA2B']
}
sc.pl.dotplot(mye, marker, groupby='unilung_ann_level2')





## SMC
smc = adata[adata.obs.unilung_ann_level1 == "Smooth muscle cell", :]
sc.pp.neighbors(smc, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(smc, min_dist=0.3)
sc.tl.leiden(smc, resolution=0.3, key_added='second_cluster')
sc.tl.rank_genes_groups(smc, 'second_cluster', method='t-test')

## add cell annotation
second_anno = {
    "0": "Bronchial smooth muscle cell", "1": "Bronchial smooth muscle cell", "2": "Bronchial smooth muscle cell",
    "3": "Vascular smooth muscle cell", "4": "Bronchial smooth muscle cell",
    "5": "Vascular smooth muscle cell", "6": "Bronchial smooth muscle cell", "7": "Erythrocyte"
}
smc.obs["unilung_ann_level2"] = smc.obs.second_cluster.map(second_anno)
sc.pl.umap(
    smc, color=['second_cluster', 'unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['second_cluster', 'unilung_ann_level2', 'map_ann_level2']
)

marker = {
    'Vascular smooth muscle cell': ['TAGLN', 'ACTA2'],
    'Bronchial smooth muscle cell': ['ACTA2'],
}
sc.pl.dotplot(smc, marker, groupby='unilung_ann_level2')




## SMG
smg = adata[adata.obs.unilung_ann_level1 == "Submucosal gland cell", :]
sc.pp.neighbors(smg, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(smg, min_dist=0.3)
sc.tl.leiden(smg, resolution=0.3, key_added='second_cluster')
sc.tl.rank_genes_groups(smg, 'second_cluster', method='t-test')

## add cell annotation
second_anno = {
    "0": "SMG serous cell", "1": "SMG serous cell", "2": "SMG duct cell"
}
smg.obs["unilung_ann_level2"] = smg.obs.second_cluster.map(second_anno)
sc.pl.umap(
    smg, color=['second_cluster', 'unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['second_cluster', 'unilung_ann_level2', 'map_ann_level2']
)

marker = {
    'SMG serous cell': ['BPIFA1', 'DMBT1'],
    'SMG duct cell': ['MGST1', 'KRT19']
}
sc.pl.dotplot(smg, marker, groupby='unilung_ann_level2')


other = adata[adata.obs.unilung_ann_level1.isin(["Neuroendocrine cell", "Cancer cell"])]
sc.pp.neighbors(other, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(other, min_dist=0.3)
sc.tl.leiden(other, resolution=0.3, key_added='second_cluster')
other.obs["unilung_ann_level2"] = 'Unclassified'
other.obs.unilung_ann_level2 = other.obs.unilung_ann_level2.astype('category')




## merge second level datasets and add level2 anno to adata
def add_batch(adata,key='BATCH'):
    new_adata = adata.copy()
    new_adata.obs[key] = new_adata.obs.batch
    return new_adata

endonew = add_batch(endo)
fibnew = add_batch(fib)
epinew = add_batch(epi)
myenew = add_batch(mye)
lymnew = add_batch(lym)
smcnew = add_batch(smc)
smgnew = add_batch(smg)
othernew = add_batch(other)

adata1 = sc.AnnData.concatenate(endonew, fibnew, epinew, myenew, lymnew, smcnew, smgnew, othernew, batch_key='BATCH')
adata1.obs.drop(columns='BATCH', inplace=True)
adata1.obs.unilung_ann_level2 = adata1.obs.unilung_ann_level2.astype('category')
sc.pp.neighbors(adata1, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(adata1, min_dist=0.3)
sc.pl.umap(
    adata1, color=['unilung_ann_level1', 'map_ann_level1'],
    ncols=1, title=['unilung_ann_level1', 'map_ann_level1']
)
plt.savefig(dir+'/second_anno_lv1.pdf', dpi=300, format='pdf', bbox_inches='tight')
sc.pl.umap(
    adata1, color=['unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['unilung_ann_level2', 'map_ann_level2']
)
plt.savefig(dir+'/second_anno_lv2.pdf', dpi=300, format='pdf', bbox_inches='tight')
sc.pl.umap(
    adata1, color=['unilung_ann_level1', 'unilung_ann_level2'],
    ncols=1, title=['unilung_ann_level1', 'unilung_ann_level2']
)
plt.savefig(dir+'/second_anno.pdf', dpi=300, format='pdf', bbox_inches='tight')




## third annotation
## Alveolar cell
ale = adata1[adata1.obs.unilung_ann_level2 == "Alveolar cell", :]
sc.pp.neighbors(ale, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(ale, min_dist=0.3)
sc.tl.leiden(ale, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(ale, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "Type II alveolar cell", "1": "Type I alveolar cell", "2": "Type II alveolar cell",
    "3": "Type II alveolar cell", "4": "Type II alveolar cell", "5": "Type II alveolar cell",
    "6": "Type I alveolar cell", "7": "Type II alveolar cell", "8": "Type II alveolar cell"
}
ale.obs["unilung_ann_level3"] = ale.obs.third_cluster.map(third_anno)
sc.pl.umap(
    ale, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'Type I alveolar cell': ['AGER', 'CAV1'],
    'Type II alveolar cell': ['SFTPA1', 'SFTPC', 'SFTPD']
}
sc.pl.dotplot(ale, marker, groupby='unilung_ann_level3')




## Basal cell
bas = adata1[adata1.obs.unilung_ann_level2 == "Basal cell", :]
sc.pp.neighbors(bas, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(bas, min_dist=0.3)
sc.tl.leiden(bas, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(bas, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "Suprabasal cell", "1": "Basal resting cell", "2": "Basal resting cell",
    "3": "Basal resting cell", "4": "Basal resting cell", "5": "Suprabasal cell",
    "6": "Basal resting cell"
}
bas.obs["unilung_ann_level3"] = bas.obs.third_cluster.map(third_anno)
sc.pl.umap(
    bas, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'Suprabasal cell': ['KRT19', 'TP63'],
    'Basal resting cell': ['KRT15', 'KRT17']
}
sc.pl.dotplot(bas, marker, groupby='unilung_ann_level3')




## Dendritic cell
den = adata1[adata1.obs.unilung_ann_level2 == "Dendritic cell", :]
sc.pp.neighbors(den, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(den, min_dist=0.3)
sc.tl.leiden(den, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(den, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "Conventional dendritic cell", "1": "Conventional dendritic cell", "2": "Conventional dendritic cell",
    "3": "Conventional dendritic cell", "4": "Plasmacytoid dendritic cell", "5": "Migratory dendritic cell",
    "6": "Conventional dendritic cell", "7": "Plasmacytoid dendritic cell", "8": "Conventional dendritic cell"
}
den.obs["unilung_ann_level3"] = den.obs.third_cluster.map(third_anno)
sc.pl.umap(
    den, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'Conventional dendritic cell': ['CD1C', 'THBD'],
    'Plasmacytoid dendritic cell': ['CSF3R', 'S100A8'],
    'Migratory dendritic cell': ['TMSB10', 'CCR7']
}
sc.pl.dotplot(den, marker, groupby='unilung_ann_level3')




## Monocyte
mon = adata1[adata1.obs.unilung_ann_level2 == "Monocyte", :]
sc.pp.neighbors(mon, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(mon, min_dist=0.3)
sc.tl.leiden(mon, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(mon, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "Classical monocyte", "1": "Classical monocyte", "2": "Classical monocyte",
    "3": "Non-classical monocyte", "4": "Promonocyte", "5": "Promonocyte"
}
mon.obs["unilung_ann_level3"] = mon.obs.third_cluster.map(third_anno)
sc.pl.umap(
    mon, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'Classical monocyte': ['CD14'],
    'Non-classical monocyte': ['CD14'],
    'Promonocyte': ['VCAN']
}
sc.pl.dotplot(mon, marker, groupby='unilung_ann_level3')




## Macrophage
mac = adata1[adata1.obs.unilung_ann_level2 == "Macrophage", :]
sc.pp.neighbors(mac, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(mac, min_dist=0.3)
sc.tl.leiden(mac, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(mac, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "M2 macrophage", "1": "M2 macrophage", "2": "M2 macrophage",
    "3": "M1 macrophage", "4": "M1 macrophage", "5": "M2 macrophage", "6": "M2 macrophage"
}
mac.obs["unilung_ann_level3"] = mac.obs.third_cluster.map(third_anno).astype("category")
sc.pl.umap(
    mac, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'M1 macrophage': ['CD86'],
    'M2 macrophage': ['CD163']
}
sc.pl.dotplot(mac, marker, groupby='unilung_ann_level3')




## b cell
bcell = adata1[adata1.obs.unilung_ann_level2 == "B cell", :]
sc.pp.neighbors(bcell, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(bcell, min_dist=0.3)
sc.tl.leiden(bcell, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(bcell, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "Plasma cell", "1": "Naive B cell", "2": "Memory B cell", "3": "Naive B cell",
    "4": "Naive B cell", "5": "Plasma cell", "6": "Plasma cell", "7": "Follicular B cell"
}
bcell.obs["unilung_ann_level3"] = bcell.obs.third_cluster.map(third_anno)
sc.pl.umap(
    bcell, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'Plasma cell': ['MZB1', 'IGHG1'],
    'Plasmablast cell': ['MKI67', 'CD38'],
    'Naive B cell': ['MS4A1', 'CD37'],
    'Memory B cell': ['CD27', 'CD19'],
    'Mature B cell': ['CD19', 'CD5'],
    'Follicular B cell': ['CD22']
}
sc.pl.dotplot(bcell, marker, groupby='unilung_ann_level3')




## t cell
tcell = adata1[adata1.obs.unilung_ann_level2 == "T cell", :]
sc.pp.neighbors(tcell, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(tcell, min_dist=0.3)
sc.tl.leiden(tcell, resolution=0.3, key_added='third_cluster')
sc.tl.rank_genes_groups(tcell, 'third_cluster', method='t-test')

## add cell annotation
third_anno = {
    "0": "CD8 T cell", "1": "CD4 T cell", "2": "CD4 T cell", "3": "NKT cell",
    "4": "CD8 T cell", "5": "NKT cell", "6": "CD8 T cell", "7": "Cytotoxic T cell"
}
tcell.obs["unilung_ann_level3"] = tcell.obs.third_cluster.map(third_anno)
sc.pl.umap(
    tcell, color=['third_cluster', 'unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['third_cluster', 'unilung_ann_level3', 'map_ann_level3']
)

marker = {
    'CD4 T cell': ['CD4', 'IL7R', 'CCR7'],
    'NKT cell': ['PITPNC1', 'MBNL1'],
    'Treg cell': ['IL2RA', 'FOXP3'],
    'Cycling T cell': ['MKI67', 'CDK1'],
    'CD8 T cell': ['CD8A', 'CD8B'],
    'Cytotoxic T cell': ['PRF1'],
}
sc.pl.dotplot(tcell, marker, groupby='unilung_ann_level3')




escape = adata1.obs.unilung_ann_level2.isin(["B cell", "T cell", "Alveolar cell", "Basal cell", "Dendritic cell", "Monocyte", "Macrophage"])
other = adata1[~escape]
sc.pp.neighbors(other, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(other, min_dist=0.3)
sc.tl.leiden(other, resolution=0.3, key_added='third_cluster')
other.obs["unilung_ann_level3"] = 'Unclassified'
other.obs.unilung_ann_level3 = other.obs.unilung_ann_level3.astype('category')


## merge third level datasets and add level2 anno to adata
def add_batch(adata,key='BATCH'):
    new_adata = adata.copy()
    new_adata.obs[key] = new_adata.obs.batch
    return new_adata

alenew = add_batch(ale)
basnew = add_batch(bas)
dennew = add_batch(den)
monnew = add_batch(mon)
macnew = add_batch(mac)
bcellnew = add_batch(bcell)
tcellnew = add_batch(tcell)
othernew = add_batch(other)

adata2 = sc.AnnData.concatenate(alenew, basnew, dennew, monnew, macnew, bcellnew, tcellnew, othernew, batch_key='BATCH')
adata2.obs.drop(columns='BATCH', inplace=True)
adata2.obs.unilung_ann_level3 = adata2.obs.unilung_ann_level3.astype('category')
sc.pp.neighbors(adata2, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(adata2, min_dist=0.3)
sc.pl.umap(
    adata2, color=['unilung_ann_level1', 'map_ann_level1'],
    ncols=1, title=['unilung_ann_level1', 'map_ann_level1']
)
plt.savefig(dir+'/third_anno_lv1.pdf', dpi=300, format='pdf', bbox_inches='tight')
sc.pl.umap(
    adata2, color=['unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['unilung_ann_level2', 'map_ann_level2']
)
plt.savefig(dir+'/third_anno_lv2.pdf', dpi=300, format='pdf', bbox_inches='tight')
sc.pl.umap(
    adata2, color=['unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['unilung_ann_level3', 'map_ann_level3']
)
plt.savefig(dir+'/third_anno_lv3.pdf', dpi=300, format='pdf', bbox_inches='tight')
sc.pl.umap(
    adata2, color=['unilung_ann_level1', 'unilung_ann_level2', 'unilung_ann_level3'],
    ncols=1, title=['unilung_ann_level1', 'unilung_ann_level2', 'unilung_ann_level3']
)
plt.savefig(dir+'/third_anno.pdf', dpi=300, format='pdf', bbox_inches='tight')




## fourth annotation
## cd4t cell
cd4 = adata2[adata2.obs.unilung_ann_level3 == "CD4 T cell", :]
sc.pp.neighbors(cd4, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(cd4, min_dist=0.3)
sc.tl.leiden(cd4, resolution=0.3, key_added='fourth_cluster')
sc.tl.rank_genes_groups(cd4, 'fourth_cluster', method='t-test')

## add cell annotation
fourth_anno = {
    "0": "Memory CD4 T cell", "1": "Memory CD4 T cell", "2": "Naive CD4 T cell", "3": "Naive CD4 T cell"
}
cd4.obs["unilung_ann_level4"] = cd4.obs.fourth_cluster.map(fourth_anno)
sc.pl.umap(
    cd4, color=['fourth_cluster', 'unilung_ann_level4', 'map_ann_level4'],
    ncols=1, title=['fourth_cluster', 'unilung_ann_level4', 'map_ann_level4']
)

marker = {
    'Naive CD4 T cell': ['CCR7', 'LEF1', 'TMSB4X'],
    'Memory CD4 T cell': ['S100A11', 'TRAC', 'CD2']
}
sc.pl.dotplot(cd4, marker, groupby='unilung_ann_level4')




## cd8t cell
cd8 = adata2[adata2.obs.unilung_ann_level3 == "CD8 T cell", :]
sc.pp.neighbors(cd8, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(cd8, min_dist=0.3)
sc.tl.leiden(cd8, resolution=0.3, key_added='fourth_cluster')
sc.tl.rank_genes_groups(cd8, 'fourth_cluster', method='t-test')

## add cell annotation
fourth_anno = {
    "0": "Naive CD8 T cell", "1": "Naive CD8 T cell", "2": "Memory CD8 T cell",
    "3": "Memory CD8 T cell", "4": "Naive CD8 T cell"
}
cd8.obs["unilung_ann_level4"] = cd8.obs.fourth_cluster.map(fourth_anno)
sc.pl.umap(
    cd8, color=['fourth_cluster', 'unilung_ann_level4', 'map_ann_level4'],
    ncols=1, title=['fourth_cluster', 'unilung_ann_level4', 'map_ann_level4']
)

marker = {
    'Naive CD8 T cell': ['TMSB4X'],
    'Memory CD8 T cell': ['GZMB', 'TRAC', 'CD2'],
}
sc.pl.dotplot(cd8, marker, groupby='unilung_ann_level4')




escape = adata2.obs.unilung_ann_level3.isin(["CD8 T cell", "CD4 T cell"])
other = adata2[~escape]
sc.pp.neighbors(other, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(other, min_dist=0.3)
sc.tl.leiden(other, resolution=0.3, key_added='fourth_cluster')
other.obs["unilung_ann_level4"] = 'Unclassified'
other.obs.unilung_ann_level4 = other.obs.unilung_ann_level4.astype('category')



## merge third level datasets and add level2 anno to adata
def add_batch(adata,key='BATCH'):
    new_adata = adata.copy()
    new_adata.obs[key] = new_adata.obs.batch
    return new_adata

cd4new = add_batch(cd4)
cd8new = add_batch(cd8)
othernew = add_batch(other)

adata3 = sc.AnnData.concatenate(cd4new, cd8new, othernew, batch_key='BATCH')
adata3.obs.drop(columns='BATCH', inplace=True)
adata3.obs.unilung_ann_level4 = adata3.obs.unilung_ann_level4.astype('category')
sc.pp.neighbors(adata3, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(adata3, min_dist=0.3)
sc.pl.umap(
    adata3, color=['unilung_ann_level1', 'map_ann_level1'],
    ncols=1, title=['unilung_ann_level1', 'map_ann_level1']
)
plt.savefig(dir+'/fourth_anno_lv1.pdf', dpi=300, format='pdf', bbox_inches='tight')
sc.pl.umap(
    adata3, color=['unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['unilung_ann_level2', 'map_ann_level2']
)
plt.savefig(dir+'/fourth_anno_lv2.pdf', dpi=300, format='pdf', bbox_inches='tight')
sc.pl.umap(
    adata3, color=['unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['unilung_ann_level3', 'map_ann_level3']
)
plt.savefig(dir+'/fourth_anno_lv3.pdf', dpi=300, format='pdf', bbox_inches='tight')
sc.pl.umap(
    adata3, color=['unilung_ann_level4', 'map_ann_level4'],
    ncols=1, title=['unilung_ann_level4', 'map_ann_level4']
)
plt.savefig(dir+'/fourth_anno_lv4.pdf', dpi=300, format='pdf', bbox_inches='tight')
sc.pl.umap(
    adata3, color=['unilung_ann_level1', 'unilung_ann_level2', 'unilung_ann_level3', 'unilung_ann_level4'],
    ncols=1, title=['unilung_ann_level1', 'unilung_ann_level2', 'unilung_ann_level3', 'unilung_ann_level4']
)
plt.savefig(dir+'/fourth_anno.pdf', dpi=300, format='pdf', bbox_inches='tight')




def caculate_same(ad, lv1, lv2, lv3):
    num = ((ad.obs['unilung_ann_level' + lv1].astype('str') == ad.obs['map_ann_level' + lv1].astype('str'))|
           (ad.obs['unilung_ann_level' + lv2].astype('str') == ad.obs['map_ann_level' + lv2].astype('str'))|
           (ad.obs['unilung_ann_level' + lv3].astype('str') == ad.obs['map_ann_level' + lv3].astype('str')))
    return num

def get_percent(ad, lv, out_csv):
    output = pd.DataFrame()
    for i in ad.obs['unilung_ann_level'+lv].cat.categories:
        escape = ad.obs['unilung_ann_level'+lv].isin(["Unknown"])
        escape2 = ad[~escape]
        kk = escape2[escape2.obs['unilung_ann_level'+lv].isin([i])]
        samecell = kk.obs['unilung_ann_level'+lv].astype('str') == kk.obs['map_ann_level'+lv].astype('str')
        correct_labeled = (samecell.sum() / len(kk.obs)) * 100
        subcell = kk[~samecell]

        if lv == '1':
            covercell = caculate_same(subcell, '2','3','4').sum()
            cover_labeled = (covercell / len(kk.obs)) * 100
        if lv == '2':
            covercell = caculate_same(subcell, '1','3','4').sum()
            cover_labeled = (covercell / len(kk.obs)) * 100
        if lv == '3':
            covercell = caculate_same(subcell, '1', '2', '4').sum()
            cover_labeled = (covercell / len(kk.obs)) * 100
        else:
            covercell = caculate_same(subcell, '1', '2', '3').sum()
            cover_labeled = (covercell / len(kk.obs)) * 100
        mislabeled = 100 - correct_labeled - cover_labeled
        output.loc['correct_labeled', i] = correct_labeled
        output.loc['cover_labeled', i] = cover_labeled
        output.loc['mislabeled', i] = mislabeled

        print('the percent_match of {} is {}'.format(i, correct_labeled))
        print('the percent_cover of {} is {}'.format(i, cover_labeled))
        print('the percent_mistake of {} is {}'.format(i, mislabeled))
        output.T.to_csv(out_csv)
        gc.collect()

get_percent(adata3, '1', dir+'/anno_stat_lv1.csv')
get_percent(adata3, '2', dir+'/anno_stat_lv2.csv')
get_percent(adata3, '3', dir+'/anno_stat_lv3.csv')
get_percent(adata3, '4', dir+'/anno_stat_lv4.csv')


## add cell annotation to final celltype
adata4 = adata3.copy()
for i in adata4.obs_names:
    if adata4.obs.loc[i, 'unilung_ann_level2'] == 'Unclassified':
        adata4.obs.loc[i, 'final_celltype'] = adata4.obs.loc[i, 'unilung_ann_level1']
    elif adata4.obs.loc[i, 'unilung_ann_level3'] == 'Unclassified':
        adata4.obs.loc[i, 'final_celltype'] = adata4.obs.loc[i, 'unilung_ann_level2']
    elif adata4.obs.loc[i, 'unilung_ann_level4'] == 'Unclassified':
        adata4.obs.loc[i, 'final_celltype'] = adata4.obs.loc[i, 'unilung_ann_level3']
    else:
        adata4.obs.loc[i, 'final_celltype'] = adata4.obs.loc[i, 'unilung_ann_level4']
    adata4.obs.final_celltype = adata4.obs.final_celltype.astype('category')

sc.pl.umap(
    adata4, color=['unilung_ann_level1', 'unilung_ann_level2', 'unilung_ann_level3', 'unilung_ann_level4', 'final_celltype'],
    ncols=1, title=['unilung_ann_level1', 'unilung_ann_level2', 'unilung_ann_level3', 'unilung_ann_level4', 'final_celltype']
)
plt.savefig(dir+'/final_anno.pdf', dpi=300, format='pdf', bbox_inches='tight')
adata4.write_h5ad('final_anno.h5ad')
adata4.obs.to_csv(dir+'/meta.csv', sep=",")





## high percent_mistake cell type
mislv1 = adata4[adata4.obs['unilung_ann_level1'].isin(['Fibroblast', 'Submucosal gland cell'])]

sc.pp.neighbors(mislv1, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(mislv1, min_dist=0.3)
sc.tl.leiden(mislv1, resolution=0.1, key_added='mis_percent')
sc.tl.rank_genes_groups(mislv1, 'mis_percent', method='t-test')

marker = {
    'Fibroblast': ['DCN', 'LUM', 'COL1A2', 'COL1A1'],
    'Submucosal gland cell': ['LTF', 'SLPI', 'LYZ']
}
sc.pl.umap(
    mislv1, color=['unilung_ann_level1', 'map_ann_level1'],
    ncols=1, title=['unilung_celltype', 'original_celltype']
)
plt.savefig(dir+'/mislv1_umap.pdf', dpi=300, format='pdf', bbox_inches='tight')
sc.pl.dotplot(mislv1, marker, groupby='unilung_ann_level1')



mislv2 = adata4[adata4.obs['unilung_ann_level2'].isin(['Adventitial fibroblast', 'Alveolar fibroblast', 'Deuterosomal cell', 'SMG duct cell', 'SMG serous cell', 'Secretory cell'])]
sc.pp.neighbors(mislv2, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(mislv2, min_dist=0.3)
sc.tl.leiden(mislv2, resolution=0.2, key_added='mis_percent')
sc.tl.rank_genes_groups(mislv2, 'mis_percent', method='t-test')

marker = {
    'Adventitial fibroblast': ['SFRP2', 'SERPINF1', 'TAGLN'],
    'Alveolar fibroblast': ['GPC3', 'FGFR4'],
    'Deuterosomal cell': ['HES6'],
    'SMG duct cell': ['MGST1', 'KRT19'],
    'SMG serous cell': ['BPIFA1', 'DMBT1'],
    'Secretory cell': ['SERPINB3', 'S100A4', 'SCGB3A1']
}

sc.pl.umap(
    mislv2, color=['unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['unilung_celltype', 'original_celltype']
)
plt.savefig(dir+'/mislv2_umap.pdf', dpi=300, format='pdf', bbox_inches='tight')
sc.pl.dotplot(mislv2, marker, groupby='unilung_ann_level2')




mislv3 = adata4[adata4.obs['unilung_ann_level3'].isin(['Suprabasal cell'])]
sc.pp.neighbors(mislv3, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(mislv3, min_dist=0.3)
sc.tl.leiden(mislv3, resolution=0.2, key_added='mis_percent')

marker = {
    'Suprabasal cell': ['KRT19', 'KRT5', 'TP63']
}

sc.pl.umap(
    mislv3, color=['unilung_ann_level3', 'map_ann_level3'],
    ncols=1, title=['unilung_celltype', 'original_celltype']
)
plt.savefig(dir+'/mislv3_umap.pdf', dpi=300, format='pdf', bbox_inches='tight')

sc.pl.dotplot(mislv3, marker, groupby='unilung_ann_level3')
plt.savefig(dir+'/mislv3_dotplot.pdf', dpi=300, format='pdf', bbox_inches='tight')



## rare celltype (<3000 cells)
rarelv1 = adata4[adata4.obs['unilung_ann_level1'].isin(['Neuroendocrine cell', 'Submucosal gland cell'])]
rarelv2 = adata4[adata4.obs['unilung_ann_level2'].isin(['Deuterosomal cell', 'SMG duct cell', 'Megakaryocyte', 'Squamous cell'])]
rarelv3 = adata4[adata4.obs['unilung_ann_level3'].isin(['T helper cell', 'Promonocyte', 'Plasmablast cell', 'Mature B cell', 'Follicular B cell'])]
rarelv4 = adata4[adata4.obs['unilung_ann_level4'].isin(['Th1 cell'])]

rare = sc.concat([rarelv1, rarelv2, rarelv3])
sc.pp.neighbors(rare, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(rare, min_dist=0.3)
sc.pl.umap(
    rare, color=['final_celltype', 'map_celltype'],
    ncols=1, title=['unilung_celltype', 'original_celltype']
)
plt.savefig('rare_umap.pdf', dpi=300, format='pdf', bbox_inches='tight')


sc.pp.neighbors(rarelv1, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(rarelv1, min_dist=0.3)
sc.tl.rank_genes_groups(rarelv1, 'unilung_ann_level1', method='t-test')
marker = {
    'Neuroendocrine cell': ['CHGB', 'COL1A2'],
    'Submucosal gland cell': ['LTF', 'SLPI', 'LYZ']
}
sc.pl.dotplot(rarelv1, marker, groupby='unilung_ann_level1')



sc.pp.neighbors(rarelv2, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(rarelv2, min_dist=0.3)
sc.tl.rank_genes_groups(rarelv2, 'unilung_ann_level2', method='t-test')
sc.pl.umap(
    rarelv2, color=['unilung_ann_level2', 'map_ann_level2'],
    ncols=1, title=['unilung_celltype', 'original_celltype']
)
plt.savefig(dir+'/rarelv2_umap.pdf', dpi=300, format='pdf', bbox_inches='tight')

marker = {
    'Deuterosomal cell': ['HES6'],
    'Megakaryocyte': ['PPBP', 'PF4', 'ITGA2B'],
    'SMG duct cell': ['MGST1', 'KRT19', 'LTF'],
    'Squamous cell': ['SCEL', 'SPRR1B']
}
sc.pl.dotplot(rarelv2, marker, groupby='unilung_ann_level2')



sc.pp.neighbors(rarelv3, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(rarelv3, min_dist=0.3)
sc.tl.rank_genes_groups(rarelv3, 'unilung_ann_level3', method='t-test')
marker = {
    'T helper cell': ['CD3E', 'CD3D', 'IL7R'],
    'Promonocyte': ['MPO', 'VCAN', 'S100A8'],
    'Plasmablast cell': ['MKI67', 'CD38'],
    'Mature B cell': ['CD83', 'CXCR4'],
    'Follicular B cell': ['CD22', 'CD81', 'CD74']
}
sc.pl.dotplot(rarelv3, marker, groupby='unilung_ann_level3')
plt.savefig(dir+'/rarelv3_dotplot.pdf', dpi=300, format='pdf', bbox_inches='tight')


sc.pp.neighbors(rarelv4, use_rep="X_scANVI", n_neighbors=30)
sc.tl.umap(rarelv4, min_dist=0.3)
marker = {
    'Th1 cell': ['CXCR3', 'ITGA4', 'CD4']
}
sc.pl.dotplot(rarelv4, marker, groupby='unilung_ann_level4')
plt.savefig(dir+'/rarelv4_dotplot.pdf', dpi=300, format='pdf', bbox_inches='tight')

