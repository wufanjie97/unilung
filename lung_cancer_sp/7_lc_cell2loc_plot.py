#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2024/9/27 18:42

# @Author : Fanjie
"""
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

def select_slide(adata, s, batch_key="sample"):
    slide = adata[adata.obs[batch_key].isin([s]), :].copy()
    s_keys = list(slide.uns["spatial"].keys())
    s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]
    slide.uns["spatial"] = {s_spatial: slide.uns["spatial"][s_spatial]}
    return slide

def cell2loc_slide_plot(spatial_adata, adata_ref, sample, cancertype, choose_celltype=None, TLS='cell2loc'):
    # Select one slide
    from cell2location.utils import select_slide
    slide = select_slide(spatial_adata, sample)
    # Determine cell type
    if choose_celltype is not None:
        celltype = choose_celltype
    else:
        celltype = adata_ref.obs.cell_type_major.astype('str').unique().tolist()
    # Plot in spatial coordinates
    with mpl.rc_context({'axes.facecolor': 'black',
                         'figure.figsize': [4.5, 5]}):
        sc.pl.spatial(slide, cmap='magma',
                      # Show cell types
                      color=celltype,
                      ncols=4, size=1.3,
                      img_key='hires',
                      # Limit color scale at 99.2% quantile of cell abundance
                      vmin=0, vmax='p99.2'
                      )
    plt.savefig(f'cell2loc/{cancertype}_{sample}_{TLS}.pdf', dpi=300, format='pdf', bbox_inches='tight')
    plt.close()

# load luad adata_vis and model
adata_vis = sc.read_h5ad("cell2loc/map_luad/sp.h5ad")
adata_ref = sc.read_h5ad("cell2loc/map_luad/sc.h5ad")
mod = cell2location.models.Cell2location.load("cell2loc/map_luad", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# all celltypes
for sample in adata_vis.obs['sample'].unique().tolist():
    cell2loc_slide_plot(spatial_adata=adata_vis,
                        adata_ref=adata_ref,
                        cancertype='luad',
                        sample=sample,
                        choose_celltype=None,
                        TLS='cell2loc')

# select TLS celltypes (the diff between high and low NSCLC-like are TILs). show in P10_T1 and T4
choose_celltype = ['cDC1', 'cDC2', 'pDC', 'DC mature', 'T cell CD4', 'T cell CD8', 'B cell', 'Plasma cell', 'Tumor cells']
for sample in adata_vis.obs['sample'].unique().tolist():
    cell2loc_slide_plot(spatial_adata=adata_vis,
                        adata_ref=adata_ref,
                        cancertype='luad',
                        sample=sample,
                        choose_celltype=choose_celltype,
                        TLS='TLS')




# load lusc adata_vis and model
adata_vis = sc.read_h5ad("cell2loc/map_lusc/sp.h5ad")
adata_ref = sc.read_h5ad("cell2loc/map_lusc/sc.h5ad")
mod = cell2location.models.Cell2location.load("cell2loc/map_lusc", adata_vis)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']

# all celltypes
for sample in adata_vis.obs['sample'].unique().tolist():
    cell2loc_slide_plot(spatial_adata=adata_vis,
                        adata_ref=adata_ref,
                        cancertype='lusc',
                        sample=sample,
                        choose_celltype=None,
                        TLS='cell2loc')

# select TLS celltypes (the diff between high and low NSCLC-like are 'cDC1', 'cDC2', 'pDC', 'transitional club/AT2'). show in P17
choose_celltype = ['cDC1', 'cDC2', 'pDC', 'DC mature', 'T cell CD4', 'T cell CD8', 'B cell', 'Plasma cell', 'transitional club/AT2', 'Tumor cells']
for sample in adata_vis.obs['sample'].unique().tolist():
    cell2loc_slide_plot(spatial_adata=adata_vis,
                        adata_ref=adata_ref,
                        cancertype='lusc',
                        sample=sample,
                        choose_celltype=choose_celltype,
                        TLS='TLS')




# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial
from cell2location.utils import select_slide

# select up to 6 clusters
clust_labels = ['cDC1', 'cDC2', 'pDC', 'DC mature', 'T cell CD4',
                'T cell CD8', 'B cell']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = select_slide(adata_vis, 'P10_T1')

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=1,
        # size of locations (adjust depending on figure size)
        circle_diameter=6,
        colorbar_position='right'
    )

plt.show()
plt.savefig('kk.pdf', dpi=300, bbox_inches='tight')

