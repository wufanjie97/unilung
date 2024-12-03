#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2024/9/25 20:52

# @Author : Fanjie
"""
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import os
import warnings
warnings.filterwarnings('ignore')

sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3
dir = 'lung_cancer_sp'
os.chdir(dir)


nsclc_sclc = sc.read_h5ad('addmodulescore/nsclc_sclc.h5ad')
sc.tl.rank_genes_groups(nsclc_sclc, 'group2', method='t-test')
result = nsclc_sclc.uns['rank_genes_groups']
groups = result['names'].dtype.names
deg = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']})

marker = deg['NSCLC-like SCLC_n'].head(50).tolist()

luad = sc.read_h5ad('luad_latest_all.h5ad')
sc.tl.score_genes(luad, marker)

lusc = sc.read_h5ad('lusc_latest_all.h5ad')
sc.tl.score_genes(lusc, marker)


def spatial_featureplot(adata, patient, cancertype, marker):
    samples = list(set(s for s in adata.obs['sample'] if patient in s))
    fig, axs = plt.subplots(1, len(samples), figsize=(len(samples) * 12, 10))
    for i, sample in enumerate(samples):
        ad = adata[adata.obs['sample'] == sample, :].copy()
        sc.pl.spatial(
            ad,
            img_key="hires",
            library_id=sample,
            color=marker,
            size=1.5,
            legend_fontsize="x-large",
            ax=axs.flatten()[i],
            show=False
        )
        axs.flatten()[i].set_title(samples[i])
    plt.tight_layout()
    plt.show()
    plt.savefig(cancertype+'_'+patient+'_'+marker+'.pdf', dpi=300, bbox_inches='tight')


for patient in ['P10', 'P15', 'P16', 'P24', 'P25']:
    spatial_featureplot(adata=luad, patient=patient, cancertype='luad', marker='score')


for patient in ['P11', 'P17', 'P19']:
    spatial_featureplot(adata=lusc, patient=patient, cancertype='lusc', marker='score')



## DC-LAMP marker --- TLS
patient = 'P15'
marker = 'LAMP3'
samples = list(set(s for s in luad.obs['sample'] if patient in s))
fig, axs = plt.subplots(1, len(samples), figsize=(len(samples) * 12, 10))
for i, sample in enumerate(samples):
    ad = luad[luad.obs['sample'] == sample, :].copy()
    sc.pl.spatial(
        ad,
        img_key="hires",
        library_id=sample,
        color=marker,
        size=1.5,
        legend_fontsize="x-large",
        ax=axs.flatten()[i],
        show=False
    )
    axs.flatten()[i].set_title(samples[i])
plt.tight_layout()
plt.show()