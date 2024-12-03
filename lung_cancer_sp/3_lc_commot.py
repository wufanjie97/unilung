#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Created on : 2024/10/4 15:10

# @Author : Fanjie
"""
import matplotlib.pyplot as plt
import os
import commot as ct
import scanpy as sc
import numpy as np

dir = 'lung_cancer_sp'
os.chdir(dir)


def select_slide(adata, s, batch_key="sample"):
    slide = adata[adata.obs[batch_key].isin([s]), :].copy()
    s_keys = list(slide.uns["spatial"].keys())
    s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]
    slide.uns["spatial"] = {s_spatial: slide.uns["spatial"][s_spatial]}
    return slide


def run_commot(adata, sample, cancertype):
    ad = select_slide(adata=adata, s=sample)
    sc.pp.normalize_total(ad, inplace=True)
    sc.pp.log1p(ad)
    species = 'human'
    df_cellchat = ct.pp.ligand_receptor_database(species=species, signaling_type='Secreted Signaling', database='CellChat')
    df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, ad, min_cell_pct=0.05)
    df_cellchat_filtered.to_csv(f'commot/{cancertype}_{sample}_cellchat.csv', header=True, index=False)
    ct.tl.spatial_communication(ad, database_name='cellchat', df_ligrec=df_cellchat_filtered, dis_thr=500, heteromeric=True, pathway_sum=True)
    ad.write(f'commot/{cancertype}_{sample}_commot.h5ad')
    return ad, df_cellchat_filtered


def choose_pathwy(adata, pathway, cancertype, sample, save_or_show):
    ct.tl.communication_direction(adata, database_name='cellchat', pathway_name=pathway, k=5)
    ct.pl.plot_cell_communication(adata, database_name='cellchat', pathway_name=pathway,
                                  plot_method='grid', background_legend=True,
                                  scale=0.000008, ndsize=8, grid_density=0.4, summary='sender', background='image',
                                  clustering='leiden', cmap='Alphabet',
                                  normalize_v=True, normalize_v_quantile=0.995)
    if save_or_show == 'save':
        plt.savefig(f'commot/{cancertype}_{sample}_{pathway}.pdf', dpi=300, bbox_inches='tight')
    elif save_or_show == 'show':
        plt.show()



luad_raw = sc.read_h5ad('commot/luad_latest_all.h5ad')
adata_dis = luad_raw.copy()
adata_dis.var_names_make_unique()
luad_p10_t1, cellchat_p10_t1 = run_commot(adata=adata_dis, sample='P10_T1', cancertype='luad')
luad_p10_t4, cellchat_p10_t4 = run_commot(adata=adata_dis, sample='P10_T4', cancertype='luad')
luad_p24_t2, cellchat_p24_t2 = run_commot(adata=adata_dis, sample='P24_T2', cancertype='luad')


choose_pathwy(luad_p10_t1, 'MIF', 'luad', 'P10_T1', 'save')
choose_pathwy(luad_p10_t4, 'MIF', 'luad', 'P10_T4', 'save')
choose_pathwy(luad_p24_t2, 'MIF', 'luad', 'P24_T2', 'save')






lusc_raw = sc.read_h5ad('commot/lusc_latest_all.h5ad')
adata_dis = lusc_raw.copy()
adata_dis.var_names_make_unique()
lusc_p11_t1, cellchat_p11_t1 = run_commot(adata=adata_dis, sample='P11_T1', cancertype='lusc')
lusc_p17_t2, cellchat_p17_t2 = run_commot(adata=adata_dis, sample='P17_T2', cancertype='lusc')
lusc_p19_t1, cellchat_p19_t1 = run_commot(adata=adata_dis, sample='P19_T1', cancertype='lusc')


choose_pathwy(lusc_p11_t1, 'MK', 'lusc', 'P11_T1', 'save')
choose_pathwy(lusc_p17_t2, 'MK', 'lusc', 'P17_T2', 'save')
choose_pathwy(lusc_p19_t1, 'MK', 'lusc', 'P19_T1', 'save')



def choose_pathwy(adata, pathway, cancertype, sample):
    ct.tl.communication_direction(adata, database_name='cellchat', pathway_name=pathway, k=5)
    ct.pl.plot_cell_communication(adata, database_name='cellchat', pathway_name=pathway,
                                  plot_method='grid', background_legend=True,
                                  scale=0.000008, ndsize=8, grid_density=0.4, summary='sender', background='image',
                                  clustering='leiden', cmap='Alphabet',
                                  normalize_v=True, normalize_v_quantile=0.995)
    plt.show()



choose_pathwy(lusc_p11_t1, 'MK', 'lusc', 'P11_T1')
choose_pathwy(lusc_p17_t2, 'MK', 'lusc', 'P17_T2')
choose_pathwy(lusc_p19_t1, 'MK', 'lusc', 'P19_T1')




