#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import scanpy as sc
import pandas as pd
pd.set_option('display.max_columns', 100)
import numpy as np
import cell2location
import scvi
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
import seaborn as sns
import decoupler as dc
import pickle

dir = 'lung_cancer_sp'
os.chdir(dir)


# Get MSigDB resource
msigdb = dc.get_resource('MSigDB')
#Select WikiPathways databse
msigdb = msigdb[msigdb['collection']=='wikipathways']
msigdb = msigdb[~msigdb.duplicated(['geneset', 'genesymbol'])]
#Remove prefix
msigdb.loc[:, 'geneset'] = [name.split('WP_')[1] for name in msigdb['geneset']]
msigdb.to_csv("ncem/msigdb_WP.csv", index=True)
msigdb_WP = pd.read_csv("ncem/msigdb_WP.csv", index_col = 0)
# transfer the SYMBOL id fgene id of target data
gene_mapping = pd.read_csv('ncem/total_gene_list_43878.txt', sep='\t', header=None, names=['GeneSymbol', 'LocusGroup' ,'GeneID'])
gene_symbol_to_id = dict(zip(gene_mapping['GeneSymbol'], gene_mapping['GeneID']))
msigdb_WP_genes2 = msigdb_WP.genesymbol.map(gene_symbol_to_id).tolist()


def subset_slide(adata, samples, batch_key="sample"):
    if not isinstance(samples, list):
        samples = [samples]
    slide = adata[adata.obs[batch_key].isin(samples), :].copy()
    s_keys = list(slide.uns["spatial"].keys())
    matching_keys = [k for k in s_keys if any(sample in k for sample in samples)]
    if matching_keys:
        slide.uns["spatial"] = {k: slide.uns["spatial"][k] for k in matching_keys}
    else:
        raise ValueError("No matching spatial keys found for the provided samples.")
    return slide


def ncem_prepare(cancertype):
    adata_vis = sc.read_h5ad(f"ncem/map_{cancertype}/sp.h5ad")
    mod = cell2location.models.Cell2location.load(f"ncem/map_{cancertype}", adata_vis)
    adata_vis2 = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
    )

    ### Estimate cell-type specific expression of every gene in the spatial data (needed for NCEM)
    # Compute expected expression per cell type
    expected_dict = mod.module.model.compute_expected_per_cell_type(
        mod.samples[f"post_sample_q05"], mod.adata_manager  # mod.samples
    )

    # Add to anndata layers
    for i, n in enumerate(mod.factor_names_):
        adata_vis2.layers[n] = expected_dict['mu'][i]

    # Save anndata object with results
    pickle.dump(adata_vis2, open(f"ncem/{cancertype}_ncem_raw.pkl", 'wb'))

    # Select msigdb Wiki Pathways (WP) genes (6716 genes in total)
    adata_vis.var["keep"] = adata_vis.var_names.isin(msigdb_WP_genes2)
    # Select WikiPathways gene sets
    adata_vis = adata_vis[:, adata_vis.var.keep]
    pickle.dump(adata_vis, open(f"ncem/{cancertype}_ncem_raw_wiki.pkl", 'wb'))

    ## select high sample and low sample
    if cancertype == "luad":
        subset_samples = ["P10_T1", "P24_T2", "P24_T1", "P16_T1", "P16_T2", "P10_T4", "P15_T1", "P15_T2", "P25_T1", "P25_T2"]
    elif cancertype == "lusc":
        subset_samples = ["P17_T2", "P19_T1", "P17_T1", "P19_T2"]
    adata_vis = subset_slide(adata_vis, samples=subset_samples)
    pickle.dump(adata_vis, open(f"ncem/{cancertype}_ncem_raw_sub.pkl", 'wb'))
    adata_vis.X = adata_vis.layers["count"].copy()
    # Select all cell types
    cell_types = np.unique(list(adata_vis.layers.keys())[:25])
    cell_types = np.delete(cell_types, 21)  # Remove count finest
    # Extract cell type abundances: using 5% quantile (representing confident cell abundance)
    adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
    # Extract cell type proportions
    prop = adata_vis.obs[cell_types]
    # Set samples
    samples = adata_vis.obs["sample"]

    ### collect NCEM data
    # Prepare data for NCEM
    cell_expression = []
    node_types = []
    proportions = []
    spatial = []
    sample = []
    for i, ct in enumerate(cell_types):
        proportions.append(prop)
        cell_expression.append(adata_vis.layers[ct].toarray())
        nt = np.zeros((prop.shape[0], len(cell_types)))
        nt[:, i] = 1
        node_types.append(nt)
        spatial.append(adata_vis.obsm['spatial'])
        sample.append(adata_vis.obs['sample'])

    proportions = pd.DataFrame(np.concatenate(proportions), columns=cell_types)
    cell_expression = pd.DataFrame(np.concatenate(cell_expression), columns=adata_vis.var_names)
    node_types = pd.DataFrame(np.concatenate(node_types), columns=cell_types)
    spatial = pd.DataFrame(np.concatenate(spatial))
    sample = pd.DataFrame(np.concatenate(sample))

    from anndata import AnnData
    adata_ncem = AnnData(cell_expression)
    adata_ncem.obsm['proportions'] = np.array(proportions)
    adata_ncem.obsm['node_types'] = np.array(node_types)
    adata_ncem.obsm['spatial'] = np.array(spatial)
    if np.array(sample).ndim == 2:
        new_sample = np.array(sample).flatten()
    else:
        new_sample = np.array(sample)
    adata_ncem.obs['sample'] = new_sample
    adata_ncem.uns["node_type_names"] = {x: x for x in cell_types}

    # process NCEM data
    sc.pp.filter_genes(adata_ncem, min_cells=0)
    adata_ncem.obsm['sample'] = np.array(sample)
    adata_ncem.layers["Cell_expression"] = adata_ncem.X
    sc.pp.log1p(adata_ncem)
    h_0 = pd.DataFrame(adata_ncem.obsm['node_types'], columns=list(adata_ncem.uns['node_type_names'].values()))
    target_type = pd.DataFrame(np.array(h_0.idxmax(axis=1)), columns=["target_cell"]).reset_index()
    adata_ncem.obs = target_type
    adata_ncem.obs['sample'] = new_sample
    # save NCEM data
    adata_ncem.write(f"ncem/{cancertype}_ncem.h5ad")


ncem_prepare('luad')
ncem_prepare('lusc')
