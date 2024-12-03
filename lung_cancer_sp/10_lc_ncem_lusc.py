#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import scanpy as sc
import pandas as pd
pd.set_option('display.max_columns', 100)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
import scvi
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
import decoupler as dc
import pickle

dir = 'lung_cancer_sp'
os.chdir(dir)

adata_vis = sc.read_h5ad("ncem/map_lusc/sp.h5ad")
adata_ref = sc.read_h5ad("ncem/map_lusc/sc.h5ad")
mod = cell2location.models.Cell2location.load("ncem/map_lusc", adata_vis)
adata_vis2 = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
)

### Estimate cell-type specific expression of every gene in the spatial data (needed for NCEM)
# Compute expected expression per cell type
expected_dict = mod.module.model.compute_expected_per_cell_type(
    mod.samples[f"post_sample_q05"], mod.adata_manager #mod.samples
)

# Add to anndata layers
for i, n in enumerate(mod.factor_names_):
    adata_vis2.layers[n] = expected_dict['mu'][i]

# Save anndata object with results
pickle.dump(adata_vis2, open('ncem/lusc_ncem_raw.pkl', 'wb'))


### subset WIKIpathway gene
adata_vis = pickle.load(open('ncem/lusc_ncem_raw.pkl', 'rb'))
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

#Select msigdb Wiki Pathways (WP) genes (6716 genes in total)
adata_vis.var["keep"] = adata_vis.var_names.isin(msigdb_WP_genes2)
adata_vis = adata_vis[:, adata_vis.var.keep]
pickle.dump(adata_vis, open('ncem/lusc_ncem_raw_wiki.pkl', 'wb'))


### subset high and low nsclc-like sclc sample
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

## select high sample: P17_T2, P19_T1; and low sample: P17_T1, P19_T2
adata_vis = pickle.load(open('ncem/lusc_ncem_raw_wiki.pkl', 'rb'))
adata_vis = subset_slide(adata_vis, samples=["P17_T2", "P19_T1", "P17_T1", "P19_T2"])
pickle.dump(adata_vis, open('ncem/lusc_ncem_raw_sub.pkl', 'wb'))
adata_vis.X = adata_vis.layers["count"].copy()
#Select all cell types
cell_types = np.unique(list(adata_vis.layers.keys())[:25])
cell_types = np.delete(cell_types,21) #Remove count finest
# Extract cell type abundances: using 5% quantile (representing confident cell abundance)
adata_vis.obs[adata_vis.uns['mod']['factor_names']] = adata_vis.obsm['q05_cell_abundance_w_sf']
#Extract cell type proportions
prop = adata_vis.obs[cell_types]
#Set samples
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
adata_ncem.write("ncem/lusc_ncem.h5ad")





import os
import ncem as nc
import scanpy as sc
import squidpy as sq
import pickle
from ncem.interpretation import InterpreterDeconvolution
from ncem.train import TrainModelLinearDeconvolution
from ncem.data import get_data_custom, customLoaderDeconvolution

dir = 'lung_cancer_sp'
os.chdir(dir)

high = ["P17_T2", "P19_T1"]
low = ["P17_T1", "P19_T2"]

adata = sc.read_h5ad("ncem/lusc_ncem.h5ad")
adata_high = adata[adata.obs['sample'].isin(high)]
adata_high.write("ncem/lusc_ncem_high.h5ad")
adata_low = adata[adata.obs['sample'].isin(low)]
adata_low.write("ncem/lusc_ncem_low.h5ad")

def run_ncem_couple(group, cancertype):
    import scanpy as sc
    import pickle
    import os
    from ncem.interpretation import InterpreterDeconvolution
    from ncem.train import TrainModelLinearDeconvolution
    from ncem.data import get_data_custom, customLoaderDeconvolution
    dir = '/public8/lilab/student/fjwu/btit1/case_ca/spatial'
    os.chdir(dir)
    adata = sc.read_h5ad(f"ncem/{cancertype}_ncem_{group}.h5ad")
    ### Initialize ncem model for deconvoluted Visium
    ncem_ip = InterpreterDeconvolution()
    ncem_ip.data = customLoaderDeconvolution(
        adata=adata, patient='sample', library_id='sample', radius=None,
    )
    get_data_custom(interpreter=ncem_ip, deconvolution=True)
    ### Type coupling analysis
    ncem_ip.get_sender_receiver_effects()
    ncem_ip.cell_names = ['Alveolar cell type 1', 'Alveolar cell type 2', 'B cell', 'Ciliated', 'Club', 'DC mature',
                          'Endothelial cell', 'Macrophage', 'Macrophage alveolar', 'Mast cell', 'Monocyte', 'NK cell',
                          'Neutrophils', 'Plasma cell', 'Stromal', 'T cell CD4', 'T cell CD8', 'T cell regulatory',
                          'Tumor cells', 'cDC1', 'cDC2', 'other', 'pDC', 'transitional club/AT2']
    pickle.dump(ncem_ip, open(f"ncem/{cancertype}_ncem_ip_{group}.pkl", 'wb'))
    # Type coupling analysis with at least 500 differentially expressed genes
    type_coupling = ncem_ip.type_coupling_analysis_circular(
        edge_attr='magnitude', figsize=(11, 11), edge_width_scale=0.6, text_space=1.3, de_genes_threshold=300,
        save="ncem", suffix=f"{cancertype}_{group}_couple.pdf"
    )



run_ncem_couple(group='high', cancertype='lusc')
run_ncem_couple(group='low', cancertype='lusc')





### lusc high score group
### Set sender and receiver cell types
ncem_ip = pickle.load(open('ncem/lusc_ncem_ip_high.pkl', 'rb'))
# the gene name in ncem_ip is ensemble, change to symbol
adata = sc.read_h5ad("ncem/map_lusc/sp.h5ad")
adata_sub = adata[:, ncem_ip.node_feature_names]
ncem_ip.data.celldata.var.index = adata_sub.var['Symbol'].tolist()
ncem_ip.node_feature_names = adata_sub.var['Symbol'].tolist()

sender = "Macrophage"
receiver = "cDC2"
effect_df = ncem_ip.sender_receiver_values(sender=f'{sender}', receiver=f'{receiver}')
#Write results
effect_df.to_csv("ncem/lusc_send_high_effect.csv", index=True)
#Select statistically and biologically meaningful genes
gene_subset = effect_df[effect_df['qvalue'] < 0.05]
gene_subset = gene_subset[~gene_subset.index.str.startswith(('RPL', 'RPS'))]
gene_subset = gene_subset[abs(gene_subset['fold change']) > 1.5].index

### Plot sender and receiver effect analysis results
ncem_ip.sender_effect(receiver=f'{receiver}', gene_subset=list(gene_subset), figsize=(12,12), save = "ncem/sender_receiver", suffix= "lusc_high_sender.pdf")

### Sender similarity analysis for cDC2 cells
ncem_ip.sender_similarity_analysis(receiver='cDC2', figsize=(12,12), save = "ncem/sender_receiver", suffix= "cDC2_lusc_high_similarity.pdf")





### lusc low score group
### Set sender and receiver cell types
ncem_ip = pickle.load(open('ncem/lusc_ncem_ip_low.pkl', 'rb'))
# the gene name in ncem_ip is ensemble, change to symbol
adata = sc.read_h5ad("ncem/map_lusc/sp.h5ad")
adata_sub = adata[:, ncem_ip.node_feature_names]
ncem_ip.data.celldata.var.index = adata_sub.var['Symbol'].tolist()
ncem_ip.node_feature_names = adata_sub.var['Symbol'].tolist()

sender = "Tumor cells"
receiver = "pDC"
effect_df = ncem_ip.sender_receiver_values(sender=f'{sender}', receiver=f'{receiver}')
#Write results
effect_df.to_csv("ncem/lusc_send_low_effect.csv", index=True)
#Select statistically and biologically meaningful genes
gene_subset = effect_df[effect_df['qvalue'] < 0.01]
gene_subset = gene_subset[~gene_subset.index.str.startswith(('RPL', 'RPS'))]
gene_subset = gene_subset[abs(gene_subset['fold change']) > 2].index

### Plot sender and receiver effect analysis results
ncem_ip.sender_effect(receiver=f'{receiver}', gene_subset=list(gene_subset), figsize=(12,12), save = "ncem/sender_receiver", suffix= "lusc_low_sender.pdf")

### Sender similarity analysis for cDC2
ncem_ip.sender_similarity_analysis(receiver='cDC2', figsize=(12,12), save = "ncem/sender_receiver", suffix= "cDC2_lusc_low_similarity.pdf")
