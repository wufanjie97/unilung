#!/usr/bin/env python
# -*- coding: utf-8 -*-
import scanpy as sc
import pandas as pd
import numpy as np
import torch
import matplotlib.pyplot as plt
import os
import SEDR

random_seed = 2024
SEDR.fix_seed(random_seed)
dir = 'lung_cancer_sp'
os.chdir(dir)


def load_data_qc(cancertype, save_count=False):
    folders = [f for f in os.listdir(cancertype + '/E-MTAB-13530') if
               os.path.isdir(os.path.join(cancertype + '/E-MTAB-13530', f))]
    meta = pd.read_csv(cancertype + '/meta.csv')
    # set adata list and graph_dict
    adata_list = []
    graph_dict_combined = None
    for sample in folders:
        folder_path = os.path.join(cancertype + '/E-MTAB-13530', sample)
        adata_tmp = sc.read_visium(folder_path)
        adata_tmp.var_names_make_unique()
        ## transfer the adata.uns['spatial'].keys to sample name
        adata_tmp.uns['spatial'][sample] = adata_tmp.uns['spatial'].pop(", ".join(adata_tmp.uns['spatial'].keys()))
        sample_meta = meta[meta['sample'] == sample]
        if sample_meta.empty:
            raise ValueError(f"No matching sample found in meta.csv for sample: {sample}")
        adata_tmp.obs['sample'] = sample_meta.iloc[0, 0]
        adata_tmp.obs['age'] = sample_meta.iloc[0, 1]
        adata_tmp.obs['gender'] = sample_meta.iloc[0, 2]
        adata_tmp.obs['stage'] = sample_meta.iloc[0, 3]
        adata_tmp.obs['patient'] = sample_meta.iloc[0, 5]
        adata_tmp.obs['smoke'] = sample_meta.iloc[0, 6]
        # qc
        adata_tmp.layers['count'] = adata_tmp.X.toarray()
        adata_tmp.var["mt"] = adata_tmp.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata_tmp, qc_vars=["mt"], inplace=True)
        sc.pp.filter_cells(adata_tmp, max_counts=50000)
        sc.pp.filter_genes(adata_tmp, min_cells=50)
        adata_tmp = adata_tmp[adata_tmp.obs["pct_counts_mt"] < 20].copy()
        # build graph_dict
        graph_dict_tmp = SEDR.graph_construction(adata_tmp, 12)
        if graph_dict_combined is None:
            graph_dict_combined = graph_dict_tmp
        else:
            graph_dict_combined = SEDR.combine_graph_dict(graph_dict_combined, graph_dict_tmp)
        # combine adata list
        adata_list.append(adata_tmp)
    adata_combined = sc.concat(adata_list, join='outer', label="library_id", index_unique='-', uns_merge='unique')
    for library in adata_combined.uns['spatial']:
        adata_combined.uns['spatial'][library]['library_id'] = library
    adata_combined.write(cancertype + "/filter_sedr.h5ad")
    if save_count:
        # save count.csv for cancer_finder
        mat = adata_combined.to_df()
        mat = np.transpose(mat)
        mat_reset = mat.reset_index()
        mat_reset.columns = ['SYMBOL'] + mat.columns[0:].tolist()
        mat_reset.to_csv(cancertype + '_count.csv', index=False)
    return adata_combined, graph_dict_combined



def mclust_R(adata, n_clusters, use_rep, key_added='SEDR', random_seed=123456, modelNames='EEE'):
    np.random.seed(random_seed)
    import rpy2.robjects as robjects
    robjects.r.library("mclust")
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    r_random_seed = robjects.r['set.seed']
    r_random_seed(random_seed)
    rmclust = robjects.r['Mclust']
    # change adata.obsm[use_rep] from dataframe to numpy
    use_rep_data = adata.obsm[use_rep]
    if isinstance(use_rep_data, pd.DataFrame):
        use_rep_data = use_rep_data.to_numpy()
    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(use_rep_data), n_clusters, modelNames)
    mclust_res = np.array(res[-2])
    adata.obs[key_added] = mclust_res
    adata.obs[key_added] = adata.obs[key_added].astype('int')
    adata.obs[key_added] = adata.obs[key_added].astype('category')
    return adata




def run_all(cancertype, ncluster, save_count=False):
    import warnings
    warnings.filterwarnings('ignore')
    adata, graph_dict = load_data_qc(cancertype=cancertype,save_count=save_count)
    # run SEDR
    print("Start SEDR")
    device = 'cuda:2' if torch.cuda.is_available() else 'cpu'
    sedr_net = SEDR.Sedr(adata.obsm['X_pca'], graph_dict, mode='clustering', device=device)
    # train
    using_dec = False
    train_fn = sedr_net.train_with_dec if using_dec else sedr_net.train_without_dec
    train_fn()
    sedr_feat, *_ = sedr_net.process()
    adata.obsm['SEDR'] = sedr_feat
    # harmony + SEDR
    import harmonypy as hm
    ho = hm.run_harmony(adata.obsm['SEDR'], adata.obs[['patient']], ['patient'])
    res = pd.DataFrame(ho.Z_corr).T
    res_df = pd.DataFrame(data=res.values, columns=['X{}'.format(i + 1) for i in range(res.shape[1])],
                          index=adata.obs.index)
    adata.obsm['SEDR_harmony'] = res_df
    sc.pp.neighbors(adata, use_rep='SEDR_harmony')
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="leiden", resolution=0.5)
    mclust_R(adata, n_clusters=ncluster, use_rep='SEDR_harmony', key_added='SEDR')
    sc.pl.umap(
        adata, color=['patient', 'leiden', 'SEDR'],
        ncols=1, title=['Patients', 'Leiden', 'SEDR']
    )
    plt.savefig(cancertype + '/SEDR.pdf', dpi=300, bbox_inches='tight')
    adata.write_h5ad(cancertype + '/inte_SEDR.h5ad')
    return adata

luad = run_all(cancertype='luad', ncluster=5, save_count=True)
lusc = run_all(cancertype='lusc', ncluster=3, save_count=True)



def spatial_plot(adata, patient, cancertype, group):
    samples = list(set(s for s in adata.obs['sample'] if patient in s))
    fig, axs = plt.subplots(1, len(samples), figsize=(len(samples) * 12, 10))
    group_keys = {"SEDR": "SEDR", "cancer_finder": "cafind", "region": "Region"}
    if group == 'SEDR':
        colors = {1: '#FDE089', 2: '#D1615B', 3: '#CEE7E7', 4: '#B29ABE', 5: '#F59F9F'}
    elif group == 'region':
        colors = {'Normal': '#BC3C29FF', 'Tumor_like': '#0072B5FF', 'Tumor': '#FFDC91FF'}
    else:
        colors = {'F': '#7E66A4', 'T': '#D6C881'}
    for i, sample in enumerate(samples):
        ad = adata[adata.obs['sample'] == sample, :].copy()
        sc.pl.spatial(
            ad,
            img_key="hires",
            library_id=sample,
            color=[group],
            palette=colors,
            size=1.5,
            legend_fontsize="x-large",
            ax=axs.flatten()[i],
            show=False
        )
        axs.flatten()[i].set_title(samples[i])
    plt.tight_layout()
    plt.show()
    plt.savefig(f'{cancertype}_{patient}_{group_keys[group]}.pdf', dpi=300, bbox_inches='tight')


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
    plt.savefig(f'{cancertype}_{patient}_{marker}.pdf', dpi=300, bbox_inches='tight')


adata = sc.read_h5ad('luad/inte_SEDR.h5ad')
# python -u CancerFinder/infer.py --ckp=CancerFinder/st_pretrain_article.pkl --matrix=luad_count.csv --out=luad_ca_finder.csv --threshold=0.5
cafinder = pd.read_csv('luad_ca_finder.csv',index_col=0)
cafinder['predict_transformed'] = cafinder['predict'].map({1:'T',0:'F'})
adata.obs['cancer_finder'] = cafinder['predict_transformed']
# region define
adata.obs['region'] = adata.obs['cancer_finder']
adata.obs['region'] = adata.obs['region'].astype('str')
adata.obs.loc[adata.obs['SEDR'] == 1, 'region'] = 'Normal'
adata.obs.loc[adata.obs['region'] == 'F', 'region'] = 'Tumor_like'
adata.obs.loc[adata.obs['region'] == 'T', 'region'] = 'Tumor'
adata.write_h5ad('luad_latest.h5ad')
sc.pl.umap(
    adata, color=['patient', 'SEDR', 'cancer_finder','region'],frameon=False,
    ncols=4, title=['Patient_LUAD', 'SEDR_cluster', 'Cafinder', 'Region']
)
plt.savefig('luad_malignant.pdf', dpi=300, format='pdf', bbox_inches='tight')





luad = sc.read_h5ad('luad_latest.h5ad')
luad_raw = sc.read_h5ad('luad/filter_sedr.h5ad')
luad_raw.obs = luad.obs
luad_raw.obsm = luad.obsm
luad_raw.uns = luad.uns
luad_raw.obsp = luad.obsp
luad_raw.write('luad_latest_all.h5ad')
for patient in ['P10', 'P15', 'P16', 'P24', 'P25']:
    spatial_plot(adata=luad, patient=patient, cancertype='luad', group='SEDR')
    spatial_plot(adata=luad, patient=patient, cancertype='luad', group='region')
    spatial_plot(adata=luad, patient=patient, cancertype='luad', group='cancer_finder')
    spatial_featureplot(adata=luad_raw, patient=patient, cancertype='luad', marker='EPCAM')

luad_raw.obs = luad.obs
luad_raw.obsm = luad.obsm
sc.pl.umap(
    luad, color=['EPCAM'],frameon=False,
    ncols=1, title=['EPCAM'],cmap='magma_r'
)
plt.savefig('luad_malignant_EPCAM.pdf', dpi=300, format='pdf', bbox_inches='tight')




adata = sc.read_h5ad('lusc/inte_SEDR.h5ad')
# python -u CancerFinder/infer.py --ckp=CancerFinder/st_pretrain_article.pkl --matrix=lusc_count.csv --out=lusc_ca_finder.csv --threshold=0.5
cafinder = pd.read_csv('lusc_ca_finder.csv',index_col=0)
cafinder['predict_transformed'] = cafinder['predict'].map({1:'T',0:'F'})
adata.obs['cancer_finder'] = cafinder['predict_transformed']

# region define
adata.obs['region'] = adata.obs['cancer_finder']
adata.obs['region'] = adata.obs['region'].astype('str')
adata.obs.loc[adata.obs['SEDR'] == 3, 'region'] = 'Normal'
adata.obs.loc[adata.obs['SEDR'] == 2, 'region'] = 'Tumor'
adata.obs.loc[adata.obs['region'] == 'F', 'region'] = 'Normal'
adata.obs.loc[adata.obs['region'] == 'T', 'region'] = 'Tumor_like'
adata.write_h5ad('lusc_latest_all.h5ad')
sc.pl.umap(
    adata, color=['patient', 'SEDR', 'cancer_finder', 'region'],frameon=False,
    ncols=4, title=['Patient_LUSC', 'SEDR_cluster', 'Cafinder', 'Region']
)
plt.savefig('lusc_malignant.pdf', dpi=300, format='pdf', bbox_inches='tight')


lusc = sc.read_h5ad('lusc_latest.h5ad')
lusc_raw = sc.read_h5ad('lusc/filter_sedr.h5ad')
lusc_raw.obs = lusc.obs
lusc_raw.obsm = lusc.obsm
lusc_raw.uns = lusc.uns
lusc_raw.obsp = lusc.obsp
lusc_raw.write('lusc_latest_all.h5ad')
for patient in ['P11', 'P17', 'P19']:
    spatial_plot(adata=lusc, patient=patient, cancertype='lusc', group='SEDR')
    spatial_plot(adata=lusc, patient=patient, cancertype='lusc', group='region')
    spatial_plot(adata=lusc, patient=patient, cancertype='lusc', group='cancer_finder')
    spatial_featureplot(adata=lusc_raw, patient=patient, cancertype='lusc', marker='EPCAM')

lusc_raw.obs = lusc.obs
lusc_raw.obsm = lusc.obsm
sc.pl.umap(
    lusc_raw, color=['EPCAM'],frameon=False,
    ncols=1, title=['EPCAM'],cmap='magma_r'
)
plt.savefig('lusc_malignant_EPCAM.pdf', dpi=300, format='pdf', bbox_inches='tight')
