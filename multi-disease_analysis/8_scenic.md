```python
vim build_loom.py
import os, sys
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp;
import numpy as np;
import scanpy as sc;
celltype = sys.argv[1]
counttype = sys.argv[2]

x = sc.read_csv(celltype + "_" + counttype + ".csv");
row_attrs = {"Gene": np.array(x.obs_names),};
col_attrs = {"CellID": np.array(x.var_names)};
lp.create(celltype + "_" + counttype + ".loom", x.X, row_attrs,col_attrs);
```

```shell
vim run_scenic.sh

python build_loom.py $1 $2

dir=/scenic_data/
tfs=$dir/hs_hgnc_tfs.txt
feather=$dir/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
tbl=$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
celltype=$1
counttype=$2

## pyscenic GRN
pyscenic grn \
--num_workers 15 \
--output ${celltype}_${counttype}_grn.tsv \
--method grnboost2 \
${celltype}_${counttype}.loom \
$tfs

pyscenic ctx \
${celltype}_${counttype}_grn.tsv $feather \
--annotations_fname $tbl \
--expression_mtx_fname ${celltype}_${counttype}.loom \
--mode "dask_multiprocessing" \
--output ${celltype}_${counttype}_reg.csv \
--num_workers 10 \
--mask_dropouts

pyscenic aucell \
${celltype}_${counttype}.loom \
${celltype}_${counttype}_reg.csv \
--output ${celltype}_${counttype}_SCENIC.loom \
--num_workers 10


bash run_scenic.sh bcell sub_latest
bash run_scenic.sh mono sub_latest
```

```R
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)


celltype <- 'mono'

loom <- open_loom(paste0(celltype,'_sub_latest_SCENIC.loom'))
gse <- readRDS(paste0(celltype,'_sub_latest.rds'))

# Read information from loom file: 
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC") 
regulonAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)  
close_loom(loom)

# intergrate the seurat object and scenic result
sub_regulonAUC <- regulonAUC
colnames(sub_regulonAUC) <- colnames(gse)

anno <- paste0(celltype,"_anno") 
cellinfo <- gse@meta.data[,c('donor_status',anno)]
cellTypes <-  as.data.frame(subset(cellinfo,select = anno))
status <-  as.data.frame(subset(cellinfo,select = 'donor_status'))

save(sub_regulonAUC, cellTypes, gse, file = paste0(celltype,'_tf_scenic.rdata'))
```

