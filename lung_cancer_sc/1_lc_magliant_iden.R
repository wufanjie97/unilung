library(reticulate)
library(Seurat)
library(tidyverse)
library(harmony)
library(cowplot)
library(ggsci)
library(MuDataSeurat)
library(infercnv)
library(AnnoProbe)
library(future)

options(future.globals.maxSize = 5*1024^3)
options(scipen = 100)
options("Seurat.object.assay.version" = "v3")

dir <- "lung_cancer_sc"
cnv_dir <- file.path(dir, "infercnv")
setwd(dir)

use_condaenv("/envs/scvi")
scanpy <- import("scanpy")
pd <- import("pandas")
np <- import("numpy")


process_sc_data <- function() {
  if (!file.exists(file.path(dir, "epi_raw.rds"))) {
    in_ad <- scanpy$read_h5ad(file.path(dir, "epi_raw.h5ad"))
    
    count_matrix <- np$transpose(in_ad$X)
    dimnames(count_matrix) <- list(rownames(in_ad$var), rownames(in_ad$obs))
    
    seu <- CreateSeuratObject(
      counts = count_matrix,
      meta.data = as.data.frame(in_ad$obs)
    )
    saveRDS(seu, file.path(dir, "epi_raw.rds"))
  } else {
    seu <- readRDS(file.path(dir, "epi_raw.rds"))
  }
  
  processing_pipeline <- function(obj) {
    obj %>%
      NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
      ScaleData() %>%
      RunPCA(npcs = 30, verbose = FALSE) %>%
      RunHarmony(group.by.vars = "donor_ID", max.iter.harmony = 10) %>%
      RunUMAP(reduction = "harmony", dims = 1:30) %>%
      FindNeighbors(reduction = "harmony", dims = 1:30) %>%
      FindClusters(resolution = 0.1)
  }
  
  sce <- processing_pipeline(seu)
  
  cluster_mapping <- c(
    "0" = "Type II alveolar cell", "1" = "Basal cell",
    "2" = "Malignant cell", "3" = "Malignant cell",
    "4" = "Malignant cell", "5" = "Malignant cell",
    "6" = "Ciliated cell", "7" = "Type I alveolar cell"
  )
  sce <- RenameIdents(sce, cluster_mapping)
  sce$annotation <- Idents(sce)
  
  generate_plots(sce)
  saveRDS(sce, file.path(dir, "epi.rds"))
  WriteH5AD(sce, 'epi.h5ad', assay = "RNA")
}


generate_plots <- function(sce) {
  # DotPlot
  marker_genes <- c(
    "KRT13", "KRT17", "S100A2",    # Basal
    "CAV1", "AGER",                # AT1
    "SFTPC", "SFTPA1", "SFTPA2",   # AT2
    "SCGB1A1", "SCGB3A1",          # Club
    "TPPP3", "FOXJ1", "PIFO",      # Ciliated
    "EPCAM"                        # Cancer
  )
  
  dot_theme <- theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_blank()
  )
  
  dot_plot <- DotPlot(sce, features = marker_genes, group.by = "annotation") +
    dot_theme +
    guides(size = guide_legend("Percent Expression")) +
    scale_color_gradientn(colors = c('#330066','#336699','#66CC66','#FFCC33'))
  
  # UMAP
  umap_theme <- theme(plot.title = element_text(hjust = 0.5))
  
  p1 <- DimPlot(sce, reduction = "umap", cols = pal_nejm("default")(6),
                group.by = "annotation") +
    ggtitle("Cell type") +
    umap_theme
  
  p2 <- DimPlot(sce, reduction = "umap", 
                cols = c("#B24745FF","#79AF97FF","#7B4173FF"),
                group.by = "sub_atlas") +
    ggtitle("Cancer type") +
    umap_theme
  
  plot_grid(p1, p2, ncol = 2)
}







## infercnv
prepare_cnv_data <- function() {
  if (!file.exists(file.path(cnv_dir, "epidata.RData"))) {
    ref_epi <- scanpy$read_h5ad(file.path(dir, "ref_epi10k.h5ad"))
    epi_dat <- np$transpose(ref_epi$X)
    dimnames(epi_dat) <- list(rownames(ref_epi$var), rownames(ref_epi$obs))
    save(sce, epi_dat, file = file.path(cnv_dir, "epidata.RData"))
  }
}

run_cnv <- function(times, clusters) {
  if (missing(times) | missing(clusters)) {
    stop("Both times and clusters parameters are required")
  }
  
  load(file.path(cnv_dir, "epidata.RData"))
  gse <- sce[, sce$seurat_clusters %in% clusters]
  
  common_genes <- intersect(rownames(gse), rownames(epi_dat))
  combined_data <- cbind(
    GetAssayData(gse, "counts")[common_genes, ],
    epi_dat[common_genes, ]
  )
  
  annotations <- data.frame(
    cell_id = colnames(combined_data),
    group = c(gse$seurat_clusters, rep('ref-epithelial', ncol(epi_dat)))
  )

  gene_info <- annoGene(common_genes, "SYMBOL", 'human') %>%
    arrange(chr, start) %>%
    distinct(symbol, .keep_all = TRUE)
  
  output_prefix <- file.path(cnv_dir, paste0("cnv", times))
  dir.create(output_prefix, showWarnings = FALSE, recursive = TRUE)
  
  write.table(combined_data, file.path(output_prefix, "exp_matrix.tsv"), 
              sep = "\t", quote = FALSE)
  write.table(annotations, file.path(output_prefix, "annotations.tsv"),
              sep = "\t", quote = FALSE, col.names = FALSE)
  write.table(gene_info, file.path(output_prefix, "gene_positions.tsv"),
              sep = "\t", quote = FALSE, col.names = FALSE)
  

  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = file.path(output_prefix, "exp_matrix.tsv"),
    annotations_file = file.path(output_prefix, "annotations.tsv"),
    gene_order_file = file.path(output_prefix, "gene_positions.tsv"),
    delim = "\t",
    ref_group_names = 'ref-epithelial'
  )
  
  plan("multicore", workers = 8)
  infercnv_results <- run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = output_prefix,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = FALSE,
    output_format = "pdf"
  )
  
  saveRDS(infercnv_results, file.path(output_prefix, "results.rds"))
}



process_sc_data()
prepare_cnv_data()
## Due to the scale of data, We split it into three subsets according to the clustering results and run them separately
run_cnv(times = 1, cluster = c(0,16:34))
run_cnv(times = 2, cluster = c(1:3,13:15))
run_cnv(times = 3, cluster = c(4:12))





## CNV score
library(infercnv)
library(Seurat)
library(ggplot2)
library(data.table)
library(tidyverse)
library(ggsci)
library(CytoTRACE)
library(Matrix)


dir <- "lung_cancer_sc"
cnv_dir <- file.path(dir, "infercnv")
output_dir <- file.path(dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

future::plan("multicore", workers = 8)
options(future.globals.maxSize = 10 * 1024^3)


process_cnv_analysis <- function() {
  cnv_data_path <- file.path(cnv_dir, "cnv_all/run.final.infercnv_obj")
  group_file <- file.path(cnv_dir, "groupfile.txt")
  
  if (!file.exists(cnv_data_path)) stop("CNV data file not found")
  if (!file.exists(group_file)) stop("Group file not found")
  
  cnv_data <- readRDS(cnv_data_path)
  dat <- cnv_data@expr.data
  
  cnv_score <- function(mat) {
    scaled_mat <- t(scale(Matrix::t(mat)))
    scaled_mat <- scales::rescale(scaled_mat, to = c(-1, 1))
    data.frame(
      cellid = colnames(mat),
      score = colSums(scaled_mat^2, na.rm = TRUE),
      row.names = NULL
    )
  }
  
  scores <- cnv_score(dat)
  
  groups <- fread(group_file, header = FALSE, col.names = c("cellid", "cluster"))

  combined_data <- scores %>%
    left_join(groups, by = "cellid") %>%
    filter(!is.na(cluster)) 
  
  generate_cnv_plot <- function(data) {
    ggplot(data, aes(x = cluster, y = score, fill = cluster)) +
      geom_violin(alpha = 0.8, show.legend = FALSE) +
      scale_fill_igv() +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm")
      ) +
      labs(x = "Cell Type", y = "CNV Score")
  }
  
  cnv_plot <- generate_cnv_plot(combined_data)
  ggsave(file.path(output_dir, "cnv_score.pdf"), cnv_plot,
         width = 15, height = 10, units = "cm", dpi = 300)
}


cnvScore <- function(data){
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    rescale(to=c(-1, 1)) %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}



# CytoTRACE
run_cytotrace_analysis <- function() {
  sce_path <- file.path(dir, "epi.rds")
  if (!file.exists(sce_path)) stop("Seurat object not found")
  sce <- readRDS(sce_path)
  
  sce <- FindVariableFeatures(
    sce,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  
  sketch_data <- function(obj) {
    SketchData(
      object = obj,
      ncells = min(1e5, ncol(obj)), 
      method = "LeverageScore",
      sketched.assay = "bpcell"
    )
  }
  
  sce_sketch <- sketch_data(sce)
  
  process_expression <- function(obj) {
    counts <- GetAssayData(obj, "counts", "RNA")
    Matrix::Matrix(counts, sparse = TRUE)
  }
  
  mat <- process_expression(sce_sketch)
  
  run_cytotrace <- function(expr_mat) {
    CytoTRACE(
      mat = expr_mat,
      ncores = 8,
      subsample = FALSE
    )
  }
  
  cytotrace_res <- run_cytotrace(mat)
  
  generate_cytotrace_plots <- function(res, metadata, reduction = "umap") {
    emb <- Embeddings(metadata, reduction)
    
    plot_features <- list(
      list(var = "annotation", prefix = "anno"),
      list(var = "sub_atlas", prefix = "atlas")
    )
    
    purrr::walk(plot_features, function(feat) {
      phenotype <- setNames(as.character(metadata[[feat$var]]), colnames(metadata))
      plotCytoTRACE(
        res,
        phenotype = phenotype,
        emb = emb,
        outputDir = file.path(output_dir, "cytotrace"),
        prefix = feat$prefix
      )
    })
  }
  
  generate_cytotrace_plots(cytotrace_res, sce_sketch)
}



process_cnv_analysis()
run_cytotrace_analysis()
