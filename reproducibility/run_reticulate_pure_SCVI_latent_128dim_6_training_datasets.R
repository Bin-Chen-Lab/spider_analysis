library(Seurat)
library(dplyr)
library(anndata)

load("training_combined_6_datasets_RNA_count_file2_20230115.RData")

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Continue with SCANVI embedding:
#

SCVI_embed_ref <- function(anndata_filename = 'PBMC_GBM_BM_preprocessed_anndata.h5ad',
                             save_reference_latent_file = 'PBMC_GBM_BM_ref_RNA_latent_representations.csv',
                             ref_model_path = 'demo_scanvi/PBMC_GBM_BM/',
                             save_top_1000_HVGs_file = 'top_1000_HVGs_PBMC_GBM_BM.csv'
){
  RNA <- NormalizeData(RNA, normalization.method = "LogNormalize", scale.factor = 10000)
  RNA <- FindVariableFeatures(RNA, selection.method = "vst", assay = "RNA", nfeatures = 1000)
  top1000 <- head(VariableFeatures(RNA, assay = "RNA"), 1000)
  RNA <- RNA[top1000]
  print(RNA)
  
  library(reticulate)
  sc <- import('scanpy', convert = FALSE)
  scvi <- import('scvi', convert = FALSE)
  torch <- import('torch', convert = FALSE)
  #remove_sparsity <- import('scarches.dataset.trvae.data_handling', convert = FALSE)
  plt <- import('matplotlib.pyplot', convert = FALSE)
  np <- import('numpy', convert = FALSE)
  os <- import('os', convert = FALSE)
  scvi$settings$progress_bar_style = 'tqdm'
  adata <- sc$AnnData(
    X   = t(as.matrix(GetAssayData(RNA,slot='counts',assay = "RNA"))), #scVI requires raw counts
    obs = RNA[[]],
    var = GetAssay(RNA)[[]]
  )
  print(adata)
  
  #SCANVI model setting:
  sc$settings$set_figure_params(dpi=as.integer(200), frameon=FALSE)
  sc$set_figure_params(dpi=as.integer(200))
  sc$set_figure_params(figsize=c(4, 4))
  torch$set_printoptions(precision=as.integer(3), sci_mode=FALSE, edgeitems=as.integer(7))
  
  condition_key = 'study'
  cell_type_key = 'cell_type'
  
  #Create SCANVI model and train it on fully labelled reference dataset:
  scvi$data$setup_anndata(adata, batch_key=condition_key, labels_key=cell_type_key)
  
  vae = scvi$model$SCVI(
    adata,
    n_layers=as.integer(2),
    n_latent=as.integer(128),
    encode_covariates=TRUE,
    deeply_inject_covariates=FALSE,
    use_layer_norm="both",
    use_batch_norm="none")
  
  vae$train()
  
  #scanvae = sca$models$SCANVI$from_scvi_model(vae, "Unknown")
  #print(paste0("Labelled Indices: ", length(scanvae$`_labeled_indices`)))
  #print(paste0("Unlabelled Indices: ", length(scanvae$`_unlabeled_indices`)))
  
  #scanvae$train(max_epochs=as.integer(20))
  
  #Getting the latent representation and visualization
  reference_latent = vae$get_latent_representation()
  
  reference_latent = as.matrix(reference_latent)
  rownames(reference_latent) = colnames(RNA)
  
  RNA[['scvi']] <- CreateDimReducObject(embeddings = reference_latent, key = "scvi_", assay = DefaultAssay(RNA))
  # Find clusters, then run UMAP, and visualize
  RNA <- FindNeighbors(RNA, dims = 1:10, reduction = 'scvi')
  RNA <- FindClusters(RNA, resolution =1)
  
  RNA <- RunUMAP(RNA, dims = 1:10, reduction = 'scvi', n.components = 2)
  
  
  RNA$study[grep('GSM5008738', RNA$study)] = 'GSM5008738'
  
  
  vae$save(ref_model_path, overwrite=TRUE)
  
  write.csv(reference_latent, save_reference_latent_file)
  
  write.csv(as.character(py_to_r(adata$var$index)), save_top_1000_HVGs_file)
  
  write_h5ad(
    adata,
    anndata_filename,
    compression = NULL,
    compression_opts = NULL,
    as_dense = list()
  )
}

#Run function to embed ref data sets:
SCVI_embed_ref(anndata_filename = 'training_combined_6_datasets_RNA_preprocessed_anndata_pure_SCVI_latent_128dim_20240412.h5ad',
                 save_reference_latent_file = 'training_combined_6_datasets_RNA_ref_pure_SCVI_latent_128dim_representations_20240412.csv',
                  ref_model_path = 'demo_scvi/training_combined_6_datasets_RNA_pure_SCVI_latent_128dim_20240412/',
                 save_top_1000_HVGs_file = 'top_1000_HVGs_training_combined_6_datasets_RNA_pure_SCVI_latent_128dim_20240412.csv')



