I0X_Genomics.load <- function(dataname, title = "", download = FALSE){
  #source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")
  library("cellrangerRkit")
  pipestance_path <- paste("./data/",dataname,sep = "")
  if (download){
    download_sample(sample_name=dataname,sample_dir=pipestance_path,
                    host="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/")
  }
  gbm <- load_cellranger_matrix(pipestance_path)
  
  #data <- exprs(gbm) 
  # expression matrix #raw UMI counts
  #genes <- fData(gbm)
  # data frame of genes
  #pData(gbm)
  # data frame of cell barcodes
  
  #analysis_results <- load_cellranger_analysis_results(pipestance_path)
  #tsne_proj <- analysis_results$tsne
  #visualize_umi_counts(gbm,tsne_proj[c("TSNE.1","TSNE.2")],limits=c(3,4),marker_size=0.05)
  
  use_genes <- get_nonzero_genes(gbm)
  gbm_bcnorm<- normalize_barcode_sums_to_median(gbm[use_genes,])
  data <- exprs(gbm_bcnorm) 
  
  return(data)
}