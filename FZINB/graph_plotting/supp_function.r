library("ggplot2")
library("reshape2")
library("gplots")
Heatmap.matrix <- function(X , title = "", subtitle = "", clustering = FALSE, SAVE = FALSE, savepath = NULL){
  X <- data.matrix(X)
  if (SAVE == TRUE){
    png(savepath,    # create PNG for the heat map        
        width = 3*660,        # 5 x 300 pixels
        height = 3*660,
        res = 300,            # 300 pixels per inch
        pointsize = 8)
  }
  
  
  if (clustering == FALSE){
    colnames(X) <- 1:dim(X)[2]
    rownames(X) <- 1:dim(X)[1]
    
    melted_X <- melt(X)
    rangeF <- range(as.vector(X))
    col_breaks = c(seq(rangeF[1],rangeF[2],length=300))
    
    
    p <- ggplot(data = melted_X, aes(x=Var1, y=Var2, fill=value)) + 
      scale_fill_gradient(low = "blue", high = "yellow", space = "Lab" ,guide = "colourbar") +
      geom_tile() +
      labs(title = title,
           subtitle = subtitle,
          x =  "row of matrix",
           y ="col of matrix") + 
      coord_flip() +scale_x_reverse()
    
    return(p)
  }
  else{
    my_palette <- colorRampPalette(c("blue", "white" , "yellow"))(n = 299)
    rangeF <- range(as.vector(X))
    col_breaks = c(seq(rangeF[1],rangeF[2],length=300))
    

   heatmap.2(X,  
              main = title,
              #density.info="none",  # turns off density plot inside color legend
              dendrogram="row",
              trace="none",         # turns off trace lines inside the heat map
              margins =c(5,5),     # widens margins around plot
              col=my_palette,       # use on color palette defined earlier
              breaks = col_breaks)    # only draw a row dendrogram
  }
  if (SAVE == TRUE){
    dev.off()  
  }
}

#colMap <- colorRampPalette(c("blue","white","yellow" ))(256)
#image(Sig, axes = FALSE, asp=1,  col = colMap)


I0X_Genomics.load <- function(dataname, title = "", download = FALSE){
  #source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")
  library("cellrangerRkit")
  pipestance_path <- paste("./data/",dataname,sep = "")
  if (download){
    download_sample(sample_name=dataname,sample_dir=pipestance_path,
                    host="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/")
  }
  gbm <- load_cellranger_matrix(pipestance_path)
  
  use_genes <- get_nonzero_genes(gbm)
  gbm_bcnorm<- normalize_barcode_sums_to_median(gbm[use_genes,])
  data <- exprs(gbm_bcnorm) 
  
  return(data)
}

