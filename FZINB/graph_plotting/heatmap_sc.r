source("packageimport.r")
#source("supp_function.r")
library(gplots)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(ggpubr)

### scRNA-seq ######
load(file="../data/Buettner.RData")
load(file="../data/Kolodziejczyk.RData")
load(file="../data/Pollen.RData")
load(file="../data/Usoskin.RData")

Test <- list(Buettner , Kolodziejczyk, Pollen, Usoskin)
n_dataset <- length(Test)

#Load the result from RDATA

FK <- list()
GK <- list()
maxK <- 51
#compute FZINB kernel
# kernel scale 
for (n in 1:n_dataset){
  print(n)
  data <- Test[[n]]$in_X
  n_gene <- nrow(data)
  n_cell <- ncol(data)
  X10 <-round(10^data-1.0)
  FK[[n]] <-FZINB.matrix(X10,cores.ratio = 0.8)
  GK[[n]] <- Gaussian.Euclidean.kernel(t(data))[[maxK]]
}

#kernel matrix scale 
Gaussian <- list()
FZINB <- list()
for (n in 1:n_dataset){
  #result <- SIMLR(X=Test[[n]]$in_X,c=Test[[n]]$n_clust)
  #maxK <- max.alphaK(result$alphaK)
  #Gaussian[[n]] <- multiple.kernel.standard(x = t(Test[[n]]$in_X))[[maxK]]
  D_Gaussian[[n]] <- multiple.kernel(x = t(Test[[n]]$in_X))[[maxK]]
  D_FZINB[[n]] <- Gaussian.FZINB.kernel(x = t(Test[[n]]$in_X), Fisher_matrix = FK[[n]])[[maxK]]
}

# kernel scale 
Gaussian <- GK
FZINB <- FK



name <- c("Buettner" , "Kolodziejczyk", "Pollen", "Usoskin")
p <- list()
nmi <- list()
for (n in 1:n_dataset){
  n_clusters <- Test[[n]]$n_clust
  Kmeans_FZINB <- kmeans(FZINB[[n]], n_clusters, nstart = 30)$cluster
  Kmeans_Gaussian <- kmeans(Gaussian[[n]], n_clusters, nstart = 30)$cluster
  true_clusters <- Test[[n]]$true_labs$V1
  
  nmi[[n]] <- compare(Kmeans_Gaussian,true_clusters,method = "nmi")
  nmi[[n+4]] <- compare(Kmeans_FZINB,true_clusters,method = "nmi")
  
  clustered <- sort(true_clusters,index.return = TRUE)$ix
  
  melted_FZINB <- melt(as.matrix(FZINB[[n]])[clustered,clustered])
  melted_Gaussian <- melt(as.matrix(Gaussian[[n]])[clustered,clustered])
  pp<-0
  for (melted in list(melted_Gaussian, melted_FZINB)){
    limits <- range(quantile(melted$value,c(0,0.99)))
    
    p[[n+pp]] <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile() +
      scale_fill_gradient(low = "blue",high = "yellow", space = "Lab" ,guide = "colourbar", limits=limits)+
      labs(y ="",
           caption = paste("NMI ",signif(nmi[[n+pp]],3) ) )+
      #theme_void()+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            aspect.ratio=1,
            legend.title=element_blank())+
      scale_y_continuous(trans = "reverse")
    
    pp <-4
  }
}
p[[1]]<- p[[1]] +  labs(y = "Gaussian Kernel")
p[[5]]<- p[[5]] +  labs(y = "FZINB Kernel")

plot(p[[5]])

png(paste("../ZINB/Pics/kernel_compare/kernel_compare.png"),    # create PNG for the heat map        
    width = 6*660/2,        # 5 x 300 pixels
    height = 2.3*660/2,
    res = 300/2.5,            # 300 pixels per inch
    pointsize = 8)
grid.arrange(
  arrangeGrob(p[[1]],p[[5]],top = textGrob(name[1], gp=gpar(fontsize=12,font=8))),
  arrangeGrob(p[[2]],p[[6]],top = textGrob(name[2], gp=gpar(fontsize=12,font=8))),
  arrangeGrob(p[[3]],p[[7]],top = textGrob(name[3], gp=gpar(fontsize=12,font=8))),
  arrangeGrob(p[[4]],p[[8]],top = textGrob(name[4], gp=gpar(fontsize=12,font=8))),
  ncol=4)
dev.off()


### FZINB application ####