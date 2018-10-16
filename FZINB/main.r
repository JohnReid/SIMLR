source("packageimport.r")
dyn.load("../src/projsplx_R.so")

load(file="../data/Buettner.RData")
load(file="../data/Kolodziejczyk.RData")
load(file="../data/Pollen.RData")
load(file="../data/Usoskin.RData")

Test <- list(Buettner , Kolodziejczyk, Pollen, Usoskin)
n_dataset <- length(Test)

RESULT.full <- list()
RESULT.trun <- list()

# To compute the result of modified SIMLR INDIVIDUALLY by involving various of kernels
result <- SIMLR(X = Test[[1]]$in_X, c = Test[[1]]$n_clust, cores.ratio = 0.8,
                      IncludeGaussian = TRUE, IncludeFZINBvariants = "BOTH")
# Note that IncludeFZINBvariants = c("FZINB_kernel","Gaussian_kernel.FZINB_distance","BOTH","NONE")


# To compute the result of SIMLR BULKLY for all datasets and all methods
for (n in 1:n_dataset){
  RESULT.full[[n]] <-  SIMLR.COMBINE(X=Test[[n]]$in_X,c=Test[[n]]$n_clust, cores.ratio = 0.8) 
}

for (n in 1:n_dataset){
  X_counts <-round(10^Test[[n]]$in_X-1.0)
  Theta0 <- prior.zinb(X_counts)$Theta0
  THETA <- MLE.zinb(X_counts,Theta0)
  extract <- which(THETA[,3]>0.3 & THETA[,3]<0.9)
  RESULT.trun[[n]] <- SIMLR.COMBINE(X=Test[[n]]$in_X[extract,],c=Test[[n]]$n_clust)
}

# The analysis of the results based on the BULK computation by SIMLR.COMBINE
RESULT <- RESULT.trun
RESULT <- RESULT.full
# The result run by SIMLR.COMBINE is stored in the folder FZINB/results/RESULT_main.RDATA


###  RESULT ANALYSIS  ####
library(gplots)
library(reshape2)
library(gridExtra)
library(grid)
library(ggpubr)

  #### Heatmap for similarity matrix ###
p_S <- list()
dataname <- c("Buettner" , "Kolodziejczyk", "Pollen", "Usoskin")
nmi <-list()

for (n in 1:n_dataset){
  true_clusters <- Test[[n]]$true_labs$V1
  clustered <- sort(true_clusters,index.return = TRUE)$ix

  pp<-0
  for (result in list(RESULT[[n]]$results$SIMLR,RESULT[[n]]$results$FZINB)){
    S <- result$S
    melted <- melt(S[clustered,clustered])
    
    limits <- as.vector(quantile(result$S,c(0,.95)))
    nmi <- compare(true_clusters,result$y$cluster, method = "nmi")
    
    p_S[[n+pp]] <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile() +
      scale_fill_gradient(low = "blue",high = "yellow", space = "Lab" ,guide = "colourbar", limits=limits)+
      labs(y ="",
           caption = paste("NMI ",signif(nmi,3) ) )+
      #theme_void()+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            aspect.ratio=1,
            legend.title=element_blank())+
      scale_y_continuous(trans = "reverse")
    
    pp <-pp+n_dataset
  }
}
  
p_S[[1]]<- p_S[[1]] +  labs(y =expression(paste("Gaussian kernel ", (P[1]))))
p_S[[5]]<- p_S[[5]] +  labs(y =expression(paste("FZINB kernel ", (P[2]))))

plot(p_S[[1]])

png(paste("../FZINB/Pics/S/similarity_compare.png"),    # create PNG for the heat map        
    width = 6*660/2,        # 5 x 300 pixels
    height = 2.3*660/2,
    res = 300/2.5,            # 300 pixels per inch
    pointsize = 10)
grid.arrange(
  arrangeGrob(p_S[[1]],p_S[[5]],top = textGrob(dataname[[1]], gp=gpar(fontsize=12,font=8))),
  arrangeGrob(p_S[[2]],p_S[[6]],top = textGrob(dataname[[2]], gp=gpar(fontsize=12,font=8))),
  arrangeGrob(p_S[[3]],p_S[[7]],top = textGrob(dataname[[3]], gp=gpar(fontsize=12,font=8))),
  arrangeGrob(p_S[[4]],p_S[[8]],top = textGrob(dataname[[4]], gp=gpar(fontsize=12,font=8))),
  ncol=4)
dev.off()

#plot weights
p_w <-list()
p_w_zoom <- list()
for (n in 1:n_dataset){
  K<-data.frame("index"= 1:111,"weights"=RESULT[[n]]$results$Threekind$alphaK ,"candidates" = c(rep("Gaussian kernel based on Euclidean distance",55),rep("Gaussian kernel based on FZINB distance",55),"FZINB kernel"))
  p_w[[n]] <- ggplot(K, aes(x=index, y=weights,shape = candidates,  color = candidates)) + geom_point() + 
    labs(title = dataname[[n]],
         x = "kernel") +
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          plot.title = element_text(hjust = 0.5))+
    scale_shape_discrete(breaks=c("Gaussian kernel based on Euclidean distance","Gaussian kernel based on FZINB distance","FZINB kernel"))+
    scale_color_discrete(breaks=c("Gaussian kernel based on Euclidean distance","Gaussian kernel based on FZINB distance","FZINB kernel"))
}
for (n in 1:n_dataset){
  K<-data.frame("index"= 1:55,"weights"=RESULT[[n]]$results$Threekind$alphaK[1:55])
  p_w_zoom[[n]] <- ggplot(K,aes(x=index, y=weights) ) + geom_point(shape = 17, colour = "green") + 
    labs( x = "kernel", title = "zoom-in") +
    theme(legend.title=element_blank(),
          axis.title = element_blank())+
    ylim(c(0.0178,0.0184))
}
p_w_zoom[[1]]
png(paste("../FZINB/Pics/weights/threekind.png"),    # create PNG for the heat map        
    width = 8*660,        # 5 x 300 pixels
    height = 4*660,
    res = 500,            # 300 pixels per inch
    pointsize = 8)
ggarrange(
  ggarrange(p_w[[1]],p_w[[2]], p_w[[3]],p_w[[4]],nrow=1,ncol=4,  common.legend = TRUE, legend="bottom"),
  ggarrange(p_w_zoom[[1]],p_w_zoom[[2]], p_w_zoom[[3]],p_w_zoom[[4]],nrow=1,ncol=4,  common.legend = TRUE, legend="bottom"),
  nrow = 2, heights = c(4,2)
  )

dev.off()
ggarrange(
  ggarrange(p_w[[1]], p_w_zoom[[1]],nrow=2, heights = c(4,1)),
  ggarrange(p_w[[2]], p_w_zoom[[2]],nrow=2, heights = c(4,1)),
  ggarrange(p_w[[3]], p_w_zoom[[3]],nrow=2, heights = c(4,1)),
  ggarrange(p_w[[4]], p_w_zoom[[4]],nrow=2, heights = c(4,1)),
  nrow=1,ncol=4,  common.legend = TRUE, legend="bottom")


ggarrange(p_w[[1]], p_w_zoom[[1]],nrow=2, heights = c(4,1))

#### heatmap the distance metric : Gaussian kernel based on Euclidean distance / Gaussian kernel based on FINB kernel distance

p_m <- list()
nmi <- list()
for (n in 1:4){
  D_kernel <- as.matrix(RESULT[[n]]$matrix$Gaussian.Distance)
  DF_kernel <- FZINB.Distance.matrix(RESULT[[n]]$matrix$FZINB.kernel)
  n_clusters <- Test[[n]]$n_clust
  Kmeans_FZINB <- kmeans(DF_kernel, n_clusters, nstart = 30)$cluster
  Kmeans_Gaussian <- kmeans(D_kernel, n_clusters, nstart = 30)$cluster
  true_clusters <- Test[[n]]$true_labs$V1
  
 
  nmi[[n]] <- compare(Kmeans_Gaussian,true_clusters,method = "nmi")
  nmi[[n+4]] <- compare(Kmeans_FZINB,true_clusters,method = "nmi")
  
  clustered <- sort(true_clusters,index.return = TRUE)$ix
  
  melted_FZINB <- melt(DF_kernel[clustered,clustered])
  melted_Gaussian <- melt(as.matrix(D_kernel)[clustered,clustered])
  pp<-0
  for (melted in list( melted_Gaussian, melted_FZINB )){
    colnames(melted)[3] <- "D"
    limits <- range(melted$D[melted$D!=0])
    
    p_m[[n+pp]] <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=D)) + 
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
    
    pp <- pp+4
  }
}

p_m[[1]]<- p_m[[1]] +  labs(y = "Gaussian Kernel Distance (SIMLR)")
p_m[[5]]<- p_m[[5]] +  labs(y =  "FZINB Kernel Distance")

plot(p_m[[5]])

png(paste("../FZINB/Pics/heatmap/distance_compare.trun.png"),    # create PNG for the heat map        
    width = 6*660/2,        # 5 x 300 pixels
    height = 2.3*660/2,
    res = 300/2.5,            # 300 pixels per inch
    pointsize = 8)
grid.arrange(
  arrangeGrob(p_m[[1]],p_m[[5]],top = textGrob(dataname[[1]], gp=gpar(fontsize=12,font=8))),
  arrangeGrob(p_m[[2]],p_m[[6]],top = textGrob(dataname[[2]], gp=gpar(fontsize=12,font=8))),
  arrangeGrob(p_m[[3]],p_m[[7]],top = textGrob(dataname[[3]], gp=gpar(fontsize=12,font=8))),
  arrangeGrob(p_m[[4]],p_m[[8]],top = textGrob(dataname[[4]], gp=gpar(fontsize=12,font=8))),
  ncol=4)
dev.off()




### nmi for all
methodname<- names(RESULT[[1]]$results)
M <- length(methodname)
NMI <- array(0,c(n_dataset,M))
rownames(NMI) <- dataname
colnames(NMI) <- methodname
for (n in 1:n_dataset){
  for (m in 1:M){
    NMI[n,m] <- compare(Test[[n]]$true_labs$V1,RESULT[[n]]$results[[m]]$y$cluster, method = "nmi")
  }
}
NMI



### time taken to calculate kernel distance metric
metricname<- c("Gaussian based FZINB", "Gaussian based Euclidean")
M <- length(metricname)
TIME <- array(0,c(n_dataset,M))
rownames(TIME) <- dataname
colnames(TIME) <- metricname
for (n in 1:n_dataset){
  for (m in 1:M){
    TIME[n,m] <- RESULT[[n]]$execution.time[[m]]
  }
  TIME[n,1]<-TIME[n,1] + system.time(na<-Gaussian.FZINB.kernel(t(Test[[n]]$in_X), FZINB = RESULT[[n]]$matrix$FZINB.kernel))[3]
}
TIME


