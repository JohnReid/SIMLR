source("heatmap.r")
library(gplots)
library(reshape2)
library(gridExtra)
library(ggplot2)
library(ggpubr)


l = 1


## plot FZINB kernel and its characteristics
p_heatmap <- ggarrange(p[[l]],p_trun[[l]], nrow=2, heights = c(3,1))
p_supp <- ggarrange(ggarrange(p_0[[l]],p_diag[[l]],p_c[[l]], ncol=3, widths  = c(1,1,1), labels = c("A","B","C")), ggarrange(p_scatter, labels = "D"), nrow=2 , heights = c(3,1))

png(paste("./Pics/FZINB/FZINB_matrix",l,"_theta",theta,".png",sep = ""),    # create PNG for the heat map        
    width = 5*660*1.2,        # 5 x 300 pixels
    height = 3*660*1.2,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
ggarrange(p_heatmap,p_supp, ncol = 2)
dev.off()
## plot normalized FZINB kernel and its characteristics
p_heatmap <- ggarrange(pN[[l]],pN_trun[[l]], nrow=2, heights = c(3,1))
p_supp <- ggarrange(ggarrange(pN_0[[l]],pN_diag[[l]],pN_c[[l]], ncol=3, widths  = c(1,1,1), labels = c("A","B","C")), ggarrange(p_scatter, labels = "D"), nrow=2 , heights = c(3,1))

png(paste("./Pics/FZINB_normalised/FZINB_norm",l,"_theta",theta,".png",sep = ""),    # create PNG for the heat map        
    width = 5*660*1.5,        # 5 x 300 pixels
    height = 3*660*1.5,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
ggarrange(p_heatmap,p_supp, ncol = 2)
dev.off()




## Comparison normalised FZINB kernel or not
png(paste("./Pics/FZINB_normalised/FZINB_norm",l,"_theta",theta,".png",sep = ""),    # create PNG for the heat map        
    width = 5*660,        # 5 x 300 pixels
    height = 3*660,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
ggarrange(ggarrange(p[[l]],p_trun[[l]], nrow=2, heights = c(3,1), labels = c("A","B")) ,ggarrange(pN[[l]],pN_trun[[l]], nrow=2, heights = c(3,1), labels = c("C","D")) )
dev.off()
png(paste("./Pics/FZINB_normalised/FZINB_norm_sup",l,"_theta",theta,".png",sep = ""),    # create PNG for the heat map        
    width = 5*660,        # 5 x 300 pixels
    height = 3*660,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
ggarrange(ggarrange(pN_0[[l]],pN_diag[[l]],pN_c[[l]], ncol=3, widths  = c(1,1,1), labels = c("A","B","C")), ggarrange(p_scatter, labels = "D"), nrow=2 , heights = c(3,1))
dev.off()



l = 2

## Comparison FZINB variation
png(paste("./Pics/FZINB_variation/FZINB_variant",l,"_theta",theta,".png",sep = ""),    # create PNG for the heat map        
    width = 5*660,        # 5 x 300 pixels
    height = 3*660,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
ggarrange(ggarrange(p[[l]],p_trun[[l]],nrow = 2,  labels = c("A", "B"), heights = c(2,1)), ggarrange(p_diag[[l]], p_0[[l]], p_scatter, nrow=3, labels = c("C", "D","E")) , ncol = 2)
dev.off()

png(paste("./Pics/FZINB_variation/FZINB_variant_norm",l,"_theta",theta,".png",sep = ""),    # create PNG for the heat map        
    width = 4*660,        # 5 x 300 pixels
    height = 3*660,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
ggarrange(ggarrange(pN[[l]],pN_trun[[l]], nrow=2, heights = c(3,1), labels = c("A","B")) ,ggarrange(pN_0[[l]],labels = "C"), ncol=2, widths = c(2,1))
dev.off()




## compare with RBF kernel 
png(paste("./Pics/FZINB_compare/FZINB_compare",l,"_theta",theta,".png",sep = ""),    # create PNG for the heat map        
    width = 7*660,        # 5 x 300 pixels
    height = 3*660,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
grid.arrange(ggarrange(p[[l]],p_trun[[l]], nrow=2, heights = c(3,1)),ggarrange(p_RBF,p_RBF_trun, nrow=2, heights = c(3,1)),ggarrange(p_poly,p_poly_trun, nrow=2, heights = c(3,1)) , ncol=3)
dev.off()

### additional ###

## Compare distance

png(paste("./Pics/FZINB_distance/FD",l,"_theta",theta,".png",sep = ""),   # create PNG for the heat map        
    width = 4*660,        # 5 x 300 pixels
    height = 3*660,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
ggarrange(ggarrange(p_D[[l]],p_D_trun[[l]], nrow=2, heights = c(3,1), labels = c("A","B")) ,ggarrange(p_D0[[l]],labels = "C"), ncol=2, widths = c(2,1))
dev.off()



png(paste("./Pics/FZINB_distance/E_theta",theta,".png",sep = ""),   # create PNG for the heat map        
    width = 4*660,        # 5 x 300 pixels
    height = 2*660,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
grid.arrange(p_D[[1]],p_D[[3]], ncol=2)
dev.off()