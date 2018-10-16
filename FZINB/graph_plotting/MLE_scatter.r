source("functions.r")
load(file="../data/Test_1_mECS.RData")
load(file="../data/Test_2_Kolod.RData")
load(file="../data/Test_3_Pollen.RData")
load(file="../data/Test_4_Usoskin.RData")


library(ggplot2)
library(gridExtra)

name <- c("Buettner" , "Kolodziejczyk", "Pollen", "Usoskin")
p <- list()
Test <- list(Test_1_mECS , Test_2_Kolod, Test_3_Pollen, Test_4_Usoskin)
n_dataset <- length(Test)
MLE <- list()

for (l in 1:n_dataset){
  data0<-round(10^(Test[[l]]$in_X)-1.0)
  #sorted <- names(sort(apply(data, 1, var), decreasing = TRUE)) #1000 highly variable genes
  #data0 <- data[sorted[1:1000],]
  #data0 <- Test[[l]]$in_X*100
  Theta0 <- prior.zinb(data0)$Theta0
  print(Theta0)
  thet <- MLE.zinb(data0,Theta0,invisible =1)
  
  MLE[[l]]  <- data.frame(mu=thet[,1],theta=thet[,2],pi=thet[,3])
}

 #plot mu against pi

for (l in 1:n_dataset){
  xmax <- quantile(MLE[[l]]$mu,.95)
  ymax <- 1
  nn <- length(which(MLE[[l]]$mu < xmax))
  size <- nrow(MLE[[l]])
  
  p[[l]] <- ggplot(MLE[[l]],aes(x=mu,y = pi)) +
    geom_point(shape=18, color="blue") + 
    xlim(c(0,xmax)) +
    ylim(c(0,ymax)) +
    labs(title = name[l],
         x = expression(hat(mu)),
         y = expression(hat(pi)),
         caption = paste(signif(nn/size*100,3),"% shown")) +
    theme(plot.title = element_text(hjust = 0.5))
}

# plot mu against theta coloured by pi

for (l in 1:n_dataset){
  
  #remove genes with low pi
  xmax <- quantile(MLE[[l]]$mu,.95)
  ymax <- quantile(MLE[[l]]$theta,.95)
  nn <- length(intersect(which(MLE[[l]]$mu < xmax),which(MLE[[l]]$theta < ymax)))
  size <- nrow(MLE[[l]])
  p[[l]] <- ggplot(MLE[[l]],aes(x=mu,y=theta,colour = pi)) +
    geom_point() + 
    scale_colour_gradient(low = "blue",high = "yellow", space = "Lab" ,guide = "colourbar",limits=c(0,1)) +
    xlim(c(0,xmax)) +
    ylim(c(0,ymax)) +
    labs(title = name[l],
         x = expression(hat(mu)),
         y = expression(hat(theta)),
         caption = paste(signif(nn/size*100,3),"% shown")) +
    theme(plot.title = element_text(hjust = 0.5))
}
grid.arrange(p[[1]], p[[2]] , p[[3]], p[[4]], nrow = 2)


for (l in 1:n_dataset){
  #remove genes with low pi
  highpi <- intersect(which(MLE[[l]]$pi>0.3), which(MLE[[l]]$pi<0.9))
  xmax <- quantile(MLE[[l]]$mu[highpi],1)
  ymax <- quantile(MLE[[l]]$theta[highpi],1)
  xmax <- quantile(MLE[[l]]$mu,.95)
  ymax <- quantile(MLE[[l]]$theta,.95)
  nn <- length(intersect(which(MLE[[l]][highpi,]$mu < xmax),which(MLE[[l]][highpi,]$theta < ymax)))
  #nn <- length(intersect(which(MLE[[l]]$mu < xmax),which(MLE[[l]]$theta < ymax)))
  size <- nrow(MLE[[l]])
  #p[[l]] <- ggplot(MLE[[l]],aes(x=mu,y=theta,colour = pi)) +
  p[[l]] <- ggplot(MLE[[l]][highpi,],aes(x=mu,y=theta,colour = pi)) +
    geom_point() + 
    scale_colour_gradient(low = "blue",high = "yellow", space = "Lab" ,guide = "colourbar",limits=c(0,1)) +
  xlim(c(0,xmax)) +
    ylim(c(0,ymax)) +
    labs(title = name[l],
         x = expression(hat(mu)),
         y = expression(hat(theta)),
         caption = paste(signif(nn/size*100,3),"% shown")) +
    theme(plot.title = element_text(hjust = 0.5))
}

grid.arrange(p[[1]], p[[2]] , p[[3]], p[[4]], nrow = 2)



na <-"MLE_truncatePI"
png(paste("./Pics/MLE/",na,".png",sep = ""),    # create PNG for the heat map        
    width = 3*800,        # 5 x 300 pixels
    height = 3*660,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
grid.arrange(p[[1]], p[[2]] , p[[3]], p[[4]], nrow = 2)
#Heatmap.matrix(Fisher,title = "Aggregate Fisher kernel matrix on mESC dataset", clustering = TRUE)
dev.off()  




Theta <- c(100,1,0.3)
y1 <- rzinb(1000, Theta[2], Theta[1] , Theta[3])
Theta <- c(100,10,0.3)
y2 <- rzinb(1000, Theta[2], Theta[1] , Theta[3])
p1 <-qplot(y1, geom="histogram",
          main = "Histogram for 1000 ZINB samples", 
          xlab = "",
          ylab = "frequency") +
  labs(subtitle = paste("dispersion =",1 ))
plot(p1)

p2 <-qplot(y2, geom="histogram",
           xlab = "samples",
           ylab = "frequency") +
  labs(subtitle = paste("dispersion =",10.0 ))
plot(p2)
grid.arrange(p1, p2 , nrow = 2)


png(paste("./Pics/ZINB_hist.png"),    # create PNG for the heat map        
    width = 6*660,        # 5 x 300 pixels
    height = 3*660,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
grid.arrange(p1, p2 , nrow = 2)
dev.off()


### find some typical samples ###
X_counts <- round(10^Test_2_Kolod$in_X-1)
Theta0 <- prior.zinb(X_counts)$Theta0
thet <- MLE.zinb(X_counts,Theta0,invisible =1)

ideal <- which(thet[,2]>4 & thet[,3] > 0.5)


X_counts <- round(10^Test_2_Kolod$in_X-1)
k <- 13109

y0 <- data.frame(weight = as.vector(X_counts[k,]))
p0<-ggplot(y0, aes(x = weight)) + 
  geom_histogram(bins = 30) +
  labs(
    title = "Histogram of 704 cells for a gene", 
    subtitle = paste("The ",k,"th gene in Kolodziejczyk scRNA-seq data",sep = ""),
    y = "number of cells"
  ) + 
  theme(axis.title.x=element_blank()) +
  xlim(NA,quantile(y0$weight,0.999)) +
  theme(plot.title = element_text(hjust = 0.5))

x = c(3,9,15,30,50)
y <-rep(-10,5)
lab <- c("'&' %<->% '&'",rep("'#' %<->% '#'",2),rep("'&' %<->% '&'",2))
for (i in 1:5){
  p0 <- p0 + annotate("text", x = x[i], y = y[i], label = lab[i], colour = "red", size = 5, parse = TRUE)
}

plot(p0)



j <- 410


y00 <- data.frame(weight = as.vector(X_counts[j,]))
p00 <-ggplot(y00, aes(x = weight)) + 
  geom_histogram(bins = 30) +
  labs(
    subtitle = paste("The ",j,"th gene in Kolodziejczyk scRNA-seq data",sep = ""),
    x = "x, counts of the gene",
    y = "number of cells"
  ) + 
  xlim(NA,quantile(y00$weight,0.999)) +
  theme(plot.title = element_text(hjust = 0.5))

x = c(100,300,1500,2080)
y <-rep(-10,4)
lab <- c(rep("'#' %<->% '#'",2),rep("'&' %<->% '&'",2))
for (i in 1:4){
  p00 <- p00 + annotate("text", x = x[i], y = y[i], label = lab[i], colour = "red", size = 5, parse = TRUE)
}
  
plot(p00)




png(paste("./Pics/scRNAseq_hist.png"),    # create PNG for the heat map        
    width = 6*660,        # 5 x 300 pixels
    height = 3*660,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
grid.arrange(p0, p00 , nrow = 2)
dev.off()

