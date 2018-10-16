source("functions.r")
source("supp_function.r")
library("kernlab")



### Kernel metric ###
theta <- 1
pi <- 0.3
mu <- 100 
xmax <- 500
### Simulated samples from ZINB ######
Theta <- c(mu,theta,pi)
mu <- Theta[1]
theta <- Theta[2] #aggregation paramter, dispersion parameter
pi <- Theta[3]
heat_x <- 0:xmax

y <- rzinb(1000, theta, mu , pi)
p_scatter <- qplot(y, geom="histogram",
                  binwidth = 10,
                  main = "Histogram for 1000 ZINB samples",
                  xlab = "samples",
                  ylab = "frequency") + xlim(-10, xmax)

Fisher_heat <- numeric(0)
Fisher_heat2 <- numeric(0)
Eucliean_heat <- numeric(0)

I <- I.zinb(Theta)
I2 <- I.zinb2(Theta)
for (ll in 1:length(heat_x)){
  Fisher_heat <- cbind(Fisher_heat,Fisher.zinb(cbind(heat_x[ll],heat_x),Theta,I=I))
  Fisher_heat2 <- cbind(Fisher_heat2,Fisher.zinb2(cbind(heat_x[ll],heat_x),Theta,I=I2))
  Eucliean_heat <- cbind(Eucliean_heat,abs(heat_x[ll]-heat_x))
}

Fisher_heat2<-data.matrix(Fisher_heat2)
Eucliean_heat<-data.matrix(Eucliean_heat)
rownames(Fisher_heat) <-  as.character(heat_x); colnames(Fisher_heat) <-  as.character(heat_x)
rownames(Fisher_heat2) <-  as.character(heat_x); colnames(Fisher_heat2) <-  as.character(heat_x)
rownames(Eucliean_heat) <-  as.character(heat_x); colnames(Eucliean_heat) <-  as.character(heat_x)


### Normalized kernel function ####

NFZINB <- list()
l<-1
for (heat in list(Fisher_heat,Fisher_heat2)){
  k <- 1/sqrt(diag(heat))
  NFZINB[[l]] <- heat * (k %*% t(k))
  l<-l+1
}




### plot RBF, polynomial kernel ###
sigma = 100
RBF_heat <- exp(-(Eucliean_heat^2/(2*sigma^2)))

melted <- melt(RBF_heat)
p_RBF <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue",mid = "white",high = "yellow", space = "Lab" ,guide = "colourbar", midpoint = sum(range(RBF_heat))/2)+
  labs(title = paste("RBF kernel"),
       subtitle = paste("length-scale =",sigma),
       x = expression(bold(x)),
       y = expression(bold(x))) +
  scale_y_continuous(trans = "reverse")+
  scale_x_continuous(position = "top") +
  coord_fixed()+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank())+
  theme(panel.border = element_blank())

melted <- melt(RBF_heat[1:100,1:25])
p_RBF_trun <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue",mid = "white",high = "yellow", space = "Lab" ,guide = "colourbar", midpoint = sum(range(RBF_heat))/2)+
  theme_bw()+
  theme(panel.border = element_blank())+
  theme(axis.title=element_blank(),
        #axis.text=element_blank(),
        axis.ticks=element_blank())+
  scale_y_continuous(trans = "reverse")+
  scale_x_continuous(position = "top") +
  labs(title = "Zoom in") +
  coord_fixed()+
  theme(legend.position="none")

q = 1/2
product_heat <- matrix(rep(heat_x,length(heat_x)),ncol = length(heat_x)) * t(matrix(rep(heat_x,length(heat_x)),ncol = length(heat_x)) )
poly_heat <- (product_heat+1)^q
k <- 1/sqrt(diag(poly_heat))

melted <- melt(poly_heat)
p_poly <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue",mid = "white",high = "yellow", space = "Lab" ,guide = "colourbar", midpoint = sum(range(poly_heat))/2)+
  labs(title = paste("polynomial kernel"),
       subtitle = paste("power = ",q),
       x = expression(bold(x)),
       y = expression(bold(x))) +
  scale_y_continuous(trans = "reverse")+
  scale_x_continuous(position = "top") +
  coord_fixed()+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.title=element_blank())+
  theme(panel.border = element_blank())

melted <- melt(poly_heat[1:100,1:25])
p_poly_trun <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue",mid = "white",high = "yellow", space = "Lab" ,guide = "colourbar", midpoint = sum(range(poly_heat))/2)+
  theme_bw()+
  theme(panel.border = element_blank())+
  theme(axis.title=element_blank(),
        #axis.text=element_blank(),
        axis.ticks=element_blank())+
  scale_y_continuous(trans = "reverse")+
  scale_x_continuous(position = "top") +
  labs(title = "Zoom in") +
  coord_fixed()+
  theme(legend.position="none")


### plot FZINB kernel ### 
p <- list()
p_diag <- list()
p_trun <- list()
p_trunn <- list()
p_0 <- list()
p_c <- list()
cc <- xmax/10
name <- c(expression(paste("FZINB kernel ",kappa,"(",bold(x),", ",bold(x),")")),expression(paste("FZINB kernel ",kappa[2],"(",bold(x),", ",bold(x),")")))

l=1
for (data in list(Fisher_heat,Fisher_heat2)){
  melted <- melt(data)
  c <- sum(range(data))/2
  p[[l]] <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue",mid = "white",high = "yellow", space = "Lab" ,guide = "colourbar", midpoint = c)+
    labs(title = name[l],
         subtitle = paste("dispersion =",theta ),
         x = expression(bold(x)),
         y = expression(bold(x))) +
    scale_y_continuous(trans = "reverse")+
    scale_x_continuous(position = "top") +
    coord_fixed()+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          legend.title=element_blank())+
    theme(panel.border = element_blank())

  F_diag <- data.frame("x" = heat_x, "diag" = diag(data))
  p_diag[[l]] <-ggplot(F_diag, aes(x, diag)) + geom_point() +
    labs(title = "Diagonal",
         x = expression(x))
  
  melted <- melt(data[1:100,1:25])
  p_trun[[l]] <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue",mid = "white",high = "yellow", space = "Lab" ,guide = "colourbar", midpoint = c)+
    theme_bw()+
    theme(panel.border = element_blank())+
    theme(axis.title=element_blank(),
          #axis.text=element_blank(),
          axis.ticks=element_blank())+
    scale_y_continuous(trans = "reverse")+
    scale_x_continuous(position = "top") +
    labs(title = "Zoom in") +
    coord_fixed()+
    theme(legend.position="none")
  
  melted <- melt(data[1:40,1:10])
  p_trunn[[l]] <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue",mid = "white",high = "yellow", space = "Lab" ,guide = "colourbar", midpoint = c)+
    theme_bw()+
    theme(panel.border = element_blank())+
    theme(axis.title=element_blank(),
          #axis.text=element_blank(),
          axis.ticks=element_blank())+
    scale_y_continuous(trans = "reverse")+
    scale_x_continuous(position = "top") +
    labs(title = "Zoom in") +
    coord_fixed()+
    theme(legend.position="none")
  p_0[[l]] <-ggplot(data.frame("x" = heat_x, "y" = data[1,]), aes(x, y)) + geom_point() +
    labs(title = "First row (zeros)",
         x = expression(x))
  p_c[[l]] <-ggplot(data.frame("x" = heat_x[1:(xmax-cc)], "y" = diag(data[(1+cc):xmax,1:(xmax-cc)])), aes(x, y)) + geom_point() +
    labs(title = paste("Distance c = ",cc, sep = ""),
         x = expression(x))
  l=l+1
}
p_0[[1]] <- p_0[[1]] + labs(   y = expression(paste(kappa,"(",0,", ",x,")")))
p_0[[2]] <- p_0[[2]] +  labs(  y = expression(paste(kappa[2],"(",0,", ",x,")")))
p_diag[[1]] <- p_diag[[1]] + labs(  y = expression(paste(kappa,"(",x,", ",x,")")))
p_diag[[2]] <- p_diag[[2]] + labs(  y = expression(paste(kappa[2],"(",x,", ",x,")")))
p_c[[1]] <- p_c[[1]] + labs(  y = expression(paste(kappa,"(",x,", ",x+c,")")))
p_c[[2]] <- p_c[[2]] + labs(  y = expression(paste(kappa[2],"(",x,", ",x+c,")")))


### plot Normlised FZINB kernel ### 

pN <- list()
pN_trun <- list()
pN_trunn <- list()
pN_0 <- list()
pN_diag <- list()
pN_c <- list()
name <- c(expression(paste("Normalised FZINB kernel ",kappa)),expression(paste("Normalised FZINB kernel ",kappa[2])))

l=1
for (Ndata in NFZINB){
  melted <- melt(Ndata)
  cd <- sum(range(Ndata))/2
  pN[[l]] <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue",mid = "white",high = "yellow", space = "Lab" ,guide = "colourbar", midpoint = cd)+
    labs(title = name[l],
         subtitle = paste("dispersion =",theta ),
         x = expression(bold(x)),
         y = expression(bold(x))) +
    scale_y_continuous(trans = "reverse")+
    scale_x_continuous(position = "top") +
    coord_fixed()+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          legend.title=element_blank())+
    theme(panel.border = element_blank())
  
  FN_diag <- data.frame("x" = heat_x, "diag" = diag(Ndata))
  pN_diag[[l]] <-ggplot(FN_diag, aes(x, diag)) + geom_point() +
    labs(title = "Diagonal",
         x = expression(x))
  
  melted <- melt(Ndata[1:100,1:25])
  pN_trun[[l]] <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue",mid = "white",high = "yellow", space = "Lab" ,guide = "colourbar", midpoint = cd)+
    theme_bw()+
    theme(panel.border = element_blank())+
    theme(axis.title=element_blank(),
          #axis.text=element_blank(),
          axis.ticks=element_blank())+
    scale_y_continuous(trans = "reverse")+
    scale_x_continuous(position = "top") +
    labs(title = "Zoom in") +
    coord_fixed()+
    theme(legend.position="none")
  
  melted <- melt(Ndata[1:40,1:10])
  pN_trunn[[l]] <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue",mid = "white",high = "yellow", space = "Lab" ,guide = "colourbar", midpoint = cd)+
    theme_bw()+
    theme(panel.border = element_blank())+
    theme(axis.title=element_blank(),
          #axis.text=element_blank(),
          axis.ticks=element_blank())+
    scale_y_continuous(trans = "reverse")+
    scale_x_continuous(position = "top") +
    labs(title = "Zoom in") +
    coord_fixed()+
    theme(legend.position="none")
  
  pN_0[[l]] <-ggplot(data.frame("x" = heat_x, "diag" = Ndata[1,]), aes(x, diag)) + geom_point() +
    labs(title = "First row (zeros)",
         x = expression(x))
  pN_c[[l]] <-ggplot(data.frame("x" = heat_x[1:(xmax-cc)], "y" = diag(Ndata[(1+cc):xmax,1:(xmax-cc)])), aes(x, y)) + geom_point() +
    labs(title = paste("Distance c = ",cc, sep = ""),
         x = expression(x))
  l=l+1
}
pN_0[[1]] <- pN_0[[1]] + labs(   y = expression(paste(kappa^N,"(",0,", ",x,")")))
pN_0[[2]] <- pN_0[[2]] +  labs(  y = expression(paste(kappa[2]^N,"(",0,", ",x,")")))
pN_c[[1]] <- pN_c[[1]] + labs(  y = expression(paste(kappa^N,"(",x,", ",x+c,")")))
pN_c[[2]] <- pN_c[[2]] + labs(  y = expression(paste(kappa[2]^N,"(",x,", ",x+c,")")))
pN_diag[[1]] <- pN_diag[[1]] + labs(  y = expression(paste(kappa^N,"(",x,", ",x,")")))
pN_diag[[2]] <- pN_diag[[2]] + labs(  y = expression(paste(kappa[2]^N,"(",x,", ",x,")")))

### Distance metric ####
# Distance
Dist <- list()

l<-1
for (heat in list(Fisher_heat,Fisher_heat2)){
  k <- 1/sqrt(diag(heat))
  G <- heat * (k %*% t(k))
  Dist[[l]] <- sqrt(abs(2 - 2*G))
  l=l+1
}

### plot distance ###
d <- c(1,1,xmax/2)
p_D <- list()
p_D_trun<- list()
p_D0 <- list()
data <- list(Dist[[1]],Dist[[2]] , Eucliean_heat)
name <- c(expression(paste("FZINB kernel distance ",D,"(",bold(x),", ",bold(x),")")),
          expression(paste("FZINB kernel distance ",D[2],"(",bold(x),", ",bold(x),")")),
          expression(paste("Euclidean distance ",E,"(",bold(x),", ",bold(x),")")))

for (l in 1:3){
  melted <- melt(data[[l]])
  p_D[[l]] <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue",mid = "white",high = "yellow", space = "Lab" ,guide = "colourbar", midpoint = d[l])+
    labs(title = name[l],
         x = expression(bold(x)),
         y = expression(bold(x)) ) +
    scale_y_continuous(trans = "reverse")+
    scale_x_continuous(position = "top") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.title=element_blank())+
    theme(aspect.ratio=1)
  #plot(p_D[[l]])
  
  melted <- melt(data[[l]][1:100,1:25])
  p_D_trun[[l]] <- ggplot(data =  melted, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    scale_fill_gradient2(low = "blue",mid = "white",high = "yellow", space = "Lab" ,guide = "colourbar", midpoint = d[l])+
    theme(legend.position="none")+
    theme(axis.title=element_blank(),
          #axis.text=element_blank(),
          axis.ticks=element_blank())+
    scale_y_continuous(trans = "reverse")+
    scale_x_continuous(position = "top") +
    labs(title = "Zoom in") +
    coord_fixed()
  p_D0[[l]] <-ggplot(data.frame("x" = heat_x, "y" = data[[l]][1,]), aes(x, y)) + geom_point() +
    labs(title = "First row of FZINB distance (zeros)",
         x = expression(x))
}
p_D0[[1]] <- p_D0[[1]] + labs(   y = expression(paste(D,"(",0,", ",x,")")))
p_D0[[2]] <- p_D0[[2]] +  labs(  y = expression(paste(D[2],"(",0,", ",x,")")))
p_D0[[3]] <- p_D0[[3]] +  labs(  y = expression(paste(E,"(",0,", ",x,")")))


l<-1


