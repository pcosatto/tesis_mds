source('aux_functions.R')
source("C:/Users/pcosa/Desktop/Librerías de R/mvtools/R/procrustes_mdscal.R")
load("C:/Users/pcosa/Desktop/Tesis/tesis_mds/Cap4 Data.RData")
n_cores <-1
options(scipen=999)
palette <- brewer.pal(8,'Set1')

#Clustering con kernel density estimation via MDS-----------------

#data credit card
{
  data <- read.csv("C:/Users/pcosa/Desktop/82.17/82.17 Presentaciones/Data multivariada/credit-card-data.csv")
  data <- na.omit(data)
  data <- data[,-c(1,2,3,4,8,12,15,16,18)]
  data <- unique(data)
  data$credit_limit <- log(data$credit_limit)
  data$purchases_trx <- sqrt(data$purchases_trx)
  data$oneoff_purchases <- sqrt(data$oneoff_purchases)
  data$install_purchases <- sqrt(data$install_purchases)
  data$cash_adv <- sqrt(data$cash_adv)
}

#fig 4-1
{
  pairs(data[,c(1,2,3,7,8)], lower.panel = panel.smooth,
        upper.panel= NULL, pch=20, lwd=2,
        diag.panel = panel.hist,
        cex.labels = 1)
}

#fig 4-1 bis
{
  par(mfrow=c(2,2))
  graph_par()
  plot(density(data$oneoff_purchases_freq), main='',
       xlab='', ylab='Densidad', lwd=2)
  text(0.6,4,'oneoff_purchases_freq', cex=0.7)
  plot(density(data$install_purchases_freq), main='',
       xlab='', ylab='Densidad', lwd=2)
  text(0.6,2.3,'install_purchases_freq', cex=0.7)
  plot(density(data$cash_adv_freq), main='',
       xlab='', ylab='Densidad', lwd=2)
  text(1,5.5,'cash_adv_freq', cex=0.7)
  plot(density(data$prc_full_payment), main='',
       xlab='', ylab='Densidad', lwd=2)
  text(0.6,11,'prc_full_payment', cex=0.7)
  graph_par()
}

#fig 4-2
{
  mds <- procrustes_mdscal(data,k=2,diss='gower',seed=16497)
  plot(mds, pch=20, xlab='', ylab='')
}

#fig 4-3
{
  plot(mds, pch=20, col=alpha('black',0.1),xlab='', ylab='')
}

# Bidimensional Clustering in reduced dimension

#Fig 4-4
{
  #opt1 <- optimize_dbscan(mds,eps_range = seq(0.015, 0.05, by=0.0001),minpts_range = seq(20, 100, by=1))
  cluster1 <- dbscan::dbscan(mds, eps=opt1[1], minPts=opt1[2])$cluster
  plot(mds[cluster1==0,],pch=20,
       col= alpha('black',0.2))
  points(mds[cluster1 > 0,], pch=20,
         col=alpha(palette[cluster1],1))

}

#Bidimensional clustering in original dimension

#Fig 4-5
{
  opt2 <- optimize_dbscan(data,eps_range = seq(0.22, 0.24, by=0.0001),minpts_range = seq(10, 15, by=1))
  cluster2 <- dbscan::dbscan(data, eps=opt2[1], minPts=opt2[2])$cluster
  plot(mds[cluster2==0,],pch=20,col= alpha('black',0.2))
  points(mds[cluster2 > 0,], pch=20,col=alpha(palette[cluster2],1))

}


#Tridimensional clustering in reduced dimension
{
  #mds_3d <- procrustes_mdscal(data,k=3,diss='gower',seed=16497)
  #opt3 <- optimize_dbscan(mds_3d,eps_range = seq(0.13, 0.15, by=0.0001),minpts_range = seq(10, 15, by=1))
}

#Fig 4-6
{
  palette2 <- c('bisque2','cadetblue','orange','orchid4')
  cluster3 <- dbscan::dbscan(mds_3d, eps=opt3[1], minPts=opt3[2])$cluster
  pairs(mds_3d[cluster3!=0,], pch=20,
        col=alpha(palette2[cluster3[cluster3!=0]],0.9),
        upper.panel = NULL, labels=c('dim1','dim2','dim3'),
        cex.labels = 1.3)
  }

#Group interpretation
#Fig 4-7
{
  par(mfrow=c(2,2))
  for(i in 1:4){
    graph_par()
    plot_group(data,cluster3, i,palette2)
    abline(v=50, lwd=2, col='black')
  }
}

#Bondad de ajuste

par(mfrow=c(2,2))
graph_par()
for(i in 1:4){
  idx <- sample(1:nrow(data),100)
  d <- dist(mds_3d[idx,])
  delta <- as.dist(StatMatch::gower.dist(data[idx,]))
  plot(delta,d, pch=20, col='grey')
  lines(ksmooth(delta,d, bandwidth = 0.11), col='steelblue',lwd=3)
}


#t-SNE
library(Rtsne)

for(perp in perp_grid)
tsne <- Rtsne(data, dims = 2, initial_dims = 8,perplexity = 50)$Y

#Fig 4-8
{
  opt4 <- optimize_dbscan(tsne,
                          eps_range = seq(0.743, 0.75, by= 0.0001),
                          minpts_range = seq(20, 100, by=1))
  cluster4 <- dbscan::dbscan(data, eps=opt4[1], minPts=opt4[2])$cluster
  graph_par()
  palette <- colorRampPalette(palette)(max(cluster4))

  plot(tsne[cluster4==0,],pch=20,col= alpha('black',0.2), xlab='',ylab='')
  points(tsne[cluster4 > 0,], pch=20,col=alpha(palette[cluster4],1))

}






