#Ejemplos del capitulo 1
source('aux_functions.R')

#----------------Ciudades de Argentina ------------
{
  rm(list = setdiff(ls(), lsf.str()))
  data_rutas <- read.csv("distarg-rutas.csv", header = FALSE)
  data_recta <- read.csv("distarg-recta.csv", header = FALSE)

  names <- read.csv("distarg names.csv", header = FALSE)
  colnames(data_rutas) <- names$V1
  rownames(data_rutas) <- names$V1
  colnames(data_recta) <- names$V1
  rownames(data_recta) <- names$V1

  #Doble centrado
  n <- 12
  H <- diag(1,n,n) - (1/n) * matrix(1,n,1) %*% matrix(1,1,n)
  delta_rutas <- as.matrix(data_rutas)
  delta_recta <- as.matrix(data_recta)

  B_delta_rutas <- -1/2 * H %*% delta_rutas^2 %*% H
  B_delta_recta <- -1/2 * H %*% delta_recta^2 %*% H

  eig_rutas <- eigen(B_delta_rutas)
  Z_rutas <- eig_rutas$vectors[,1:2] %*% diag(sqrt(eig_rutas$values[1:2]))
  B_rutas <- Z_rutas %*% t(Z_rutas)

  eig_recta <- eigen(B_delta_recta)
  Z_recta <- eig_recta$vectors[,1:2] %*% diag(sqrt(eig_recta$values[1:2]))
  B_recta <- Z_recta %*% t(Z_recta)

}

#Fig 1-2
{
  plot(-Z_recta[,2],-Z_recta[,1], xlim = c(-1200,1200), ylim=c(-3000,1500),
       pch=20, type='n',
       xlab=TeX(r'($Z_1$)'), ylab=TeX(r'($Z_2$)'))
  text(-Z_recta[,2],-Z_recta[,1], names$V1, col='firebrick', cex=0.6)
  }

#Fig 1-3
{
  plot(-Z_rutas[,2],-Z_rutas[,1], xlim = c(-1200,1200), ylim=c(-3000,1500),
       pch=20, type='n',
       xlab=TeX(r'($Z_1$)'), ylab=TeX(r'($Z_2$)'))
  text(-Z_rutas[,2],-Z_rutas[,1], names$V1, col='steelblue', cex=0.6)
}

#Fig 1-4
{
  plot(eig_recta$values/1000, type='h', lwd=6, ylab='', xlab='',
       ylim=c(-4000,20000),col='firebrick')

  plot(eig_rutas$values/1000, type='h', lwd=6, ylab='', xlab='',
                                                  ylim=c(-4000,20000),col='steelblue')
  par(mfrow=c(1,1))
  }


#Bondad de ajuste
MDS_GOF(delta_recta, 'CMDS')$sigma_strain[2]
MDS_GOF(delta_rutas, 'CMDS')$sigma_strain[2]

#Prints para el latex
upper<-data.frame(data_rutas - as.matrix(dist(Z_rutas[,1:2])))
upper[upper.tri(data.frame(data_rutas - as.matrix(dist(Z_rutas[,1:2]))))]<-""
upper<-round(as.data.frame(upper),2)
xtable(upper)
cat(array_to_LaTeX(D))

#Errores de representacion por ciudad
xtable(data.frame(data_rutas - as.matrix(dist(Z_rutas[,1:2]))))

#Ciudades de argentina con Least Squares MDS (ratio)
mds <- mds(data_rutas, 2, type='ratio')

#Fig 1-11
{
  par(mfrow=c(1,1))
  plot(-mds$conf[,2],-mds$conf[,1],
       pch=20, type='n',
       xlab=TeX(r'($Z_1$)'), ylab=TeX(r'($Z_2$)'),
       xlim=c(-0.5,0.5))
  text(-mds$conf[,2],-mds$conf[,1], names$V1, col='forestgreen', cex=0.6)
}

#Fig 1-12
{
  D_z <- dist(mds$conf)
  plot(as.dist(delta_rutas),D_z,
       pch=20,xlab=TeX(r'($\delta_{ij}$)'),
       ylab=TeX(r'($\hat{d}_{ij}/d_{ij}$)'),
       col='forestgreen')
  lines(sort(as.dist(delta_rutas)),sort(mds$dhat),
        col='forestgreen',
        lwd=2)
}

#Fig 1-13
{
  plot(as.dist(data_rutas),dist(Z_rutas),
       pch=20, col='steelblue',
       xlab=TeX(r'($\delta_{ij}$)'),
       ylab=TeX(r'($d_{ij}$)'))
  abline(0,1, col='steelblue')
}

clasico <- MDS_GOF(delta_rutas, 'CMDS')
kruskal <- MDS_GOF(delta_rutas, 'LSMDS', kmax=6)

#Fig 1-14
{
  plot(clasico$sigma_1, type='l',
       col='steelblue', lwd=2, ylim=c(0,0.2),
       xlab='dimension', ylab='stress-1')
  lines(kruskal$sigma_1, col='forestgreen', lwd=2)
}


#-------Ejemplo colores --------------
{
  rm(list = setdiff(ls(), lsf.str())); graph_par()
  data <- as.matrix(ekman) + diag(1,14,14)
  delta <- 1- data

  #Escalamiento clásico
  n <- 14
  H <- diag(1,n,n) - (1/n) * matrix(1,n,1) %*% matrix(1,1,n)
  B <- -1/2 * H %*% delta^2 %*% H
  eig <- eigen(B)
  Z <- eig$vectors[,1:2] %*% diag(sqrt(eig$values[1:2]))
  B_z <- Z %*% t(Z)
  D_z <- as.matrix(dist(Z))
}

#Fig 1-5
{
  pal <- c('#2800ff','#0028ff', '#0092ff',
           '#00b2ff','#00ffff','#00ff61','#77ff00',
           '#b3ff00','#fff200','#ffbe00','#ff9b00',
           '#ff5700','#ff0000','#ff0000')
  palette(pal)
  plot(Z, col=1:14, xlab=TeX(r'($Z_1$)'), ylab=TeX(r'($Z_2$)'),
       pch=20, lwd=6, xlim=c(-0.6, 0.6), ylim=c(-0.6, 0.6))
  text(Z[,1],Z[,2], pos=1, offset =1 ,colnames(data), cex=0.6)
}

#Fig 1-6
{
  Z_tri <- eig$vectors[,1:3] %*% diag(sqrt(eig$values[1:3]))
  B_z_tri <- Z_tri %*% t(Z_tri)
  D_z_tri <- as.matrix(dist(Z_tri))
  s3d <- scatterplot3d(Z_tri, pch = 20, color = pal, type='h',
                       grid=TRUE, box=FALSE,
                       xlab=TeX(r'($Z_1$)'), ylab=TeX(r'($Z_2$)'),
                       zlab=TeX(r'($Z_3$)'),
                       angle=110)
  text(s3d$xyz.convert(Z_tri), pos=1, offset =1 ,colnames(data), cex=0.6)
  graph_par()
}

MDS_GOF(delta, 'CMDS')$sigma_strain[2]
MDS_GOF(delta, 'CMDS')$sigma_strain[3]

#Fig 1-7
{
  sigma_strain <- MDS_GOF(delta, 'CMDS')$sigma_strain
  plot(sigma_strain, type='l',
       col='firebrick', lwd=2,
       xlab='dimension', ylab='error-strain')
}

#Usando least Squares Scaling
#Fig 1-8
{
  Z <- mds(delta, 2, 'interval')$conf
  pal <- c('#2800ff','#0028ff', '#0092ff',
           '#00b2ff','#00ffff','#00ff61','#77ff00',
           '#b3ff00','#fff200','#ffbe00','#ff9b00',
           '#ff5700','#ff0000','#ff0000')
  palette(pal)
  plot(Z, col=1:14, xlab=TeX(r'($Z_1$)'), ylab=TeX(r'($Z_2$)'),
       pch=20, lwd=6,xlim=c(-0.8, 0.8), ylim=c(-0.8, 0.8))
  text(Z[,1],Z[,2], pos=1, offset =1 ,colnames(data), cex=0.6)
}

#Plot 6
{
  mds <- mds(delta, 2, 'interval')
  D_z_mds <- as.matrix(dist(mds$conf))
  D_z_mds[upper.tri(D_z, diag=TRUE)] <- NA
  scatter <- cbind('delta'=reshape2::melt(disim)$value,
                   'dist'=reshape2::melt(D_z_mds)$value)

  plot(scatter, pch=20,ylim=c(0,1.4), xlim=c(0,1),
       xlab=TeX(r'($\delta_{ij}$)'), ylab=TeX(r'($d_{\bf{Z}_{ij}}$)'))
  lines(sort(as.dist(delta)),sort(mds$dhat), col='firebrick', lwd=2)
  abline(0,1,lty=2)
}

#Plot 7
{
  clasico <- MDS_GOF(delta, 'CMDS')
  kruskal <- MDS_GOF(delta, 'LSMDS', 'interval')
  plot(clasico$sigma_1, type='l',
       col='steelblue', lwd=2,
       xlab='dimension', ylab='stress-1')
  lines(kruskal$sigma_1, col='forestgreen', lwd=2)
}





#------Vinos (biplots)-------
{
  rm(list = setdiff(ls(), lsf.str()))
  par(family = "Verdana", cex.axis=0.7, cex.lab=0.7, mar=c(4,4,2,3) - 1.5,
      mgp=c(1.1,0.25,0), tcl=0)
  wine <- read.csv("C:/Users/pcosa/Desktop/82.17/82.17 Presentaciones/Data multivariada/winequality-red.csv")
  wine <- wine[,c(7,8,9,11)]
  data <- wine
  colnames(data) <- names(wine)
  data <- unique(data)
  n <- 100; p <- ncol(wine)
  set.seed(16497); index <- sample(1:nrow(data), n, replace=FALSE)
  data <- data[index,]
  data <- round(apply(data,2,rank,ties.method='average')/20)
  data <- scale(data, scale=FALSE)
  rownames(data) <- 1:100
  }

marca <- max(abs(data))
ejes <- 7
#Fig 1-15 - PCA biplot
{
  eig <- eigen(cov(data))
  Gamma <- eig$vectors[,1:2]
  Lambda <- diag(eig$values[1:2]^(1/2))
  G_pca <- data %*% Gamma
  H_pca <- diag(marca,p,p) %*% Gamma

  plot(G_pca, type='n',
       xlab='Z1', ylab='Z2',
       xlim=c(-ejes,ejes), ylim=c(-ejes,ejes))
  text(G_pca, cex=0.7,col= alpha('firebrick',0.7))
  arrows(0,0,H_pca[,1],H_pca[,2], length=0.02)
  text(x=H_pca[,1], y=H_pca[,2], label=names(wine),
       pos=3, cex=0.5, srt=0)
}

#Fig 1-16
{
  stress <- c()
  for(j in seq(0.1,4.9,0.1)){
    stress <- c(stress,MDS_GOF(dist(data, 'minkowski', p=j), 'CMDS',
                               kmax=2)$sigma_1[2])
  }
  plot(seq(0.1,4.9,0.1),stress, type='l', ylim=c(0,1),
       xlab= 'p')
}
p_opt <- seq(0.1,4.9,0.1)[which.min(stress)]

#Fig 1-17 Linear regression biplot


#Fig 1-18 Nonlinear biplot
{
  delta <- as.matrix(dist(data, method='minkowski', p=p_opt))
  n <- nrow(data)
  H <- diag(1,n,n) - (1/n) * matrix(1,n,1) %*% matrix(1,1,n)
  B <- -1/2 * H %*% delta^2 %*% H
  eig <- eigen(B)

  k <- sum(eig$values>0)
  #Configuracion en el espacio R de dimension k
  Y <- eig$vectors[,1:k] %*% diag(eig$values[1:k]^(1/2))

  E_orig <- list()
  for(i in 1:p){
    seq <- seq(0,marca,1/6*marca)
    E_orig[[i]] <- matrix(0,length(seq),p)
    E_orig[[i]][,i] <- seq
  }

  D <- -1/2 * delta^2

  E_R <- list()
  for(i in 1:p){
    E_R[[i]] <- t(apply(E_orig[[i]],1,locus))
  }

  Y <- cbind(Y,matrix(0,n,1))
  P <- procrustes_projection(Y,Y[,1:2])
  E_L <- lapply(E_R, function(i) i %*%P)

  #Proyeccion optima al subespacio L de los puntos en R
  G_mds <- Y %*% P

  svd <- svd(t(G_pca) %*% H %*% G_mds)
  V <- svd$v
  U <- svd$u
  Q <- (V %*% t(U))[,1:2]

  G_mds <- G_mds %*% Q
  plot(G_mds, type='n',
       xlab='Z1', ylab='Z2',
       xlim=c(-ejes,ejes), ylim=c(-ejes,ejes))
  text(G_mds, cex= 0.7, col='steelblue')
  for(i in 1:ncol(data)){
    H_mds <- E_L[[i]] %*% Q
    lines(H_mds)
    text(x=H_mds[nrow(H_mds),1], y=H_mds[nrow(H_mds),2], label=names(wine)[i],
         pos=1, cex=0.5, srt=0)
  }

}






