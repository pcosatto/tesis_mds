#Ejemplos del Capitulo 2
source('aux_functions.R')

# ----------Ejemplo metodos rapidos--------------

#Simulacion de los datos
{
  library(MASS)
  max <- 1
  seed <- 200
  set.seed(seed)
  delta <- as.matrix(dist(mvrnorm(n=10, c(2,2), matrix(c(5,1,1,5),2,2)))
                     + runif(45,0,max))
  n <- 10
  H <- diag(1,n,n) - (1/n) * matrix(1,n,1) %*% matrix(1,1,n)
  B <- -1/2 * H %*% delta^2 %*% H
  options(scipen = 999)
  round(eigen(B)$values,2)
}

array_to_LaTeX(round(delta,2))

#Division y conquista---------------
{
  set.seed(seed); index <- sample(1:10,4)
  noindex <- (1:10)[-index]

  delta_1 <- delta[c(index,noindex[c(1,2,3)]),
                   c(index,noindex[c(1,2,3)])]

  delta_2 <- delta[c(index,noindex[c(4,5,6)]),
                   c(index,noindex[c(4,5,6)])]

 }
Z_1 <- cmdscale(delta_1,2)
Z_2 <- cmdscale(delta_2,2)
Z_1_c <- Z_1[c(1,2,3,4),]
Z_2_c <- Z_2[c(1,2,3,4),]


#FIG 2-2
{
  plot(Z_1, xlim=c(-6,6), ylim=c(-6,6),xlab='coord1', ylab='coord2')
  points(Z_1_c, pch=20, col='firebrick')
  text(Z_1,label=rownames(Z_1), pos=3)
  abline(h=0, lty=2)
  abline(v=0, lty=2)
}

#FIG 2-3
{
  plot(Z_2, xlim=c(-6,6), ylim=c(-6,6),xlab='coord1', ylab='coord2')
  points(Z_2_c, pch=20, col='steelblue', xlab=NULL, ylab=NULL)
  text(Z_2,label=rownames(Z_2), pos=3)
  abline(h=0, lty=2)
  abline(v=0, lty=2)
}


#Prueba con PROCRUSTES
#FIG 2-4
proc <- procrustes(Z_1_c, Z_2_c)
Z_2_rot <- proc$s * Z_2 %*% proc$T + matrix(1,7,1) %*% t(proc$t)
{
  plot(Z_2_rot, xlim=c(-6,6), ylim=c(-6,6),xlab='coord1', ylab='coord2')
  points(Z_2_rot[c(1,2,3,4),], pch=20, col='tan1', xlab=NULL, ylab=NULL)
  text(Z_2_rot,label=rownames(Z_2_rot), pos=3)
  abline(h=0, lty=2)
  abline(v=0, lty=2)
}

Z <- rbind(Z_1,Z_2_rot[-c(1,2,3,4),])
Z_real <- cmdscale(delta[rownames(Z),rownames(Z)],2)
proc <- procrustes(Z,Z_real)
Z_real <- proc$s * Z_real %*% proc$T + matrix(1,10,1) %*% t(proc$t)
colMeans(Z)
#FIG 2-5
{
  plot(Z, xlim=c(-6,6), ylim=c(-6,6),xlab='coord1', ylab='coord2')
  text(Z,label=rownames(Z), pos=3)
  points(Z_real, pch=16, col='grey60')
  text(Z_real, label=rownames(Z),col='grey',pch=20, pos=1)
  abline(h=0, lty=2)
  abline(v=0, lty=2)
}
sum((Z-Z_real)^2)

#Prueba con QR
Z_1 <- cmdscale(delta_1,2)
Z_2 <- cmdscale(delta_2,2)
Z_1_c <- Z_1[c(1,2,3,4),]
Z_2_c <- Z_2[c(1,2,3,4),]
proc <- QR(Z_1_c, Z_2_c)
Z_2_rot <- proc$s * Z_2 %*% proc$T + matrix(1,7,1) %*% t(proc$t)
#FIG 2-6
{
  plot(Z_2_rot, xlim=c(-6,6), ylim=c(-6,6),xlab='coord1', ylab='coord2')
  points(Z_2_rot[c(1,2,3,4),], pch=20, col='tan1', xlab=NULL, ylab=NULL)
  text(Z_2_rot,label=rownames(Z_2_rot), pos=3)
  abline(h=0, lty=2)
  abline(v=0, lty=2)
}

Z <- rbind(Z_1,Z_2_rot[-c(1,2,3,4),])
Z_real <- cmdscale(delta[rownames(Z),rownames(Z)],2)
proc <- procrustes(Z,Z_real)
Z_real <- proc$s * Z_real %*% proc$T + matrix(1,10,1) %*% t(proc$t)
colMeans(Z)
#FIG 2-7
{
  plot(Z, xlim=c(-6,6), ylim=c(-6,6),xlab='coord1', ylab='coord2')
  text(Z,label=rownames(Z), pos=3)
  points(Z_real, pch=16, col='grey60')
  text(Z_real, label=rownames(Z),col='grey',pch=20, pos=1)
  abline(h=0, lty=2)
  abline(v=0, lty=2)
}
sum((Z-Z_real)^2)

#Interpolación----------------------
{
  delta_1 <- delta[index,index]
  delta_21 <- delta[noindex,index]
  Z_1 <- cmdscale(delta_1,2)
}

#Interpolacion de Gower
#FIG 2-8
{

Z_2 <- Gower(Z_1,as.matrix(dist(Z_1)),delta_21)
Z <- rbind(Z_1,Z_2)
rownames(Z) <- c(index,noindex)

Z_real <- cmdscale(delta[rownames(Z),rownames(Z)],2)
proc <- procrustes(Z,Z_real)
Z_real <- proc$s * Z_real %*% proc$T + matrix(1,10,1) %*% t(proc$t)

plot(Z, xlim=c(-6,6), ylim=c(-6,6),xlab='coord1', ylab='coord2')
text(Z,label=rownames(Z), pos=3)
points(Z_real, pch=16, col='grey60')
text(Z_real, label=rownames(Z),col='grey',pch=20, pos=1)
abline(h=0, lty=2)
abline(v=0, lty=2)
}
sum((Z-Z_real)^2)

