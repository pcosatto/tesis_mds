#Analisis de la Simulacion
source('aux_functions.R')
load("Simulacion Main Results.RData")

principal_escenario1$method <- factor(principal_escenario1$method,
                                      levels = c('proc','qr','gow'))
principal_escenario2$method <- factor(principal_escenario2$method,
                                      levels = c('proc','qr','gow'))

#pal <- brewer_pal(palette = "Dark2")(4)
pal <- c('orange2','tomato3','maroon3','steelblue3')
palette(pal)

#Tabla 3-1
round(primera,2)
primera <- cbind.data.frame('0'=rep(0,4),primera)
ax <- 0:3

#fig2-1
{
  plot(ax,primera[1,], type='b',pch=20,
       xaxt='n',xlab='tamaño de muestra (miles)',
       ylab='tiempo(seg)')
  axis(1, at=ax,labels=ax)
  coef <- polyinterp(ax,
                     primera[1,])
  curve(coef[1] + coef[2]*x + coef[3]*x^2 + coef[4]*x^3,
        add=TRUE, col=1, lwd=3, lty=15)
}

#fig3-1
{
  plot(0:3,primera[2,], type='b',pch=20,
       xaxt='n',xlab='tamaño de muestra (miles)',
       ylab='tiempo(seg)', col=2, ylim=c(0,0.55),
       xlim=c(0,3.5))
  axis(1, at=0:3, labels=0:3)
  points(0:3,primera[3,], type='b',pch=20,
         col=3)
  points(0:3,primera[4,], type='b',pch=20,
         col=4)
  text(rep(3,3),primera$'3000'[2:4],
       label=c('proc','qr','gow'),pos=4, cex=0.75)
}

rm(primera)

#Plots de metricas
cmds <- lapply(principal_metricas[[1]], colMeans)
proc <- lapply(principal_metricas[[2]], colMeans)
qr <- lapply(principal_metricas[[3]], colMeans)
gow <- lapply(principal_metricas[[4]], colMeans)
rm(principal_metricas)

#fig3-2
{
  p <- 10
  plot(2:p,cmds[[2]], type='l',
       xlim=c(1.4,10.4),ylim=c(0,1),col=1,
       xlab='dimension', ylab=TeX(r'($\sigma_{B_Z}$)'),
                            lwd=2, xaxt='n')
  axis(1, at=2:10,labels=2:10)
  lines(2:p,proc[[2]], col=2, lwd=2)
  lines(2:p,qr[[2]], col=3, lwd=2)
  lines(2:p,gow[[2]], col=4, lwd=2)
  abline(v=5,col='grey50',lty=10)
  legend(7,0.85,c('qr','proc','gow','cmds'),cex=0.8,
  col=c(3,2,4,1),lwd=1.5,bty='n')
}

#fig3-3
{
  plot(2:p,cmds[[1]], type='l',
       xlim=c(1.4,10.4),ylim=c(0,1),col=1,
       xlab='dimension', ylab=TeX(r'($\sigma_{D_Z}$)'),lwd=2, xaxt='n')
  axis(1, at=2:10,labels=2:10)
  lines(2:p,proc[[1]], col=2, lwd=2)
  lines(2:p,qr[[1]], col=3, lwd=2)
  lines(2:p,gow[[1]], col=4, lwd=2)
  abline(v=5,col='grey50',lty=10)
  legend(7,0.85,c('qr','proc','gow','cmds'),cex=0.8,
         col=c(3,2,4,1),lwd=1.5,bty='n')

}

#fig3-4
{
  plot(2:p,rep(0,9), type='l',
       xlim=c(1.4,10.4),ylim=c(0,1),col=1,
       xlab='dimension', ylab=TeX(r'($\sigma_{Z}$)'),lwd=2,xaxt='n')
  axis(1, at=2:10,labels=2:10)
  lines(2:p,proc[[3]], col=2, lwd=2)
  lines(2:p,qr[[3]], col=3, lwd=2)
  lines(2:p,gow[[3]], col=4, lwd=2)
  abline(v=5,col='grey50',lty=10)
  legend(7,0.85,c('qr','proc','gow','cmds'),cex=0.8,
         col=c(3,2,4,1),lwd=1.5,bty='n')

}

#summary table
table_1 <- rbind(c(sapply(cmds,function(x) x[4]),0),
sapply(proc,function(x) x[4]),
sapply(qr,function(x) x[4]),
sapply(gow,function(x) x[4]))
colnames(table_1) <- c('sigma_D','sigma_B','loss')
rownames(table_1) <- c('cmds','proc','qr','gow')
table_1

rm(cmds,gow,proc,qr,table)
par(mar=c(1,1,1,1))
#fig 3-5
{
  beeswarm(t ~ method+size,data=principal_escenario1,
           pch=18,col=alpha(pal[rep(c(2,3,4),3)],0.7), dlim=c(-3,17),
           xaxt='n',
           method='compactswarm',priority='random',
           corral='random',corralWidth=0.8,cex=1,
           labels=1:11, at=c(1,2,3,5,6,7,9,10,11),
           glab='tamaño de muestra',dlab='tiempo(seg)')
  abline(v=4,lty=2,col='grey50')
  abline(v=8,lty=2,col='grey50')
  text(cbind(c(2,6,10),rep(-2.5,3)),
       label=c('n = 10k','n = 50k','n = 100k'),cex=0.8)
  text(cbind(c(1,2,3,5,6,7,9,10,11),
             rep(-0.75,11)),
       label=rep(c('proc','qr','gow'),3),cex=0.7)


}

#fig 3-6
{
  medidas <- function(tabla){
    mean_eig <- tabla %>%
      select('eig 1','eig 2','eig 3','eig 4','eig 5') %>%
      apply(1,function(x) mean(x) )
    return(cbind.data.frame(tabla,mean_eig))
  }
  beeswarm(mean_eig ~ size+method,data=medidas(principal_escenario1),
           pch=16,col=alpha(pal[c(2,2,2,3,3,3,4,4,4)],0.7),
           xaxt='n',dlim=c(4.4,5.6),
           method='compactswarm',priority='random',
           corral='wrap',corralWidth=0.8,cex=0.8,
           labels=1:11, at=c(1,2,3,5,6,7,9,10,11),glab='',
           dlab=TeX(r'($\bar{\lambda}$)'))
  #abline(v=4,lty=2,col='grey50')
  #abline(v=8,lty=2,col='grey50')
  abline(h=5,lty=4, col='grey50')
  text(cbind(c(2,6,10),rep(5.59,3)),
       label=c('proc','qr','gow'),cex=0.8)
  text(cbind(c(1,2,3,5,6,7,9,10,11),rep(4.41,3)),
       label=rep(c('10k','50k','100k'),3),cex=0.7)
  text(cbind(c(2,6,10),rep(4.48,3)),
       label=c('n','n','n'),cex=0.7)

}

#fig 3-7
{
  beeswarm(mean_eig ~ size+method,data=medidas(principal_escenario2),
           pch=16,col=alpha(pal[c(2,2,2,3,3,3,4,4,4)],0.7),
           xaxt='n',dlim=c(4.4,5.6),
           method='compactswarm',priority='random',
           corral='wrap',corralWidth=0.8,cex=0.8,
           labels=1:11, at=c(1,2,3,5,6,7,9,10,11),
           dlab=TeX(r'($\bar{\lambda}$)'), glab='')
  abline(h=5,lty=4, col='grey50')
  text(cbind(c(2,6,10),rep(5.59,3)),
       label=c('proc','qr','gow'),cex=0.8)
  text(cbind(c(1,2,3,5,6,7,9,10,11),rep(4.41,3)),
       label=rep(c('10k','50k','100k'),3),cex=0.7)
  text(cbind(c(2,6,10),rep(4.48,3)),
       label=c('n','n','n'),cex=0.7)

}

#Calculo del ecm de los eigenvalues
ecm_calculation <- function(tabla,metodo,n){
  PAR <- tabla %>%
    filter(method==metodo,size==n) %>%
    select(c('eig 1','eig 2','eig 3','eig 4','eig 5'))

  U <- scale(PAR,scale=FALSE)
  v <- colMeans(apply(PAR,2,function(x) (x-5)))

  squared_bias <- round(as.numeric(t(v) %*% v),4)
  variance <- round(sum(diag(cov(U))),4)

  result <- c('metodo'=metodo,'n'=n,'ecm'=squared_bias+variance,
              'var'=variance,'bias'=squared_bias)

  return(result)
}

#fig 3-8
{
  table_2 <- c()
  metodos <-  c('proc','qr','gow')
  sizes <- c('10000','50000','1e+05')
  for(metodo in metodos){
    for(size in sizes){
      table_2 <- rbind(table_2,
                       ecm_calculation(principal_escenario1,metodo,size))
    }
  }
  table_2 <- as.data.frame(table_2)
  table_2$ecm <- as.numeric(table_2$ecm); table_2$var <- as.numeric(table_2$var);
  table_2$bias <- as.numeric(table_2$bias)
  print(table_2)
  plot(c(10,50,100),table_2$ecm[1:3], type='l', col=2, lwd=2, ylim=c(0,0.2),
       xaxt='n', xlab='tamaño de muestra', ylab='', xlim=c(0,120))
  axis(1, at=c(10,50,100), labels=c('10k','50k','100k'))
  lines(c(10,50,100),table_2$ecm[4:6], col=3, lwd=2)
  lines(c(10,50,100),table_2$ecm[7:9], col=4, lwd=2)
  text(rep(110,3),c(0.009,0.02,0.09), c('gow','proc','qr'), cex=0.8)
}

#fig 3-9
{
  table_3 <- c()
  metodos <-  c('proc','qr','gow')
  sizes <- c('10000','50000','1e+05')
  for(metodo in metodos){
    for(size in sizes){
      table_3 <- rbind(table_3,
                       ecm_calculation(principal_escenario2,metodo,size))
    }
  }
  table_3 <- as.data.frame(table_3)
  table_3$ecm <- as.numeric(table_3$ecm); table_2$var <- as.numeric(table_3$var);
  table_3$bias <- as.numeric(table_3$bias)
  print(table_3)
  plot(c(10,50,100),table_3$ecm[1:3], type='l', col=2, lwd=2, ylim=c(0,0.2),
       xaxt='n', xlab='tamaño de muestra', ylab='', xlim=c(0,120))
  axis(1, at=c(10,50,100), labels=c('10k','50k','100k'))
  lines(c(10,50,100),table_3$ecm[4:6], col=3, lwd=2)
  lines(c(10,50,100),table_3$ecm[7:9], col=4, lwd=2)
  text(rep(110,3),table_3$ecm[c(3,6,9)], c('proc','qr','gow'), cex=0.8)
}




