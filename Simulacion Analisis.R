#Analisis de la Simulacion
source('aux_functions.R')
load("Simulacion Main Results.RData")

pal <- c('orange','tomato2','maroon2',
                 'steelblue2')
palette(pal)

#PRIMER BLOQUE (Preliminar de tiempos) ------------------------

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
  legend(7,0.85,c('qr','proc','gow','cmds'),cex=0.6,
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
  legend(7,0.85,c('qr','proc','gow','cmds'),cex=0.6,
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
  legend(7,0.85,c('qr','proc','gow','cmds'),cex=0.6,
         col=c(3,2,4,1),lwd=1.5,bty='n')

}

#summary table
table <- rbind(c(sapply(cmds,function(x) x[4]),0),
sapply(proc,function(x) x[4]),
sapply(qr,function(x) x[4]),
sapply(gow,function(x) x[4]))
colnames(table) <- c('sigma_D','sigma_B','loss')
rownames(table) <- c('cmds','proc','qr','gow')
table

rm(cmds,gow,proc,qr,table)

#FIG 2-2
{
  proc <- principal_escenario1 %>%
    filter(method=='proc') %>%
    select(size,t)  %>%
    group_by(size) %>%
    summarize(average=mean(t)) %>%
    select(average) %>% unlist()

  qr <- principal_escenario1 %>%
    filter(method=='qr') %>%
    select(size,t)  %>%
    group_by(size) %>%
    summarize(average=mean(t)) %>%
    select(average) %>% unlist()

  gow <- principal_escenario1 %>%
    filter(method=='gow') %>%
    select(size,t)  %>%
    group_by(size) %>%
    summarize(average=mean(t)) %>%
    select(average) %>% unlist()

}
{
  plot(c(10,50,100),proc, type='b',pch=20,
       xaxt='n',xlab='tamaño de muestra (miles)',
       ylab='tiempo(seg)', xlim=c(0,120),
       ylim=c(0,15),col=2)
  axis(1, at=c(10,50,100), labels=c(10,50,100))
  points(c(10,50,100),qr, type='b',pch=20,
         col=3)
  points(c(10,50,100),gow, type='b',pch=20,
         col=4)
}


#SEGUNDO BLOQUE----------------

#calculo de bias y ecm
sesgo_ecm <- function(tabla){
  bias <- tabla %>%
    select('eig 1','eig 2','eig 3','eig 4','eig 5') %>%
    apply(1,function(x) x-5 ) %>% t() %>%  colMeans
  ecm <- tabla %>%
    select('eig 1','eig 2','eig 3','eig 4','eig 5') %>%
    apply(1,function(x) (x-5)^2 ) %>% t() %>% colMeans
  return(list('bias'=bias,'ecm'=ecm))
}


sesgo_ecm(filtrar_tabla(principal_escenario1,'proc',10000))




principal_escenario1 <- sesgo_y_ecm(principal_escenario1)
principal_escenario2 <- sesgo_y_ecm(principal_escenario2)




swarms <- function(tabla){
  beeswarm(ecm ~ size ,data=filtrar(tabla,'proc'),
           pch=20,col=2,method='compactswarm',priority='random',spacing=sp, ylim=c(0,1))
  beeswarm(ecm ~ size ,data=filtrar(tabla,'qr'),
           pch=20,col=3,method='compactswarm',priority='random',spacing=sp, ylim=c(0,1))
  beeswarm(ecm ~ size ,data=filtrar(tabla,'gow'),
           pch=20,col=4,method='compactswarm',priority='random',spacing=sp, ylim=c(0,1))
 }
swarms(principal_escenario2)

which(filtrar(principal_escenario1,'gow')$ecm>0.4)







#Analisis de sesgo y ECM -------
grafico_autovalores_1 <- function(met,siz, col){
  data <- segundo_bloque_1 %>%
    filter(method==met,size==siz) %>%
    select('eig 1','eig 2','eig 3','eig 4','eig 5')

  plot(1:5,data[1,],
       type='l', ylim=c(14,16)
       ,col=col, xlab='autovalor',
       ylab='NULL')
  for(i in 2:100){
    lines(1:5,data[i,], col=col)
  }
}

grafico_autovalores_1('proc',10000,'tomato3')
grafico_autovalores_1('qr',10000,'magenta3')
grafico_autovalores_1('gow',10000,'steelblue3')

grafico_autovalores_1('proc',50000,'tomato3')
grafico_autovalores_1('qr',50000,'magenta3')
grafico_autovalores_1('gow',50000,'steelblue3')

grafico_autovalores_1('proc',100000,'tomato3')
grafico_autovalores_1('qr',100000,'magenta3')
grafico_autovalores_1('gow',100000,'steelblue3')

grafico_autovalores_2 <- function(met,siz, col){
  data <- segundo_bloque_2 %>%
    filter(method==met,size==siz) %>%
    select('eig 1','eig 2','eig 3','eig 4','eig 5')

  plot(1:5,data[1,],
       type='l', ylim=c(14,18)
       ,col=col, xlab='autovalor',
       ylab='NULL')
  for(i in 2:100){
    lines(1:5,data[i,], col=col)
  }
}

grafico_autovalores_2('proc',10000,'tomato3')
grafico_autovalores_2('qr',10000,'magenta3')
grafico_autovalores_2('gow',10000,'steelblue3')

grafico_autovalores_2('proc',50000,'tomato3')
grafico_autovalores_2('qr',50000,'magenta3')
grafico_autovalores_2('gow',50000,'steelblue3')

grafico_autovalores_2('proc',100000,'tomato3')
grafico_autovalores_2('qr',100000,'magenta3')
grafico_autovalores_2('gow',100000,'steelblue3')







