#Simulacion Main Script
source('aux_functions.R'); source('cmdscaling.R')
n_cores <- parallel::detectCores()

#MDS CLASICO-----------------

resultados_cmds <- c()
Nrep <- 10
sizes <- c(1000, 2000, 3000)
dims <- 2
for(size in sizes){
  for(dim in dims){

    tiempo <- 0
    for(rep in 1:Nrep){

      set.seed(16497+rep)
      X <- mvrnorm(size,rep(0,2),diag(c(1,1)))

      t0 <- Sys.time()

      B <- X %*% t(X)
      eig <- eigen(B)
      escalado <- eig$vectors %*% diag(eig$values^(1/2))

      tf <- Sys.time()

      print(c(size,rep))

      tiempo <- tiempo + as.numeric(tf-t0,units='secs')
      rm(escalado)
    }
        resultados_cmds <- c(resultados_cmds,tiempo/Nrep)
  }
}

#plot(c(1000, 2000, 3000),resultados_cmds, type='b')
#curve((x/1500)^3, add=TRUE, col='orange', lwd=2)

#PRIMER BLOQUE-----------------------
primer_bloque <- data.frame()
Nrep <- 10
sizes <- c(1000, 2000, 3000)
dims <- c(2, 10, 25, 50)
metodos <- c('proc','qr','gow')
for(size in sizes){
  for(dim in dims){
    primer_bloque <- dplyr::bind_rows(primer_bloque,
                                      mds_simulation(size, Nrep, scenario=1,
                                                     p=dim, k=dim, h=0, metodos))
  }
}

sizes <- c(1000, 2000, 3000)
dims <- c(2, 3, 4)
metodos <- c('port')
for(size in sizes){
  for(dim in dims){
    primer_bloque <- dplyr::bind_rows(primer_bloque,
                                      mds_simulation(size, Nrep, scenario=1,
                                                     p=dim, k=dim, h=0, metodos))
  }
}



#SEGUNDO BLOQUE - Escenario 1--------------------------
segundo_bloque_1 <- data.frame()
Nrep <- 100
sizes <- c(10000,50000,100000)
scenarios <- 1
metodos <- c('proc','qr','gow')
for(size in sizes){
  for(scenario in scenarios){
    segundo_bloque_1  <- dplyr::bind_rows(segundo_bloque_1 ,
                                      mds_simulation(size, Nrep, scenario,
                                                     p=10, k=5, h=5, metodos))
  }
}

#SEGUNDO BLOQUE - Metricas
metricas <- metricas_bloque2(100)
cmds <- lapply(metricas[[1]], colMeans)
proc <- lapply(metricas[[2]], colMeans)
qr <- lapply(metricas[[3]], colMeans)
gow <- lapply(metricas[[4]], colMeans)
rm(metricas)

#PLOTS
{
  p <- 10
  plot(2:p,cmds[[1]], type='l',
       ylim=c(0,2),col='orange',
       xlab='k', ylab='stress1',lwd=2)
  lines(2:p,proc[[1]], col='tomato3', lwd=2)
  lines(2:p,qr[[1]], col='magenta3', lwd=2)
  lines(2:p,gow[[1]], col='steelblue3', lwd=2)

  plot(2:p,cmds[[2]], type='l',
       ylim=c(0,1),col='orange',
       xlab='k', ylab='strain',lwd=2)
  lines(2:p,proc[[2]], col='tomato3', lwd=2)
  lines(2:p,qr[[2]], col='magenta3', lwd=2)
  lines(2:p,gow[[2]], col='steelblue3', lwd=2)

  plot(2:p,rep(0,9), type='l',
       ylim=c(0,1),col='orange',
       xlab='k', ylab='loss',lwd=2)
  lines(2:p,proc[[3]], col='tomato3', lwd=2)
  lines(2:p,qr[[3]], col='magenta3', lwd=2)
  lines(2:p,gow[[3]], col='steelblue3', lwd=2)

}




#SEGUNDO BLOQUE - Escenario 2---------------------------
segundo_bloque_2 <- data.frame()
sizes <- c(10000,50000,100000)
scenarios <- 2
metodos <- c('proc','qr','gow')
for(size in sizes){
  for(scenario in scenarios){
    segundo_bloque_2 <- dplyr::bind_rows(segundo_bloque_2,
                                      mds_simulation(size, Nrep, scenario,
                                                     p=10, k=5, h=5, metodos))
  }
}



#Preparacion de tablas---------------------------
preparar_tabla <- function(tabla){
  library(dplyr)

  tabla <- tabla[,-1]

  to_factor <- c('scenario','size', 'p','k','h','method')

  tabla <- tabla %>%
    mutate(method = as.factor(method)) %>%
    mutate_if(is.character,as.numeric) %>%
    mutate_at(to_factor, as.factor)

  return(tabla)
}

primer_bloque <- preparar_tabla(primer_bloque)
segundo_bloque_1 <- preparar_tabla(segundo_bloque_1)
segundo_bloque_2 <- preparar_tabla(segundo_bloque_2)

#borrar todo lo que no sirve
remove(dim, dims, metodos, n_cores,
         Nrep, scenario, scenarios, size, sizes)
rm(cmdscaling_test,eigenvalues_calculation,mds_simulation,
   preparar_tabla, eig, B, X, rep, t0, tf, tiempo,
   metricas_bloque2)

#guardar environment
rm(list=setdiff(ls(),setdiff(ls(),lsf.str())))
save.image("Simulacion Main Results.RData")






