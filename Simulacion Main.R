#Simulacion Main Script
source('aux_functions.R'); source('cmdscaling.R')
n_cores <- parallel::detectCores()
n_cores <- 1

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
metricas_bloque2 <- function(Nrep){
  p <- 10; d <- 5
  sigma <- diag(c(rep(15,d)))
  sigma_err <- diag(rep(1,p-d))

  {
    stress1_cmds <- matrix(NA,Nrep,p-1)
    strain_cmds <- matrix(NA,Nrep,p-1)
    stress1_proc <- matrix(NA,Nrep,p-1)
    strain_proc <- matrix(NA,Nrep,p-1)
    loss_proc <- matrix(NA,Nrep,p-1)
    stress1_qr <- matrix(NA,Nrep,p-1)
    strain_qr <- matrix(NA,Nrep,p-1)
    loss_qr <- matrix(NA,Nrep,p-1)
    stress1_gow <- matrix(NA,Nrep,p-1)
    strain_gow <- matrix(NA,Nrep,p-1)
    loss_gow <- matrix(NA,Nrep,p-1)

  }

  for(seed in 1:Nrep){
    {
      set.seed(seed)
      X <- mvrnorm(1000,rep(0,d),sigma)
      E <- mvrnorm(1000,rep(0,p-d),sigma_err)
      Y <- cbind(X,E)
      delta_x <- dist(X)
      delta_y <- dist(Y)
      B_x <- X %*% t(X)
      B_y <- Y %*% t(Y)
    }

    #CMDS
    {
      Z_cmds <- cmdscale(delta_y,k=p)

      for(k in 2:p){
        D_z <- dist(Z_cmds[,1:k]); B_z <- Z_cmds[,1:k] %*% t(Z_cmds[,1:k])

        stress1_cmds[seed,k-1] <- frobenius(delta_x - D_z)/frobenius(D_z)
        strain_cmds[seed,k-1] <- frobenius(B_x - B_z)/frobenius(B_x)
      }
    }

    #Procrustes
    {

      for(k in 2:p){
        Z_proc <- cmdscaling_test(Y,k,l=400,c=2*k,
                                  m=400,seed=seed+1000-1,'proc',n_cores=1)$conf
        Z_proc <- procrustes(Z_cmds[,1:k],Z_proc)$conf
        D_z <- dist(Z_proc); B_z <- Z_proc %*% t(Z_proc)

        stress1_proc[seed,k-1] <- frobenius(delta_x - D_z)/frobenius(D_z)
        strain_proc[seed,k-1] <- frobenius(B_x - B_z)/frobenius(B_x)
        loss_proc[seed,k-1] <- frobenius(Z_cmds[,1:k]-Z_proc)/frobenius(Z_cmds[,1:k])
      }
    }

    #QR
    {

      for(k in 2:p){
        Z_qr <- cmdscaling_test(Y,k,l=400,c=2*k,
                                m=400,seed=seed+1000-1,'qr',n_cores=1)$conf
        Z_qr <- procrustes(Z_cmds[,1:k],Z_qr)$conf
        D_z <- dist(Z_qr); B_z <- Z_qr %*% t(Z_qr)

        stress1_qr[seed,k-1] <- frobenius(delta_x - D_z)/frobenius(D_z)
        strain_qr[seed,k-1] <- frobenius(B_x - B_z)/frobenius(B_x)
        loss_qr[seed,k-1] <- frobenius(Z_cmds[,1:k]-Z_qr)/frobenius(Z_cmds[,1:k])
      }
    }

    #GOW
    {

      for(k in 2:p){
        Z_gow <- cmdscaling_test(Y,k,l=400,c=2*k,
                                 m=400,seed=seed+1000-1,'gow',n_cores=1)$conf
        Z_gow <- procrustes(Z_cmds[,1:k],Z_gow)$conf
        D_z <- dist(Z_gow); B_z <- Z_gow %*% t(Z_gow)

        stress1_gow[seed,k-1] <- frobenius(delta_x - D_z)/frobenius(D_z)
        strain_gow[seed,k-1] <- frobenius(B_x - B_z)/frobenius(B_x)
        loss_gow[seed,k-1] <- frobenius(Z_cmds[,1:k]-Z_gow)/frobenius(Z_cmds[,1:k])
      }
    }

    cat(paste(seed,'-'))

  }

  metricas_cmds <- list(stress1_cmds,strain_cmds)
  metricas_proc <- list(stress1_proc,strain_proc,loss_proc)
  metricas_qr <- list(stress1_qr,strain_qr,loss_qr)
  metricas_gow <- list(stress1_gow,strain_gow,loss_gow)

  return(list(metricas_cmds,
              metricas_proc,
              metricas_gow))
  }
metricas <- metricas_bloque2(100)

#PLOTS
{
  plot(2:p,stress1_cmds, type='l',
       ylim=c(0,1),col='orange',
       xlab='dim', ylab='stress1',lwd=2)
  lines(2:p,stress1_proc, col='tomato3', lwd=2)
  lines(2:p,stress1_qr, col='magenta3', lwd=2)
  lines(2:p,stress1_gow, col='steelblue3', lwd=2)

  plot(2:p,strain_cmds,
       ylim=c(0,1),type='l', col='orange',
       xlab='dim', ylab='strain',lwd=2)
  lines(2:p,strain_proc, col='tomato3', lwd=2)
  lines(2:p,strain_qr, col='magenta3', lwd=2)
  lines(2:p,strain_gow, col='steelblue3', lwd=2)

  plot(2:p,rep(0,length(2:p)),
       ylim=c(0,1),type='l', col='orange', lwd=2,
       xlab='dim', ylab='perdida_Z')
  lines(2:p,loss_proc, col='tomato3', lwd=2)
  lines(2:p,loss_qr, col='magenta3', lwd=2)
  lines(2:p,loss_gow, col='steelblue3', lwd=2)

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






