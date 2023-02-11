#Aux Functions

# Batería de librerias y funciones iniciales TESIS MDS-----------
rm(list = ls())
lib <- c('extrafont','latex2exp','xtable','scales',
               'MASS','cmna','readxl','smacof','scatterplot3d',
         'dplyr','beeswarm','biotools')
lapply(lib, require, character.only = TRUE)
#remotes::install_version("Rttf2pt1", version = "1.3.8")
#extrafont::font_import()
loadfonts(device="win")

graph_par <- function(){
  par(family = "Verdana", cex.axis=0.7, cex.lab=0.7, mar=c(4,4,2,3) - 1.5,
      mgp=c(1.1,0.25,0), tcl=0)
}
graph_par()
array_to_LaTeX <- function(arr){
  rows <- apply(arr, MARGIN=1, paste, collapse = " & ")
  matrix_string <- paste(rows, collapse = " \\\\ ")
  return(paste("\\begin{pmatrix}", matrix_string, "\\end{pmatrix}"))
}


# Funciones generales MDS (cuentas)----------------------
frobenius <- function(A){
  sqrt(sum(A^2))
}
MDS_GOF <- function(delta, method=c('CMDS', 'LSMDS'), type=NULL, kmax=NULL){

  delta <- as.matrix(delta)

  sigma_strain <- c()
  sigma_1 <- c()

  n <- nrow(delta)
  H <- diag(1,n,n) - (1/n) * matrix(1,n,1) %*% matrix(1,1,n)
  B <- -1/2 * H %*% delta^2 %*% H
  s <- sum(eigen(B)$values > 0.001)

  switch(method,
         'CMDS' = {
           for(k in 1:s){
             Z <- cmdscale(delta, k)
             D_z <- as.matrix(dist(Z))
             B_z <- Z %*% t(Z)

             sigma_1[k] <- frobenius(D_z - delta)/frobenius(D_z)
             sigma_strain[k] <- frobenius(B - B_z)/frobenius(B)
           }
         },
         'LSMDS' = {
           if(is.null(kmax)){ kmax <- n-1 }
           if(is.null(type)){ type <- 'interval' }
           for(k in 1:kmax){
             mds <- smacof::mds(delta, k, type)
             Z <- mds$conf
             d_hat <- as.matrix(mds$dhat)
             B_hat <- -1/2 * H %*% d_hat^2 %*% H
             D_z <- as.matrix(dist(Z))
             B_z <- Z %*% t(Z)

             sigma_1[k] <- frobenius(D_z - d_hat)/frobenius(D_z)
             sigma_strain[k] <- frobenius(B_hat - B_z)/frobenius(B_hat)
           }
         })

  return(list('sigma_strain'= sigma_strain,
              'sigma_1' = sigma_1))
}
I <- function(n){
  matrix(1, n, 1)
}
H <- function(n,m){
  diag(rep(1,n)) - 1/n * matrix(1, n, 1) %*% matrix(1,1,n)
}
pc.plot3d <- function(X,col=NULL,id=NULL, size=6, cube=TRUE){

  library(plotly)
  library(RColorBrewer)

  n <- nrow(X)
  p <- ncol(X)

  data <- data.frame(scale(X, scale=FALSE))
  names(data) <- c('x1','x2','x3')

  if(is.null(col)==TRUE){
    data$col <- rep('black',n)
  } else {
    data$col <-col}

  if(is.null(id)==TRUE){
    data$id<-1:n
  } else {data$id <- id}

  fig <- plot_ly(data,
                 x = ~data[,1], y = ~data[,2], z = ~data[,3],
                 colors = brewer.pal(p,'Set1'), text=~id,
                 marker = list(size=size))
  fig <- fig %>% add_markers(color = ~col)
  fig <- fig %>% layout(scene = list(xaxis = list(title = colnames(X)[1],
                                                  range = c(min(data$x1),max(data$x1))),
                                     yaxis = list(title = colnames(X)[2],
                                                  range = c(min(data$x2),max(data$x2))),
                                     zaxis = list(title = colnames(X)[3],
                                                  range = c(min(data$x3),max(data$x3))),
                                     aspectmode = ifelse(cube==TRUE,'cube','auto')))
  fig

}



# Metodos rapidos --------------
procrustes <- function(OBJ,PART){

  A <- OBJ
  B <- PART

  n <- nrow(A)
  H <- diag(1,n,n) - (1/n) * matrix(1,n,1) %*% matrix(1,1,n)
  C <- t(A) %*% H %*% B
  svd <- svd(C)

  W <- svd$v
  L <- svd$u

  T <- W %*% t(L)
  s <- sum(diag(C %*% T)) / sum(diag(t(B) %*% H %*% B))
  t <- 1/n * t(A - s * B %*% T) %*% matrix(1,n,1)

  return(list('T'=T, 's'=s,'t'=t,
                'conf'= s * PART %*% T + matrix(1,n,1) %*% t(t)))

}
QR <- function(OBJ,PART){

  Act <- t(scale(OBJ,scale=FALSE))
  Bct <- t(scale(PART,scale=FALSE))

  #QR decompositions
  qa <- qr(Act); qb <- qr(Bct)
  QA <- qr.Q(qa)
  RA <- qr.R(qa); RB <- qr.R(qb)

  #change signs accordingly
  QB <- matrix(1,nrow(Act),1) %*%
    (2*(sign(diag(RA)) == sign(diag(RB)))-1) * qr.Q(qb)


  T <- QA %*% t(QB)
  s <- 1
  t <-  colMeans(OBJ) - t(T) %*% colMeans(PART)

  return(list('T'=T, 's'=s,'t'=t,
              'conf'= s * PART %*% T + matrix(1,n,1) %*% t(t)))

}
Gower <- function(EXIS, D_EXIS, DELTA_NUEVAS){
  A <- EXIS
  D2A <- D_EXIS^2
  D2BA <- DELTA_NUEVAS^2

  m <- nrow(D2BA)
  n <- nrow(D2A)

  I <- function(n){
    matrix(1, n, 1)
  }
  H <- function(n,m){
    diag(rep(1,n)) - 1/n * matrix(1, n, 1) %*% matrix(1,1,n)
  }

  C1 <- (I(m) %*% t(I(n))) %*% (1/n * D2A  -  mean(D2A))
  C2 <- D2BA %*% H(n)
  C3 <- A %*% solve(t(A) %*% A)

  return(1/2 * (C1 - C2) %*% C3)

}


# Funciones para biplots (cap1)-----------------
#Locus para biplots
locus <- function(e){
  I <- function(n){
    matrix(1, n, 1)
  }
  H <- function(n){
    diag(rep(1,n)) - 1/n * matrix(1, n, 1) %*% matrix(1,1,n)
  }

  d_plus <- -1/2 * apply(data,1, function(i) ( sum(abs(i - e)^p_opt))^(2/p_opt))

  y <-  diag(eig$values[1:k]^(-1)) %*% t(Y) %*% (d_plus - 1/n * D %*% I(n))

  y_plus <- sqrt(abs(mean(D) - 2*mean(d_plus) - t(y)%*%y))

  y <- c(y, y_plus)

  return(y)
}

#Procrustes projection para biplots
procrustes_projection <- function(X,Y){
  #Y matriz objetivo - n x q

  #X matriz de partida - n x p

  #A matrix is required - p x q projection matrix
  H <- function(n){
    diag(rep(1,n)) - 1/n * matrix(1, n, 1) %*% matrix(1,1,n)
  }
  r2 <- function(A){
    sum((Y - X %*% A)^2)
  }
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  Y <- cbind(Y,matrix(0,n,p-q))
  err <- 2

  while (err > 0.01){
    svd <- svd( t(Y) %*% H(n) %*% X )
    A_star <- (svd$v) %*% t(svd$u)
    Y[,(q+1):p] <- (X %*% A_star)[,(q+1):p]

    err <- r2(A_star)
  }

  return(A_star[,1:q])
}

# Funciones para la simulacion (cap3) -----------
primera_simulacion <- function(Nrep){
  cat('Primera simulación - Avance %: ')

  sizes <- c(1000, 2000, 3000)
  resultados <- matrix(0,nrow=4,ncol=length(sizes))

  metodos <- c('cmds','proc','qr','gow')

  w <- 0; W <- length(sizes)*Nrep*length(metodos)
  for(j in seq_along(sizes)){
    for(rep in 1:Nrep){

      set.seed(16497+rep)
      X <- mvrnorm(sizes[j],rep(0,2),diag(c(1,1)))

      for(i in seq_along(metodos)){

        if(metodos[i] == 'cmds'){
          t0 <- Sys.time()
          B <- X %*% t(X); eig <- eigen(B)
          escalado <- eig$vectors[,1:2] %*% diag(eig$values[1:2]^(1/2))
          tf <- Sys.time()
        } else {
          t0 <- Sys.time()
          escalado <- cmdscaling_test(X,k=2,l=400,c=4,
                                      m=400,
                                      method=metodos[i],seed = 16497+rep,
                                      n_cores=n_cores)$conf
          tf <- Sys.time()
        }

        resultados[i,j] <- resultados[i,j] + as.numeric(tf-t0,units='secs')
        rm(escalado); w <- w+1; cat(round(100*w/W),' - ')
      }

    }
  }


  resultados <- as.data.frame(resultados/Nrep)
  colnames(resultados) <- as.character(sizes)
  rownames(resultados) <- metodos
  return(resultados)
}
calcular_autovalores <- function(X_star){
  n <- nrow(X_star)
  val <- eigen((n-1)/n*cov(X_star))$values
  names(val) <- paste('eig',1:length(val))
  return(val)
}
simulacion_core <- function(n,Nrep,scenario,p,k,h=0,metodos){

  #Empty data frames and lists
  table <- data.frame()

  id <- 0

  #Variances
  if(h<=1){
    sigma1 <- diag(c(rep(5,p)))
    sigma2 <- diag(c(rep(5,p)))
  } else {
    sigma1 <- diag(c(rep(5,h),rep(1,p-h)))
    sigma2 <- diag(c(rep(5,h),rep(25,2),rep(1,p-h-2)))
  }

  w <- 0; W <- Nrep*length(metodos)
  for(i in 1:Nrep){
    #Simulation

    data <- switch(scenario,{
      set.seed(n+i-1)
      mvrnorm(n,rep(0,p),sigma1)},
      {
        set.seed(n+i-1)
        rbind(mvrnorm(0.95*n,rep(0,p),sigma1),
              mvrnorm(0.05*n,rep(0,p),sigma2))[sample(n),]
      }
    )

    #MDS
    for(j in seq_along(metodos)){

      id <- id + 1

      t0 <- Sys.time()
      #scaling itself - core of the core
      X_star <- cmdscaling_test(data, k, l=400, c=2*k, m=400, method=metodos[j],
                                seed = n+i-1, n_cores=n_cores)$conf
      tf <- Sys.time()

      #information for the summary
      line <- c(id, i, scenario, n, p,k, h, metodos[j],
                as.numeric(tf-t0,units='secs'))
      names(line) <- c('id','nrep','scenario','size','p','k','h','method',
                       't')

      #eigenvalues
      line <- c(line,calcular_autovalores(X_star))

      #table
      table <- dplyr::bind_rows(table,line)

      w <- w+1; cat(round(100*w/W),' - ')
    }
  }

  return(table)



}
metricas_n_fijo <- function(Nrep){
  cat('Metricas con n=1000 - Avance %: ')
  p <- 10; d <- 5
  sigma <- diag(c(rep(5,d)))
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

  w <- 0; W <- Nrep
  set.seed(16497); X <- mvrnorm(1000,rep(0,d),sigma)
  for(seed in 1:Nrep){
    {
      set.seed(seed); E <- mvrnorm(1000,rep(0,p-d),sigma_err)
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
        strain_cmds[seed,k-1] <- frobenius(B_y - B_z)/frobenius(B_y)
      }
    }

    #Procrustes
    {

      for(k in 2:p){
        Z_proc <- cmdscaling_test(Y,k,l=400,c=2*k,
                                  m=400,seed=seed+1000-1,'proc',n_cores = n_cores)$conf
        Z_proc <- procrustes(Z_cmds[,1:k],Z_proc)$conf
        D_z <- dist(Z_proc); B_z <- Z_proc %*% t(Z_proc)

        stress1_proc[seed,k-1] <-  frobenius(delta_x - D_z)/frobenius(D_z)
        strain_proc[seed,k-1] <- frobenius(B_y - B_z)/frobenius(B_y)
        loss_proc[seed,k-1] <- frobenius(Z_cmds[,1:k]-Z_proc)/frobenius(Z_cmds[,1:k])
      }
    }

    #QR
    {

      for(k in 2:p){
        Z_qr <- cmdscaling_test(Y,k,l=400,c=2*k,
                                m=400,seed=seed+1000-1,'qr',n_cores = n_cores)$conf
        Z_qr <- procrustes(Z_cmds[,1:k],Z_qr)$conf
        D_z <- dist(Z_qr); B_z <- Z_qr %*% t(Z_qr)

        stress1_qr[seed,k-1] <- frobenius(delta_x - D_z)/frobenius(D_z)
        strain_qr[seed,k-1] <- frobenius(B_y - B_z)/frobenius(B_y)
        loss_qr[seed,k-1] <- frobenius(Z_cmds[,1:k]-Z_qr)/frobenius(Z_cmds[,1:k])
      }
    }

    #GOW
    {

      for(k in 2:p){
        Z_gow <- cmdscaling_test(Y,k,l=400,c=2*k,
                                 m=400,seed=seed+1000-1,'gow',n_cores = n_cores)$conf
        Z_gow <- procrustes(Z_cmds[,1:k],Z_gow)$conf
        D_z <- dist(Z_gow); B_z <- Z_gow %*% t(Z_gow)

        stress1_gow[seed,k-1] <- frobenius(delta_x - D_z)/frobenius(D_z)
        strain_gow[seed,k-1] <- frobenius(B_y - B_z)/frobenius(B_y)
        loss_gow[seed,k-1] <- frobenius(Z_cmds[,1:k]-Z_gow)/frobenius(Z_cmds[,1:k])
      }
    }


    w <- w+1; cat(round(100*w/W),' - ')

  }

  metricas_cmds <- list('stress1-cmds'=stress1_cmds,
                        'strain_cmds'=strain_cmds)
  metricas_proc <- list('stress1_proc'=stress1_proc,
                        'strain_proc'=strain_proc,
                        'loss_proc'=loss_proc)
  metricas_qr <- list('stress1_qr'=stress1_qr,
                      'strain_qr'=strain_qr,
                      'loss_qr'=loss_qr)
  metricas_gow <- list('stress1_gow'=stress1_gow,
                       'strain_gow'=strain_gow,
                       'loss_gow'=loss_gow)

  return(list(metricas_cmds,
              metricas_proc,metricas_qr,
              metricas_gow))
}
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
simulacion_principal <- function(Nrep,scenario){
  sizes <- c(10000,50000,100000)
  p <- 10; k <- 5; h <- 5
  metodos <- c('proc','qr','gow')
  resultado <- c()

  for(size in sizes){
    cat("\n","\n",'Simulacion principal -','scen: ',scenario,
        ' - Tamaño: ', size, ' - Av %: ')
      resultado  <- dplyr::bind_rows(resultado,
                                     simulacion_core(size, Nrep, scenario,
                                                     p, k, h, metodos))
  }
  return(preparar_tabla(resultado))
}

# Funciones para el analisis y graficos (cap3)-----------
recuperar_solucion <- function(results, id, recover_scaling=FALSE){

  nrep <- results$rep[id]
  scenario <- as.numeric(results$scenario[id])
  n <- as.numeric(results$size[id])
  p <- as.numeric(results$p[id])
  k <- as.numeric(results$k[id])
  h <- as.numeric(results$h[id])
  method <- as.numeric(results$method[id])

  #Simulation of original data
  library(MASS)
  #Variances
  if(h<=1){
    sigma1 <- diag(c(rep(15,p)))
    sigma2 <- diag(c(rep(60,2),rep(15,p-2)))
  } else {
    sigma1 <- diag(c(rep(15,h),rep(1,p-h)))
    sigma2 <- diag(c(rep(60,2),rep(15,h-2),rep(1,p-h)))
  }

  X0 <- switch(scenario,{
    set.seed(n+nrep-1)
    mvrnorm(n,rep(0,p),sigma1)},
    {
      set.seed(n+nrep-1)
      rbind(mvrnorm(0.95*n,rep(0,p),sigma1),
            mvrnorm(0.05*n,rep(0,p),sigma2))[sample(n),]
    }
  )
  return(X0)

  if(recover_scaling == TRUE){
    t0 <- Sys.time()
    X <- cmdscaling_test(X0, k, l=400, c=2*k, m=400, method=method,
                         seed = n+nrep-1, n_cores=n_cores)$conf
    tf <- Sys.time()
    return(list('X0'=X0, 'X'=X, 'tiempo'=as.numeric(tf-t0,units='secs')))
  }

}


