#Aux Functions


#Batería de funciones TESIS MDS-----------
rm(list = ls())
load <- function(){

  library(extrafont); library(latex2exp); library(xtable)
  library(scales); library(MASS)
  library(readxl); library(smacof); library(scatterplot3d)
  #remotes::install_version("Rttf2pt1", version = "1.3.8")
  #extrafont::font_import()
  loadfonts(device="win")       #Register fonts for Windows bitmap output
}
load()

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

#Bondad de ajuste MDS
MDS_GOF <- function(delta, method=c('CMDS', 'LSMDS'), type=NULL, kmax=NULL){

  frobenius <- function(A){
    sqrt(sum(A^2))
  }

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

#Funciones sencillas
I <- function(n){
  matrix(1, n, 1)
}
H <- function(n,m){
  diag(rep(1,n)) - 1/n * matrix(1, n, 1) %*% matrix(1,1,n)
}

#Metodos rapidos --------------
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

  return(list('T'=T, 's'=s,'t'=t))

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

  return(list('T'=T, 's'=s,'t'=t))

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

#Funciones para biplots-----------------
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


# Funciones para la simulacion principal
eigenvalues_calculation <- function(X_star){
  n <- nrow(X_star)
  val <- eigen((n-1)/n*cov(X_star))$values
  names(val) <- paste('eig',1:length(val))
  return(val)
}
strain_loss_calculation <- function(data,X_star){
  B <- data %*% t(data)
  B_1 <- X_star %*% t(X_star)

  sum((B - B_1)^2)/sum(B^2)

}
mds_simulation <- function(n,Nrep,scenario,p,k,h=0,metodos){

  #Empty data frames and lists
  table <- data.frame()

  id <- 0

  #Variances
  if(h<=1){
    sigma1 <- diag(c(rep(15,p)))
    sigma2 <- diag(c(rep(60,2),rep(15,p-2)))
  } else {
    sigma1 <- diag(c(rep(15,h),rep(1,p-h)))
    sigma2 <- diag(c(rep(60,2),rep(15,h-2),rep(1,p-h)))
  }


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

      if(n <= 2000){
        strain_loss <- strain_loss_calculation(data,X_star)
      } else {strain_loss <- 0}


      #information for the summary
      line <- c(id, i, scenario, n, p,k, h, metodos[j],
                as.numeric(tf-t0,units='secs'), strain_loss)
      names(line) <- c('id','nrep','scenario','size','p','k','h','method',
                       't','strain')

      #eigenvalues
      line <- c(line,eigenvalues_calculation(X_star))

      #table
      table <- dplyr::bind_rows(table,line)

      print(paste(c(metodos[j], n, i, scenario)))
    }
  }

  return(table)



}
