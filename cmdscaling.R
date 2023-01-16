#This version is for testing
cmdscaling_test <- function(X0, k, l, c = 2*k, m,
                            method = c('proc', 'qr', 'gow', 'port'),
                            seed=NULL,
                            selection = c('random', 'paradis'),
                            FUN=dist, ...,
                            n_cores=1){

  X0 <- as.matrix(X0) ; n <- nrow(X0)
  method <- match.arg(method)

  if(method == 'proc' || method == 'qr'){

    #function to get partitions
    get_partitions <- function(n, l, c, k, seed) {

      if (l-c <= 0) {
        stop("\"l\" must be greater than \"c\"")
      } else if (l-c <= c) {
        stop("\"l-c\" must be greater than \"c\"")
      } else if (l-c <= k) {
        stop ("\"l-c\" must be greater than \"r\"")
      }

      if(is.null(seed) == FALSE) set.seed(seed)
      permutation <- sample(x = n, size = n, replace = FALSE)

      if (n<=l) {
        list_indexes <- list(permutation)
      } else {
        min_sample_size <- max(k+2, c)
        p <- 1 + ceiling((n-l)/(l-c))
        last_partition_sample_size <- n - (l + (p-2) * (l-c))

        if (last_partition_sample_size < min_sample_size & last_partition_sample_size > 0) {
          p <- p - 1
          last_partition_sample_size <- n - (l + (p-2) * (l-c))
        }

        first_partition <- permutation[1:l]
        last_partition <- permutation[(n-last_partition_sample_size+1):n]
        list_indexes <- split(x = permutation[(l+1):(n-last_partition_sample_size)],
                              f = 1:(p-2))
        names(list_indexes) <- NULL
        list_indexes[[p-1]] <- list_indexes[[1]]
        list_indexes[[p]] <- last_partition
        list_indexes[[1]] <- first_partition
      }

      return(list_indexes)
    }

    #function for alignment itself (both methods)
    alignment <- function(X1, X2, Xa, method){

      n<-nrow(X1)

      if(method == 'proc'){

        #Procrustes transformation
        H <- diag(rep(1,n))-(1/n)*matrix(1,n,1) %*% t(matrix(1,n,1))
        C <- t(X1) %*% H %*% X2
        svd <- svd(C)

        R <- svd$v %*% t(svd$u) #rotation matrix
        s <- sum(diag(C %*% R))/sum(diag(t(X2) %*% H %*% X2)) #dilation factor
        t <- (1/n)*t(X1 - s * X2 %*% R) %*% matrix(1,n,1) #translation

        return(s * Xa %*% R + matrix(1,nrow(Xa),1) %*% t(t))
      }

      if(method == 'qr'){

        X1.mean <-  matrix(colMeans(X1))
        X2.mean <- matrix(colMeans(X2))
        X1 <- t(X1) - X1.mean %*% matrix(1,1,nrow(X1))
        X2 <- t(X2) - X2.mean %*% matrix(1,1,nrow(X2))

        #QR decompositions
        q1 <- qr(X1); q2 <- qr(X2)
        Q1 <- qr.Q(q1)
        R1 <- qr.R(q1); R2 <- qr.R(q2)

        #change signs accordingly
        Q2 <- matrix(1,nrow(X1),1) %*%
          (2*(sign(diag(R1)) == sign(diag(R2)))-1) * qr.Q(q2)

        R <- Q1 %*% t(Q2) #rotation matrix
        t <- X1.mean - R %*% X2.mean #translation

        return(Xa %*% t(R) + matrix(1,nrow(Xa),1) %*% t(
          t))


      }

    }

    #core sequence for methods proc and qr
    X <- matrix(0,nrow = nrow(X0),ncol = k) #create 'empty' configuration

    #seed for partition selection is w * 10^6 + seed
    w <- switch(method, 'proc'= 1 , 'qr' = 2)
    idx <- get_partitions(nrow(X0),l,c,k, seed + w * 1e06) #list of indexes

    idm <- idx[[1]][1:c] #select alignment points
    X[idx[[1]],] <- cmdscale(FUN(X0[idx[[1]],]),k) #first CMDS

    results_list <- parallel::mclapply(idx[-1],
                                       function(i) {
                                         Xj <- cmdscale(FUN(X0[c(idm,i),]),k)
                                         alignment(X[idm,],Xj[1:c,],Xj[-(1:c),],
                                                   method)
                                       },
                                       mc.cores = n_cores)
    #scaling and alignment of the rest

    X[Reduce(c,idx[-1]),] <- Reduce(rbind,results_list)

  }

  if(method == 'gow' || method == 'port'){

    selection <- match.arg(selection)

    #function for selection of landmarks
    landmark_selection <- function(X0, m, selection, seed){

      idm <- c()

      if(selection == 'random'){
        if(is.null(seed) == FALSE) set.seed(seed)
        idm <- sample(n,m)
      }

      if(selection == 'paradis'){
        if(is.null(seed) == FALSE) set.seed(seed)
        i <- ceiling(runif(1, 0, n)) ## randomize the 1st observation picked
        o <- 1
        idm[o] <- i
        X0T <- t(X0)
        repeat {
          D <- sqrt(colSums((X0T-X0[i,])^2))
          D[idm[1:o]] <- NA_real_
          k <- which.max(D)

          D <- sqrt(colSums((X0T-X0[k,])^2))
          D[idm[1:o]] <- NA_real_

          i <- which.min(abs(D - median(D, na.rm = TRUE)))

          o <- o + 1
          idm[o] <- i
          if (o == m) break
        }
      }

      return(idm)
    }

    #selects landmarks
    #seed for landmark selection is w * 10^6 + seed
    w <- switch(method, 'gow'= 3 , 'port' = 4)
    idm <- landmark_selection(X0, m, selection,  seed + w * 1e06)

    #does CMDS for the landmarks
    X <- matrix(NA,n,k)
    X[idm,] <- cmdscale(FUN(X0[idm,]),k)

    if(method=='gow'){

      #optimizing stress with Gower gowolation Formula

      #calculates matrices for Gower Interpolation
      S.inv <- solve((m-1)/m * cov(X[idm,]))
      Delta_1 <- as.matrix(FUN(X[idm,])^2)
      P <- diag(1,m,m) - (1/m)* matrix(1,m,1) %*% matrix(1,1,m)
      q1 <- as.matrix(diag(-1/2* P %*% Delta_1 %*% t(P)))

      #splits the points to embed in groups of size m
      M <- n-m

      #no seed needed here
      permutation <- sample(x = (1:n)[-idm], size = M, replace = FALSE)

      p <- ceiling(M/m)
      last_partition_sample_size <- M - (m + (p-2) * m)
      first_partition <- permutation[1:m]
      idx <- split(x = permutation[(m+1):(M-last_partition_sample_size)],
                   f = 1:(p-2))

      if(last_partition_sample_size >0){
        last_partition <- permutation[(M-last_partition_sample_size+1):M]
        idx[[p]] <- last_partition
      }

      names(idx) <- NULL
      idx[[p-1]] <- idx[[1]]
      idx[[1]] <- first_partition

      #the following function does interpolation itself for a given group
      interpolation <-function(id){

        #sq eucl distances of new points to the first points in original space
        sq_A21 <- apply(X0[id,],1,function(x) {
          d <- 0
          for(i in 1:ncol(X0)) d <- d + (x[i] - X0[idm,i])^2
          return(d) #vectorizing this is not faster
        })

        #interpolation
        return(((matrix(1,length(id),1) %*% t(q1) - t(sq_A21)) %*%
                  X[idm,] %*% S.inv)/(2*m))

      }

      #core sequence for method interp
      results_list <- parallel::mclapply(idx,
                                         function(i) {
                                           interpolation(i)
                                         },
                                         mc.cores = n_cores)

      X[Reduce(c,idx),] <- Reduce(rbind,results_list)

    }

    if(method == 'port'){
      #optimizing stress with PORT routines

      Z <- X[idm,]

      projection <-function(x){

        #calculates delta vector for x
        delta <- 0
        for(i in 1:ncol(X0)) delta <- delta + (x[i]-X0[idm,i])^2
        delta <- sqrt(delta)

        #function to calculate f for a given z
        f <- function(z){
          d <- 0
          for(i in 1:k) d <- d + (z[i]-Z[,i])^2
          d <- sqrt(d)
          sum((delta - d)^2)
        }

        #optimization with PORT routine
        nlminb(rep(0,k),f)$par

      }

      #core sequence for method port
      index_for_port <- as.list(((1:n)[-idm]))
      results_list <- parallel::mclapply(index_for_port,
                                         function(i) {
                                           projection(X0[i,])
                                         },
                                         mc.cores = n_cores)

      X[Reduce(c,index_for_port),] <- Reduce(rbind,results_list)
    }


  }

  return(list('conf'=X, 'idm'=idm))
}

