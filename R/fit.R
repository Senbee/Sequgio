
## == Main fitting function == ##


fitModels <- function(iGene, design, counts, probeLen = 50L, minoverlap=5L,
                      robust = TRUE, use.joint = TRUE, verbose = FALSE, ls.start=FALSE,
                      ridge.lambda = 0, maxit = 50L, std.err = FALSE, Q1 = 0.75, DF = 3L, fix.lambda = exp(20),
                      max.exons = 100, use.trueiLen=FALSE, trueiLen, attr.iLen=FALSE, attr.df=FALSE,useC=FALSE)
{
  
  if (verbose) 
    cat(iGene, " - ")
  
  juncLength <- probeLen - 2L * minoverlap

  X <- design[[iGene]][,,1,drop=FALSE]  ## rows=Transcripts, cols=exons
  dim(X) <- dim(X)[-3]
  dimnames(X) <- dimnames(design[[iGene]])[1:2]

  matInd <- design[[iGene]][,,2,drop=FALSE]
  dim(matInd) <- dim(matInd)[-3]
  dimnames(matInd) <- dimnames(X)

  txNames <- rownames(X)
  iExons <- colnames(X)

  len <- colSums(counts)
  
  Y <- Ye <- counts[iExons, ,drop=FALSE]

  ## All counts zeros: nothing to do!
  if (all(Y == 0)) #|| is.na(design[[iGene]][["condNumber"]]))
    {
      Theta <- matrix(0, nrow = nrow(X), ncol = ncol(Y), dimnames = list(rownames(X), 
                                                           colnames(Y)))
      attr(Theta,"time-elapsed") <- 0
      return(Theta)
    }

  dupTx <- FALSE
  if (nrow(X) > 1)
    {
      ## Keep only unique design matrix rows/transcript. For the identical ones we should just copy over the estimate
      dX <- as.matrix(dist(X, "binary"))
      dX[upper.tri(dX, diag = TRUE)] <- NA

      coef.mat <- dX == 0
      mode(coef.mat) <- 'integer'
      coef.mat[is.na(coef.mat)] <- 0
      diag(coef.mat) <- 1
      coef.mat <- coef.mat[,!rowSums(dX == 0, na.rm = TRUE) >= 1,drop = FALSE]

      dupTx <- any(rowSums(dX == 0, na.rm = TRUE) >= 1)
      X <- X[!rowSums(dX == 0, na.rm = TRUE) >= 1, , drop = FALSE]
    } else ## only 1 tx
      {
        ls.start <- TRUE
      }


  exLen <- matrix(as.numeric(as.logical(X)), ncol=ncol(X), nrow=nrow(X),dimnames=dimnames(X))

  storage.mode(exLen) <- 'double'
  storage.mode(len) <- 'double'
  
    
  ##-- First run: uniform poisson

  txNames <- rownames(X)
  sampNames <- colnames(Y)

  L_j <- attributes(Design[[iGene]])[["txLen"]][txNames]
  Theta <- sapply(sampNames, .fitIt, L_j=L_j, exLen = exLen, XX = X, Y = Y,ridge.lambda=ridge.lambda,
                  len = len, std.err = std.err, Q1 = Q1, robust=robust,use.ls=ls.start)
   
  if(is.vector(Theta))
    {
      dim(Theta) <- list(1,length(Theta))
      colnames(Theta) <- sampNames
    } else
  Theta <- Theta[1:nrow(X), ,drop=FALSE]


  if(std.err)
    ThetaFit <- Theta

  rownames(Theta) <- txNames

  if (!use.joint || all(Theta == 0) || median(Y) <= 1 || ncol(X)==1)
    {
      if(dupTx)
        Theta <- coef.mat %*% Theta
      attr(Theta,"time-elapsed") <- 0
      return(Theta)
    }

  
  tx <- sapply(txNames, function(i) iExons[as.logical(X[i, ])],simplify=FALSE)

  
  ##-- Joint model

  ## 1) Not allowing std.err & resd right now
  ## 2) Fit only with robust = TRUE
  ## 3) ridge.lambda not used

  d <- 1L

  if(ncol(X) >= max.exons)
    lambda <- fix.lambda else
  lambda <- seq(0.05, 1000, len = 5)
  
  time <- proc.time()


  
  if(useC)
      Warning("C code is still under testing. Reverting to pure R computing")
      ## oTheta <- .Call("fitJoint_R",Theta, as.integer(maxit), 0.01, X,  Y, len, 0L, 0L, lambda, as.integer(d), Q1,
      ##                 as.integer(DF), exLen, tx, 1L, 0L)
#  else
      oTheta <- .fitJoint(iTheta = Theta, matInd=matInd, L_j=L_j, maxit = maxit, error.limit = 0.01, X = X,
                          lambda=lambda,
                          Y = Y, len = len, verbose = verbose, std.err = std.err, d = 1, Q1 = Q1, exLen = exLen,
                          tx = tx, DF = DF, use.trueiLen=use.trueiLen, attr.df=attr.df, attr.iLen=attr.iLen,
                          ThetaFit=ThetaFit)  
  new.time <- proc.time()

  rownames(oTheta) <- txNames
    
  if (verbose) 
    cat("DONE\n")

  if(dupTx)
    oTheta <- coef.mat %*% oTheta
  attr(oTheta,"time-elapsed") <- (new.time-time)["elapsed"]
  return(oTheta)
}


######################
## Hidden functions ##
######################

## These are workhorse functions not supposed to be used directly by the user

.fitIt <- function(iSample, L_j = L_j, exLen = exLen, std.err = std.err, ridge.lambda=0,
                   XX = X, Y = Y, len = len, Q1 = Q1, robust=TRUE, use.ls=FALSE)
{
    ##Change
    ## iX <- t(XX) * matrix(exLen * L_j * len[iSample]/1e+09, nr=ncol(XX), nc=nrow(XX), byrow=T)

    ## exLen => same dimension XX
    ## L_j => length=nrow(XX) => XX*L_j will do the proper multiplication
    ## len[iSample] => constant
    
    iX <- t(XX * exLen * L_j * len[iSample]/1e+09)
    iY <- Y[, iSample]
    
    ## 1 exon only!
    if (nrow(iX) == 1)
        return(iY/(exLen * L_j * len[iSample]/1e+09))
    
    if(use.ls)
        return(.fitter(iX,iY))
    
    if (robust)
        .cpois.fitter(x = iX, y = iY, std.err = std.err, Q1 = Q1)
    else     ## setting ridge.lambda to any x > 0 will perform ridge penalization
        .pois.fitter(x = iX, y = iY, lbd = ridge.lambda)
}



.fitter <- function(A, b) {
  sol <- nnls(A,b)
  sol$x
}

.pois.fitter <- function(x, y, maxit = 100, error.limit = 0.01, lbd = 0) {
    x <- cbind(x)
    betas = .fitter(x, y)
    for (it in 1:maxit) {
        old = betas
        xb = c(x %*% betas)
        mu = xb
        res <- y - mu
        q3 <- quantile(abs(res), 0.75, na.rm = T)
        mu[mu < 0.1] <- 0.1
        wt <- 1/mu
        nwt <- (q3/(abs(res) + 0.01))/mu  #avoid 0 in mu
        wt[abs(res) >= q3] <- nwt[abs(res) >= q3]
        wt[wt < 0.1] <- 0.1
        wt[wt > 10] <- 10
        Y = xb + res
        if (lbd > 0) 
            betas <- solve(t(x * wt) %*% x + lbd * diag(ncol(x)), t(x * wt) %*% c(Y))
        else betas <- solve(t(x * wt) %*% x, t(x * wt) %*% c(Y))
        
        secureBetas <- betas
        secureBetas[secureBetas < 0.001] <- 0.001  #Chen: if the secureBetas==0
        err <- max(abs(betas - old)/abs(secureBetas), na.rm = TRUE)
        if (err < error.limit) 
            break
    }
    betas <- ifelse(betas < 0, 0, betas)
    return(round(betas, 7))
}


.cpois.fitter <- function(x, y, maxit = 100, error.limit = 0.01, std.err = FALSE, 
                         resd = FALSE, Q1 = 0) {
  

    betas = .fitter(x, y)
    
    for (it in 1:maxit) {
        xb = c(x %*% betas)
        mu = xb
        res <- y - mu
        mu[mu < 0.1] <- 0.1
        wt <- 1/mu
        
        q3 <- quantile(abs(res), Q1, na.rm = T)
        
        nwt <- (q3/(abs(res) + 0.01))/mu  #avoid 0 in mu
        wt[abs(res) >= q3] <- nwt[abs(res) >= q3]

        wt[wt < 1e-3] <- 1e-3
        wt[wt > 10] <- 10


        old = betas
        
        betas <- .fitter(t(x * wt) %*% x, t(x * wt) %*% y)
        secureBetas <- betas
        secureBetas[secureBetas < 0.001] <- 0.001
        err <- max(abs(betas - old)/abs(secureBetas), na.rm = TRUE)

        if (err < error.limit) 
          break
      }

    
    std_error <- NULL
    if (std.err) {
        wt2 = wt^2 * c(res)^2
        J = t(x * wt2) %*% x
        if (ncol(x) >= 2) {
            Iinv = apply(diag(1, ncol(x)), 2, function(i) {
                .fitter(t(x * wt) %*% x, i)
            })
        }
        else {
            Iinv = apply(diag(1, 1), 2, function(i) {
                .fitter(t(x * wt) %*% x, i)
            })
        }
        cov_matrix = Iinv %*% J %*% Iinv
        std_error <- sqrt(diag(cov_matrix))
    }

    
    if (resd & std.err) {
        res <- y - c(x %*% betas)
        return(list(betas = betas, std.err = std_error, res = res))
      }
    if (!(resd & std.err))
        return(cbind(betas = betas, std.err = std_error))
}
