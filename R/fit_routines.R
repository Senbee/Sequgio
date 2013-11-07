.fitCr <- function(iRegion, exlen, std.err = FALSE, res = FALSE, XX,
                   YY, itheta, Len, Q1, maxit, error.limit, L_j)
  {
      ##Change
    iX <- t(XX[,iRegion] * L_j * itheta*
            matrix(Len, ncol = ncol(itheta),nrow = nrow(itheta), byrow = TRUE)/1e+09 * exlen[iRegion])
    
    iY <- as.numeric(YY[iRegion, ])

    if (res) {
      return(var(.cpois.fitter(x = iX, y = iY, std.err = std.err, out.res = res,
                                  Q1 = Q1)$res))
    }
    else return(.cpois.fitter(x = iX, y = iY, std.err = std.err, Q1=Q1))
  }


Zi2V <- function(Zi, wt, r)
{

  n <- nrow(Zi)

  asum <- colSums(Zi * Zi* matrix(wt, nrow=n, byrow=F))

  diag(asum[r],ncol=sum(r))
}


ZiVy = function(Zi, y, wt, r)
{

  n <- nrow(Zi)

  asum = colSums(Zi * matrix(y, nrow=n, byrow=F) * matrix(wt, nrow=n, byrow=F))

  asum[r]
}


Zi_b <- function(Xf, b, i,ni)
{
                                        #completment of i; b=p
  Xf <- Xf[, !(colnames(Xf) %in% i),drop=FALSE]
  b <- b[!rownames(b) %in% i, ,drop=FALSE]
  ## rownames(b) <- colnames(Xf) <- rownames(X)[!(rownames(X) %in% i)]

  abb <- Xf * rep(t(b), each = ni)
  rownames(abb) <- rep(colnames(b), each = ni)
  return(apply(abb, 1, sum))
}


Zb <- function(Xf, b, Y)
{
  abb <- Xf * rep(t(b), each = ncol(Y))
  rownames(abb) <- rep(colnames(b), each = ncol(Y))
  return(apply(abb, 1, sum))
}


brb <- function(b, r) {
  return(as.vector(b) %*% r %*% b)
}



Yfun <- function(Y, b, Q1 = 0, Xf)
{
  eta <- Zb(Xf, b, Y)
  mu <- eta

  y <- as.vector(t(Y))  # vectorized data **** by regions  ****

  err <- y - mu

  q3 <- quantile(abs(err), Q1, na.rm = T)

  mu[mu < 0.1] <- 0.1
  wt <- 1/mu
  nwt <- (q3/(abs(err) + 0.01))/mu  #avoid 0 in mu
  wt[abs(err) >= q3] <- nwt[abs(err) >= q3]

  #Change 2)
  #wt[wt < 0.1] <- 0.1
  wt[wt > 10] <- 10

  Y.out <- eta + err
  return(list(wt = c(wt), Y = c(Y.out)))
}

.getR <- function(x,d)
{
  npar <- length(x) + 1

  if (npar > 2)
    {
      d1 <- d2 <- diag(length(x)-1)
      d1 <- cbind(d1*(-1),0)
      d2 <- cbind(0,d2)
      delta <- d1+d2
                                        #degree of differencing
      if (d == 1)
        R <- t(delta) %*% delta

      if (d == 2) {
        delta2 <- delta[-1, -1] %*% delta
        R <- t(delta2) %*% delta2
      }

      colnames(R) <- rownames(R) <- x
    }

  else R <- 1

  return(R)
}

.getStart <- function(i,exLen,tx,iter,p,DF, XX)
{

  txLen <- exLen[tx[[i]]]

  if (length(txLen) > 3) {
    if (iter == 1) {
      ssp <- smooth.spline(cumsum(txLen), p[i, tx[[i]]], df = 3)
    }
    else {
      ssp <- smooth.spline(cumsum(txLen), p[i, tx[[i]]], df = DF)
    }
    p[i, tx[[i]]] <- fitted(ssp)
  }
  return(p[i, ])
                                        #}
}

.getXf <- function(i,X,iTheta,len,exLen, L_j)
{
  Mat1 <- matrix(len,byrow=TRUE,nrow=nrow(iTheta),ncol=ncol(iTheta))
    ##Change
  iX <- (X[,i] * L_j * iTheta * Mat1 * exLen[i]) / 1e+09
  return(iX)
}

.getZi <- function(i, Xf, Y)
  {
    Zi <- matrix(Xf[, i],nrow=ncol(Y),byrow=FALSE)
    colnames(Zi) <- rownames(Y)
    return(Zi)
  }


## Iterative Backfitting Algorithm

.getLambda <- function(la,Zi,Xf,p,Y,XX,R,Q1,iTheta,tx)
  {

    tx.ids <- rownames(XX)
    names(tx.ids) <- tx.ids

    ## Update beta_j
    .getTestbeta <- function(i,Xf,Zi,p,a,XX,tx,ni)
      {

        #sel <- XX[i, ] == 1
        sel <- XX[i, ] != 0

        if (all(iTheta[i, ] == 0) | length(tx[[i]]) == 1)
          out <- p[i, sel]
        else {
          wt <- a$wt
          if (nrow(XX) == 1) {
            yc <- a$Y
          }
          else {
            yc <- a$Y - Zi_b(Xf, p, i, ni) ## get lambda_c
          }

          ## do we really need nnls here?
          out <-  .fitter(Zi2V(Zi[[i]], wt, sel) + (R[[i]] * la), ZiVy(Zi[[i]], yc, wt, sel))
          ## out <-  solve(Zi2V(Zi[[i]], wt, sel) + (R[[i]] * la), ZiVy(Zi[[i]], yc, wt, sel))
          ## out[out < 0] <- 0
        }

        names(out) <- tx[[i]]
        return(out)
      }

    ## Update variance components
    .getdf <- function(i,tx,iTheta,R,wt,Zi,la,XX)
      {

        ## No need to use NA, given that result of .getLambda (dfall) is used to get si where we already limit
        ## to !all(iTheta[i, ] == 0)
        ## ans <- NA

        ans <- 0 ## We would get zero if all Theta == 0. Save computing time
        #sel <- XX[i, ] == 1
        sel <- XX[i, ] != 0

        if (!all(iTheta[i, ] == 0))
          {
            z2v <- Zi2V(Zi[[i]], wt, sel)
            imat <- z2v + (R[[i]] * la)
            imat[imat==0] <- 0.001
            smat <- solve(imat)
            ans <- sum(smat * z2v) ## get sigma^2_j expressed as df
          }

        return(ans)
      }



      for (it in 1:5)
        { ## start iteration. This will iterate to update p

          a <- Yfun(Y=Y, b=p, Q1=Q1, Xf=Xf)
          wt <- a$wt

          b <- lapply(tx.ids, .getTestbeta, Xf=Xf,Zi=Zi,p=p,a=a,XX=XX,tx=tx,ni=ncol(Y))

          ## let's update p matrix
          for(k in names(b))
            {
              if(!(all(iTheta[k, ] == 0) | length(tx[[k]]) == 1)) ## no p update needed in this case!!!
              p[k,names(b[[k]])] <- b[[k]]
            }

        } ## end iteration

    ## .getdf will use the latest version of wt so I guess we need latest wt to feed into .getdf.
    ## To achieve 5 iteration this must be done here too
    ## a <- Yfun(Y=Y, b=p, Q1=Q1, Xf=Xf)
    ## wt <- a$wt

    dfall <- sapply(rownames(XX), .getdf,tx=tx,iTheta=iTheta,R=R,wt=wt,Zi=Zi,la=la,XX=XX)
    return(dfall)
  }


.getbeta <- function(i,a,Xf,p,ni,iTheta,tx,XX,R,si,Zi)
{
  out <- p[i,]

  wt <- a$wt

  if (nrow(XX) == 1) {
    yc <- a$Y
  } else {
    yc <- a$Y - Zi_b(Xf=Xf, b=p, i=i, ni=ni)
  }
  if (!(all(iTheta[i, ] == 0) | length(tx[[i]]) == 1))
    {
      #out[XX[i, ] == 1] <- .fitter(Zi2V(Zi[[i]], wt, XX[i, ] == 1) + (R[[i]]/si[i]), ZiVy(Zi[[i]], yc, wt, XX[i, ] == 1))
      out[XX[i, ] != 0] <- .fitter(Zi2V(Zi[[i]], wt, XX[i, ] != 0) + (R[[i]]/si[i]), ZiVy(Zi[[i]],
                                                           yc, wt, XX[i, ] != 0))

    }
  ## names(out) <- tx[[i]]
  return(out)
}


.getsi <- function(i,tx,Zi,R,XX,wt,si,d,p)
{
    smat <- apply(diag(1, length(tx[[i]])), 2, function(x)
                  {
                                        #.fitter(rbind(Zi2V(Zi[[i]], wt, XX[i, ] == 1) + R[[i]]/si[i]), x)
                      .fitter(rbind(Zi2V(Zi[[i]], wt, XX[i, ] != 0) + R[[i]]/si[i]), x)
                  })
    
    ## mat[[i]] <- smat * R[[i]] ## why index?
    mat <- smat * R[[i]]
    
                                        #sel <- X[i,] == 1
    sel <- XX[i,] != 0
    ## What is this supposed to do? Where is df used afterwards?
    ## Do not use df, which stands for density of F
    ## df[i] <- sum(smat * Zi2V(Zi[[i]], wt, XX[i, ] == 1))
    
    ## if (DF == 0) { ## Useless: this condition is already used to execute getsi or not!
    ## si[i] <- 1/(length(tx[[i]]) - d) * (brb(b[[i]], R[[i]]) + sum(mat))
    isi <- (p[i,sel] %*% (R[[i]]+ sum(mat)) %*% p[i,sel])/(length(tx[[i]]) - d)
    ## }
    
    return(isi)
}

# for normalization s.t sum(pex) equals to the original effective tx length
.getPexLen <- function(i,exlen, pexlen, tx, XX) {
  txLen <- sum(exlen[tx[[i]]])

  sel <- as.logical(XX[i, ])
  pex <- pexlen[i,]

  #Change 3) ####
  #pex[pex < min(exLen] <- min(exLen)
  pex[pex < 1] <- 1

  if(all(pex[sel]==1))
    return(exlen)
  ###############

  pex[!sel] <- exlen[!sel]

  pex[sel] <- pex[sel]/sum(pex[sel]) * txLen

  return(pex)
}

.updateP <- function(Xf, tx, Zi, p, X, R, Y, si,iTheta,DF,d,iter,Q, attr.df)
  {
    for (it in 1:iter)
      {
                                        # step 1: update bi
        a <- Yfun(Y=Y, b=p, Q1=Q, Xf=Xf)
        p <- do.call("rbind",sapply(names(tx), .getbeta,a=a,Xf=Xf,p=p,ni=ncol(Y),iTheta=iTheta,tx=tx,
                                    XX=X,R=R,si=si, Zi=Zi,simplify=FALSE))

                                        # step 2: update variance components
        if (DF == 0) {
          si <- sapply(rownames(X), .getsi,tx=tx,Zi=Zi,R=R,XX=X,wt=a$wt,si=si,d=d,p=p)
          names(si) <- rownames(X)
        }
      }
    if(attr.df){
          .getdfout <- function(i,tx,iTheta,R,wt,Zi,la,XX)
      {

        ## No need to use NA, given that result of .getLambda (dfall) is used to get si where we already limit
        ## to !all(iTheta[i, ] == 0)
        ## ans <- NA
        la=la[i]
        ans <- 0 ## We would get zero if all Theta == 0. Save computing time
        #sel <- XX[i, ] == 1
        sel <- XX[i, ] != 0

        if (!all(iTheta[i, ] == 0))
          {
            z2v <- Zi2V(Zi[[i]], wt, sel)
            imat <- z2v + (R[[i]] * la)
            imat[imat==0] <- 0.001
            smat <- solve(imat)
            ans <- sum(smat * z2v) ## get sigma^2_j expressed as df
          }

        return(ans)
      }
      dfout <- sapply(rownames(X), .getdfout,tx=tx,iTheta=iTheta,R=R, wt=a$wt, Zi=Zi, la=1/si, XX=X)
      return(list(p=p,si=si, dfout=dfout))
    }

    return(list(p=p,si=si))
  }

.fitJoint <- function(iTheta, maxit = 50, error.limit = 0.01, lambda,
                      X, matInd, L_j, Y, len, verbose = FALSE, std.err = FALSE,
                      spar = NULL, d = 1,
                      Q1 = 0, DF = 3, exLen, tx, robust = TRUE, ridge.lambda = 0,
                      use.trueiLen, attr.iLen, attr.df, trueiLen,ThetaFit)
{

  ## COMMENTS
  ## r = number of regions (i.e. exons)
  ## j = number of isoforms
  ## i = number of samples

  ##Change
  #grouping count, X, exLen, to 1D
  r1d <- sapply(rownames(Y), function(x) strsplit(x, '\\.')[[1]][1]) 
  
  Y1d <- apply(Y, 2, function(x) tapply(x, r1d, sum)) 
  msort <- sort(rownames(Y1d))
  Y1d <- Y1d[msort, ]

  n <- nrow(Y1d)
  ni <- rep(ncol(Y1d), n)
  N <- length(Y1d)
   
  #X1d <- aggregate(matrix(X, nr=nrow(Y), byrow=T), by=list(r1d), FUN=sum)
  X1d <- apply(matrix(X, ncol=nrow(Y)), 1, function(x) tapply(x, r1d, sum))
  X1d <- X1d[msort,]
  X1d <- matrix(X1d, ncol=n, byrow=T, dimnames=list(rownames(X), rownames(Y1d)))

  exLen1d <- rep(1, nrow(Y1d))
  names(exLen1d) <- rownames(Y1d)
  exLen1d <- exLen1d[rownames(Y1d)]
    
  tx1d <- list()
  for(i in rownames(X)){
    y1dBytx <- colnames(X)[as.logical(X[i, ])]
    tx1d[[i]] <- unique(na.omit(sapply(y1dBytx, function(x) strsplit(x, '\\.')[[1]][1])))
  }

  R <- lapply(tx1d, .getR, d=d) ## compute first-order difference matrix R for smoothing random effects.

  effEL = apply(X1d * L_j, 2, function(x) x[as.logical(x)][1])
  iLen <- matrix(effEL, byrow=T, ncol=n, nrow=length(tx1d), dimnames=list(rownames(X), rownames(Y1d)))
      
  for (iter in 1:maxit)
    {
      ##Change
      ## we use Theta from Uniform Poisson model to get an estimate of c_rj
      p <- sapply(rownames(Y1d), .fitCr, L_j=L_j, exlen = exLen1d, XX = X1d, YY = Y1d, Len = len,
                  itheta = iTheta, Q1 = Q1)

      if(is.vector(p))
        {
          dim(p) <- list(1,length(p))
          colnames(p) <- rownames(Y1d)
        }

      rownames(p) <- rownames(X)

      ##Change
      p <- t(sapply(rownames(X), .getStart,exLen=effEL,tx=tx1d,iter=iter,p=p,DF=DF,XX=X1d))

      if(is.vector(p)| ncol(X)==1)
        {
          dim(p) <- list(1,length(p))
          colnames(p) <- rownames(Y1d)
          rownames(p) <- rownames(X)
        }

      ## Are these the same DF as in .fitJoint argument? They cannot be DF=0 in smooth.spline...so...

      ##Change
      ## design matrix Xf and other related quantities for fixed effects. This a matrix r*i x j and holds a_rj elements
      Xf <- t(do.call("cbind",lapply(rownames(Y1d),
                                     function(x) .getXf(x,X=X1d,iTheta=iTheta,len=len,exLen=exLen1d, L_j=L_j))))
      XfXf <- t(Xf) %*% Xf
      ii <- rep(1:ncol(Y1d), n) + rep((0:(n - 1) * ncol(Y1d) * (n + 1)), each = ncol(Y1d))
      Zi <- lapply(rownames(X), .getZi, Xf=Xf, Y=Y1d)
      names(Zi) <- rownames(X)

      if(length(lambda)==1)
        {
          si <- 1/rep(lambda,nrow(X))
          names(si) <- rownames(X)
        } else
        {
          dfall <- sapply(lambda,.getLambda,Zi=Zi,Y=Y1d,Xf=Xf,p=p,XX=X1d,R=R,Q1=Q1,iTheta=iTheta,tx=tx1d)
          if(is.vector(dfall))
            {
              dim(dfall) <- list(1,length(dfall))
              rownames(dfall) <- rownames(iTheta)
            }

          if(all(dfall==0)){
            si <- 1/rep(max(lambda), nrow(X))
            names(si) <- rownames(X)
          } else {
            si <- sapply(rownames(iTheta), function(i)
                       {
                         ifelse(all(iTheta[i, ] == 0), 0.05, 1/approx(dfall[i, ], lambda, xout = DF,
                                            rule = 2)$y)})
          }
        }
      
      ##Change
      ## This will update si => wt => p
      pUP <- .updateP(Xf, tx1d, Zi, p, X1d, R, Y1d, si,iTheta,DF,d,iter=10,Q1,attr.df=attr.df)

      p <- pUP$p

      if(DF==0)
        si <- pUP$si

      if(verbose && any(is.na(p)))
        warning("Missing values in p")

      p[is.na(p)] <- min(p, na.rm = T) ## When do we get this??

      ##Change
      pexLen <- p * matrix(effEL,byrow=TRUE,ncol=ncol(p),nrow=nrow(p))

      ntx <- names(tx)
      names(ntx) <- ntx

      pexLen <- do.call("rbind",lapply(ntx, .getPexLen,exlen=effEL,pexlen=pexLen,tx=tx1d,XX=X1d))

      p <- pexLen/matrix(effEL,byrow=TRUE,ncol=ncol(p),nrow=nrow(p))

      oldLen <- iLen
      iLen <- pexLen

      secureP <- iLen
      secureP[iLen == 0] <- min(iLen[iLen > 0])
      secureP[iLen < 0.001] <- 0.001
      checkError <- max(abs(iLen - oldLen)/secureP, na.rm = TRUE)


      if (checkError < error.limit)
        break
      oldTheta <- iTheta
      
      ##Change
      p2d <- apply(p, 1, function(x) rep(x, table(r1d)))
      rownames(p2d) <- names(sort(r1d))
      p2d <- p2d[names(r1d),]
      p2d[is.na(p2d)] <- 0
      
      if(use.trueiLen)
          iTheta <- sapply(colnames(Y), .fitIt,exLen = t(trueiLen), XX = X, Y = Y, len = len,
                           std.err = std.err, Q1 = Q1)
      else
          iTheta <- sapply(colnames(Y), .fitIt, exLen = t(p2d), XX = X, Y = Y, len = len, std.err = std.err,
                           Q1 = Q1, L_j=L_j)
      
      if(std.err){
          iFit <- iTheta
          if(is.vector(iFit))
              {
                  dim(iFit) <- list(1,length(iFit))
                  colnames(iFit) <- colnames(Y)
              }
          iTheta <- iTheta[1:nrow(X), ]
      }
      if(is.vector(iTheta))
          {
              dim(iTheta) <- list(1,length(iTheta))
              colnames(iTheta) <- colnames(Y)
          }
      
      
      
      ## iTheta <- rbind(iFit[1:nrow(X), ])
      rownames(iTheta) <- rownames(X)
      
      secureTheta <- iTheta
      secureTheta[secureTheta < 0.001] <- 0.001
      checkError <- max(apply(abs(iTheta - oldTheta)/secureTheta, 2, max, na.rm = TRUE),
                        na.rm = TRUE)

      if (checkError < error.limit)
        break

      ## if(iter == maxit)
      ##   return(out_list)

    } ## End iteration

  if (verbose)
    cat(iter, "\n")

  attr(iTheta, "iter") <- iter
  if(attr.iLen)
    attr(iTheta, "iLen") <- p
  if(attr.df)
    attr(iTheta, "df") <- pUP$dfout
  ## NOT IMPLEMENTED FROM HERE

  if (std.err) {
    ##iSd <- rbind(sapply(colnames(Y),.fitIt,exLen=iLen,std.err=TRUE))
    if (iter > 1) {
      iSd <- iFit[-c(1:nrow(X)), ]
    }
    else {
      iSd <- ThetaFit[-c(1:nrow(X)), ]
    }
                                        #resvar <- colSums(iLm$res^2)/(ncol(Y)-1)
                                        #R <- chol2inv(iLm$qr$qr)
                                        #iCr <- sqrt(diag(R) * resvar)

    attr(iTheta, "expSd") <- iSd
                                        #attr(iTheta,'Cr') <- iLen
                                        #attr(iTheta,'CrSd') <- iCr
  }


  return(iTheta)

}

