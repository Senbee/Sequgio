makeXmatrixSE <- function(object,use.tx.names=TRUE,exlen,probelen,mulen=probelen,mcpar,verbose=FALSE)
  {

    if(missing(mcpar))
      mcpar <- registered()[[1]]
    
    if(missing(object) || !is(object,'GRangesList'))
      stop('Input object must be of class \'GRangesList\'')
    
    if(!attributes(object)$reshaped)
      stop('Input object must be an object returned by \'reshapeTxDb\'')
    
    if(missing(exlen))
      stop("exlen argument required")

    if(missing(probelen))
      stop("probelen argument required")


    val <- getOption("stringsAsFactors")
    options(stringsAsFactors=FALSE)
    on.exit(options(stringsAsFactors=val))
    
    .local <- function(id,df,use.tx.names,exlen,probelen,mulen,verbose)
      {
        if(verbose)
          cat(id,"\n")
        
        obj <- subset(df,region_id == id)
        EX <- obj$exon_name
        EX.pos <- order(obj$start,obj$end)
        
        if(use.tx.names)
          TX <- factor(obj$tx_name)
        else
          TX <- factor(obj$tx_id)
        
        if(length(unique(TX))==1)
          {
            X <- rbind(rep(1,length(EX)))
            colnames(X) <- EX
            rownames(X) <- unique(TX)
          }
        else
          {
            X <- model.matrix(~0+TX)
            X <- aggregate(X,list(EX),sum)
            rownames(X) <- X[,1]
            X <- t(X[,-1,drop=FALSE])
            rownames(X) <- levels(TX)
          }

        
        X <- X[,unique(EX[EX.pos]),drop=FALSE]
        

        exons.position <- obj[,c('exon_name','start')]
        names(exons.position)[2] <- 'exon_pos'
        exons.position <- exons.position[!duplicated(exons.position$exon_name),,drop=FALSE]
        rownames(exons.position) <- exons.position$exon_name

        
        
        gnames <- paste(colnames(X),id,sep=".")
        exL <- exlen[[id]]
        exL <- exL[match(gnames,names(exL))]

        mEx <- t(X)*exL
        
        
        if(nrow(mEx) > 1) ## more then 1 exon
          {
            ## if(ncol(mEx) > 1)
            ##   {
            dEx <- apply(mEx[nrow(mEx):1, , drop = FALSE], 2, function(x) cumsum(x) - mulen)
            ##   }
            ## else {
            ##   dEx <- matrix(cumsum(mEx[nrow(mEx):1,1]) - mulen,ncol=1)
            ## }
          } else {

            dEx <- rbind(mEx[nrow(mEx):1,] - mulen)
            dimnames(dEx) <- dimnames(mEx)
          }
        
        selTx <- !as.logical(colMeans(dEx < 0)==1)

        if(all(!selTx)) ## all tx have length non compatible with fraglen
          {
            exons.position <- exons.position[colnames(X),]
            XL <- matrix(0,nrow=nrow(X),ncol=ncol(X),dimnames=dimnames(X))
            ans <- list(design=XL,condNumber=NA,exons.position=exons.position)
            return(ans)
          }
        
        if(!all(selTx))
          {
            X <- X[selTx,,drop=FALSE]
            selEx <- as.logical(colMeans(X))
            X <- X[,selEx,drop=FALSE]
            gnames <- gnames[selEx]
            exL <- exL[names(exL) %in% gnames]
            
          }

        dEx <- dEx[colnames(X),rownames(X),drop=FALSE]
        exons.position <- exons.position[colnames(X),]
        

        
        ## selRow <- apply(dEx >= 0,2,function(x) min(which(x)))
        ## last.exon <- nrow(dEx) - selRow + 1 ## last exon not fully in fragment length
        ## rest <- dEx[t(col(t(dEx)) == selRow)] + 1

        last.exon <- apply(dEx >= 0,2,function(x) max(which(x)))
        rest <- dEx[t(col(t(dEx)) == last.exon)] + 1
        
        fL <- (1-(probelen-1)/exL)
        
        XL <- matrix(fL,byrow=TRUE,nrow=nrow(X),ncol=ncol(X))*X

        ## N.B. last.exon is recycled row-major, that's why the following selection works!
        XL[col(X) == last.exon] <- rest/exL[last.exon]
        XL[col(X) > last.exon] <- 0
        XL[XL < 0] <- 0
        
        cN <- svd(XL)$d
        condN <- max(cN)/min(cN[cN > 0])
        
        ans <- list(design=XL,condNumber=condN,exons.position=exons.position)
        
        return(ans)
        
      }
    

    which.cols <- c("tx_id",'tx_name','exon_name','start','end','region_id')
    ## df <- IRanges:::as.data.frame(elementMetadata(object@unlistData))[,which.cols]
    df <- as.data.frame(object@unlistData)[,which.cols]
    lnames <- names(object)
    names(lnames) <- lnames

    
    out <- bpmapply(.local,lnames,
                    MoreArgs = list(df=df,use.tx.names=use.tx.names,exlen=exlen,probelen=probelen,
                        mulen=mulen,verbose=verbose),
                    USE.NAMES=TRUE,SIMPLIFY=FALSE,BPPARAM=mcpar)

    return(out)

  }
