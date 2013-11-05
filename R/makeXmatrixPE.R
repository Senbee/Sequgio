## Paired-end mapping fragments 
##   computing expected contributions of ONE transcript
#
## transcript info:
## loc = location of boundaries of *consecutive* and *non-overlapping* regions, 
##       each can be exon, junction, or both
## nr = number of regions, within this computation 
##      they will be numbered 1 to nr,
##      but they will be mapped back into the original
##      region codes in the transcription unit# 


sort_ids <- function(x)
  {
    id1 <- as.numeric(gsub("^ex_(\\d+).*","\\1",x))
    id2 <- suppressWarnings(as.numeric(gsub("^ex_\\d+-ex_(\\d+)$","\\1",x)))
    id2[is.na(id2)] <- 0
    return(x[order(id1,id2)])
  }

plen1 = function(x,y2,rlen,mulen,sdlen)
  {
    a = pnorm(y2-x+rlen,mulen,sdlen) - pnorm(2*rlen,mulen,sdlen)
    return(a)
  }

plen2 = function(x,y2,y1,rlen,mulen,sdlen)
  {
    a = pnorm(y2-x+rlen,mulen,sdlen) - pnorm(y1-x+rlen,mulen,sdlen)
    return(a)
  }


makeXmatrixPE <- function(object,probelen=50,mulen=300,sdlen=30,mcpar)
  {
    
    
    if(missing(mcpar))
      mcpar <- registered()[[1]]
    
    if(missing(object) || !is(object,'GRangesList'))
      stop('Input object must be of class \'GRangesList\'')
    
    if(!attributes(object)$reshaped)
      stop('Input object must be an object returned by \'reshapeTxDb\'')
    
    if(missing(probelen))
      stop("probelen argument required")
    
    
    val <- getOption("stringsAsFactors")
    options(stringsAsFactors=FALSE)
    on.exit(options(stringsAsFactors=val))


    
    ####################################################################################################

    .local <- function(id,df,rlen,mulen,sdlen)
      {
        
        subobj <- subset(df,region_id == id)
        ##obj <- obj[order(obj$start,obj$end),]
        ## id1 <- as.numeric(gsub("^ex_(\\d+).*","\\1",subobj$exon_name))
        ## id2 <- suppressWarnings(as.numeric(gsub("^ex_\\d+-ex_(\\d+)$","\\1",subobj$exon_name)))
        ## id2[is.na(id2)] <- 0
        ## subobj <- subobj[order(id1,id2),,drop=FALSE]
        
        if(!"exLen" %in% names(subobj))
          subobj$exLen <- cigarToQWidth(subobj$cigar)

        ## ## Let's create the all possible pairwise
        ## cNames <- unique(subobj$exon_name)
        ## all_names <- outer(cNames,cNames,FUN=paste,sep=".")
        ## all_names <- all_names[upper.tri(all_names,diag=TRUE)]
        
        df.split <- split(subobj[,c("exon_name","tx_name","exLen","rank")],subobj$tx_name)
        

        nsplt <- names(df.split)
        names(nsplt) <- nsplt

        ffun <- function(itx,inobj,rlen,mulen,sdlen,id)
          {
            ans <- tryCatch(.computeX_R(itx,inobj=df.split,rlen=rlen,mulen=mulen,sdlen=sdlen,id=id),
                            error=function(e) e)

            ## ans <- tryCatch(.Call("computeX",itx,df.split,rlen,mulen,
            ##                       sdlen,cnames),error=function(e) e)

            if('error' %in% class(ans))
              {
                msg = gsub("'","",ans$message)
                txt = paste("Error: ",msg,". Transcript=",itx,". Region=",id,". Process=",Sys.getpid(),sep="")
                cat(txt,"\n")
                return(NA)
              }
            else
              return(ans)
            
          }

        ## every element is a 2-x-n_exons matrix. First row are the weights, second row the indicators
        res <- lapply(nsplt,ffun,inobj=df.split,rlen=rlen,mulen=mulen,sdlen=sdlen,id=id)

        res_names <- unique(unlist(lapply(res,colnames)))

        res_names <- sort_ids(res_names)


        
        dmat <- matrix(0,nrow=length(res),ncol=length(res_names),dimnames=list(names(res),res_names))
        ans <- array(dim=c(nrow(dmat),ncol(dmat),2),dimnames=list(rownames(dmat),colnames(dmat),c("Design","Indicator")))
        ans[,,1] <- dmat

        for(i in 1:length(res))
          {
            ans[i,colnames(res[[i]]),1] <- res[[i]][1,]
            ans[i,colnames(res[[i]]),2] <- res[[i]][2,]
          }
        
        ## Using C++ we already iterate within region
        ## res <- ffun(itx=nsplt,inobj=df.split,rlen=rlen,mulen=mulen,sdlen=sdlen,cnames=cNames)

        txLen <- rowSums(ans[,,1,drop=FALSE])

        if(all(txLen == 0))
          return(NA)

        if(any(txLen==0))
          {
            ans <- ans[txLen > 0,,,drop=FALSE]
            txLen <- txLen[txLen > 0]
          }
        
        ans[,,1] <- ans[,,1]/txLen
        attr(ans,"txLen") <- txLen

        return(ans)
        
      }
    
    
####################################################################################################
    
    which.cols <- c("tx_id",'tx_name','exon_name','start','end','region_id','cigar','rank','exLen')
    df <- GenomicRanges::as.data.frame(object@unlistData)[,which.cols]
    lnames <- names(object)
    names(lnames) <- lnames

    out <- bplapply(lnames,FUN=.local,df=df,
                    rlen=probelen,mulen=mulen,sdlen=sdlen,
                    BPPARAM=mcpar)

    if(any(is.na(out)))
      out <- out[!is.na(out)]
    
    return(out)

  }


.computeX_R <- function(itx,inobj,rlen,mulen,sdlen,cnames,id)
  {
    obj <- inobj[[itx]]
    exLen <- obj$exLen
    
    regs = cumsum(c(0,exLen-rlen+1))
    nr <- length(exLen) ## number of regions
    
    
    frac <- rev((cumsum(rev(exLen)) - 2*rlen + 1)/rev(exLen))
    frac[frac < 0] <- 0
    frac <- pmin(frac,(exLen-rlen+1)/exLen)
    
    loc <- unique(cumsum(exLen*frac))
    
    ## nloc <= nr. If last region(s) is(are) smaller than 2*rlen, it doesn'count for integration
    nloc <- length(loc) 
    
    mloc <- loc-rlen ## point WITHIN regions rlen from 'upper boundary'
    
    ## for regions as long as rlen we actually get the boundary. Remove from mloc
    mloc <- unique(mloc[!mloc %in% loc])
    
    
    ## x on rows, y on cols => the marginal (sums along the columns in Fig 1) is computed summing along ROWS!
    res = matrix(0,nr,nr,dimnames=list(obj$exon_name,obj$exon_name))

    xloc <- sort(c(0,loc,mloc))
    yloc <- xloc+rlen
    
    
    if(max(regs) < max(yloc))
      {
        yloc <- yloc[yloc <= max(regs)]
        yloc[length(yloc)] <- max(regs)
      }
    
    ## Loop along i = x-index (regions)
    for (i in 1:nloc)
      {
        ## lc1 & lc2: boundaries of region #1
        if (i==1)
          lc1=0
        else
          lc1=loc[i-1]
        
        lc2 = loc[i]
        pick = (xloc >=lc1) & (xloc<= lc2)
        xloci = xloc[pick]  # x locations of the same index
        
        ## x intervals for integration within region nr. Note: xloci be longer then 3,
        ## if effective length of the second region is smaller than rlen
        ## xloci => loc_inf - mloc ( - possibly another mloc) - loc_sup
        for (ii in seq(along=head(xloci,-1)))
          {
            xlim1 = xloci[ii]
            xlim2 = xloci[ii+1]
            
            ## iterate along the column, i.e. y-coordinates
            
            ## first y-coordinate correspond to triangle regions => plen1
            yloci <- yloc[yloc >= xlim1+rlen]
            
            for(idy in seq(along=head(yloci,-1)))
              {
                ylim1 <- yloci[idy]
                ylim2 <- yloci[idy+1]
                yind <- findInterval(ylim1,regs,all.inside=TRUE)
                
                if(idy == 1)
                  prop =integrate(plen1,xlim1,xlim2,ylim2,rlen,mulen,sdlen)$val
                else
                  prop =integrate(plen2,xlim1,xlim2,ylim2,ylim1,rlen,mulen,sdlen)$val
                
                
                res[i,yind] = res[i,yind] + prop
              }
          }
      }



    
    res_names <- outer(colnames(res),colnames(res),FUN=paste,sep=".")
    res_names <- res_names[upper.tri(res_names,diag=TRUE)]

    ind_res <- row(res)
    
    res <- res[upper.tri(res,diag=TRUE)]
    ind_res <- ind_res[upper.tri(ind_res,diag=TRUE)]
    names(res) <- names(ind_res) <- paste(res_names,id,sep="__")
    
    return(rbind(res,ind_res))
    
    ## This is the marginal 'mimicking' a single end situation.
    ## This would deliver a single row (tx) in the single end design matrix

    ## res <- rowSums(res)/exLen ## return marginal for the moment
    
    ## ans <- rep(0,length(cnames))
    ## names(ans) <- cnames
    ## ans[match(obj$exon_name,cnames)] <- res
    
    ## return(ans)
  }
