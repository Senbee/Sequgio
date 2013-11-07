
reshapeTxDb <- function(txdb,disjoin=TRUE,with.junctions=TRUE,probelen,ignore.strand=TRUE,junction.overlap=5L,
                        exclude.non.std=TRUE,exclude.mir=TRUE,what=c('exon','cds'),mcpar,verbose=FALSE,test.genes)
  {
    
    if(missing(mcpar))
      mcpar <- registered()[[1]]

    if(!ignore.strand)
      {
        warning("'ignore.strand = FALSE' not yet implemented. Setting it to 'TRUE'")
        ignire.strand = TRUE
      }
    
    val <- getOption("stringsAsFactors")
    options(stringsAsFactors=FALSE)
    on.exit(options(stringsAsFactors=val))

    if(missing(probelen))
      stop("Missing 'probelen' is not allowed.")
    
    ## Adapted from: GenomicFeatures:::.featuresBy

    if (!is(txdb, "TranscriptDb")) 
      stop("'txdb' must be a TranscriptDb object")

    ## to avoid problems with findOverlaps(...,'within'): not working with circular chr
    isActiveSeq(txdb)[isCircular(txdb)] <- FALSE
    
    if(exclude.non.std) ## exclude non standard chromosome: usually have "_" in their seqname
      isActiveSeq(txdb)[grep('_',seqlevels(txdb))] <- FALSE

    ## if(!missing(what.genome))
    ##   genome(txdb) <- what.genome
    
    short <- match.arg(what)

    selectClause <- 'SELECT short._short_id AS short_id, short_name, short_chrom, short_strand, short_start, short_end, gene_id, tx_name, splicing._tx_id AS tx_id, exon_rank'
    
    fromClause <- 'FROM short INNER JOIN splicing ON (short._short_id=splicing._short_id)  INNER JOIN gene ON (splicing._tx_id=gene._tx_id) INNER JOIN transcript ON (splicing._tx_id=transcript._tx_id)'
    
    whereClause <- 'WHERE tx_id IS NOT NULL AND gene_id IS NOT NULL'

    if(ignore.strand)
      orderByClause <- 'ORDER BY short_end, short_start, gene_id'
    else
      orderByClause <- 'ORDER BY exon_rank, gene_id'
    
    whereSeqsClause <- paste("AND short_chrom IN ('", paste(GenomicFeatures:::.getOnlyActiveSeqs(txdb), 
                                                           collapse = "','"), "')", sep = "")
    
    sql <- paste(selectClause, fromClause, whereClause, whereSeqsClause, 
                 orderByClause)

    sql <- gsub('short',short,sql)

    if(verbose)
      {
        cat("Reading in data from TranscriptDb: ")
        time <- proc.time()
      }
    
    geneGRs <- GenomicFeatures:::dbEasyQuery(AnnotationDbi:::dbConn(txdb), sql)

    if(short == 'cds')
      names(geneGRs) <- gsub('cds_','exon_',names(geneGRs))
    
    names(geneGRs)[c(4:6)] <- gsub("exon_","",names(geneGRs)[c(4:6)])
    names(geneGRs)[3] <- 'seqnames'
    geneGRs$width <- geneGRs$end - geneGRs$start + 1
    geneGRs <- geneGRs[,c("seqnames","start","end","width","strand","exon_id","exon_name",
                          "exon_rank","tx_id","tx_name","gene_id")]
    geneGRs$cigar <- paste(geneGRs$width,'M',sep='')
    geneGRs$ngaps <- 0
    geneGRs$exon_pos <- geneGRs$start ## again...we need to account for '-' strand    

    
    
    isActSeq <- isActiveSeq(txdb) ## here we should remove the 'random' chr
    seqlev <- names(isActSeq)[isActSeq]
    seqinfo <- seqinfo(txdb)[seqlev]
    strandlev <- c('+','-','*')

    txdf <- transcripts(txdb)

    ## This is very fragile and also relies on gene_ids being gene symbols. What if entrez codes?
    if(exclude.mir)
      {
        geneGRs <- subset(geneGRs,!grepl('^Mir',gene_id,ignore.case=TRUE))
        txdf <- txdf[values(txdf)$tx_name %in% geneGRs$tx_name]
      }


    redtx <- reduce(txdf,ignore.strand=TRUE)

    OL <- findOverlaps(query=redtx,subject=txdf,ignore.strand=TRUE)
    Group <- queryHits(OL)
    names(Group) <- values(txdf)$tx_name[subjectHits(OL)]
    
    Group.Ex <- Group[match(geneGRs$tx_name,names(Group))]

    geneGRs$region_id <- Group.Ex
    
    geneGRs <- mcsplit(geneGRs,Group.Ex,mcpar=mcpar)

    if(verbose)
      {
        delta <- proc.time() - time
        cat("DONE (",delta[3L],"secs)\n\n")
      }
    
    if(!missing(test.genes))
      geneGRs <- geneGRs[test.genes]

    
    if(disjoin)
      {
        if(verbose)
          {
            cat('Starting disjoin: ')
            time <- proc.time()
          }
        
        geneGRs <- mcDisjoin(genedb=geneGRs,ignore.strand = ignore.strand,
                             probelen=probelen,overlap.exons=1,mcpar=mcpar)
        if(verbose)
          {
            delta <- proc.time()-time
            cat("DONE (",delta[3L],"secs)\n\n")
          }
      }

    ##################################
    ## Now get biological junctions. #
    ##################################

    lnames <- names(geneGRs)
    names(lnames) <- lnames
    
    
    if(with.junctions)
      {

        if(verbose)
          {
            cat('Starting junction creation: ')
            time <- proc.time()
          }
        
        ## This returns a list of data.frames

        geneGRs <- bplapply(lnames,.getBioJuncs,inputGR=geneGRs,junction.overlap = junction.overlap,
                            probelen=probelen,BPPARAM=mcpar)

        if(verbose)
          {
            delta <- proc.time()-time
            cat("DONE (",delta[3L],"secs)\n\n")
          }
      }

    
    .dfToGRanges <- function(x)
      {

        subobj <- geneGRs[[x]]

        subobj$exLen <- cigarToQWidth(subobj$cigar)
        subobj <- subset(subobj,exLen >= 50)

        if(nrow(subobj) == 0)
          return(NA)
        
        id1 <- as.numeric(gsub("^ex_(\\d+).*","\\1",subobj$exon_name))
        id2 <- suppressWarnings(as.numeric(gsub("^ex_\\d+-ex_(\\d+)$","\\1",subobj$exon_name)))
        id2[is.na(id2)] <- 0
        subobj <- subobj[order(id1,id2),,drop=FALSE]
        
        subobj$rank <- as.integer(factor(subobj$exon_name,levels=unique(subobj$exon_name)))
        
        return(GRanges(seqnames=factor(subobj$seqnames,levels=seqlev),
                       ranges=IRanges(subobj$start,subobj$end),
                       strand=factor(subobj$strand,levels=strandlev),subobj[,-c(1:5)]))
      }

    if(verbose)
      {
        cat('Create GRanges objects: ')
        time <- proc.time()
      }
    
    ans.geneGRs <- bplapply(lnames,.dfToGRanges,BPPARAM=mcpar)
    ans.geneGRs <- ans.geneGRs[!is.na(ans.geneGRs)]
    
    if(verbose)
      {
        delta <- proc.time()-time
        cat("DONE (",delta[3L],"secs)\n\n")

        cat('Create GRangesList object: ')
        time <- proc.time()
      }

    ## This is helpful in order to get unlistedData afterwards
    ansGRs <- GRangesList(ans.geneGRs)

    if(verbose)
      {
        delta <- proc.time()-time
        cat("DONE (",delta[3L],"secs)\n\n") 
      }

    seqinfo(ansGRs) <- seqinfo
    
    attr(ansGRs,'reshaped') <- TRUE

    if(verbose)
      cat("Processing Complete\n\n")

    return(ansGRs)

  }

## TODO: disjoining has to be done within gene & seqnames, because some gene have exon mapped to different positions
## in the genome => see GFF/GTF files
mcDisjoin <- function(genedb,ignore.strand=FALSE,probelen,overlap.exons,mcpar)
  {
    geneId <- names(genedb)
    names(geneId) <- geneId

    .local <- function(x)
      {
        out.rng <- .Disjoin(genedb[[x]],ignore.strand=ignore.strand,probelen=probelen,
                            overlap.exons=overlap.exons)
        out.rng$width <- out.rng$end - out.rng$start + 1
        out.rng$cigar <- paste(out.rng$width,'M',sep='')
        out.rng$ngaps <- 0
        out.rng$exon_pos <- out.rng$start

        ####################################################################################
        ## This is a temporary remedy to small 'region' size. Must be fixed in a better way
        
        ## out.rng <- subset(out.rng,width > probelen)
        
        ####################################################################################

        out.rng
      }

    ans <- bplapply(geneId,.local,BPPARAM=mcpar)
    ln <- unlist(bplapply(ans,nrow,BPPARAM=mcpar))
    ans <- ans[ln >= 1]
    return(ans)
    
  }


## Actually ignore.strand doesn't do anything here
.Disjoin <- function (x,ignore.strand = FALSE, probelen,overlap.exons, ...) 
{

  ## From: getMethod(disjoin,signature=c(x="Ranges"))
  .local <- function (x) 
    {
      ffun <- function(x,...)
        {
   
          
          thisIR <- IRanges(start=x$start,end=x$end)
          starts <- unique(x$start)
          ends <- unique(x$end)
          adj_start <- sort(unique(c(starts, ends + 1L)))
          adj_end <- sort(unique(c(ends, starts - 1L)))
          adj <- IRanges(start = head(adj_start, -1L), end = tail(adj_end,-1L))
          UL <- findOverlaps(adj,thisIR,maxgap = 0L, minoverlap = 1L,type="any")
          iNames=NULL

          if(length(subjectHits(UL)) > 0)
            {

              ansIR <- x[subjectHits(UL),]
              ansIR$exon_name <- paste("ex",as.integer(factor(queryHits(UL))),sep="_")
              ansIR$start <- start(adj[queryHits(UL)])
              ansIR$end <- end(adj[queryHits(UL)])


              tabExTx <- table(paste(ansIR$exon_id,ansIR$tx_id,sep="_"))

              ## get artificial junctions (after disjoin)
              if(any(tabExTx > 1))
                {
                  nn <- names(tabExTx)[tabExTx > 1]
                  obj <- ansIR[paste(ansIR$exon_id,ansIR$tx_id,sep='_') %in% nn,]
                  addIR <- lapply(nn,getJunctions,inobj=obj,probelen=probelen,
                                  overlap.exons=overlap.exons)
                  ansIR <- do.call('rbind',c(list(ansIR),addIR))
                }

              ansIR <- ansIR[order(ansIR$start,ansIR$end),]
              return(ansIR)
                
            } else return(x)
        }
      
      return(lapply(x,ffun))
      
    }
  
  
  if(ignore.strand)
    xIRangesList <- split(x,x$seqnames)
  else
    xIRangesList <- split(x,list(x$seqnames,x$strand))
  

  xIRangesList <- xIRangesList[sapply(xIRangesList,nrow) > 0]

  ansIRangesList <- .local(xIRangesList, ...)

  out <- do.call('rbind',ansIRangesList)
  rownames(out) <- NULL

  return(out)
  
}


getJunctions <- function(iname,inobj,probelen,overlap.exons)
  {
    idx <- unlist(strsplit(iname,'_'))
    object <- inobj[inobj$exon_id == idx[1] & inobj$tx_id == idx[2],]
    rng <- IRanges(start=object$start,end=object$end)
    ovl <- probelen - overlap.exons + 1
    ll <- length(rng) ## how many exons
    howmuch <- probelen-1
    new.rng <- resize(flank(rng[-1],howmuch),width=howmuch*2)

    ## ISSUE: What if three exons are consecutive with the central one (much) smaller than probelen?
    ## |----------| ex1
    ##            |---| ex2 
    ##                |--------------| ex3
    ##          |.......| a read

    ## To account for starting/ending exons smaller than 24bp (e.g. gene "100019")
    new.rng <- restrict(new.rng,start=start(rng)[-length(rng)],end=end(rng)[-1])

    out <- object[rep(1,length(new.rng)),]
    out$start <- start(new.rng)
    out$end <- end(new.rng)
    out$width <- width(new.rng)
    
    out$exon_name <- paste(object$exon_name[1:(ll-1)],object$exon_name[2:ll],sep="-")

    rownames(out) <- NULL

    return(out)
  }


.getBioJuncs <- function(id,inputGR,junction.overlap,probelen)
  {
    
    .local <- function(x,junction.overlap,probelen)
      {

        if(length(x) == 1 || length(unique(x$exon_id))==1) ## only 1 exon
          return(NULL)

        rlength <- probelen - junction.overlap

        jNomi <- x$exon_name[grepl("ex_\\d+-ex_\\d+",x$exon_name)]
        subtmp <- subset(x,!grepl("ex_\\d+-ex_\\d+",x$exon_name)) ## remove junctions
        subtmp <- subtmp[order(subtmp$start,subtmp$end,subtmp$exon_name),]

        thisIR <- IRanges(subtmp$start,subtmp$end)
        ll <- nrow(subtmp)
        jNames <- paste(subtmp$exon_name[-ll],'-',subtmp$exon_name[-1],sep="")
        
        ans <- subtmp[-1,]

        ans$start <- subtmp$end[-ll] - (rlength - 1)
        ans$end <- subtmp$start[-1] + (rlength - 1)
        
        iWidths <- width(thisIR)

        ## Why this? CHECK!!!!
        pos <- ans$start + iWidths[-ll] - rlength + 1

        
        iWidths[iWidths > rlength] <- rlength
        
        iLeft <- iWidths[-ll]
        iRight <- iWidths[-1]
        
        thisCigar <- paste(iLeft,'M', width(gaps(thisIR)),'N',iRight,'M',sep='')

        ans$exon_name <- jNames

        ans$cigar <- thisCigar
        ans$ngaps <- 1

        ans$exon_pos <- pos

        ans <- ans[!jNames %in% jNomi,]
        
        return(ans)
      }


    if(nrow(inputGR[[id]])==1) ## only 1 exon
      return(inputGR[[id]])

    
    ansGRs <- split(inputGR[[id]],list(inputGR[[id]]$seqnames,inputGR[[id]]$tx_id))
    ansGRs <- ansGRs[sapply(ansGRs,nrow) > 0]
    addGRs <- lapply(ansGRs,.local,junction.overlap=junction.overlap,probelen=probelen)

    if(length(addGRs))
      {
        ansGR <- do.call('rbind',c(inputGR[id],addGRs))
        rownames(ansGR) <- NULL
        return(ansGR)
      }

    return(inputGR[[id]])
    
  }


## Parallel version of split.data.frame
mcsplit <- function (x, f, drop = FALSE,mcpar,...)
  {
    bplapply(split(seq_len(nrow(x)), f, drop = drop,...), function(ind) x[ind,, drop = FALSE],BPPARAM=mcpar)
  }
