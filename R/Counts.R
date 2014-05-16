## Time-stamp: <15-05-2014 17:07:25 on Masklin.med.unibs.it>

## 1) input is BamFile: e.g. BamFile(fl,asMate=TRUE,yieldsize=10^5) for pair-end
## 2) Will use parallel to work along samples
## TODO:
## - use BatchJobs


## Create the matrix of counts.

setClass(
    Class = "seqCounts",
    slots = c(counts = "ANY", exonNames = 'integer', sampleNames = 'character',
        files = "ANY", Exons = "GAlignments"),
    validity = function(object)
    {
        ## class(object@counts) %in% c("ff_matrix","big.matrix") &&
        ##     class(object@files) %in% c("BamFile","BamFileList")


        is.character(object@counts) &&
            class(object@files) %in% c("BamFile","BamFileList")
    }
    )


setMethod(
    f = "initialize",
    signature = "seqCounts",
    definition = function(.Object,counts,exonNames,sampleNames,files,Exons,template)
    {
        .Object@counts = counts
        .Object@files = files
        .Object@Exons = Exons
        .Object@exonNames = exonNames
        .Object@sampleNames = sampleNames
        
        validObject(.Object)

        return(.Object)
    }

    )


setMethod(
    f = "show",
    signature = "seqCounts",
    definition = function(object)
    {
        cat("BAM file(s):\n")
        show(object@files)
        cat("\nCounts object files:\n")
        f2 <- object@counts[1]
        f2 <- gsub('desc$','bin',f2)
        cat('Description file:',object@counts[1],"\n")
        cat('Binary file:',f2,"\n")
        ln <- length(object@sampleNames)
        cat("\nSample name(s) (N=",ln,"): ",sep="")
        if(ln<5)
            cat(paste(object@sampleNames,collapse=", "),"\n")
        else
            cat(paste(head(object@sampleNames,2),collapse=", "),',..., ',tail(object@sampleNames,1),"\n",sep="")

        ln <- length(object@exonNames)
        cat("\nExon/region name(s) (N=",ln,"): ",sep="")
        cat(paste(head(names(object@exonNames),2),collapse=", "),',..., ',tail(names(object@exonNames),1),"\n\n",sep="")
        
        
    })

setGeneric("resetCounts",
           def=function(Object,...){standardGeneric("resetCounts")})

setMethod("resetCounts",signature(Object="seqCounts"),
          definition =
          function(Object)
          {
              counts <- attach.big.matrix(Object@counts)
              counts[,] <- 0L
              flush(counts)
          })

setCounts <- function(bfl,txdb,sample.names=NULL,fileName=NULL,rootpath=getwd())
    {

        if(!(is(bfl,'BamFile') || is(bfl,'BamFileList')))
            stop("Argument 'bfl' must be either a 'BamFile' or a 'BamFileList' object")
        
        n.samples <- length(bfl) ## bfl=BamFile/BamFileList
        
        if(!is.null(names(path(bfl))))
            {
                if(!is.null(sample.names))
                    warning("Using BamFile/BamFileList file names")
                
                sample.names <- names(path(bfl))
            }
        else
            {
                if(is.null(sample.names))
                    sample.names <- gsub("\\.bam$","",basename(path(bfl)),ignore=TRUE)
                
                if(is(bfl,'BamFile'))
                    names(bfl$path) <- sample.names
                else
                    names(bfl@listData) <- sample.names
            }
        
        if(is.null(fileName))
            fileName <- basename(tempfile(tmpdir=getwd()))
        
        Exons <- txdb@unlistData
        
        with.junctions <- any(values(Exons)$ngaps > 0)
        
        
        ## Create vector with exon names combinations within tx.
        ##   - Exons already sorted according to rank
        ##   - Each exon combined all the following (in pairs)
        
        ex_list <- split(values(Exons)$exon_name,values(Exons)$tx_name)
        reg_vec <- sapply(split(values(Exons)$region_id,values(Exons)$tx_name),function(x) x[1])
        
        exons.names <- .Call("makeExNames",ex_list,reg_vec)
        
        n.exons <- length(exons.names)
        exons.vec <- seq.int(n.exons) - 1 ## zero based indexing in C++
        storage.mode(exons.vec) <- 'integer'
        names(exons.vec) <- exons.names
        
        
        counts.file <- file.path(rootpath,paste(fileName,'.desc',sep=''))
        
        options(bigmemory.allow.dimnames=FALSE)
        counts <- big.matrix(n.exons,n.samples,type="integer", init=0,
                             backingfile=paste(fileName,'.bin',sep=''),
                             descriptorfile=paste(fileName,'.desc',sep=''),
                             binarydescriptor=TRUE,shared=TRUE,
                             backingpath=rootpath)
        
        
        Exons <- txdb@unlistData
        
        Exons.gap <- GAlignments(names=as.character(1:length(Exons)),
                                 seqnames=seqnames(Exons),
                                 strand=strand(Exons),
                                 cigar=values(Exons)$cigar,pos=start(Exons),seqlengths=seqlengths(Exons),
                                 elementMetadata(Exons)[,c("exon_name","exon_id","exon_rank","tx_id",
                                                           "tx_name","gene_id",'region_id')])
        
        ans <- new("seqCounts", counts = counts.file, exonNames=exons.vec,
                   sampleNames = sample.names,
                   files = bfl, Exons = Exons.gap)
        
        return(ans)
        
    }




setMethod('yieldSize',signature(object='seqCounts'),
          definition=function(object)
          {
              yieldSize(object@files)
          }
          )

setMethod('yieldSize<-',signature(object='seqCounts'),
          definition=function(object,...,value)
          {
              if (1L != length(value)) 
                  stop("'value' must be length 1")
              object@files$yieldSize <- value
              object
          }
          )


setGeneric("getCounts",
           def=function(Object,...){standardGeneric("getCounts")})


setMethod("getCounts",signature(Object="seqCounts"),
          definition =
          function(Object)
          {

              BigMat <- attach.big.matrix(Object@counts)
              counts <- BigMat[,,drop=FALSE]
              dimnames(counts) <- list(names(Object@exonNames),Object@sampleNames)
              return(counts)
          })

setGeneric("doCounts",
           def=function(Object,...){standardGeneric("doCounts")})


setMethod("doCounts",signature(Object="seqCounts"),
          definition =
          function(Object,bam.params=NULL,minoverlap=5L,ignore.strand=TRUE,
                   which.sample=NULL,mapq.filter=NULL,unique.only=FALSE,mcpar,verbose=FALSE)
          {
              

              if(!is.null(which.sample) && !(is.integer(which.sample) ||  is.character(which.sample)))
                  stop("Argument \"which.sample\" must be a vector of either integers (sample #) or characters (sample names)")
              
              if(missing(mcpar))
                  mcpar <- registered()[[1]]
              
              if(!ignore.strand)
                  {
                      warning("'ignore.strand = FALSE' not yet implemented. Setting it to 'TRUE'")
                      ignore.strand = TRUE
                  }
              

              if(!is.null(mapq.filter))
                  mapq.filter <- as.integer(mapq.filter)


              
              Exons <- Object@Exons
              bfl <- Object@files

              exonNames <- Object@exonNames

              pathFiles <- path(bfl)
              
              colIDs <- seq.int(length(pathFiles))
              colIDs <- setNames(colIDs,names(pathFiles))
              
              if(length(pathFiles) > 1 && !is.null(which.sample))
                  {
                      if(length(pathFiles) < length(which.sample))
                          stop("seqCounts object length < 'which.sample' length: please check!")

                      if(is.character(which.sample) && !all(which.sample %in% names(pathFiles)))
                          stop("Sample names in argument 'which.sample' do not match sample names in 'seqCounts' object: please check!")

                      
                      colIDs <- colIDs[which.sample]
                      
                  }

              open(bfl)
              on.exit(close(bfl))




              
              if(is.null(bam.params) || !is(bam.params,'ScanBamParam'))
                  {
                      cat('bam.params NULL or not a ScanBamParam object. A default one will be used\n')
                      
                      ## This works for PE ends. Make sure it works for SE too
                      
                      bam.params <- ScanBamParam(simpleCigar = FALSE, reverseComplement = FALSE,
                                                 what=c('qname',"qwidth",'mapq'),
                                                 flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE,
                                                     isNotPassingQualityControls=FALSE))

                      if(verbose)
                          {
                              ll <- bamFlag(bam.params)
                              print(ll[!is.na(ll)])
                              rm(ll)
                          }
                  }
              
              


              if(is(bfl,"BamFile"))
                  invisible(.getGapped(colIDs,bfl=bfl,bam.params=bam.params,
                                       Exons=Exons,counts=counts,ignore.strand=ignore.strand,
                                       mapq.filter=mapq.filter,unique.only=unique.only,
                                       seqObj=Object@counts,exonNames=exonNames
                                       ))
              else
                  invisible(bplapply(colIDs,.getGapped,
                                     bfl=bfl,bam.params=bam.params,
                                     Exons=Exons,counts=counts,ignore.strand=ignore.strand,
                                     mapq.filter=mapq.filter,unique.only=unique.only,
                                     seqObj=Object@counts,exonNames=exonNames,
                                     BPPARAM=mcpar))
          })
    

.left_right <- function(x,...)
{
  x_first <- x@first
  x_last <- x@last
  left_is_last <- which(start(x_first) > start(x_last))
  idx <- seq_len(length(x))
  idx[left_is_last] <- idx[left_is_last] + length(x)

  ans_l <- c(x_first, x_last)[idx]
  ans_r <- c(x_last, x_first)[idx]

  setNames(ans_l, names(x))
  setNames(ans_r, names(x))

  return(list(left=ans_l,right=ans_r))

}



.getGapped <- function(i,bfl,bam.params,Exons,counts,ignore.strand,mapq.filter,unique.only=FALSE,
                       seqObj,exonNames=exonNames)
    {
        ## which column in big.matrix are we going to consider? 0-base index
        sampID <- as.integer(i-1)
        
        if(is(bfl,'BamFileList'))
            iBfl <- bfl[[i]]
        else
            iBfl <- bfl
        
        counts <- attach.big.matrix(seqObj)
        
        .local <- function(object)
            {

                sbv <- suppressWarnings(readGAlignmentPairsFromBam(iBfl,param=bam.params))

                flushDumpedAlignments() ## just in case!
                
                if(!is.null(mapq.filter))
                    {
                        if(is(sbv,"GAlignments"))
                            sbv <- sbv[values(sbv)$mapq >= mapq.filter]
                        else ## CHECK are the first & last mapq values always the same? If so no need to test both
                            sbv <- sbv[values(sbv@first)$mapq >= mapq.filter & values(sbv@last)$mapq >= mapq.filter]
                    }
                
                ## get leftmost & rightmost mates base on + strand position, i.e. the one reported in annotation
                ## - makes sure that right start > left start
                sbv_lr <- .left_right(sbv)
                
                OL.l <- .getOverlaps(sbv_lr[["left"]],Exons,ignore.strand=ignore.strand,unique.only=unique.only)
                OL.r <- .getOverlaps(sbv_lr[["right"]],Exons,ignore.strand=ignore.strand,unique.only=unique.only)
                
                
                
                ## Matrix holding the indeces of matches. The indeces refer to the rows id of Exons object
                mat.idx <- matrix(NA,ncol=2,nrow=length(sbv_lr[["right"]]),dimnames=list(NULL,c('left','right')))
                rm(sbv_lr)
                
                ## assign at every row both the right and left matches. Row numbers identify the read id
                mat.idx[queryHits(OL.l),'left'] <- subjectHits(OL.l)
                mat.idx[queryHits(OL.r),'right'] <- subjectHits(OL.r)
                
                rownames(mat.idx) <- seq_len(nrow(mat.idx)) ## so we can keep track of the reads ids
                
                mat.idx <- na.omit(mat.idx)
                
                
                ## We do consider only reads mapping within the same region. We should check these though.
                regid <- values(Exons)$region_id[mat.idx[,1]]
                regid2 <- values(Exons)$region_id[mat.idx[,2]]
                
                if(any(regid != regid2))
                    {
                        rsel <- which(regid==regid2)
                        mat.idx <- mat.idx[rsel,,drop=FALSE]
                        regid <- regid[rsel]
                    }
                
                
                mat.ex <- cbind(values(Exons)$exon_name[mat.idx[,1]],values(Exons)$exon_name[mat.idx[,2]])
                
                mode(mat.ex) <- "character"
                
                ## For Dhany's debugging
                ## debmat <- cbind(values(sbv@first)$qname[as.integer(rownames(mat.idx))],paste(mat.ex[,1],mat.ex[,2],sep="."),
                ##                 regid)
                ## colnames(debmat) <- c("qname","exon","region")
                ## assign("debugMe",debmat,env=.GlobalEnv)
                
                rm(mat.idx)

                ## Count reads
                totReads = countExBM(mat.ex,regid,exonNames,sampID,object@address)
                
            }

        while(isIncomplete(iBfl))
            .local(counts)


    }


.getOverlaps <- function(sbv,Exons,ignore.strand,unique.only)
  {
    
    OL <- findOverlaps(sbv,Exons,type='within',ignore.strand=ignore.strand)
    
    if(unique.only)
      OL <- OL[!duplicated(queryHits(OL))]
    
    ww <- njunc(sbv) == 0

    ## When data are not gapped we just stop here
    if(all(ww))
        return(OL)

    
    ok1 <- which(ww)
    OL2 <- OL[queryHits(OL) %in% ok1] ## Hits mapping reads without gap

    
    ok2 <- which(!ww) ## gapped reads in the annotattion (Exons)
    OL3 <- OL[queryHits(OL) %in% ok2]
    
    sel1 <- queryHits(OL3)
    
    rd1 <- sbv[sel1]
    
    ref <- Exons[subjectHits(OL3)]
    
    ref_gr1 = extractAlignmentRangesOnReference(cigar(ref),start(ref))
    
    ## we avoid mapping gapped <-> non gapped (long ref exons)
    sel <- width(ref_gr1@partitioning) > 1
    OL3 <- OL3[sel]
    ref_gr1 <- ref_gr1[sel]
    rd1 <- rd1[sel]
    
    gr1 <- extractAlignmentRangesOnReference(cigar(rd1),start(rd1))
    en <- start(gr1@unlistData)[end(gr1@partitioning)]-1
    st <- end(gr1@unlistData)[start(gr1@partitioning)]+1
    gp1 <- IRanges(start=st,end=en) ## reads gaps
    
    ref_en = start(ref_gr1@unlistData)[end(ref_gr1@partitioning)]-1
    ref_st = end(ref_gr1@unlistData)[start(ref_gr1@partitioning)]+1
    ref_gp1 = IRanges(start=ref_st,end=ref_en)
    sel = which(start(gp1)==start(ref_gp1) & end(gp1)==end(ref_gp1))
    
    OL3 <- OL3[sel]

    ansH <- new("Hits")
    ansH@subjectHits <- c(subjectHits(OL2),subjectHits(OL3))
    ansH@queryHits <- c(queryHits(OL2),queryHits(OL3))
    
    ## subH <- rbind(cbind(subjectHits(OL2),queryHits(OL2)),
    ##               cbind(subjectHits(OL3),queryHits(OL3)))
    ## colnames(subH) <- c("subjectHits","queryHits")
    return(ansH)
    
  }

getWidth <- function(x)
  {
    ww <- cigarWidthAlongQuerySpace(values(x@unlistData)$cigar)
    gnomi <- paste(values(x@unlistData)$exon_name,values(x@unlistData)$region_id,sep=".")
    names(ww) <- gnomi

    ww <- ww[!duplicated(gnomi)]
    iGene <- values(x@unlistData)$region_id
    iGene <- iGene[!duplicated(gnomi)]
    ww <- ww[order(iGene)]

    tab <- table(iGene)
    
    thisPart = PartitioningByEnd(cumsum(tab))
    thisPart@NAMES <- names(tab)
    
    new2("CompressedIntegerList", unlistData = ww,
         partitioning = thisPart, check = FALSE)
    
  }



jointOverlaps <- function (features, reads,minoverlap = 5L,ignore.strand=FALSE)
{
  
  .local <- function (features, reads,  minoverlap = 5L, ignore.strand = FALSE) 
    {

      featRng <- ranges(features)
      featGaps <- unlist(gaps(featRng),use.names = FALSE)


      ## TODO: what if a reads is something like: 1M1220N22M6375N2M?
      ## see sbv.gaps[7549] first simulated sample e.g.
      ## for the moment discard them!
      
      allL <- elementLengths(reads)
      reads <- reads[allL == 2]

      minSide <- min(width(reads))
      reads <- reads[minSide >= minoverlap]
      
      readRng <- ranges(reads)
      readGaps <- unlist(gaps(readRng),use.names = FALSE)

      readStrand <- unlist(runValue(strand(reads)))
      featStrand <- unlist(runValue(strand(features)))

      readSeq <- unlist(runValue(seqnames(reads)))
      featSeq <- unlist(runValue(seqnames(features)))

      readGR <- GRanges(readGaps,seqnames=readSeq,strand=readStrand)
      featGR <- GRanges(featGaps,seqnames=featSeq,strand=featStrand)

      findOverlaps(featGR,readGR,type="equal",ignore.strand=ignore.strand)
      
    }
  
  .local(grglist(features), grglist(reads), minoverlap)
}
