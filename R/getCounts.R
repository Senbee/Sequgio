## We need these to deal with the problem of QNAMEs which differ between mates e.g. NAME1/1 & NAME1/2

readGAlignmentPairsFromBam_char <- function (file, index = file, use.names = FALSE, param = NULL,regesp=NULL) 
{
  if (missing(index) && (is.null(param) || 0L == length(bamWhich(param)))) 
    index <- character(0)
  bam <- open(BamFile(file, index), "rb")
  on.exit(close(bam))

  if(is.null(regesp))
    readGAlignmentPairsFromBam(bam, character(), use.names = use.names,param = param)
  else
    readGAlignmentPairsFromBam_bam(bam, character(), use.names = use.names,param = param, regesp=regesp)
}

readGAlignmentPairsFromBam_bam <- function (file, index = file, use.names = FALSE, param = NULL, regesp=NULL) 
{
  if (!isTRUEorFALSE(use.names)) 
    stop("'use.names' must be TRUE or FALSE")
  if (!isTRUE(obeyQname(file)) && !is.na(yieldSize(file))) {
    warning("'yieldSize' set to 'NA'")
    yieldSize(file) <- NA_integer_
  }
  flag0 <- scanBamFlag(isPaired = TRUE, hasUnmappedMate = FALSE)
  what0 <- c("flag", "mrnm", "mpos")
  param2 <- Rsamtools:::.normargParam(param, flag0, what0)
  galn <- readGAlignmentsFromBam(file, use.names = TRUE, param = param2)

  names(galn) <- gsub(regesp[1],regesp[2],names(galn))
  
  if (is.null(param)) {
    use.mcols <- FALSE
  }
  else {
    use.mcols <- c(bamWhat(param), bamTag(param))
  }
  makeGAlignmentPairs(galn, use.names = use.names, use.mcols = use.mcols)
}

################################################################################################


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

## .getBAM <- function(i,target,bam.params,Exons,counts,ignore.strand)
##   {
##     bamfile <- target[i,"filenames"]
##     baifile <- gsub(".bam$","",bamfile)
    
##     sampName <- target[i,"samplenames"]
    
##     sbv <- scanBam(bamfile,index=baifile,param=bam.params)
    
##     sqNames <- sbv[[1]][["rname"]]
    
##     if(with.junctions)
##       {
##         juncs <- grep("\\d+_\\d+$",sqNames)
        
##         iSR <- GRanges(seqnames=Rle(factor(sqNames[-juncs])),ranges=IRanges(start=sbv[[1]][["pos"]][-juncs],
##                                                                width=sbv[[1]][["qwidth"]][-juncs]),
##                        strand=strand(sbv[[1]][["strand"]][-juncs]))
##         iSR.j <- as.character(sqNames[juncs])
##       }
##     else
##       iSR <- GRanges(seqnames=Rle(sqNames),ranges=IRanges(start=sbv[[1]][["pos"]],
##                                              width=sbv[[1]][["qwidth"]]),
##                      strand=strand(sbv[[1]][["strand"]]))
    
##     rm(sbv)
    
    
##     ## subject = sample exons reads, query = reference exons Db
##     if(ignore.strand)
##       strand(Exons) <- "*"

##     ## i.e. => within! If reads have all the same length we get the same result.
##     ## NOTE: If 'type' is 'within', the QUERY interval must be wholly contained within the SUBJECT interval
##     ## That is SUBJECT => reference genome, QUERY = reads
    
##     OL <- findOverlaps(subject=Exons, query=iSR,type="within")

  
   
##     exonCounts <- tabulate(subH,length(Exons)) ## vector with counts for each exon
    
##     ## 6May11 - with disjoint exons we have repeated exons_id depending on how many Tx they map to
##     ## This makes it easier to compute design matrix but returns multiple counts here. Just get the first count
##     N.names <- as.character(paste(values(Exons)[["exon_name"]],values(Exons)[["gene_id"]],sep="."))
##     names(exonCounts) <- unique(N.names)
    
##     ## ## CHECK: why do we need this? Is it because of exons on _random_ chromosomes?
##     ## exonCounts <- exonCounts[names(exonCounts) %in% rownames(counts)]
    
##     counts[names(exonCounts),sampName] <- exonCounts
    
    
##   }


.getOverlaps <- function(sbv,Exons,ignore.strand,unique.only)
  {
    
    OL <- findOverlaps(sbv,Exons,type='within',ignore.strand=ignore.strand)
    
    if(unique.only)
      OL <- OL[!duplicated(queryHits(OL))]
    
    ww <- ngap(sbv) == 0

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
    
    ref_gr1 = cigarToIRangesListByAlignment(cigar(ref),start(ref))
    
    ## we avoid mapping gapped <-> non gapped (long ref exons)
    sel <- width(ref_gr1@partitioning) > 1
    OL3 <- OL3[sel]
    ref_gr1 <- ref_gr1[sel]
    rd1 <- rd1[sel]
    
    gr1 <- cigarToIRangesListByAlignment(cigar(rd1),start(rd1))
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


.getGapped <- function(i,target,bam.params,Exons,counts,ignore.strand,
                       mapq.filter,use.samtools,unique.only=FALSE,regesp=NULL,
                       only.proper=FALSE)
  {

    bamfile <- target[i,"filenames"]
    baifile <- target[i,"index"]
    
    sampName <- target[i,"samplenames"]



    if(use.samtools && (!is.null(mapq.filter) || only.proper ))
      {

        tmp <- tempfile(tmpdir=".",fileext=".bam")
        
        if(!is.null(mapq.filter) && only.proper)
            cmd <- paste('samtools view -f 2 -bq',mapq.filter,bamfile,'>',tmp)

        if(!is.null(mapq.filter) && !only.proper)
          cmd <- paste('samtools view -bq',mapq.filter,bamfile,'>',tmp)

        if(only.proper && is.null(mapq.filter))
          cmd <- paste('samtools view -f 2 -b',bamfile,'>',tmp)

        system(cmd)


        cmd <- paste('samtools index ',tmp)
        system(cmd)
        
        baifile <- bamfile <- tmp
        
      }

    
    if(any(is.na(baifile)))
        stop("BAM Index file missing. You need to provide them (see samtools)")
    
    sbv <- readGAlignmentPairsFromBam_char(file=bamfile,param=bam.params,index=baifile,regesp=regesp)

    
    

    if(exists(deparse(quote(tmp))))
      {
        ## remove temp files
        unlink(paste(baifile,'bai',sep='.'))
        unlink(bamfile)
      }


    if(!use.samtools && !is.null(mapq.filter))
      sbv <- sbv[values(sbv)$mapq >= mapq.filter]
    
    if(!use.samtools && only.proper)
      sbv <- sbv[isProperPair(sbv)]

    ## get leftmost & rightmost mates base on + strand position, i.e. the one reported in annotation
    ## - makes sure that right start > left start
    sbv_lr <- .left_right(sbv)
      
    OL.l <- .getOverlaps(sbv_lr[["left"]],Exons,ignore.strand=ignore.strand,unique.only=unique.only)
    OL.r <- .getOverlaps(sbv_lr[["right"]],Exons,ignore.strand=ignore.strand,unique.only=unique.only)



    ## Matrix holding the indeces of matches. The indeces refer to the rows id of Exons object
    mat.idx <- matrix(NA,ncol=2,nrow=length(sbv_lr[["right"]]),dimnames=list(NULL,c('left','right')))
    rm(sbv_lr)

    mat.idx[queryHits(OL.l),'left'] <- subjectHits(OL.l)
    mat.idx[queryHits(OL.r),'right'] <- subjectHits(OL.r)

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

    rm(mat.idx)

    mode(mat.ex) <- "character"
    exonCounts <- .Call("countEx",mat.ex,regid)
    exonCounts <- exonCounts[names(exonCounts) %in% rownames(counts)]
    
    counts[names(exonCounts),sampName] <- exonCounts
    
    ## ##########################################################################################################
    ## Exons is a list, with nogap and possibly gap. nogap goes on normally, gap need gapped alignment counting
    

    ## readlen <- min(width(sbv))
    ## ## -- Exons with no gaps
    
    ## OL <- findOverlaps(query=sbv,subject=Exons[['NoGaps']],type="within")
    ## ## OL <- findOverlaps(query=sbv,subject=Exons[['NoGaps']],minoverlap=readlen)
    ## exonCounts <- tabulate(subjectHits(OL),length(Exons[['NoGaps']]))
    ## N.names <- as.character(paste(values(Exons[['NoGaps']])[["exon_name"]],
    ##                               values(Exons[['NoGaps']])[["gene_id"]],sep="."))
    ## names(exonCounts) <- N.names
    ## exonCounts <- exonCounts[names(exonCounts) %in% rownames(counts)]
    ## counts[names(exonCounts),sampName] <- exonCounts

    ## ## -- Exons with gaps. Let's also consider only gapped reads. These can be the only one spanning a junction

    ## sbv.gaps <- sbv[ngap(sbv)>0]
    
    ## OL <- jointOverlaps(Exons[['Gaps']],sbv.gaps,ignore.strand=ignore.strand,minoverlap=minoverlap)
    ## exonCounts <- tabulate(queryHits(OL),length(Exons[['Gaps']]))
    ## N.names <- paste(values(Exons[['Gaps']])[["exon_name"]],values(Exons[['Gaps']])[["gene_id"]],sep=".")
    ## names(exonCounts) <- N.names
    ## counts[names(exonCounts),sampName] <- exonCounts

    ## ###############################################################################################################

  }

getCounts <- function(target,txdb,type="BAM",bam.params=NULL,minoverlap=5L,ignore.strand=TRUE,
                      mapq.filter=NULL,use.samtools=FALSE,unique.only=FALSE,regesp=NULL,mcpar,only.proper=FALSE)
  {

    if(!ignore.strand)
      {
        warning("'ignore.strand = FALSE' not yet implemented. Setting it to 'TRUE'")
        ignire.strand = TRUE
      }
    
    
    if(type != "BAM")
      stop("Only reading BAM files implemented")

    if(!all(file.exists(target$filenames)))
      stop('Some or all files specified do not exist in your filesystem. Please check the \'target\' content')

    if(missing(mcpar))
      mcpar <- registered()[[1]]

    if(use.samtools && !is.null(mapq.filter))
      {
        if(Sys.which('samtools')=="")
          {
            cat("samtools executable not found. Filter will be performed in R.")
            use.samtools <- FALSE
          }
        mapq.filter <- as.integer(mapq.filter)
      }
    
    Exons <- txdb@unlistData
    
    ## exons.names <- paste(values(Exons)[,"exon_name"],values(Exons)[,"region_id"],sep=".")
    
    with.junctions <- any(values(Exons)$ngaps > 0)


    ## Create vector with exon names combinations within tx.
    ##   - Exons already sorted according to rank
    ##   - Each exon combined all the following (in pairs)
    
    ex_list <- split(values(Exons)$exon_name,values(Exons)$tx_name)
    reg_vec <- sapply(split(values(Exons)$region_id,values(Exons)$tx_name),function(x) x[1])
    sizes_ex_list <- sapply(ex_list,length)
    n.exons <- sum((sizes_ex_list^2+sizes_ex_list)/2) ## total number of different exons names expected

    exons.names <- .Call("makeExNames",ex_list,reg_vec,as.integer(n.exons))

    ## The same pair in two or more Tx same region would be repeated.
    ## This is not efficient: can we identify them before makeExNames?
    exons.names <- unique(exons.names)
    n.exons <- length(exons.names)
    
    ## --
    
    n.samples <- nrow(target)
    sample.names <- target[,"samplenames"]
    
    options(bigmemory.allow.dimnames=TRUE)
    counts <- big.matrix(n.exons,n.samples,type="integer", init=0,
                         dimnames=list(exons.names,sample.names))

    if(is.null(bam.params) || !is(bam.params,'ScanBamParam'))
    {
      cat('bam.params NULL or not a ScanBamParam object. A default one will be used\n')

      ## This works for PE ends. Make sure it works for SE too
      
      bam.params <- ScanBamParam(simpleCigar = FALSE, reverseComplement = FALSE,
                                 what=c('qname',"qwidth",'mapq'),
                                 flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE,#isNotPrimaryRead=FALSE,
                                     isNotPassingQualityControls=FALSE))

      ## bam.params <- ScanBamParam(simpleCigar = FALSE, reverseComplement = FALSE,
      ##                            what=c("rname","strand","pos","qwidth",'mapq'),
      ##                            flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE,
      ##                                isNotPassingQualityControls=FALSE,
      ##                                isNotPrimaryRead=FALSE))
    }

    ## if(with.junctions)
    ##   {
    
    Exons.gap <- GAlignments(names=as.character(1:length(Exons)),
                             seqnames=seqnames(Exons),
                             strand=strand(Exons),
                             cigar=values(Exons)$cigar,pos=start(Exons),seqlengths=seqlengths(Exons),
                             elementMetadata(Exons)[,c("exon_name","exon_id","exon_rank","tx_id",
                                                       "tx_name","gene_id",'region_id')])
    
    
    invisible(bplapply(1:nrow(target),.getGapped,target,bam.params,Exons.gap,counts,ignore.strand,
                       mapq.filter=mapq.filter,use.samtools=use.samtools,unique.only=unique.only,regesp=regesp,
                       only.proper=only.proper,BPPARAM=mcpar))
    
    ##   }
    ## else
    ##   {
    ##   invisible(bplapply(1:nrow(target),.getGapped,target,bam.params,Exons,counts,ignore.strand,
    ##                      mapq.filter=mapq.filter,use.samtools=use.samtools,unique.only=unique.only,regesp=regesp,
    ##                      only.proper=only.proper,BPPARAM=mcpar))
    ## }
    
    return(counts)
    
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
