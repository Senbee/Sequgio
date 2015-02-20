#########################################
##                                     ##
## Get summarized expression by "gene" ##
##                                     ##
## Created: 18 Feb 2015                ##
## Last modified:                      ##
#########################################

setGeneric("sumByGene",
           def=function(Object,TxDb,...){standardGeneric("sumByGene")})

setMethod("sumByGene",signature(Object="list",TxDb="data.frame"),
          definition =
              function(Object,TxDb)
                  {
                      
                      TxDb <- TxDb[!duplicated(TxDb$tx_name),]
                      rownames(TxDb) <- TxDb$tx_name
                      
                      bmat <- do.call('rbind',Object)
                      TxDb <- TxDb[rownames(bmat),]
                      
                      bmat <- data.table(bmat,Gene=TxDb$gene_id)
                      sumExpr <- bmat[,lapply(.SD,sum),by=Gene]
                      
                      return(sumExpr)
                      
                  }
          )



setMethod("sumByGene",signature(Object="list",TxDb="GRangesList"),
          definition =
              function(Object,TxDb)
                  {
                      
                      df <- elementMetadata(txdb@unlistData)[,c('tx_name','gene_id')]
                      df <- df[!duplicated(df$tx_name),]
                      rownames(df) <- df$tx_name
                      
                      bmat <- do.call('rbind',Object)
                      df <- df[rownames(bmat),]
                      
                      bmat <- data.table(bmat,Gene=df$gene_id)
                      sumExpr <- bmat[,lapply(.SD,sum),by=Gene]
                      
                      return(sumExpr)
                      
                  }
          )

