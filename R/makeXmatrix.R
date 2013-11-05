## Time-stamp: <28-10-2013 17:49:11 on Masklin.med.unibs.it>


makeXmatrix <- function(object,method=c("SE","PE"),probelen,mulen,mcpar,...)
    {
        method <- match.arg(method)

        if(missing(mcpar))
          mcpar <- registered()[[1]]
        
        if(missing(object) || !is(object,'GRangesList'))
          stop('Input object must be of class \'GRangesList\'')
        
        if(!attributes(object)$reshaped)
          stop('Input object must be an object returned by \'reshapeTxDb\'')
        
        
        if(missing(probelen))
          stop("probelen argument required")
        
        thisFun <- match.call(expand.dots = TRUE)

        thisFun[[1]] <- as.name(paste(thisFun[[1]],method,sep=""))
        thisFun$method <- NULL
        
        if(method == "SE")
            {
                thisFun$sdlen <- NULL
                if(is.null(thisFun$exlen))
                    stop("exlen argument required")
            }
        else
            {
             thisFun$exlen <- thisFun$use.tx.names <- thisFun$verbose <- NULL
            }


        out <- eval(thisFun)
        return(out)
    }
