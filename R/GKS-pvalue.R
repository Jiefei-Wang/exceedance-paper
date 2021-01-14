#' Compute the pvalue of generalized Kolmogorov-Smirnov test statistics
#' 
#' Compute the pvalue of generalized Kolmogorov-Smirnov test statistics,
#' see details on how to use the function.
#' 
#' @param x Numeric, the sample data.
#' @param stat a `generalKSStat` object or numeric, the statistic that the p-value
#' is computed for. If `stat` is a `generalKSStat` object, all other parameters will 
#' be ignored. Otherwise, this parameter will be suppressed if the parameter `x` 
#' is not null.
#' If a numeric value is provided to the parameter `stat`, you must at least 
#' specify the sample size `n`.
#' @param n Integer, the sample size of the data.
#' 
#' @examples 
#' ## Generate samples
#' x <- rbeta(10, 1, 2)
#' 
#' 
#' ## Perform KS test
#' ks_res <- GKSStat(x = x, statName = "KS")
#' 
#' ## Compute the pvalue for the KS test
#' GKSPvalue(stat = ks_res)
#' 
#' ## For any observed statistic
#' GKSPvalue(stat = 0.2, n = 10, statName = "KS")
#' 
#' ## Change the detection range of the KS test
#' ## to test only the first 3 ordered samples
#' ## All gives the same result
#' GKSPvalue(stat = 0.2, n = 10, alpha0 = 0.3, statName = "KS")
#' GKSPvalue(stat = 0.2, n = 10, index = 1:3, statName = "KS")
#' GKSPvalue(stat = 0.2, n = 10, indexL = 1:3, indexU = 1:3, statName = "KS")
#' 
#' 
#' 
#' @return A numeric value representing the pvalue
#' @inheritParams GKSStat
#' @inherit GKSStat details
#' @rdname pvalue
#' @export
GKSPvalue<-function(stat=NULL , n=NULL, 
                    index = NULL,
                    indexL=NULL,indexU=NULL,
                    x=NULL, statName = NULL){
    
    if(is(stat,"GKSStat")){
        statName <- getStatName(stat)
        statValue <- getStatValue((stat))
        n <- getSampleSize(stat)
        indexL <- stat$indexL
        indexU <- stat$indexU
    }else{
        statName <- match.arg(statName,
                              c("KS","KS+","KS-","BJ","BJ+","BJ-",
                                "HC","HC+","HC-","Simes")
                              )
        if(!is.null(x)){
            n <- length(x)
            stat <- GKSStat(x=x,
                            index= index, indexL=indexL,indexU=indexU,
                            statName = statName, pvalue = FALSE)
            statValue <- getStatValue(stat)
        }else{
            stopifnot(!is.null(n))
            statValue <- stat
        }
        
    }
    stopifnot(!is.null(n))
    idx <- getGKSIndex(statName = statName, n = n, 
                       index= index, indexL = indexL, indexU = indexU)
    indexL <- idx$indexL
    indexU <- idx$indexU
    statName <- idx$statName
    pvalue <- call_func(root = "Pvalue", prefix = statName,
                        statValue=statValue,n=n,indexL=indexL,indexU=indexU)
    
    return(pvalue)
}




uniformProbability<-function(LocalCriticalFunc,statValue ,n,indexL,indexU){
    localCritical <- LocalCriticalFunc(statValue=statValue,n=n)
    l=localCritical$l
    h=localCritical$h
    if(length(indexL)!=0){
        l[-indexL] <- 0
    }else{
        l=rep(0,length(l))
    }
    if(length(indexU)!=0){
        h[-indexU] <- 1
    }else{
        h=rep(1,length(h))
    }
    orderedProb(l,h)
}


HCPvalue<-function(statValue,n,indexL=seq_len(n),indexU=seq_len(n)){
    1-uniformProbability(LocalCriticalFunc=HCLocalCritical,
                        statValue=statValue,n=n,indexL=indexL,indexU =indexU)
}


BJPvalue<-function(statValue,n,indexL=seq_len(n),indexU=seq_len(n)){
    1-uniformProbability(LocalCriticalFunc=BJLocalCritical,
                        statValue=statValue,n=n,indexL=indexL,indexU =indexU)
}

KSPvalue<-function(statValue,n,indexL=seq_len(n),indexU=seq_len(n)){
    1-uniformProbability(LocalCriticalFunc=KSLocalCritical,
                        statValue=statValue,n=n,indexL=indexL,indexU =indexU)
}
SimesPvalue<-function(statValue,n,indexL=seq_len(n),indexU=seq_len(n)){
    statValue
}
