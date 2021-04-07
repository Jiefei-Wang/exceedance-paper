all.null<-function(...){
    args<-list(...)
    res <- vapply(args,is.null,logical(1))
    all(res)
}

get_index_from_proportion<-function(n,param){
    if(is.null(param)) return(param)
    param <- ceiling(param*n)
    param[param==0] <- 1
    seq_len(param[2]-param[1]+1) + as.integer(param[1] -1L)
}

## Get the ordered index for the ordered x
## such that x[index] == ordered_x[ordered_index]
## ri: rejected i
## sri: sorted rejected i
## rx: rejected x
get_ordered_index<-function(x_rank,ri=NULL,sri = NULL,rx=NULL){
    if(all.null(ri,sri,rx))
        stop("The rejection set must be provided")
    if(!is.null(sri)){
        sorted_i <- sri
    }else{
        if(!is.null(rx)){
            ri <- which(x%in%rx)
        }
        sorted_i <- x_rank[ri]
    }
    sorted_i
}
######################################################
## GW: general JW algorithm
######################################################
get_local_critical <-function(statName, n, critical, indexL,indexU){
    level <- do.call(paste0(statName,"LocalCritical"),args = list(
        stat= critical,
        n=n
    ))
    level2 <- process_local_critical(level, indexL,indexU)
    level2
}

## modify the local criticals according to the index
process_local_critical <- function(bound,indexL,indexU){
    l <- bound$l
    h <- bound$h
    n <- length(l)
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
    for(i in seq_len(n-1)){
        if(l[i]>l[i+1]) l[i+1] <- l[i]
        j <- n-i
        if(h[j] > h[j+1]) h[j] <- h[j+1]
    }
    bound$h <-h
    bound$l <- l
    bound
}
#################################
##   critical cache
#################################

get_index_key <- function(n,index){
    if(length(index)==0){
        return("")
    }
    if(length(index)==1){
        return(as.character(index))
    }
    if(length(index)==n){
        return("*")
    }
    key <- NULL
    if(length(index)==index[length(index)]-index[1]+1){
        key <- paste0(index[1],"-",index[length(index)])
    }else{
        key <- paste0(length(index), ":",
                      digest::digest(index, algo = "crc32")
        )
    }
    key
}

get_cache_key <- function(statName, n, alpha, indexL, indexU){
    keyL <- get_index_key(n,indexL)
    keyU <- get_index_key(n,indexU)
    cache_key <- paste0(statName,",", n, ",", alpha, "L", keyL,"U",keyU)
    cache_key
}
exist_cache_value <- function(cache_key){
    load_criticals()
    exists(cache_key,envir=pkg_data$criticals)
}
get_cache_value <- function(cache_key){
    load_criticals()
    pkg_data$criticals[[cache_key]]
}
set_cache_value <- function(cache_key, value){
    load_criticals()
    pkg_data$criticals[[cache_key]] <- value
}

compute_key_critical <- function(statName, n, alpha, indexL, indexU){
    if(length(indexL)!=0&&max(indexL)>n){
        indexL <- indexL[indexL<=n]
    }
    if(length(indexU)!=0&&max(indexU)>n){
        indexU <- indexU[indexU<=n]
    }
    if(length(indexU)==0&&length(indexL)==0){
        return(NULL)
    }
    cache_key <- get_cache_key(statName = statName, n=n, 
                               alpha=alpha, indexL=indexL, indexU=indexU)
    if(use_cache()&&exist_cache_value(cache_key)){
        if(verbose()){
            message("use cache value")
        }
        critical <- get_cache_value(cache_key)
    }else{
        if(verbose()){
            message("compute criticals")
        }
        critical <- GKSCritical(n=n,alpha=alpha,indexL=indexL,indexU=indexU,statName=statName)
        set_cache_value(cache_key, critical)
    }
    list(key = cache_key, critical = critical)
}


## Get critical value for BJ, KS, HC...
## If the critical has been cached, we will get it from cache
get_critical <- function(statName, n, alpha, indexL, indexU){
    compute_key_critical(statName, n, alpha, indexL, indexU)$critical
    
}