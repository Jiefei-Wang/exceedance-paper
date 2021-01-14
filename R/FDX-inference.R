## dispatch rule: if postfit is not null, 
## dispatch according to the postfix
## Otherwise, dispatch to the general algorithm

#' Perform multiple hypothesis testing while controlling the exceedance rate
#' 
#' Perform multiple hypothesis testing while controlling the exceedance rate. 
#' The exceedance rate is Prob(FDP > bound), where FDP is the false discover 
#' proportion defined by the number of false positives devided by the number 
#' of total rejections. Small FDP with large number of reject is favorable 
#' in practice.
#' 
#' @inheritParams exceedance_confidence
#' @param bound The upper bound of the false discover proportion
#' @examples 
#' ## The 3rd pvalue statistic
#' param <- param_fast_GW(statistic = "kth_p", param1 = 3)
#' 
#' ## generate p-values
#' x <- rbeta(10, 1, 10)
#' 
#' ## profile the data
#' profile <- profile_pvalue(x,param)
#' 
#' 
#' ## compute the 95% confidence envolop
#' alpha <- 0.05
#' 
#' ## reject the first three hypotheses
#' exceedance_confidence(profile, alpha, ri = 3)
#' 
#' ## reject the hypothese which pvalues are equal to
#' ## the first three samples.
#' ## In other word, this is equivalent to reject the first three hypotheses
#' exceedance_confidence(profile, alpha, rx = x[1:3])
#' 
#' ## reject the hypotheses which have the lowest 3 p-values
#' exceedance_confidence(profile, alpha, sri = 3)
#' 
#' 
#' ## Determine which hypotheses can be rejected while controlling the
#' ## exceedance rate: P(FDP > bound) < alpha
#' alpha <- 0.05
#' bound <- 0.2
#' exceedance_inference(profile, alpha, bound)
#' 
#' @return Indices of the hypotheses that are rejected in the procedure.
#' @export
exceedance_inference<-function(profiled_data, alpha, bound){
    inference_func <- get(profiled_data$param$inference_func)
    result <- inference_func(
        profiled_data = profiled_data, alpha = alpha,
        bound=bound)
    result
}

inference_general<-function(profiled_data, alpha, bound){
    profile <- profiled_data$profile
    param <- profiled_data$param
    m <- profile$m
    x_sort_index <- profile$x_sort_index
    
    reject <- integer(0)
    for(j in rev(seq_len(m))){
        gammabar = exceedance_confidence(profiled_data,alpha,sri = 1L:j)
        if(gammabar<=bound){
            reject<-1L:j
            break
        }
    }
    x_sort_index[reject]
}

inference_fast_GW_kth_p_index<-function(profiled_data, alpha, bound){
    profile <- profiled_data$profile
    param <- profiled_data$param
    max_alpha <- profile$max_alpha
    x_sort_index <- profile$x_sort_index
    k<- param$param1
    m <- profile$m
    
    if(alpha>max_alpha){
        smallest_FDR <- (k-1)/m
        if(smallest_FDR > bound){
            return(integer(0))
        }else{
            return(x_sort_index)
        }
    }
    
    inference_general(profiled_data, alpha, bound)
}


## Given the samples and subset size n, compute
## the CI for all possible rejection sets
## R is a list of reject sets
get_fast_GW_subset_all_FDX <- function(statName, sx, n, indexL, indexU, alpha){
    m <- length(sx)
    FDP <- rep(0, m)
    critical <- get_critical(statName= statName, n=n, alpha=alpha, 
                             indexL=indexL, indexU=indexU)
    bound <- get_local_critical(statName = statName, n= n, critical=critical,
                                indexL=indexL,indexU=indexU)
    
    # x_range <- get_range_by_bound(sx=x_sort,bound=bound)
    x_range <- C_get_range_by_bound2(R_sx=sx,R_l=bound$l,R_h=bound$h)
    if(!is.null(x_range)){
        P <- x_range$P
        for(i in seq_len(m)){
            FDP[i] <- sum(P<=i)/i
        }
    }
    FDP
}


inference_fast_GW_order_general<-function(profiled_data, alpha, bound){
    profile <- profiled_data$profile
    param <- profiled_data$param
    range_type <- param$range_type
    param1 <- param$param1
    param2 <- param$param2
    m <- profile$m
    x_sort_index <- profile$x_sort_index
    x_sort<- profile$x_sort
    statistic <- param$statistic
    
    FDR <- rep(0, m)
    for(n in seq_len(m)){
        if(range_type=="proportion"){
            indexL <- get_index_from_proportion(n=n,param=param1)
            indexU <- get_index_from_proportion(n=n,param=param2)
        }else{
            indexL <- param1[param1<=n]
            indexU <- param2[param2<=n]
            ## Check if the current sample size satiefies
            ## The minimum requirement of the test.
            if(length(indexL)==0&&length(indexU)==0){
                FDR <- pmax(FDR, pmin(n,seq_len(m)) / seq_len(m))
                next
            }
        }
        cur_FDR <- get_fast_GW_subset_all_FDX(statName = statistic, sx = x_sort, n=n,
                                             indexL = indexL, indexU = indexU, alpha = alpha)
        FDR <- pmax(FDR, cur_FDR)
    }
    idx <- which(FDR<=bound)
    if(length(idx)!=0){
        x_sort_index[seq_len(max(idx))]
    }else{
        integer(0)
    }
}

