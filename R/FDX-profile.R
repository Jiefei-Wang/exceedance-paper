#' Profile the data
#' 
#' Preprocessing the data and collecting information before performing 
#' the exceedance control. This step is developed for reducing the 
#' computational burder only. It is not related to statistics.
#' 
#' @param x numeric, p-values from some tests
#' (E.g T-test for high-thoughput gene data.), with each element representing a 
#' hypothesis.
#' @param param exceedance_parameters object
#' 
#' @inherit exceedance_inference examples
#' @return an exceedance_profile object
#' @export
exceedance_profile<-function(x, param){
    m <- length(x)
    sx <- sort(x,index.return = TRUE)
    x_rank <- rep(0,m)
    x_rank[sx$ix] <- seq_len(m)
    basic_info <- list(
        m = m,
        x = x,
        x_sort = as.numeric(sx$x),
        x_sort_index = sx$ix,
        x_rank = x_rank
    )
    
    profile_func <- get(param$profile_func)
    
    profile <- profile_func(
        x=x,
        param = param,
        profile= basic_info
    )
    
    profiled_data <- .exceedance_profile(
        param = param, profile = profile
    )
    profiled_data
}


#pvalue_func <- function(x)ks.test(x,punif)$p.value
profile_general_GW_general<-function(x,param,profile){
    pvalue_func <- param$pvalue_func
    m <-length(x)
    x_sort <- profile$x_sort
    ## The set that contains all possible combination
    ## of the data at c++ level
    ## print_subset_list(unreject_set)
    unreject_set <- general_GW_construct_subset(pvalue_func,parent.frame(),x_sort)
    profile <- c(profile,
                 list(unreject_set=unreject_set))
    
    profile
}


profile_general_GW_JW<-function(x,param,profile){
    cache <- list()
    cache$search_path <- new.env()
    cache$pvalues <- new.env()
    
    profile$cache <- cache
    profile
}



profile_fast_GW_kth_p_index <- function(x, param,profile){
    k <- param$param1
    m <- profile$m
    x_sort <- profile$x_sort
    
    if(k>length(x)){
        stop("the order k(param1) cannot be greater than the sample size")
    }
    
    local_level <- pbeta(x_sort,k,m-seq_len(m)+1)
    max_alpha <- max(local_level[k:m])
    
    profile <- c(profile,
               list(
                   local_level=local_level,
                   max_alpha=max_alpha
               ))
    profile
}

profile_fast_GW_kth_p_proportion <- function(x, param,profile){
    k <- param$param1
    m <- profile$m
    x_sort <- profile$x_sort
    local_level <- pbeta(x_sort,k,m-seq_len(m)+1)
    
    profile<-c(profile,
               list(
                   local_level=local_level
               ))
    profile
}





profile_fast_GW_order_general<-function(x,param,profile){
    range_type <- param$range_type
    param1 <- param$param1
    param2 <- param$param2
    statistic <- param$statistic
    
    profile$params_key <- digest::digest(list(range_type,param1,param2,statistic))
    profile
}

profile_combine_GW<-function(x,params,profile){
    test_params <- params$test_params
    profile$profiles <- lapply(test_params,function(param)exceedance_profile(x=x,param=param))
    profile
}


