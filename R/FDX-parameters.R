#######################
## Dispatching rule:
## if postfix_XX is not null, postfix_XX will first be used to dispatch
## otherwise, dispatch based on postfix
#######################

#' Fast Exceedance control using the GW method
#' 
#' Performing a fast exceedance control of the false discovery proportion(FDP)
#' under the framework proposed by Genovese, C., & Wasserman, L. (2004),
#' where FDP is defined by the number of false positives devided by the 
#' number of total rejections.
#' The GW method requires an uniform(0,1) distributional test statistic
#' and each value in the input data represents a p-value from a hypothesis.
#' The sources of the data can be derived from any testing procedure
#' (e.g. pvalues from testing high-throughput gene data).
#' This function only supports specific uniform(0,1) distributional tests.
#' please consider `param_general_GW` if you want a generic version of the GW method, 
#' 
#' @param statistic character, the name of the statistic for the test, see details.
#' @param param1 integer, the first parameter for the statistic, see details.
#' @param param2 integer, the second parameter for the statistic, see details. 
#' @param range_type character, the type of the parameters, see details.
#' 
#' @details 
#' \emph{Note: We will use the term `samples` and `pvalues` interchangebly to refer
#' the data gathered by an inference procedure.}
#' 
#' This function perform the fast GW algorithm, currently it supports 
#' the following test statistics:
#' \itemize{
#' \item kth_p: The kth pvalue statistic
#' \item KS: The Kolmogorov-Smirnov statistic
#' \item HC: The higher criticism statistic
#' \item BJ: The Berk-Jones statistic
#' }
#' 
#' For each statistic, you can specify the index of the ascending ordered samples
#' to control which data will be considered in the test statistics. 
#' For example, the index of the kth pvalue statistic determines the value of `k`
#' and use the kth smallest sample as its statistic.
#' similarly, the index of KS, HC and BJ determines which ordered samples 
#' will be used to compute the test statistic. 
#' 
#' By default, `range_type = "index"`, which means `param1` and `param2` represent 
#' the index. However, the index can also depend on the sample size. Therefore,  
#' if `range_type = "proportion"`, The index is determined by the formula 
#' `max(floor(n*param),1)`, where `n` is the sample size. 
#' 
#' For the kth pvalue statistic, `param1` determine the value of `k` and must be 
#' a single integer or a 0-1 value depending on `range_type`.
#' 
#' For KS, HC and BJ, the formula to decide the index is a little bit complicated.
#' If `range_type = "index"`, `param1` determines which 
#' small sample(s) will be considered as the evidence of significance. 
#' For example, if `param1 = 2` and the second smallest sample is significantly small, 
#' it can lead to a significant result. Conversely,
#' `param2` determines which large sample(s) can be treated as significance.
#' Both `param1` and `param2` can be a vector of integer. By default,
#' if both `param1` and `param2` are null, it is equivalent to `param1 = c(0, 1)`, 
#' `param2 = NULL` and `range_type = "proportion"`.
#' 
#' If `range_type = "proportion"`, `param1` and `param2` can be length 1 vectors, 
#' which will be explained as the index from `1` to `max(floor(n*param),1)`, or they can
#' be length 2 vectors, where the index ranges from `max(floor(n*param[1]),1)` to 
#' `max(floor(n*param[2]),1)`.
#' 
#' @examples 
#' ## The 3rd pvalue statistic
#' param_fast_GW(statistic = "kth_p", param1 = 3)
#' 
#' ## One-sided KS statistic
#' param_fast_GW(statistic = "KS", param1 = c(0,1), range_type = "proportion")
#' 
#' ## One-sided KS statistic, 
#' ## Test first 10 smallest pvalues only
#' param_fast_GW(statistic = "KS", param1 = 1:10, range_type = "index")
#' 
#' @return an exceedance_parameters object
#' @export
param_fast_GW<-function(statistic = c("kth_p",
                                      "KS","HC","BJ","Simes"), 
                        param1=NULL,
                        param2=NULL,
                        range_type=c("index","proportion")){
    statistic <- match.arg(statistic)
    range_type <- match.arg(range_type,c("index","proportion"))
    
    if(statistic%in%c("KS","HC","BJ")){
        if(range_type == "proportion"){
            stopifnot(length(param1)<=2)
            stopifnot(length(param1)<=2)
            if(length(param1)==1){
                param1 <- c(0,param1)
            }
            if(length(param2)==1){
                param2 <- c(0,param2)
            }
        }else{
             param1 <- as.numeric(sort(param1))
             param2 <- as.numeric(sort(param2))
        }
        if(length(param1)==0&&length(param2)==0){
            param1 <- c(0,1)
            param2 <- c(0,1)
            range_type <- "proportion"
        }
        profile_func <- "profile_fast_GW_order_general"
        confidence_func <- "confidence_fast_GW_order_general"
        inference_func <- "inference_fast_GW_order_general"
    }
    if(statistic =="kth_p"){
        if(is.null(param1)){
            range_type <- "index"
            param1 <- 1L
        }else{
            stopifnot(length(param1)==1)
            stopifnot(is.null(param2))
        }
        if(range_type=="index"){
            profile_func <- "profile_fast_GW_kth_p_index"
            confidence_func <- "confidence_fast_GW_kth_p_index"
            inference_func <- "inference_fast_GW_kth_p_index"
        }else{
            profile_func <- "profile_fast_GW_kth_p_proportion"
            confidence_func <- "confidence_fast_GW_kth_p_proportion"
            inference_func <- "inference_general"
        }
    }
    if(statistic =="Simes"){
        param1 <- c(0,1)
        param2 <- c(0,1)
        range_type <- "proportion"
        profile_func <- "profile_fast_GW_order_general"
        confidence_func <- "confidence_fast_GW_order_general"
        inference_func <- "inference_fast_GW_order_general"
    }
    param <- .exceedance_parameter(method = "fast_GW",
                  statistic = statistic,
                  param1 = param1,
                  param2 = param2,
                  range_type=range_type,
                  profile_func = profile_func,
                  confidence_func = confidence_func,
                  inference_func = inference_func)
    param
}

#' General exceedance control using the GW method
#' 
#' Performing a general exceedance control of the false discovery proportion(FDP)
#' under the framework proposed by Genovese, C., & Wasserman, L. (2004),
#' where FDP is defined by the number of false positives devided by the 
#' number of total rejections.
#' The GW method requires an uniform(0,1) distributional test statistic
#' and each value in the input data represents a p-value from a hypothesis.
#' The sources of the data can be derived from any testing procedure
#' (e.g. pvalues from testing high-throughput gene data).
#' This function works for all distributional test statistics.
#' There are fast algorithms exist for some statistics. 
#' Please see the function `param_fast_GW` for more details.
#' 
#' @param pvalue_func function, The pvalue function for a uniform(0,1)
#' distributional test. The inpute samples(pvalues) will be ascending sorted.
#' @param algorithm character, The searching algorithm that will be used to find
#' the upper bound of the FDR, see details.
#' 
#' @details 
#' \emph{Note: We will use the term `samples` and `pvalues` interchangebly to refer
#' the data gathered by a large-scale inference procedure.}
#' 
#' Computing the upper bound of the FDP using the GW method is computationally 
#' challenging.
#' Given the number of hypotheses `n`, the naive algorithm requires approximately `2^n` 
#' operations, which make it only available for small inference. Therefore, 
#' JW algorithm is developed for reducing the computational burden, but it also
#' add a restriction on the test statistic.
#' 
#' If `general = "general"`, it means the naive method proposed by GW.
#' There is no additional assumptions on the pvalue function except the 
#' regular definition of the pvalue.
#' 
#' If `general = "JW"`, the pvalue function must be a non-decreasing function
#' for its input samples. That is, if input1 and input2 are two set of 
#' ascending sorted samples and input1 is not smaller than the input2 for all 
#' elements, then the resulting pvalue for input1 should not be smaller than the 
#' pvalue for input2.
#' 
#' With an additional restriction on the pvalue function, the JW algorithm has 
#' a polynomial computational complexity.
#' 
#' @examples 
#' ## kth pvalue test
#' k_test <- function(x,k){
#' n <- length(x)
#' if(length(x) >= k)
#'     pbeta(x[k], k, n-k+1)
#' else 
#'     1
#' }
#' 
#' ## general algorithm, very slow
#' param_general_GW(pvalue_func = k_test, algorithm = "general")
#' 
#' 
#' ## JW algorithm, relatively qucik
#' ## You need to verify if the pvalue function
#' ## satisfies the requirement of the JW algorithm
#' ## before use it.
#' param_general_GW(pvalue_func = k_test, algorithm = "JW")
#' 
#' @return an exceedance_parameters object
#' @export
param_general_GW<-function(
    pvalue_func,
    algorithm = c("general","JW")){
    algorithm <- match.arg(algorithm)
    if(algorithm == "general"){
        profile_func = "profile_general_GW_general"
        confidence_func = "confidence_general_GW_general"
    }else{
        profile_func = "profile_general_GW_JW"
        confidence_func = "confidence_general_GW_JW"
    }
    inference_func <- "inference_general"
    
    
    param <- .exceedance_parameter(method = "general_GW",
                  pvalue_func = pvalue_func,
                  statistic = algorithm,
                  profile_func = profile_func,
                  confidence_func = confidence_func,
                  inference_func = inference_func)
    param
}

#' @export
param_combine<-function(...,param_list = NULL, alpha_weight = NULL){
    if(!is.null(param_list)){
        test_params <- param_list
    }else{
        test_params <- list(...)
    }
    if(is.null(alpha_weight)){
        alpha_weight <- rep(1,length(test_params))
    }
    param <- .exceedance_parameter(method = "combine_GW",
                  test_params = test_params,
                  statistic = "combine_GW",
                  alpha_weight=alpha_weight,
                  profile_func = "profile_combine_GW",
                  confidence_func = "confidence_combine_GW",
                  inference_func = "inference_general"
                  )
    param
}

