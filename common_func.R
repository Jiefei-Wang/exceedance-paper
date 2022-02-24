
test_summary <- function(rejected_idx, hypo_num, H0_num, FDP_bound){
    n_reject <- max(length(rejected_idx),1)
    
    H1_num <- hypo_num-H0_num
    FDP <- sum(rejected_idx<=H0_num)/n_reject
    FDX <- FDP>FDP_bound
    power <- sum(rejected_idx>H0_num)/max(H1_num, 1)
    setNames(c(FDP,FDX,power), c("FDP", "FDX", "power"))
}



## Get combined kth order statistic parameter
get_combined_param<-function(n){
    param_individual <- list()
    for(j in 1:n){
        param_individual[[j]] <- 
            param_fast_GW(statistic = "kth_p", param1 = j, 
                          range_type = "index")
    }
    param_combine(param_list = param_individual)
}


sample_pvalue <- function(nrep, n, mu, n_null, sigma) {
    if(!all(sigma==1)){
        n_alt <- n - n_null
        1-pnorm(rmvnorm(nrep, c(rep(0,n_null),rep(mu,n_alt)),sigma, method= "chol"))
    } else {
        data <- matrix(0, nrep, n)
        for (i in 1:nrow(data))
            data[i,] <- runif(1)
        if(n_null!=n){
            data[,(n_null+1):n] <- data[,(n_null+1):n] + mu
        }
        data
    }
}

#compute the covariance matrix
getCovMat <- function(n, params, type = c("ind", "cs", "ar1", "toeplitz")){
    type <- match.arg(type)
    covMat <- matrix(0, n, n)
    if (type == "ind") {
        for(i in 1:n) {
            covMat[i,i] <- 1
        }
    }
    if (type == "cs") {
        for (i in 1:n) {
            for (j in 1:n) {
                if (i == j) 
                    covMat[i,j] <- 1
                else 
                    covMat[i,j] <- params
            }
        }
    }
    if (type == "ar1") {
        for (i in 1:n) {
            for (j in 1:n) {
                covMat[i,j] <- params ^ abs(i-j)
            }
        }
    }
    if (type == "toeplitz") {
        for (i in 1:n) {
            for (j in 1:n) {
                if (i == j) 
                    covMat[i,j] <- 1
                else if (abs(i-j) == 1)
                    covMat[i,j] <- params
            }
        }
    }
    covMat
}


collapse_list <- function(x){
    do.call(rbind, x)
}


KR_fdp <- function(x_sort, alpha){
    pmin(1, sapply(seq_along(x_sort), function(k){
        floor(log(1/alpha)/log(1+log(1/alpha))*(1+length(x_sort)*x_sort[k]))
    })/seq_along(x_sort))
}
## return the index of rejections
KR_inference <- function(x, alpha, bound){
    x_sort <- sort(x)
    fdp_bar<- KR_fdp(x_sort, alpha)
    i <- which(fdp_bar<bound)
    if(length(i)){
        i <- i[length(i)]
        which(x<=x_sort[i])
    }else{
        integer()
    }
}
