## test the samples where the first `n_true` samples
## are from true null hypothesis
test_summary<-function(param, x, alpha, gamma, n_true){
    profile <- exceedance_profile(x, param)
    rejected_idx <- exceedance_inference(profile, alpha, gamma)
    n_reject <- max(length(rejected_idx),1)
    
    FDP <- sum(rejected_idx<=H0_num)/n_reject
    FDX <- FDP>gamma
    power <- sum(rejected_idx>H0_num)/H1_num
    c(FDP,FDX,power)
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

#compute the covariance matrix
getCovMat <- function(n, params, type = c("ind", "ec", "ar1", "toeplitz")){
    type <- match.arg(type)
    covMat <- matrix(0, n, n)
    if (type == "ind") {
        for(i in 1:n) {
            covMat[i,j] <- params[[1]]
        }
    }
    if (type == "ec") {
        for (i in 1:n) {
            for (j in 1:n) {
                if (i == j) 
                    covMat[i,j] <- params[[1]]
                else 
                    covMat[i,j] <- params[[1]] * params[[2]]
            }
        }
    }
    if (type == "ar1") {
        for (i in 1:n) {
            for (j in 1:n) {
                covMat[i,j] <- params[[1]] * params[[2]] ^ abs(i-j)
            }
        }
    }
    if (type == "toeplitz") {
        for (i in 1:n) {
            for (j in 1:n) {
                if (i == j) 
                    covMat[i,j] <- params[[1]]
                else if (abs(i-j) == 1)
                    covMat[i,j] <- params[[1]] * params[[2]]
            }
        }
    }
    covMat
}