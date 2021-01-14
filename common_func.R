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