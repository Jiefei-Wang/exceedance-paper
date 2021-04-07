## Test for performance
library(exceedance)

gamma <- 0.1
alpha <- 0.1
nSim <- 100
m_list <- c(50,100,200,500,1000,2000,5000,10000)
final <- c()
for(j in seq_along(m_list)){
    m <- m_list[j]
    message(j,":",m)
    x <- runif(m)
    
    ## combined stat
    param_list <- list()
    for(i in 1:m){
        param_list[[i]] <-
            param_fast_GW(statistic = "kth_p", param1 = i,
                          range_type = "index")
    }
    param_combined <- param_combine(param_list = param_list)
    profile_combined <- exceedance_profile(x, param_combined)

    combined_time <- system.time(
        for(i in 1:nSim){
            message(i)
            exceedance_inference(profile_combined,alpha,gamma)
        }
    )
    message("combined")


    ## BJ stat
    param_BJ <- param_fast_GW(statistic = "BJ", param1 = c(0,1),
                              range_type = "proportion")
    profile_BJ <- exceedance_profile(x, param_BJ)

    BJ_time <- system.time(
        for(i in 1:nSim){
            message(i)
            exceedance_inference(profile_BJ,alpha,gamma)
        }
    )
    message("BJ")

    
    collections <- c(m, 
                     combined_time[3]/nSim,
                     BJ_time[3]/nSim)
    final <- rbind(final, collections)
}

colnames(final) <- c("m", "FDP_combined", "power_combined")
final

