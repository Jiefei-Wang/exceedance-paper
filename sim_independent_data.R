## Please make sure your working dir is the path of this file
source("common_func.R")
devtools::load_all()
H1_p <- function(n, mu){
    pnorm(rnorm(n,mu,1), lower.tail = FALSE)
}

nSim <- 100
alpha <- 0.1
gamma <- 0.1
a_mu <- expand.grid(a=c(0,0.2,0.5),mu=3)
m_list <- c(rep(100,nrow(a_mu)),rep(500,nrow(a_mu)),
            rep(1000,nrow(a_mu)),rep(2000,nrow(a_mu)),
            rep(5000,nrow(a_mu)))
a_list <- rep(a_mu$a,length(m_list)/length(a_mu$a))
mu_list <- rep(a_mu$mu,length(m_list)/length(a_mu$mu))

set.seed(1)
FDP_result<-c()
FDX_result<-c()
power_result<-c()
for(j in seq_along(a_list)){
    a <- a_list[j]
    mu <- mu_list[j]
    m <- m_list[j]
    
    H1_num <- m * a
    H0_num <- m - H1_num
    
    mydata <- t(sapply(1:nSim, function(x) c(runif(H0_num),H1_p(H1_num,mu))))
    
    ## combined stat: all region
    param_combined_all <- get_combined_param(m)
    result_combined_all<-matrix(NA,nSim,3)
    for(i in 1:nSim){
        x <- mydata[i,]
        result_combined_all[i,]<-
            test_summary(param_combined_all, x, alpha, gamma, H0_num)
    }
    
    ## BJ+ stat: all region
    param_BJ_plus_all <- param_fast_GW(statistic = "BJ", param1 = c(0,1),
                                       range_type = "proportion")
    result_BJ_plus_all<-matrix(NA,nSim,3)
    for(i in 1:nSim){
        x <- mydata[i,]
        result_BJ_plus_all[i,]<-
            test_summary(param_BJ_plus_all, x, alpha, gamma, H0_num)
    }
    
    ## BJ stat: all region
    param_BJ_all <- param_fast_GW(statistic = "BJ", 
                                  param1 = c(0,1), param2 = c(0,1),
                                  range_type = "proportion")
    result_BJ_all<-matrix(NA,nSim,3)
    for(i in 1:nSim){
        x <- mydata[i,]
        result_BJ_all[i,]<-
            test_summary(param_BJ_all, x, alpha, gamma, H0_num)
    }
    
    ## combined stat: 1-10
    param_combined_1_10 <- get_combined_param(10)
    result_combined_1_10<-matrix(NA,nSim,3)
    for(i in 1:nSim){
        x <- mydata[i,]
        result_combined_1_10[i,]<-
            test_summary(param_combined_1_10, x, alpha, gamma, H0_num)
    }
    
    ## BJ+ stat: 1-10
    param_BJ_plus_1_10 <- param_fast_GW(statistic = "BJ", param1 = 1:10,
                                        range_type = "index")
    result_BJ_plus_1_10<-matrix(NA,nSim,3)
    for(i in 1:nSim){
        x <- mydata[i,]
        result_BJ_plus_1_10[i,]<-
            test_summary(param_BJ_plus_1_10, x, alpha, gamma, H0_num)
    }
    
    ## BJ stat: 1-10
    param_BJ_1_10 <- param_fast_GW(statistic = "BJ", 
                                   param1 = 1:10, param2 = 1:10,
                                   range_type = "index")
    result_BJ_1_10<-matrix(NA,nSim,3)
    for(i in 1:nSim){
        x <- mydata[i,]
        result_BJ_1_10[i,]<-
            test_summary(param_BJ_1_10, x, alpha, gamma, H0_num)
    }
    
    
    collections <- c(
        colMeans(result_combined_all),
        colMeans(result_BJ_plus_all),
        colMeans(result_BJ_all),
        colMeans(result_combined_1_10),
        colMeans(result_BJ_plus_1_10),
        colMeans(result_BJ_1_10)
    )
    FDP_collections <- collections[(0:((length(collections)-1)/3))*3+1]
    FDX_collections <- collections[(0:((length(collections)-1)/3))*3+2]
    power_collections <- collections[(0:((length(collections)-1)/3))*3+3]
    
    FDP_result<-rbind(FDP_result,
                      c(m, a, mu, FDP_collections))
    FDX_result<-rbind(FDX_result,
                      c(m, a, mu,FDX_collections))
    power_result<-rbind(power_result,
                        c(m, a, mu,power_collections))
}

FDP_result<-cbind(alpha,gamma,FDP_result)
FDX_result<-cbind(alpha,gamma,FDX_result)
power_result<-cbind(alpha,gamma,power_result)

stat_names <- c("alpha","c","m", "a", "theta","CB all","BJ+ all","BJ all",
                "CB 1-10","BJ+ 1-10","BJ 1-10")
colnames(FDP_result) <- stat_names
colnames(FDX_result) <- stat_names
colnames(power_result) <- stat_names

FDP_result
FDX_result
power_result