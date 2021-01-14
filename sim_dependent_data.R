## Please make sure your working dir is the path of this file
source("common_func.R")
library(mvtnorm)
devtools::load_all()
sample_pvalue<-function(n,mu,n_null,sigma){
    n_alt <- n - n_null
    1-pnorm(rmvnorm(1, c(rep(0,n_null),rep(mu,n_alt)),sigma))
}
#compute the covariance matrix
getCovMat<-function(n,cor_parm){
    covMat=matrix(0,n,n)
    for(i in 1:n){
        for(j in 1:n){
            if(i==j){
                covMat[i,j]=1
            }else{
                covMat[i,j]=cor_parm
            }
        }
    }
    covMat
}


nSim <- 100
alpha <- 0.1
gamma <- 0.1
m <- 100
mu <- 3
a <- 0.5
H1_num <- m * a
H0_num <- m - H1_num
cor_list <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

set.seed(1)
FDP_result<-c()
FDX_result<-c()
power_result<-c()
for(j in seq_along(cor_list)){
    cor_parm<-cor_list[j]
    sigma <- getCovMat(m,cor_parm)
    
    mydata <- t(sapply(1:nSim, function(x) sample_pvalue(m,mu,H0_num,sigma)))
    
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
                      c(m, a, mu,cor_parm, FDP_collections))
    FDX_result<-rbind(FDX_result,
                      c(m, a, mu,cor_parm,FDX_collections))
    power_result<-rbind(power_result,
                        c(m, a, mu,cor_parm,power_collections))
}

FDP_result<-cbind(alpha,gamma,FDP_result)
FDX_result<-cbind(alpha,gamma,FDX_result)
power_result<-cbind(alpha,gamma,power_result)

stat_names <- c("alpha","c","m", "a", "theta","cor","CB all","BJ+ all","BJ all",
                "CB 1-10","BJ+ 1-10","BJ 1-10")
colnames(FDP_result) <- stat_names
colnames(FDX_result) <- stat_names
colnames(power_result) <- stat_names

FDP_result
FDX_result
power_result