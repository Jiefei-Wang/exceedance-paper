## Please make sure your working dir is the path of this file
source("common_func.R")
library(exceedance)
library(mvtnorm)
## We use BiocParallel + RedisParam to do the parallel computing
library(BiocParallel)
library(RedisParam)

sample_pvalue <- function(n, mu, n_null, sigma) {
    n_alt <- n - n_null
    1-pnorm(rmvnorm(1, c(rep(0,n_null),rep(mu,n_alt)),sigma))
}


nSim <- 100
alpha <- 0.1
gamma <- 0.1
m <- 100
a <- 0.5
mu <- 3
cor_list <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

FDP_result<-c()
FDX_result<-c()
power_result<-c()
for(j in seq_along(cor_list)){
    message("j:",j)
    cor_parm<-cor_list[j]
    sigma <- getCovMat(m,cor_parm)
    
    H1_num <- m * a
    H0_num <- m - H1_num
    
    mydata <- t(sapply(1:nSim, function(x) sample_pvalue(m,mu,H0_num,sigma)))
    
    ## combined stat: all region
    message("combined stat: all region")
    param_combined_all <- get_combined_param(m)
    result_combined_all<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        test_summary(param_combined_all, x, alpha, gamma, H0_num)
    }
    ## BJ stat: all region
    message("BJ stat: all region")
    param_BJ_all <- param_fast_GW(statistic = "BJ", 
                                  param1 = c(0,1), param2 = c(0,1),
                                  range_type = "proportion")
    result_BJ_all<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        test_summary(param_BJ_all, x, alpha, gamma, H0_num)
    }
    ## BJ+ stat: all region
    message("BJ+ stat: all region")
    param_BJ_plus_all <- param_fast_GW(statistic = "BJ", param1 = c(0,1),
                                       range_type = "proportion")
    result_BJ_plus_all<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        test_summary(param_BJ_plus_all, x, alpha, gamma, H0_num)
    }
    
    
    ## combined stat: 1-k
    message("combined stat: 1-k")
    result_combined_auto<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        alpha_hat <- 2*(mean(x<0.5)-0.5)
        k <- max(round(alpha_hat*gamma*m),1)
        param_combined_auto <- get_combined_param(k)
        test_summary(param_combined_auto, x, alpha, gamma, H0_num)
    }
    
    ## BJ+ stat: 1-k
    message("BJ+ stat: 1-k")
    result_BJ_plus_auto<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        alpha_hat <- 2*(mean(x<0.5)-0.5)
        k <- max(round(alpha_hat*gamma*m),1)
        param_BJ_plus_auto <- param_fast_GW(statistic = "BJ", param1 = 1:k,
                                            range_type = "index")
        test_summary(param_BJ_plus_auto, x, alpha, gamma, H0_num)
    }
    
    
    
    collections <- c(
        colMeans(result_BJ_all),
        colMeans(result_BJ_plus_all),
        colMeans(result_combined_all),
        colMeans(result_BJ_plus_auto),
        colMeans(result_combined_auto)
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

FDP_result1<-cbind(alpha,gamma,FDP_result)
FDX_result1<-cbind(alpha,gamma,FDX_result)
power_result1<-cbind(alpha,gamma,power_result)

stat_names <- c("alpha","c","m", "a", "theta","sigma","BJ all","BJ+ all","CB all",
                "BJ+ 1-k","CB 1-k")
colnames(FDP_result1) <- stat_names
colnames(FDX_result1) <- stat_names
colnames(power_result1) <- stat_names

FDP_result1
FDX_result1
power_result1