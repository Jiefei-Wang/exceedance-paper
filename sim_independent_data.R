## Please make sure your working dir is the path of this file
source("common_func.R")
library(exceedance)
## We use foreach to do the parallel computing
library(foreach)
library(doParallel)
registerDoParallel(cores=6)

H1_p <- function(n, mu){
    pnorm(rnorm(n,mu,1), lower.tail = FALSE)
}

nSim <- 100
alpha <- 0.1
gamma <- 0.1
a_mu <- expand.grid(a=c(0.2,0.5),mu=3)
m_list <- c(rep(100,nrow(a_mu)),rep(500,nrow(a_mu)),
            rep(1000,nrow(a_mu)),rep(2000,nrow(a_mu)),
            rep(5000,nrow(a_mu)))
a_list <- rep(a_mu$a,length(m_list)/length(a_mu$a))
mu_list <- rep(a_mu$mu,length(m_list)/length(a_mu$mu))

FDP_result<-c()
FDX_result<-c()
power_result<-c()
for(j in seq_along(a_list)){
    a <- a_list[j]
    mu <- mu_list[j]
    m <- m_list[j]
    message(j,":",a,",",mu,",",m)
    
    H1_num <- m * a
    H0_num <- m - H1_num
    
    mydata <- t(sapply(1:nSim, function(x) c(runif(H0_num),H1_p(H1_num,mu))))
    
    message("BJ stat: all region")
    param_BJ_all <- param_fast_GW(statistic = "BJ", 
                                  param1 = c(0,1), param2 = c(0,1),
                                  range_type = "proportion")
    # result_BJ_all<-matrix(NA,nSim,3)
    result_BJ_all<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        test_summary(param_BJ_all, x, alpha, gamma, H0_num)
    }
    
    message("BJ+ stat: all region")
    param_BJ_plus_all <- param_fast_GW(statistic = "BJ", param1 = c(0,1),
                                       range_type = "proportion")
    # result_BJ_plus_all<-matrix(NA, nSim,3)
    result_BJ_plus_all<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        test_summary(param_BJ_plus_all, x, alpha, gamma, H0_num)
    }
    
    message("combined stat: all region")
    param_combined_all <- get_combined_param(m)
    # result_combined_all<-matrix(NA, nSim,3)
    result_combined_all<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        test_summary(param_combined_all, x, alpha, gamma, H0_num)
    }
    
    message("BJ+ stat: auto scale")
    result_BJ_plus_auto<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        alpha_hat <- 2*(mean(x<0.5)-0.5)
        k <- max(round(alpha_hat*gamma*m),1)
        param_BJ_plus_auto <- param_fast_GW(statistic = "BJ", param1 = 1:k,
                                            range_type = "index")
        test_summary(param_BJ_plus_auto, x, alpha, gamma, H0_num)
    }
    
    message("combined stat: auto scale")
    # result_combined_auto<-matrix(NA, nSim,3)
    result_combined_auto<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        alpha_hat <- 2*(mean(x<0.5)-0.5)
        k <- max(round(alpha_hat*gamma*m),1)
        param_combined_auto <- get_combined_param(k)
        test_summary(param_combined_auto, x, alpha, gamma, H0_num)
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
                      c(m, a, mu, FDP_collections))
    FDX_result<-rbind(FDX_result,
                      c(m, a, mu,FDX_collections))
    power_result<-rbind(power_result,
                        c(m, a, mu,power_collections))
}

FDP_result1<-cbind(alpha,gamma,FDP_result)
FDX_result1<-cbind(alpha,gamma,FDX_result)
power_result1<-cbind(alpha,gamma,power_result)

stat_names <- c("alpha","c","m", "a", "theta","BJ all","BJ+ all","CB all",
                "BJ+ 1-k","CB 1-k")
colnames(FDP_result1) <- stat_names
colnames(FDX_result1) <- stat_names
colnames(power_result1) <- stat_names

FDP_result1
FDX_result1
power_result1
