## Please make sure your working dir is the path of this file
source("common_func.R")
library(exceedance)
library(mvtnorm)
## We use foreach to do the parallel computing
library(foreach)
library(doParallel)
library(parallel)

cl <- makeCluster(12)
registerDoParallel(cl)


nSim <- 1000
alpha <- 0.1
gamma <- 0.1
m <- 100
corList <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
corTypeList <- c("ind", "ar1", "cs")

FDP_result<-c()
FDX_result<-c()
power_result<-c()
for(k in seq_along(corTypeList)){
    for(j in seq_along(corList)){
        message("j:",j)
        cor_parm <- corList[j]
        
        corType <- corTypeList[k]
        if(corType == "ind" && cor_parm!=0)
            next
        
        sigma <- getCovMat(m,cor_parm, corType)
        
        mydata <- t(sapply(1:nSim, function(x) sample_pvalue(m,0,m,sigma)))
        
        ## combined stat: 0 ~ k hat
        message("combined stat: 0 ~ k hat ")
        result_combined_auto<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
            x <- mydata[i,]
            alpha_hat <- 2*(mean(x<0.5)-0.5)
            k <- max(round(alpha_hat*gamma*m),1)
            param_combined_auto <- get_combined_param(k)
            test_summary(param_combined_auto, x, alpha, gamma, m)
        }
        
        ## BJ+ stat: 0~0.5 region
        message("BJ+ stat: 0~0.5 region")
        param_BJ_plus_all <- param_fast_GW(statistic = "BJ", param1 = c(0,0.5),
                                           range_type = "proportion")
        # result_BJ_plus_all<-matrix(NA, nSim,3)
        result_BJ_plus_all<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
            x <- mydata[i,]
            test_summary(param_BJ_plus_all, x, alpha, gamma, m)
        }
        
        collections <- c(
            colMeans(result_combined_auto),
            colMeans(result_BJ_plus_all)
        )
        FDP_collections <- collections[(0:((length(collections)-1)/3))*3+1]
        FDX_collections <- collections[(0:((length(collections)-1)/3))*3+2]
        power_collections <- collections[(0:((length(collections)-1)/3))*3+3]
        
        FDP_result<-rbind(FDP_result,
                          c(corType,cor_parm, FDP_collections))
        FDX_result<-rbind(FDX_result,
                          c(corType,cor_parm,FDX_collections))
        power_result<-rbind(power_result,
                            c(corType,cor_parm,power_collections))
    }
}

stat_names <- c("cor_type","rho", "combined", "BJ+")
colnames(FDP_result) <- stat_names
colnames(FDX_result) <- stat_names
colnames(power_result) <- stat_names


# save(FDX_result, file = "sim1")

cor_type <- FDX_result[,1]
rho <- as.numeric(FDX_result[,2])
FDX_bj <- as.numeric(FDX_result[,4])

plot(rho[cor_type == "ar1"], FDX_bj[cor_type == "ar1"], type = 'b',
     ylim = c(0,0.2),xlim = c(0, 1),  xlab = expression(rho), ylab = "Exceedance Rate")
lines(rho[cor_type == "cs"], FDX_bj[cor_type == "cs"], type = 'b', lty=2)
abline(0.1,0, lty=2)

FDX_result[,3]

