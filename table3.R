## Please make sure your working dir is the path of this file
source("common_func.R")
library(exceedance)
library(mvtnorm)
## We use foreach to do the parallel computing
library(foreach)
library(doParallel)
registerDoParallel(cores=12)


nSim <- 100
alpha <- 0.1
FDP_bound <- 0.1
hypo_num <- 1000
H0_num <- 500
H1_num <- 500
theta <- 1.5
corList <- 0:10/10
corType <- "ar"

FDP_result<-c()
FDX_result<-c()
power_result<-c()
for(k in seq_along(corList)){
    message("k:",k)
    cor_parm <- corList[k]
    
    ## sample data
    sigma <- getCovMat(hypo_num, cor_parm, corType)
    mydata <- sample_pvalue(nSim, hypo_num, theta, H0_num, sigma)
    
    ## BJ+ stat: 0~1 region
    message("BJ+ all")
    param_BJ_plus_all <- param_fast_GW(statistic = "BJ", param1 = c(0,1),
                                       range_type = "proportion")
    result_BJ_plus_all<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        profile <- exceedance_profile(x, param_BJ_plus_all)
        rejected_idx <- exceedance_inference(profile, alpha, FDP_bound)
        test_summary(rejected_idx, hypo_num, H0_num, FDP_bound)
    }
    
    message("CB 1~k")
    result_combined_auto<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        alpha_hat <- 2*(mean(x<0.5)-0.5)
        k <- max(round(alpha_hat*FDP_bound*hypo_num),1)
        param_combined_auto <- get_combined_param(k)
        
        profile <- exceedance_profile(x, param_combined_auto)
        rejected_idx <- exceedance_inference(profile, alpha, FDP_bound)
        test_summary(rejected_idx, hypo_num, H0_num, FDP_bound)
    }
    
    message("KR")
    result_KR<-foreach(i = 1:nSim,.combine=rbind,.packages="exceedance")%dopar%{
        x <- mydata[i,]
        rejected_idx <- KR_inference(x, alpha, FDP_bound)
        test_summary(rejected_idx, hypo_num, H0_num, FDP_bound)
    }
    
    collections <- rbind(
        colMeans(result_BJ_plus_all),
        colMeans(result_combined_auto),
        colMeans(result_KR)
    )
    
    FDP_collections <- collections[,1]
    FDX_collections <- collections[,2]
    power_collections <- collections[,3]
    
    FDP_result<-rbind(FDP_result,
                      c(theta, corType, cor_parm, H1_num, FDP_collections))
    FDX_result<-rbind(FDX_result,
                      c(theta, corType, cor_parm, H1_num, FDX_collections))
    power_result<-rbind(power_result,
                        c(theta, corType, cor_parm, H1_num, power_collections))
}

statName <- c("BJ+", "CB", "KR")
tblName <- c("theta","cor_type","rho", "H1_num",  statName)
colnames(FDP_result) <- tblName
colnames(FDX_result) <- tblName
colnames(power_result) <- tblName

save(FDP_result, FDX_result, power_result, file = paste0("table3_a0.5_theta",theta))


load(file = paste0("table3_a0.5_theta",theta))
results <- cbind(FDX_result[,c("rho", statName)], power_result[,statName])
results <-round(matrix(as.numeric(results), nrow=nrow(results)),3)
colnames(results) <- 
    c(
        "$\\rho$", 
        "BJ+ all", "CB $(1-\\widehat{k})$", "KR", 
        "BJ+ all", "CB $(1-\\widehat{k})$", "KR"
    )



write.csv(results, file = paste0("table3_a0.5_theta",theta, ".csv"))
