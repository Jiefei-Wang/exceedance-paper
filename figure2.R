## Please make sure your working dir is the path of this file
source("common_func.R")
library(exceedance)
library(mvtnorm)
## We use foreach to do the parallel computing
library(foreach)
library(doParallel)
library(parallel)

registerDoParallel(cores=12)


nSim <- 1000
alpha <- 0.1
FDP_bound <- 0.1
hypo_num <- 1000
theta <- 1.5
corType <- "cs"
cor <- 0.2
aList <- (0:10)/10


FDP_result<-c()
FDX_result<-c()
power_result<-c()
start <- Sys.time()
for(k in seq_along(aList)){
    message("k:",k)
    H1Prop <- aList[k]
    H1_num <- H1Prop * hypo_num
    H0_num <- hypo_num - H1_num
    
    ## sample data
    sigma <- getCovMat(hypo_num, cor, corType)
    mydata <- sample_pvalue(nSim, hypo_num, theta, H0_num, sigma)
    
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
                      c(H1Prop,FDP_collections))
    FDX_result<-rbind(FDX_result,
                      c(H1Prop,FDX_collections))
    power_result<-rbind(power_result,
                        c(H1Prop,power_collections))
}
end <- Sys.time()
end-start

statName <- c("BJ+", "CB", "KR")
tblName <- c("H1Prop", statName)
colnames(FDP_result) <- tblName
colnames(FDX_result) <- tblName
colnames(power_result) <- tblName

save(FDP_result, FDX_result, power_result, file = paste0("figure2_m",hypo_num,"_theta",theta, "_cor", cor))

load(file = paste0("figure2_m",hypo_num,"_theta",theta, "_cor", cor))

H1 <- as.numeric(power_result[,"H1Prop"])
BJ_power <- as.numeric(power_result[,"BJ+"])
CB_power <- as.numeric(power_result[,"CB"])
KR_power <- as.numeric(power_result[,"KR"])

png(file=paste0("figure2_m",hypo_num,"_theta", "_cor", cor, ".png"),width=600, height=500)
plot(H1, BJ_power, type = 'b',
     ylim = c(0,1),xlim = c(0, 1), xlab = "prop of H1", ylab = "Power", lty=2)
lines(H1, CB_power, type = 'b', lty=2,pch  = 2)
lines(H1, KR_power, type = 'b', lty=2, pch  = 3)
legend("topleft", 
       legend=c("BJ+ all", expression(paste("CB(1-",hat(k),")")), "KR"), 
       lty=c(2,2,2), pch=c(1,2,3))
dev.off()

