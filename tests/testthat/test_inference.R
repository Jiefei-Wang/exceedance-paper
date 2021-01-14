context("Test inference function")

sample_size <- 10

##############################
## statistic: k order
## algorithm: general JW vs fast
##############################
k_test <-function(y,k){
    n <- length(y)
    if(length(y)>=k)
        pbeta(y[k],k,n-k+1)
    else 
        1
}
test_that("combined k order",{
    # set.seed(123)
    for(i in 1:100){
        alpha <- 0.05
        bound <- 0.4
        k <- 2
        y<-rbeta(sample_size,1,10)
        sy <- sort(y,index.return = TRUE)
        
        param1 <- param_general_GW(function(y)k_test(y,k),algorithm = "JW")
        param2 <- param_fast_GW(statistic = "kth_p",param1 = k,range_type = "index")
        
        profile1 <- exceedance_profile(y,param1)
        profile2 <- exceedance_profile(y,param2)
        
        result1 <- exceedance_inference(profile1,alpha=alpha,bound=bound)
        result2 <- exceedance_inference(profile2,alpha=alpha,bound=bound)
        
        gammabar2 = c()
        for(j in seq_along(y)){
            gammabar2[j] = exceedance_confidence(profile2,alpha,sri = 1:j)
        }
        ind <- which(gammabar2<=bound)
        if(length(ind)!=0){
            result3 <- sy$ix[seq_len(max(ind))]
        }else{
            result3 <- integer(0)
        }
        
        expect_equal(result1,result2)
        expect_equal(result1,result3)
    }
})


##############################
## statistic: BJ
## algorithm: general JW vs fast
##############################
BJ_test <-function(y){
    n <- length(y)
    GKSStat(y,statName = "BJ")$pvalue
}
test_that("BJ: general JW vs fast",{
    # set.seed(123)
    for(i in 1:100){
        alpha <- 0.05
        bound <- 0.4
        y<-rbeta(sample_size,1,10)
        sy <- sort(y,index.return = TRUE)
        
        param1 <- param_general_GW(BJ_test,algorithm = "JW")
        param2 <- param_fast_GW(statistic = "BJ")
        
        profile1 <- exceedance_profile(y,param1)
        profile2 <- exceedance_profile(y,param2)
        
        result1 <- exceedance_inference(profile1,alpha=alpha,bound=bound)
        result2 <- exceedance_inference(profile2,alpha=alpha,bound=bound)
        
        gammabar2 = c()
        for(j in seq_along(y)){
            gammabar2[j] = exceedance_confidence(profile2,alpha,sri = 1:j)
        }
        ind <- which(gammabar2<=bound)
        if(length(ind)!=0){
            result3 <- sy$ix[seq_len(max(ind))]
        }else{
            result3 <- integer(0)
        }
        
        expect_equal(result1,result2)
        expect_equal(result1,result3)
    }
})

##############################
## statistic: BJ
## algorithm: fast vs fast
## Inference: general vs fast
##############################
test_that("BJ Inference: general vs fastr",{
    # set.seed(123)
    n <- 50
    for(i in 1:10){
        alpha <- 0.05
        bound <- 0.4
        y<-rbeta(n,1,10)
        sy <- sort(y,index.return = TRUE)
        
        param1 <- param_fast_GW(statistic = "BJ", param1 = seq_len(n))
        profile1 <- exceedance_profile(y,param1)
        result1 <- exceedance_inference(profile1,alpha=alpha,bound=bound)
        
        gammabar2 = c()
        for(j in seq_along(y)){
            gammabar2[j] = exceedance_confidence(profile1,alpha,sri = 1:j)
        }
        ind <- which(gammabar2<=bound)
        if(length(ind)!=0){
            result2 <- sy$ix[seq_len(max(ind))]
        }else{
            result2 <- integer(0)
        }
        expect_equal(result1,result2)
    }
})




##############################
## statistic: BJ partial
## algorithm: general JW vs fast
##############################
BJ_test <-function(y,k){
    n <- length(y)
    k_max <- k[length(k)]
    if(length(y)<k_max){
        k <- k[k<=length(y)]
    }
    if(length(k)!=0)
        GKSStat(y,indexL=k,statName = "BJ")$pvalue
    else 
        1
}
test_that("combined k order",{
    # set.seed(123)
    for(i in 1:100){
        alpha <- 0.05
        bound <- 0.4
        k <- c(2,4,6)
        y<-rbeta(sample_size,1,10)
        sy <- sort(y,index.return = TRUE)
        
        param1 <- param_general_GW(function(y)BJ_test(y,k),algorithm = "JW")
        param2 <- param_fast_GW(statistic = "BJ",param1 = k,range_type = "index")
        
        profile1 <- exceedance_profile(y,param1)
        profile2 <- exceedance_profile(y,param2)
        
        result1 <- exceedance_inference(profile1,alpha=alpha,bound=bound)
        result2 <- exceedance_inference(profile2,alpha=alpha,bound=bound)
        
        gammabar2 = c()
        for(j in seq_along(y)){
            gammabar2[j] = exceedance_confidence(profile2,alpha,sri = 1:j)
        }
        ind <- which(gammabar2<=bound)
        if(length(ind)!=0){
            result3 <- sy$ix[seq_len(max(ind))]
        }else{
            result3 <- integer(0)
        }
        
        expect_equal(result1,result2)
        expect_equal(result1,result3)
    }
})



##############################
## statistic: combined k order
## algorithm: manually do k order vs combined
##############################
test_that("combined k order",{
    # set.seed(123)
    for(i in 1:100){
        alpha <- 0.05
        bound <- 0.4
        y<-rbeta(sample_size,1,10)
        
        param1 <- param_fast_GW(statistic = "kth_p",param1 = 2,range_type = "index")
        param2 <- param_fast_GW(statistic = "kth_p",param1 = 3,range_type = "index")
        params3 <- param_combine(param1,param2)
        
        profile1 <- exceedance_profile(y,param1)
        profile2 <- exceedance_profile(y,param2)
        profile3 <- exceedance_profile(y,params3)
        
        result1 <- exceedance_inference(profile1,alpha=alpha,bound=bound)
        result2 <- exceedance_inference(profile2,alpha=alpha,bound=bound)
        result3 <- exceedance_inference(profile3,alpha=alpha*2,bound=bound)
        
        expect_equal(sort(result3),sort(unique(c(result1,result2))))
    }
})
