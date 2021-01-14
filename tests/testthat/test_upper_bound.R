context("Test upper bound function")

##############################
## statistic: kth pvalue
## algorithm: general vs general JW vs k order
##############################
sample_size <- 10
k_test <-function(y,k){
    n <- length(y)
    if(length(y)>=k)
        pbeta(y[k],k,n-k+1)
    else 
        1
}
test_that("General vs JW vs k order",{
    # set.seed(123)
    for(i in 1:100){
        k<-5
        alpha <- 0.05
        y<-rbeta(sample_size,1,3)
        
        ## General: General
        params1 <- param_general_GW(function(y)k_test(y,k),algorithm = "general")
        profile1 <- exceedance_profile(y,params1)
        gammabar1 = c()
        for(j in seq_along(y)){
            gammabar1[j] = exceedance_confidence(profile1,alpha,sri = 1:j)
        }
        gammabar1
        
        ## General: JW
        params2 <- param_general_GW(function(y)k_test(y,k),
                                 algorithm = "JW")
        profile2 <- exceedance_profile(y,params2)
        gammabar2 = c()
        for(j in seq_along(y)){
            gammabar2[j] = exceedance_confidence(profile2,alpha,sri = 1:j)
        }
        gammabar2
        
        
        ## GW:k order
        params3 <- param_fast_GW(statistic = "kth_p",param1 = k)
        profile3 <- exceedance_profile(y,params3)
        gammabar3 = c()
        for(j in seq_along(y)){
            gammabar3[j] = exceedance_confidence(profile3,alpha,sri = 1:j)
        }
        gammabar3
        
        
        expect_equal(gammabar1,gammabar2)
        expect_equal(gammabar1,gammabar3)
        # message(i)
    }
})


##############################
## statistic: BJ+
## algorithm: general vs general JW vs fast
##############################
# sample_size <- 10
BJ_test <-function(y,k){
    n <- length(y)
    k_max <- k[length(k)]
    if(length(y)<k_max){
        k <- k[k<=length(y)]
    }
    
    if(length(k)!=0)
        GKSStat(y,index=k,statName = "BJ+")$pvalue
    else 
        1
}
test_that("General vs JW vs k order",{
    # set.seed(123)
    for(i in 1:20){
        k<-5
        alpha <- 0.05
        y<-rbeta(sample_size,1,3)
        
        ## General: General
        params1 <- param_general_GW(function(y)BJ_test(y,k),algorithm = "general")
        profile1 <- exceedance_profile(y,params1)
        gammabar1 = c()
        for(j in seq_along(y)){
            gammabar1[j] = exceedance_confidence(profile1,alpha,sri = 1:j)
        }
        gammabar1
        
        ## General: JW
        params2 <- param_general_GW(function(y)BJ_test(y,k),
                                    algorithm = "JW")
        profile2 <- exceedance_profile(y,params2)
        gammabar2 = c()
        for(j in seq_along(y)){
            gammabar2[j] = exceedance_confidence(profile2,alpha,sri = 1:j)
        }
        gammabar2
        
        
        ## GW:k order
        params3 <- param_fast_GW(statistic = "BJ",param1 = k)
        profile3 <- exceedance_profile(y,params3)
        gammabar3 = c()
        for(j in seq_along(y)){
            gammabar3[j] = exceedance_confidence(profile3,alpha,sri = 1:j)
        }
        gammabar3
        
        
        expect_equal(gammabar1,gammabar2)
        expect_equal(gammabar1,gammabar3)
        # message(i)
    }
})




##############################
## statistic: KS
## algorithm: general vs fast
##############################
k_test2 <-function(y,k){
    n <- length(y)
    k_max <- k[length(k)]
    if(length(y)<k_max){
        k <- k[k<=length(y)]
    }
    
    if(length(k)!=0)
        GKSStat(y,indexL=k,indexU=k,statName = "KS")$pvalue
    else 
        1
}

test_that("General KS vs KS",{
     # set.seed(123)
    for(i in 1:40){
        k<-c(2,4,5)
        alpha <- 0.05
        y<-rbeta(sample_size,1,10)
        
        params1 <- param_general_GW(function(y)k_test2(y,k),
                                 algorithm = "general")
        profile1 <- exceedance_profile(y,params1)
        gammabar1 = c()
        for(j in seq_along(y)){
            gammabar1[j] = exceedance_confidence(profile1,alpha,sri = 1:j)
        }
        gammabar1
        
        params2 <- param_fast_GW(statistic = "KS",param1 = k,param2 = k)
        profile2 <- exceedance_profile(y,params2)
        gammabar2 = c()
        for(j in seq_along(y)){
            gammabar2[j] = exceedance_confidence(profile2,alpha,sri = 1:j)
        }
        gammabar2
        
        
        expect_equal(gammabar1,gammabar2)
        # message(i)
    }
})
##############################
## statistic: KS+
## algorithm: general vs fast
##############################
k_test3 <-function(y,k){
    n <- length(y)
    k_max <- k[length(k)]
    if(length(y)<k_max){
        k <- k[k<=length(y)]
    }
    
    if(length(k)!=0)
        GKSStat(y,indexL=k,statName = "KS+")$pvalue
    else 
        1
}
test_that("General KS+ vs KS+",{
    # set.seed(123)
    for(i in 1:40){
        k<-c(2,4,5)
        alpha <- 0.05
        y<-rbeta(sample_size,1,10)
        
        params1 <- param_general_GW(function(y)k_test3(y,k),
                                 algorithm = "general")
        profile1 <- exceedance_profile(y,params1)
        gammabar1 = c()
        for(j in seq_along(y)){
            gammabar1[j] = exceedance_confidence(profile1,alpha,sri = 1:j)
        }
        gammabar1
        
        params2 <- param_fast_GW(statistic = "KS",param1 = k)
        profile2 <- exceedance_profile(y,params2)
        gammabar2 = c()
        for(j in seq_along(y)){
            gammabar2[j] = exceedance_confidence(profile2,alpha,sri = 1:j)
        }
        gammabar2
        
        expect_equal(gammabar1,gammabar2)
        # message(i)
    }
})

##############################
## statistic: BJ+
## algorithm: general vs fast
##############################
BJ_test <-function(y,k){
    n <- length(y)
    k_max <- k[length(k)]
    if(length(y)<k_max){
        k <- k[k<=length(y)]
    }
    
    if(length(k)!=0)
        GKSStat(y,index=k,statName = "BJ+")$pvalue
    else 
        1
}

test_that("General BJ+ vs BJ+",{
    # set.seed(123)
    for(i in 1:40){
        k<-c(2,4,5)
        alpha <- 0.05
        y<-rbeta(sample_size,1,10)
        
        params1 <- param_general_GW(function(y)BJ_test(y,k),
                                    algorithm = "general")
        profile1 <- exceedance_profile(y,params1)
        gammabar1 = c()
        for(j in seq_along(y)){
            gammabar1[j] = exceedance_confidence(profile1,alpha,sri = 1:j)
        }
        gammabar1
        
        params2 <- param_fast_GW(statistic = "BJ",param1 = k)
        profile2 <- exceedance_profile(y,params2)
        gammabar2 = c()
        for(j in seq_along(y)){
            gammabar2[j] = exceedance_confidence(profile2,alpha,sri = 1:j)
        }
        gammabar2
        
        
        expect_equal(gammabar1,gammabar2)
        # message(i)
    }
})



##############################
## statistic: q quantile
## algorithm: general JW vs q quantile
##############################
q_quantile <-function(y,q){
    n <- length(y)
    k<- max(ceiling(seq_len(n)*q),1L)
    pbeta(y[k],k,n-k+1)
}
test_that("General JW vs GW q quantile",{
    # set.seed(123)
    for(i in 1:100){
        q<-0.1
        alpha <- 0.05
        y<-rbeta(sample_size,1,10)
        
        params1 <- param_general_GW(function(y)q_quantile(y,q),
                                 algorithm = "general")
        profile1 <- exceedance_profile(y,params1)
        gammabar1 = c()
        for(j in seq_along(y)){
            gammabar1[j] = exceedance_confidence(profile1,alpha,sri = 1:j)
        }
        gammabar1
        
        
        params2 <- param_fast_GW(statistic = "kth_p",param1 = q,
                            range_type = "proportion")
        profile2 <- exceedance_profile(y,params2)
        gammabar2 = c()
        for(j in seq_along(y)){
            gammabar2[j] = exceedance_confidence(profile2,alpha,sri = 1:j)
        }
        gammabar2
        
        expect_equal(gammabar1,gammabar2)
    }
})

##############################
## statistic: combined k order
## algorithm: manually do k order vs combined
##############################


test_that("combined k order",{
    # set.seed(123)
    for(i in 1:100){
        q<-0.1
        alpha <- 0.05
        y<-rbeta(sample_size,1,10)
        
        params1 <- param_fast_GW(statistic = "kth_p",param1 = 2,range_type = "index")
        params2 <- param_fast_GW(statistic = "kth_p",param1 = 3,range_type = "index")
        params3 <- param_combine(params1,params2)
        
        
        
        profile1 <- exceedance_profile(y,params1)
        gammabar1 = c()
        for(j in seq_along(y)){
            gammabar1[j] = exceedance_confidence(profile1,alpha,sri = 1:j)
        }
        gammabar1
        
        profile2 <- exceedance_profile(y,params2)
        gammabar2 = c()
        for(j in seq_along(y)){
            gammabar2[j] = exceedance_confidence(profile2,alpha,sri = 1:j)
        }
        gammabar2
        
        profile3 <- exceedance_profile(y,params3)
        gammabar3 = c()
        for(j in seq_along(y)){
            gammabar3[j] = exceedance_confidence(profile3,alpha*2,sri = 1:j)
        }
        gammabar3
        
        expect_equal(gammabar3,pmin(gammabar1,gammabar2))
    }
})







