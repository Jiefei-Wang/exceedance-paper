context("Simulation")
set.seed(123)
nSim <- 1000
n <- 10L
alpha <- 0.4
index <- NULL
statSign <- c(1,1,1,-1,-1,-1,1,1,1)
statNameList <-c("KS","KS+","KS-","BJ","BJ+","BJ-","HC","HC+","HC-")
for(k in seq_along(statNameList)){
    statName <- statNameList[k]
    test_that(statName,{
        critical <- GKSCritical(alpha=alpha,n=n,index = index, statName=statName)
        record <- rep(0,nSim)
        for(i in seq_len(nSim)){
            x <- runif(n)
            stat <- GKSStat(x =x ,index = index, statName = statName,pvalue = FALSE)
            record[i] <- stat$statValue
        }
        stdVar <-  sqrt(alpha*(1-alpha)/nSim)
        if(statSign[k]==1){
            typeI <- mean(record>critical)
        }else{
            typeI <- mean(record<critical)
        }
        expect_true(abs(typeI-alpha)<stdVar * qnorm(0.9999))
    })
}

index <- 2:4
for(k in seq_along(statNameList)){
    statName <- statNameList[k]
    test_that(statName,{
        critical <- GKSCritical(alpha=alpha,n=n,index = index, statName=statName)
        record <- rep(0,nSim)
        for(i in seq_len(nSim)){
            x <- runif(n)
            stat <- GKSStat(x =x ,index = index, statName = statName,pvalue = FALSE)
            record[i] <- stat$statValue
        }
        stdVar <-  sqrt(alpha*(1-alpha)/nSim)
        if(statSign[k]==1){
            typeI <- mean(record>critical)
        }else{
            typeI <- mean(record<critical)
        }
        
        expect_true(abs(typeI-alpha)<stdVar * qnorm(0.9999))
    })
}

