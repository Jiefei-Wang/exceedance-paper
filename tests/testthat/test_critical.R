context("critical")
n <- 10L
alpha=0.4
index <- NULL
statNameList <-c("KS","KS+","KS-","BJ","BJ+","BJ-","HC","HC+","HC-")

for(k in seq_along(statNameList)){
    statName <- statNameList[k]
    test_that(statName,{
        critical <- GKSCritical(alpha=alpha,n=n,index=index, statName=statName)
        pvalue <- GKSPvalue(stat = critical,n=n,index=index,statName=statName)
        expect_true(abs(alpha-pvalue)<0.001)
    })
}


index <- 2:3
for(k in seq_along(statNameList)){
    statName <- statNameList[k]
    test_that(statName,{
        critical <- GKSCritical(alpha=alpha,n=n,index = index, statName=statName)
        pvalue <- GKSPvalue(stat = critical,n=n,index=index,statName=statName)
        expect_true(abs(alpha-pvalue)<0.001)
    })
}


