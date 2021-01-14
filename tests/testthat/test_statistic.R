context("Test statistics")
x <- 1:9/10
index <- c(2,3)

test_that("empty data",{
    stat <- GKSStat(c(),statName = "KS")
    expect_true(is.na(stat$pvalue))
    expect_equal(stat$n,0)
    expect_equal(stat$statValue,0)
    
    stat <- GKSStat(c(),statName = "BJ")
    expect_true(is.na(stat$pvalue))
    expect_equal(stat$n,0)
    expect_equal(stat$statValue,1)
    
    stat <- GKSStat(c(),statName = "HC")
    expect_true(is.na(stat$pvalue))
    expect_equal(stat$n,0)
    expect_equal(stat$statValue,0)
})

test_that("KS",{
    stat <- GKSStat(x=x,index=index,statName = "KS",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.089)
    expect_equal(round(stat$pvalue,3),0.946)
    
    stat <- GKSStat(x=x,statName = "KS",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.1)
    expect_equal(round(stat$pvalue,3),1)
    
})
test_that("KS+",{
    stat <- GKSStat(x=x,index=index,statName = "KS+",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.033)
    expect_equal(round(stat$pvalue,3),0.643)
    
    stat <- GKSStat(x=x,statName = "KS+",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.1)
    expect_equal(round(stat$pvalue,3),0.786)
})
test_that("KS-",{
    stat <- GKSStat(x=x,index=index,statName = "KS-",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.089)
    expect_equal(round(stat$pvalue,3),0.542)
    
    stat <- GKSStat(x=x,statName = "KS-",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.1)
    expect_equal(round(stat$pvalue,3),0.786)
})


test_that("BJ",{
    stat <- GKSStat(x=x,index=index,statName = "BJ",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.436)
    expect_equal(round(stat$pvalue,3),0.978)
    
    stat <- GKSStat(x=x,statName = "BJ",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.387)
    expect_equal(round(stat$pvalue,3),1)
})

test_that("BJ+",{
    stat <- GKSStat(x=x,index=index,statName = "BJ+",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.537)
    expect_equal(round(stat$pvalue,3),0.646)
    
    stat <- GKSStat(x=x,statName = "BJ+",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.387)
    expect_equal(round(stat$pvalue,3),0.826)
})

test_that("BJ-",{
    stat <- GKSStat(x=x,index=index,statName = "BJ-",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.436)
    expect_equal(round(stat$pvalue,3),0.543)
    
    stat <- GKSStat(x=x,statName = "BJ-",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.387)
    expect_equal(round(stat$pvalue,3),0.826)
})

test_that("HC",{
    stat <- GKSStat(x=x,index=index,statName = "HC",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.667)
    expect_equal(round(stat$pvalue,3),0.941)
    
    stat <- GKSStat(x=x,statName = "HC",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),1)
    expect_equal(round(stat$pvalue,3),0.996)
})
test_that("HC+",{
    stat <- GKSStat(x=x,index=index,statName = "HC+",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.218)
    expect_equal(round(stat$pvalue,3),0.648)
    
    stat <- GKSStat(x=x,statName = "HC+",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),1)
    expect_equal(round(stat$pvalue,3),0.746)
})

test_that("HC-",{
    stat <- GKSStat(x=x,index=index,statName = "HC-",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),0.667)
    expect_equal(round(stat$pvalue,3),0.527)
    
    stat <- GKSStat(x=x,statName = "HC-",pvalue=TRUE)
    expect_equal(round(stat$statValue,3),1)
    expect_equal(round(stat$pvalue,3),0.746)
})




statNameList <-c("KS","BJ","HC")
x <- runif(10L)
for(k in seq_along(statNameList)){
    statName <- statNameList[k]
    test_that(paste0("oneside ",statName),{
        ## + side
        stat1 <- GKSStat(x=x,indexL=index,statName = statName,pvalue=TRUE)
        stat2 <- GKSStat(x=x,indexL=index,statName = paste0(statName,"+"),pvalue=TRUE)
        stat3 <- GKSStat(x=x,index=index,statName = paste0(statName,"+"),pvalue=TRUE)
        expect_equal(stat1,stat2)
        expect_equal(stat1,stat3)
        
        ## - side
        stat4 <- GKSStat(x=x,indexU=index,statName = statName,pvalue=TRUE)
        stat5 <- GKSStat(x=x,indexU=index,statName = paste0(statName,"-"),pvalue=TRUE)
        stat6 <- GKSStat(x=x,index=index,statName = paste0(statName,"-"),pvalue=TRUE)
        expect_equal(stat4,stat5)
        expect_equal(stat4,stat6)
        
        ## error on purpose
        stat7 <- GKSStat(x=x,indexL=index,statName = paste0(statName,"-"), pvalue=TRUE)
        stat8 <- GKSStat(x=x,statName = paste0(statName,"-"), pvalue=TRUE)
        expect_false(stat7$statValue==stat8$statValue)
    })
}

