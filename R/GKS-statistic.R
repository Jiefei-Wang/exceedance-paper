#' Compute generalized Kolmogorov-Smirnov test statistics
#' 
#' Compute the Kolmogorov-Smirnov, Berk-Jones or the higher criticism statistics
#' to test whether the data is from an uniform(0,1) distribution. 
#' The function `GKSStat` provides an uniform way to computes different 
#' test statistics.
#' To be consistent with the other statistics, the traditional higher criticism 
#' statistic is named `HC+` and the statistic `HCStat` computes the 
#' two-sided higher criticism statistic.
#' 
#' @param x Numeric, the samples that the test statistics will be based on.
#' @param alpha0 Numeric, controlling which ordered samples will be used in the
#' statistics, the default value is `1`. see details.
#' @param index Integer, controlling which ordered samples will be used in the
#' statistics, see details.
#' @param indexL Integer, controlling which ordered samples will be used in the
#' statistics, see details.
#' @param indexU Integer, controlling which ordered samples will be used in the
#' statistics, see details.
#' @param statName Character, the name of the statistic that will be computed.
#' The default is "KS".
#' @param pvalue Logical, whether to compute the p-value of the statistic. 
#' The default is `TRUE`
#' 
#' @details 
#' \bold{statistics definitions}
#' 
#' The function compute the test statistics which aggregate the significant signal
#' from the order statistics of the samples, that is, if `T` is a statistic 
#' and `X_1`,`X_2`,...,`X_n` are the samples, the value of `T` is purely 
#' based on the value of `X_(1)`,`X_(2)`,...,`X_(n)`, 
#' where `X_(i)` is the ith ascending sorted samples of `X1`,`X2`,...,`Xn`. 
#' Moreover, the rejection region of the statistic `T` can be written as 
#' a set of rejection regions of the ordered samples `X_(1)`,`X_(2)`,...,`X_(n)`.
#'  In other words, there exist two sequences `{l_i}` and `{u_i}` for `i=1,...,n` 
#' and the statistic `T` is rejected if and only if there exist
#' one `i` such that `X_(i) < l_i` or `X_(i) > u_i`.
#' 
#' The most well-known statistic which takes this form is the Kolmogorov-Smirnov 
#' statistic. Other statistics like Berk-Jones or the higher criticism also have
#' similar formulas but define different sets of `{l_i}` and `{u_i}`. 
#' 
#' \bold{alpha0, index, indexL and indexU}
#' 
#' As mentioned previouly, the rejection of a test can be determined by the 
#' sequences of `{l_i}` and `{u_i}`. Therefore, the parameter `alpha0`, `index`
#' `indexL` and `indexU`. provide a way to control which `l_i` and `u_i` 
#' will be considered in the test procedure. If no argument is provided, all `l_i`s
#' and `u_i`s will be compared with their corresponding sorted sample `X_(i)`. 
#' This yields the traditional test statistics. If `alpha0` is used, only
#' the data `X_(1),...X_(k)` will be used in the test where `k` is the nearest 
#' integer of `alpha0*n`. If `index` is provided, only `X_(i)` for `i` in `index` 
#' will be considered in the test. If `indexL` and/or `indexU` is not `NULL`, 
#' only `l_i` for `i` in `indexL` and `u_i` for `i` in `indexU` will be used as the 
#' rejection boundary for the test. These can be used to generate an one-sided version 
#' of the test statistic. For example, if `indexL` is from `1` to the length of `x` and
#' `indexU` is `NULL`, this will yield a test specifically sensitive to smaller samples.
#' The test statistics like `KS+`, `HC+` and `BJ+` are implemented by calling
#' `GKSStat(..., indexU = NULL)`, where `indexU` is always `NULL`.
#' 
#' @examples 
#' ## Generate samples
#' x <- rbeta(10, 1, 2)
#' 
#' ## Perform KS test
#' GKSStat(x = x, statName = "KS")
#' 
#' ## Perform one-sided KS test
#' GKSStat(x = x, statName = "KS+")
#' GKSStat(x = x, statName = "KS-")
#' 
#' @return a `GKSStat` S3 object
#' @rdname statistics
#' @export
GKSStat <- function(
    x, index = NULL, indexL = NULL, indexU = NULL,
    statName = c("KS","KS+","KS-","BJ","BJ+","BJ-","HC","HC+","HC-","Simes"),
    pvalue = TRUE){
    statName <- match.arg(statName)
    idx <- getGKSIndex(statName = statName, n = length(x), 
                       index = index, indexL = indexL, indexU = indexU)
    indexL <- idx$indexL
    indexU <- idx$indexU
    statName <- idx$statName
    
    statValue <- call_func(root = "Stat",
                           prefix = statName,
                           x=x, 
                           indexL=indexL,
                           indexU=indexU)
    
    stat <- .GKSStat(statName = statName,
                           statValue= statValue,
                           n=length(x),
                           indexL=indexL,
                           indexU=indexU,
                           data = x)
    if(pvalue)
        stat$pvalue <- GKSPvalue(stat=stat)
    stat
}

genericStatFunc <- function(statFunc, x, indexL, indexU) {
    n <- length(x)
    stopifnot(length(indexL)<=n)
    stopifnot(length(indexU)<=n)
    # stopifnot(all(indexU>=1&indexU<=n))
    # stopifnot(all(indexL>=1&indexL<=n))
    sx <- sort(x)
    sx[sx == 0] <- min(10 ^ -6, sx[sx != 0])
    sx[sx == 1] <- max(1 - 10 ^ -6, sx[sx != 1])
    statFunc(
        n = n,
        x = x,
        sx = sx,
        indexU = indexU ,
        indexL = indexL
    )
}

HCStatFunc <- function(n, x, sx, indexU, indexL) {
    HCPlus <- max(HCPlusLevel(n, x, sx,indexL),0)
    HCMinus <- max(HCMinusLevel(n, x, sx,indexU),0)
    max(c(HCPlus, HCMinus))
}
BJStatFunc <- function(n, x, sx, indexL, indexU) {
    BJPlus <- min(BJPlusLevel(n, x, sx,indexL),1)
    BJMinus <- min(BJMinusLevel(n, x, sx,indexU),1)
    min(BJPlus, BJMinus)
}
KSStatFunc <- function(n, x, sx, indexL, indexU) {
    KSPlus <- max(KSPlusLevel(n, x, sx,indexL),0)
    KSMinus <- max(KSMinusLevel(n, x, sx,indexU),0)
    max(KSPlus, KSMinus)
}

HCStat<-function(x,indexL = seq_along(x),indexU = seq_along(x)){
    genericStatFunc(statFunc = HCStatFunc,
                     x = x,
                     indexL = indexL,
                     indexU = indexU
    )
}

BJStat<-function(x,indexL=NULL,indexU=NULL){
    genericStatFunc(statFunc = BJStatFunc,
                     x = x,
                     indexL = indexL,
                     indexU = indexU
    )
}

KSStat<-function(x,indexL=NULL,indexU=NULL){
    genericStatFunc(statFunc = KSStatFunc,
                     x = x,
                     indexL = indexL,
                     indexU = indexU
    )
}
SimesStat<-function(x,indexL=NULL,indexU=NULL){
    min(sort(x)*length(x)/seq_along(x))
}

