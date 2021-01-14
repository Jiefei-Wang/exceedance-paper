## Indivial level function
HCPlusLevel <- function(n, x, sx, index){
    if (length(index) == 0)
        return(numeric(0))
    sqrt(n) * (index / n - sx[index]) / sqrt(sx[index] * (1 - sx[index]))
}
HCMinusLevel<- function(n, x, sx, index){
    if (length(index) == 0)
        return(numeric(0))
    sqrt(n) * (sx[index]-(index - 1) / n) / sqrt(sx[index] * (1 - sx[index]))
}
BJPlusLevel <- function(n, x, sx, index) {
    if (length(index) == 0)
        return(numeric(0))
    vapply(index, function(i)
        pbeta(sx[i], i, n - i + 1),numeric(1))
}
BJMinusLevel <- function(n, x, sx, index) {
    if (length(index) == 0)
        return(numeric(0))
    1 - BJPlusLevel(n, x, sx, index)
}
KSPlusLevel <- function(n, x, sx, index) {
    if (length(index) == 0)
        return(numeric(0))
    index / n - sx[index]
}
KSMinusLevel <- function(n, x, sx, index) {
    if (length(index) == 0)
        return(numeric(0))
    sx[index] - (index - 1) / n
}


## get local critical value
HCLocalCritical<-function(statValue,n){
    statValue <- statValue/sqrt(n)
    a<-1+statValue^2
    ## lower
    const<-seq_len(n)/n
    b<--2*const-statValue^2
    l <- (-b-sqrt(b^2-4*a*const^2))/2/a
    ## upper
    const<-(seq_len(n)-1)/n
    b<--2*const-statValue^2
    h <- (-b+sqrt(b^2-4*a*const^2))/2/a
    list(l =l,h= h)
}

BJLocalCritical<-function(statValue,n){
    l=vapply(seq_len(n),function(x)qbeta(statValue,x,n-x+1),numeric(1))
    h=vapply(seq_len(n),function(x)qbeta(1 - statValue,x,n-x+1),numeric(1))
    list(l =l,h= h)
}

KSLocalCritical<-function(statValue,n){
    l <- seq_len(n)/n - statValue
    h <- statValue + seq_len(n)/n-1/n
    l[l<0]=0
    h[h>1]=1
    list(l =l,h= h)
}
SimesLocalCritical<-function(statValue,n){
    l <- statValue/n*seq_len(n)
    h <- rep(1, n)
    list(l =l,h= h)
}

