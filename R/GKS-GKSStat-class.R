.GKSStat<-function(statName,statValue,n,indexL=1:n,indexU=1:n, data){
    stat <- list(statName = statName,
                 statValue = statValue,
                 n=n,
                 indexL=indexL,
                 indexU=indexU,
                 data = data)
    structure(stat,class = "GKSStat")
}
#' Class methods for the GKSStat S3 object
#' 
#' 
#' @param x the GKSStat S3 object
#' @param ... Ignored.
#' @examples 
#' 
#' ## Generate samples
#' x <- rbeta(10, 1, 2)
#' 
#' ## Perform KS test
#' GKSStat(x = x, statName = "KS")
#' 
#' @return  
#' print: invisible `x`
#' other fucntions: a numeric value
#' 
#' @rdname classMethod
#' @export
print.GKSStat <- function(x,...){
    #print(x$statValue)
    # class(x)=NULL
    # print(x)
    cat("The", getStatName(x), "test statistics\n")
    cat("Sample size:",getSampleSize(x),"\n")
    cat("Stat value:", getStatValue(x),"\n")
    if(!is.null(getPvalue(x)))
        cat("P-value:", getPvalue(x))
    invisible(x)
}

#' @rdname classMethod
#' @export
getStatName<-function(x){
    x$statName
}
#' @rdname classMethod
#' @export
getStatValue<-function(x){
    x$statValue
}
#' @rdname classMethod
#' @export
getSampleSize <- function(x){
    x$n
}
#' @rdname classMethod
#' @export
getPvalue<-function(x){
    x$pvalue
}
#' @rdname classMethod
#' @export
getData<-function(x){
    x$data
}