getGKSIndex <- function(statName, n, index, indexL, indexU){
    if(!is.null(index)){
        stopifnot(is.null(indexL))
        stopifnot(is.null(indexU))
        indexL <- index
        indexU <- index
    }else{
        if(is.null(indexL)&&is.null(indexU)){
            indexL <- seq_len(n)
            indexU <- seq_len(n)
        }
    }
    
    side <- substring(statName,nchar(statName))
    sideIndicator <- which(side==c("+","-"))
    if(length(sideIndicator)!=0){
        statName <- substr(statName,1,nchar(statName)-1)
        if(sideIndicator == 1){
            indexU <- NULL
        }else{
            indexL <- NULL
        }
    }
    list(indexL = indexL, indexU = indexU, statName = statName)
}


call_func <- function(root, prefix=NULL, postfix=NULL, ...){
    func_name <- paste0(c(prefix,root,postfix),collapse="")
    do.call(func_name,args = list(...))
}
