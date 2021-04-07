combineEnv <- function(dest, source) {
    source_names <- ls(source)
    dest_names <- ls(dest)
    for(v in source_names) {
        if(v %in% dest_names) {
            if(dest[[v]]==source[[v]]){
                next
            }else{
                stop(paste0("data does not match: ",v))
            }
        }
        dest[[v]] <- source[[v]]
    }
}

evalText<-function(x, env){
    if(!is.character(x))
        stop("x is not a character")
    eval(parse(text=x),envir=env)
}


compute_critical<-function(package_cache, cl, statName, alpha, n_list, indexLTxt, indexUTxt){
    # clusterExport(cl, c("alpha","statName","indexLTxt","indexUTxt"),envir = environment())
    
    cur_criticals <- foreach(n = n_list, 
                             .combine=c, .multicombine = TRUE, .inorder= FALSE,
                             .export = c("evalText"))%dopar%{
        env <- environment()
        indexL <- evalText(indexLTxt,env)
        indexU <- evalText(indexUTxt,env)
        
        key_critical <- 
            exceedance:::compute_key_critical(
                statName=statName, n=n, alpha=alpha,
                indexL=indexL,indexU=indexU)
        if(is.null(key_critical)){
            return(NULL)
        }
        result <- list()
        result[[key_critical$key]] <- key_critical$critical
        result
                             }
    message("Finish loop")
    # return(cur_criticals)
    package_cache[names(cur_criticals)] <- cur_criticals
    package_cache
}
