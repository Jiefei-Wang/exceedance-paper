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


compute_critical<-function(package_cache, cl, statName, alpha, n_list, indexL, indexU){
    clusterExport(cl, c("alpha","statName","indexL","indexU"),envir = environment())
    cur_criticals <- parSapplyLB(cl = cl, n_list, function(n){
        env <- environment()
        critical <- 
            exceedance:::get_critical(
                statName=statName, n=n, alpha=alpha,
                indexL=evalText(indexL,env),indexU=evalText(indexU,env))
    })
    cache <- clusterEvalQ(cl, exceedance:::pkg_data$criticals)
    for(e in cache){
        combineEnv(package_cache,e)
    }
}