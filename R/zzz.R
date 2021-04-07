#' @useDynLib exceedance, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom digest digest
NULL

pkg_data <- new.env()
pkg_data$criticals <- NULL
pkg_data$use_cache <- TRUE
pkg_data$verbose <- FALSE

load_criticals <- function(){
    if(is.null(pkg_data$criticals)){
        if(exists("package_cached_critical")){
            pkg_data$criticals <- package_cached_critical
        }else{
            pkg_data$criticals <- new.env(parent = emptyenv())
        }
    }
}




.onLoad <- function(libname, pkgname){
    load_criticals()
}

.onUnload <- function(libpath){
    # env = new.env()
    # capture.output(tryCatch({
    #     env$mutex <- synchronicity::boost.mutex(packageCacheName, create = FALSE)
    # },
    # error = function(e) {
    #     env$mutex <- synchronicity::boost.mutex(packageCacheName, create = TRUE)
    # }), type = "message")
    # on.exit(synchronicity::unlock(env$mutex))
    # synchronicity::lock(env$mutex)
    # old_cache <- R.cache::loadCache(key = list(packageCacheName))
    # new_cache <- combine_env(old_cache,pkg_data$criticals)
    # R.cache::saveCache(new_cache, key = list(packageCacheName))
}