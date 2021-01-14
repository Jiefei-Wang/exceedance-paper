use_cache <- function(x){
    if(!missing(x))
        pkg_data$use_cache <- x
    else
        pkg_data$use_cache
}

verbose <- function(x){
    if(!missing(x))
        pkg_data$verbose <- x
    else
        pkg_data$verbose
}

