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

add_criticals <- function(x){
    vars <- names(x)
    vars <- vars[!vars%in%names(pkg_data$criticals)]
    for(i in vars){
        pkg_data$criticals[[i]] <- x[[i]]
    }
}