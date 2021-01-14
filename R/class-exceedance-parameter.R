.exceedance_parameter <- function(...){
    x <- structure(list(...), class = "exceedance_parameter")
    .valid_exceedance_object(x)
    x
}
.exceedance_profile <- function(...){
    structure(list(...), class = "exceedance_profile")
}

.valid_exceedance_object <- function(x){
    stopifnot(is.character(x$profile_func))
    stopifnot(is.character(x$confidence_func))
    stopifnot(is.character(x$inference_func))
}

title <-function(x){
    is.null(x$title)||x$title
}

`title<-` <- function(x, value){
    x$title <- value
    x
}

truncated_print <- function(x, len = 5L){
    if(length(x)==0){
        return("NULL")
    }
    x_len <- length(x)
    x_print_len <- min(x_len, len)
    x_diff <- x_len - x_print_len
    result <- paste0(x[seq_len(x_print_len)], collapse = ",")
    if(x_diff != 0){
        result <- paste0(result, ",...(skip ",x_diff," elements)")
    }
    result
}

show_method_name <- function(name){
    name["fast_GW"%in%name] <- "fast GW"
    name["general_GW"%in%name] <- "general GW"
    name["combine_GW"%in%name] <- "combined GW"
    name
}


#' @export
print.exceedance_parameter<-function(x,...){
    if(title(x)){
        cat("An S3 `exceedance_parameter` object:\n") 
    }
    method <- x$method
    if(method == "fast_GW"){
        show_params_fast_GW(x)
    }
    if(method == "general_GW"){
        show_params_general_GW(x)
    }
    if(method == "combine_GW"){
        show_params_combine_GW(x)
    }
    invisible(x)
}

# parms <- list(method = "fast_GW",
#               postfix = postfix,
#               postfix_profile = postfix_profile,
#               postfix_bound = postfix_bound,
#               postfix_inference = postfix_inference,
#               statistic = statistic,
#               param1 = param1,
#               param2 = param2,
#               range_type=range_type)

show_params_fast_GW <- function(x){
    param1 <- x$param1
    param2 <- x$param2
    statistic <- x$statistic
    range_type <- x$range_type
    method <- x$method
    
    
    if(is.null(param1))
        param1 <- "NULL"
    if(is.null(param2))
        param2 <- "NULL"
    
    cat("Method:", show_method_name(method),"\n")
    cat("Statistic:", statistic,"\n")
    cat("param1:", truncated_print(param1), "\n")
    cat("param2:", truncated_print(param2), "\n")
    cat("range type:",range_type,"\n")
}

show_params_general_GW<-function(x){
    statistic <- x$statistic
    method <- x$method
    
    cat("Method:", show_method_name(method),"\n")
    cat("Algorithm:", statistic,"\n")
}

show_params_combine_GW<-function(x){
    weight <- x$alpha_weight
    method <- x$method
    test_params <- x$test_params
    test_methods <- vapply(test_params,function(x)show_method_name(x$method),character(1))
    test_algorithms <- vapply(test_params,function(x)x$statistic,character(1))
    
    
    cat("Method:", show_method_name(method),"\n")
    cat("Contained methods:", truncated_print(test_methods),"\n")
    cat("Contained algorithms:", truncated_print(test_algorithms),"\n")
    cat("Weight:", truncated_print(weight),"\n")
}

