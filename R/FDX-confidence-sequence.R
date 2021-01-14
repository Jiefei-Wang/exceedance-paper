get_confidence_sequence<-function(x, param, alpha){
    profile <- exceedance_profile(x, param)
    CS <- rep(0, length(x))
    for(i in seq_along(x)){
        CS[i] <- exceedance_confidence(profile,alpha,sri=seq_len(i))
    }
    CS
}