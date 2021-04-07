## code to prepare `criticals` dataset goes here
library(exceedance)

source("data-raw/functions.R")
package_cache <- as.list(exceedance:::pkg_data$criticals)

n<-10000
#####################################
## two-sided
## BJ, KS, HC
#####################################
n_list <- seq_len(n)
statName <- "BJ"
alpha_list <- c(0.1)
for(alpha in alpha_list){
    package_cache <- compute_critical(package_cache,cl,statName, 
                                      alpha, n_list,
                                      indexL="seq_len(n)", indexU="seq_len(n)")
}

#####################################
## BJ+, KS+, HC+
#####################################
n_list <- seq_len(n)
statName <- "BJ"
alpha_list <- c(0.1)
for(alpha in alpha_list){
    package_cache <- compute_critical(package_cache,cl,statName, 
                                      alpha, n_list,
                                      indexL="seq_len(n)", indexU="NULL")
}

#####################################
## paper
#####################################
n_list <- seq_len(2000)
statName <- "BJ"
alpha <- 0.1
k_list <- 40:120
for(k in k_list){
    message(k)
    package_cache <- compute_critical(package_cache,cl,statName, 
                                      alpha, n_list,
                                      indexL=paste0("seq_len(",k,")"), indexU="NULL")
}

n_list <- seq_len(5000)
statName <- "BJ"
alpha <- 0.1
k_list <- 70:130
for(k in k_list){
    message(k)
    package_cache <- compute_critical(package_cache,cl,statName, 
                                      alpha, n_list,
                                      indexL=paste0("seq_len(",k,")"), indexU="NULL")
    message(k)
}

#####################################
## paper timing
#####################################
n_list <- seq_len(1000)
statName <- "BJ"
alpha <- 0.1
k <- 4
package_cache <- compute_critical(package_cache,cl,statName, 
                                  alpha, n_list,
                                  indexL=paste0("seq_len(",k,")"), indexU="NULL")




save_criticals <- function(){
    package_cached_critical <- as.environment(package_cache)
    usethis::use_data(package_cached_critical, internal = TRUE, overwrite = TRUE)
}

save_criticals()


removeQueue(queue)
