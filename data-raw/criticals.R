## code to prepare `criticals` dataset goes here
library(exceedance)
library(BiocParallel)
library(RedisParam)
#p <- SerialParam(progressbar = T)
p <- RedisParam(4)
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
    package_cache <- compute_critical(package_cache,p,statName, 
                                      alpha, n_list,
                                      indexLTxt="seq_len(n)", indexUTxt="seq_len(n)")
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
                                      indexLTxt="seq_len(n)", indexUTxt="NULL")
}

#####################################
## paper
#####################################
n_list <- seq_len(200)
statName <- "BJ"
alpha <- 0.1
package_cache <- compute_critical(package_cache,p,statName, 
                                  alpha, n_list,
                                  indexLTxt="seq_len(ceiling(n/2))", indexUTxt="NULL")

n_list <- seq_len(5000)
statName <- "BJ"
alpha <- 0.1
k_list <- 70:130
for(k in k_list){
    message(k)
    package_cache <- compute_critical(package_cache,cl,statName, 
                                      alpha, n_list,
                                      indexLTxt=paste0("seq_len(",k,")"), indexUTxt="NULL")
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
                                  indexLTxt=paste0("seq_len(",k,")"), indexUTxt="NULL")




save_criticals <- function(){
    package_cached_critical <- as.environment(package_cache)
    usethis::use_data(package_cached_critical, internal = TRUE, overwrite = TRUE)
}

save_criticals()

