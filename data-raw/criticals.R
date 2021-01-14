## code to prepare `criticals` dataset goes here
## devtools::load_all()
library(exceedance)
library(parallel)
source("data-raw/functions.R")
cl <- makeCluster(detectCores(logical = TRUE))
invisible(clusterEvalQ(cl,library(exceedance)))
package_cache <- exceedance:::pkg_data$criticals
clusterExport(cl, "evalText")


n<-100
#####################################
## two-sided
## BJ, KS, HC
#####################################
n_list <- seq_len(n)
statName <- "BJ"
alpha_list <- c(0.05, 0.1)
for(alpha in alpha_list){
    compute_critical(package_cache,cl,statName, 
                     alpha, n_list,
                     indexL="seq_len(n)", indexU="seq_len(n)")
}

#####################################
## BJ+, KS+, HC+
#####################################
n_list <- seq_len(n)
statName <- "BJ"
alpha_list <- c(0.05, 0.1)
for(alpha in alpha_list){
    compute_critical(package_cache,cl,statName, 
                     alpha, n_list,
                     indexL="seq_len(n)", indexU="NULL")
}

#####################################
## BJ 1-10
#####################################
n_list <- seq_len(n)
statName <- "BJ"
alpha_list <- c(0.05, 0.1)
for(alpha in alpha_list){
    compute_critical(package_cache,cl,statName, 
                     alpha, n_list,
                     indexL="seq_len(10)", indexU="NULL")
}
for(alpha in alpha_list){
    compute_critical(package_cache,cl,statName, 
                     alpha, n_list,
                     indexL="seq_len(10)", indexU="seq_len(10)")
}



package_cached_critical <- package_cache
usethis::use_data(package_cached_critical, internal = TRUE, overwrite = TRUE)
