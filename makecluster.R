## This is my secret file to create a parallel cluster
## please DO NOT run them

library(parallel)
library(doParallel)
library(foreach)

arg <- c(
    rep(list(list(
        host = "localhost",
        user = "jiefei"
    )),12)
)

arg <- c(arg,
         rep(list(list(
             host = "192.168.1.186",
             user = "jiefei"
         )),8))

# arg <- c(arg,
#          rep(list(list(
#              host = "192.168.1.250",
#              user = "lenovo"
#          )),8))



cl <- makeCluster(arg)
registerDoParallel(cl)


stopCluster(cl)
# 
# parLapply(cl, 1, function(x) .libPaths())
# 
# foreach(i=1:2)%dopar%{
#     i
# }
