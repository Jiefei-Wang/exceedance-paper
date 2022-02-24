library(parallel)
library(doParallel)
library(foreach)

arg <- c(
    rep(list(list(
        host = "localhost",
        user = "jiefei"
    )),12),
    rep(list(list(
        host = "192.168.1.186",
        user = "jiefei"
    )),8),
    rep(list(list(
        host = "192.168.1.250",
        user = "lenovo"
    )),8),
)

cl <- makeCluster(arg)
