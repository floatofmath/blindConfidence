simulate_maximum <- function(G,runs){
    maxsim <- rbindlist(mclapply(1:nrow(G),
                                 function(i) {c(G[i,],
                                                simVBIA(delta=G[i,]$delta,
                                                        sigma=G[i,]$sigma,
                                                        d=G[i,]$d,
                                                        s=G[i,]$s,
                                                        cf=1,
                                                        runs=runs,
                                                        alpha=.025,
                                                        beta=.2))}))
    ## compute relative bias
    maxsim[,var.rbias:= variance.bias/sigma^2]

    ## compute the theoretical lower bound for the variance bias
    maxsim[,bound:=lowerBound(n1,d)]

    ## runs <- 10^2
    ## bsim <- mclapply(1:nrow(gsim),function(i) {c(gsim[i,],
    ##                                              simMBIA(delta=gsim[i,]$delta,
    ##                                                      sigma=gsim[i,]$sigma,
    ##                                                      n1=gsim[i,]$n1,
    ##                                                      runs=runs))})

    ## bsim <- do.call('cbind',bsim)
    ## bsim <- as.data.frame(apply(bsim,1,unlist),row.names=NA)


    maxsim[,brannath:=.4*sqrt(2)*sigma/sqrt(n1)]
    return(maxsim)
}


select_mmeanbias <- function(maxsim){
    rbind(maxsim[,.SD[which.min(mean.bias)],by=n1,.SDcols=c('delta','sigma','d','n1','s','mean.bias')][,mean.bias:=NULL],
          maxsim[,.SD[which.min(uc.mean.bias)],by=n1,.SDcols=c('delta','sigma','d','n1','s','uc.mean.bias')][,uc.mean.bias:=NULL])
}


add_epsilon <- function(G,epsilon,by){
    rbindlist(apply(G,1,function(x) expand.grid(n1=x[1],delta=x[2]+seq(-epsilon,epsilon,by),sigma=x[3]+seq(-epsilon,epsilon,by),d=x[4],s=x[6])))
}


twostage_sim <- function(G,runs){
    simulate_maximum(
             add_epsilon(
                 select_mmeanbias(
                     simulate_maximum(G,runs)),epsilon=1,by=.1)
       ,runs)

}

library(parallel)
library(blindConfidence)
library(ggplot2)
library(reshape2)
library(plyr)
options(mc.cores=min(16,detectCores()-1))

load('maxsim_node2_150620.Rd')

G2s <- select_mmeanbias(maxsim)
G2e <- add_epsilon(G2s,.1,.01)
dim(G2e)
resim1 <- simulate_maximum(G2e,10^7)
save(resim1,file="resim_150639_node1.Rd")
