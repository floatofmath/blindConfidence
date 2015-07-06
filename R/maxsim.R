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
    rbind(maxsim[,.SD[which.min(mean.bias)],by=n1,.SDcols=c('delta','sigma','d','n1','s','mean.bias')][,c("what.max","mean.bias"):=list("mean.bias",NULL)],
          maxsim[,.SD[which.min(uc.mean.bias)],by=n1,.SDcols=c('delta','sigma','d','n1','s','uc.mean.bias')][,c("what.max","uc.mean.bias"):=list("uc.mean.bias",NULL)])
}

select_mvarbias <- function(maxsim){
    rbind(maxsim[,.SD[which.min(variance.bias)],by=n1,.SDcols=c('delta','sigma','d','n1','s','mean.bias')][,c("what.max","variance.bias"):=list("variance.bias",NULL)],
          maxsim[,.SD[which.min(uc.variance.bias)],by=n1,.SDcols=c('delta','sigma','d','n1','s','uc.mean.bias')][,c("what.max","uc.variance.bias"):=list("uc.variance.bias",NULL)])
}

select_mncoverage <- function(maxsim){
    rbind(maxsim[,.SD[which.max(total.prob)],by=n1,.SDcols=c('delta','sigma','d','n1','s','mean.bias')][,c("what.max","total.prob"):=list("total.prob",,NULL)],
          maxsim[,.SD[which.max(uc.total.prob)],by=n1,.SDcols=c('delta','sigma','d','n1','s','uc.mean.bias')][,c("what.max","uc.total.prob"):=list("uc.total.prob",NULL)])
    rbind(maxsim[,.SD[which.max(upper.prob)],by=n1,.SDcols=c('delta','sigma','d','n1','s','mean.bias')][,c("what.max","upper.prob"):=list("upper.prob",NULL)],
          maxsim[,.SD[which.max(uc.upper.prob)],by=n1,.SDcols=c('delta','sigma','d','n1','s','uc.mean.bias')][,c("what.max","uc.upper.prob"):=list("uc.upper.prob",NULL)])
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
options(mc.cores=min(30,detectCores()-1))

load('maxsim_node2_150620.Rd')
set.seed(980920)

G2s <- select_mmeanbias(maxsim)
G2e <- add_epsilon(G2s,.1,.01)
dim(G2e)
resim1 <- simulate_maximum(G2e,4*10^6)
print(fname <- paste('resim1_',Sys.info()['nodename'],'_',format(Sys.time(),"%y%m%d"),'.Rd',sep=''))
save(resim1,file=fname)




G3s <- select_mmeanbias(maxsim)
G3e <- add_epsilon(G3s,.1,.002)
dim(G3e)
resim2 <- simulate_maximum(G3e,4*10^6)
save(resim2,file="resim_150703_node1.Rd")
