##' run coverage sims 
##'
##' @title coverage sim
##' @param base initial simulation across a grid to base the search on
##' @param search number of simulation runs for the wide search
##' @param resim number of simulation runs for the resimulation of found maxima
##' @param use_mclapply2 use mclapply2 in simulations
##' @param mem_factor number of cores that can handle 10^7 simulation runs, will be used to automatically determine the number of cores 
##' @param epsilon width of parameter grid around initial values
##' @param by fineness of grid around initial values
##' @return dataframe with simulation results
##' @author float
##' @export
coverage_sim <- function(base=maxsim,search=10^4,resim=10^6,use_mclapply2=TRUE,mem_factor=16,epsilon=.18,by=.06){
## Seperate file to run coverage simulations as they tend to crash the rstudio server
    library(parallel)
    ## devtools::install_github('floatofmath/blindConfidence')
    library(bt88.03.704)
    library(ggplot2)
    library(reshape2)
    library(plyr)

    data(maxsim)
    maxsim$ylow <- 0

    maxcoverageopt1 <- select_results(base,c('total.prob','upper.prob','uc.total.prob','uc.upper.prob'),base_columns=c('n1','tn1','delta','sigma','d','s','total.prob','upper.prob'),functional=which.max)

    G2coverage <- add_epsilon(maxcoverageopt1,epsilon,by,c('delta','sigma'))
##plot_grid(maxsim,maxcoverageopt1,G2coverage,'sigma')
    cat("Run search\n")
    s.cores <- floor(mem_factor* 10^7/search)
    r.cores <- floor(mem_factor* 10^7/resim)
    options(mc.cores=min(s.cores,detectCores()-1))
    maxcoveragesim2 <- simulate_batch(G2coverage,search,use_mclapply2)
    gc()

## reshape and resim
    maxcoverageopt2 <- select_results(maxcoveragesim2,c('total.prob','upper.prob','uc.total.prob','uc.upper.prob'),base_columns=c('n1','tn1','delta','sigma','d','s'),functional=which.max)
    options(mc.cores=min(r.cores,detectCores()-1))
    cat("Run resimulation\n")
    maxcoveragesim3 <- simulate_batch(maxcoverageopt2,resim,use_mclapply2)
    gc()


    fname <- paste('resim_coverage',Sys.info()['nodename'],'_',format(Sys.time(),"%y%m%d-%H%M"),'.Rd',sep='')
    save(maxcoveragesim3,file=fname)
    cat("Saved results to:",fname,"\n")
    maxcoveragesim3
}

