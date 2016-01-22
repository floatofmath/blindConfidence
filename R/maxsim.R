##' Simulate bias of blinded sample size review for a range of parameter settings
##'
##' To use mclapply2 from bt88.03.704 you have to install this library (e.g. \code{devtools::install_github('floatofmath/bt88.03.704')}) unfortunately the statusbar implementation does not work with Rstudio.
##' 
##' @title Batch simulate blinded sample size review bias
##' @param G Object with paramter settings for \code{delta}, \code{sigma} and \code{n1} in each row 
##' @param runs Number of simulation runs 
##' @param use_mclapply2 use mclapply2 from package bt88.03.704 which implements a status bar (see details)
##' @param multicore should multicore parallellization be used
##' @return Object with parameters and simulation results in each row
##' @author Florian Klinglmueller
##'
##' @examples
##' data(maxsim)
##' G2s <- select_results(maxsim,c('mean.bias'))
##' G2e <- add_epsilon(G2s,.05,.05,subset=c('delta','sigma'))
##' \dontrun{plot_grid(maxsim,G2s,G2e,what='delta')
##' resim1 <- simulate_batch(G2e,10^3)}
##' 
##' @export
simulate_batch <- function(G,runs,use_mclapply2=FALSE,multicore=TRUE){
    if(multicore) {
        require(parallel)
        mcla <- mclapply
        if(use_mclapply2) {
            require(bt88.03.704)
            mcla <- mclapply2
        }
    } 
    if(!multicore) {
        mcla  <-  lapply
    }
    args <- names(formals(simVBIA))
    grid <- G[,names(G) %in% args,with=F]
    maxsim <- rbindlist(mcla(1:nrow(grid),
                                 function(i) {c(G[i,],
                                                do.call('simVBIA',c(grid[i,],
                                                        cf=1,
                                                        runs=runs,
                                                        alpha=.025,
                                                        beta=.2)))}))
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

##' Select something for each first stages sample size from a \code{\link{data.table}} of simulation results 
##'
##' @title Select simulation results
##' @param maxsim \code{\link{data.table}} of simulation results
##' @param what character vector of column names according to which the selection should be performed
##' @param base_columns columns to keep from original data.table use \code{NULL} to keep all columns
##' @param functional function applied columns \code{what} whithin each unique first stage sample size that returns an index of the selected rows 
##' @return \code{data.table} with selected rows, plus one row giving the column name the selection was based on
##' @author Florian Klinglmueller
##' @examples
##'
##' data(maxsim)
##' ## select largest (absolute) mean.bias for each first stage sample size
##' select_results(maxsim,'mean.bias')
##'
##' @export
select_results <- function(maxsim,what,base_columns=c('delta','sigma','d','n1','s'),functional=which.min){
    if(is.null(base_columns)){
        base_columns <- names(maxsim)
    }
    select <- function(what){
        .c <- substitute(.I[functional(column)],list(column=as.name(what)))
        return(maxsim[maxsim[,eval(.c),by=n1]$V1,][,base_columns,with=F][,what.max:=what])
    }
    rbindlist(lapply(what,select))
}


##' Expand and refine parameter grid for simulation
##'
##' @title Expand grid
##' @param G matrix with one parameter setting in each line
##' @param epsilon width of interval to be added around each points in G
##' @param by distance between points to be added within intervals around points in G
##' @param subset subset of parameters for which grid points should be added
##' @return object with the same columns as G 
##' @author Florian Klinglmueller
##'
##' @export
add_epsilon <- function(G,epsilon,by,subset=colnames(G)){
    set_class <- FALSE
    if(!'data.table' %in% class(G)){
        set_class <- TRUE
        G_class <- class(G)
    }
    if(length(epsilon) != 1 | length(by) != 1){
        if(length(epsilon) != length(subset) | length(by) != length(subset))
            stop("Vector of interval widths has to have same length as number of columns to refine")
        if(length(by) != length(epsilon)){
            ## replicate the shorter one
            by <- by+epsilon*0
            epsilon <- epsilon+by*0
        }
    } else {
        epsilon <- rep(epsilon,length(subset))
        by <- rep(by,length(subset))
    }
    int <- lapply(1:length(by),function(i) seq(-epsilon[[i]],+epsilon[[i]],by[[i]]))
    for(j in 1:length(subset)){
        .c <- substitute(`:=`(what,as.numeric(outer(int[[j]],G[,what],'+'))),list(what=as.name(subset[j])))
        G <- G[rep(1:.N,each=length(int[[j]]))][,eval(.c)]
    }
    if(set_class)
        G <- as(G,G_class)
    return(G)
}


##' Plot to compare different grids of parameter values
##'
##' @title Plot parameter grids
##' @param maxsim Results of previous simulation
##' @param G2s Parameter values at which optimum occured
##' @param G2e Refined parameter grid
##' @param what For which parameter to plot the grid
##' @author Florian Klinglmueller
##'
##' @export
plot_grid <- function(maxsim,G2s,G2e,what){
  ggplot(maxsim,aes_string(x='n1',y=what))+geom_point()+geom_point(aes_string(x='n1',y=what),data=G2e,col='green')+geom_point(aes_string(x='n1',y=what),data=G2s,col='red')
}

    


smooth_criminal <- function(nominees,what){
    nominees <- subset(nominees,what.max==what)
    candidates <- list()
    candidates$nominees$n1 <- nominees$n1
    candidates$nominees$delta <- nominees$delta
    candidates$nominees$sigma <- nominees$sigma
    candidates$models <- list()
    candidates$models$delta <- loess(delta~n1,data=nominees,subset={what.max=="total.prob"})
    candidates$models$sigma <- loess(sigma~n1,data=nominees,subset={what.max=="total.prob"})
    candidates$plx$delta <- predict(candidates$models$delta,se=T)
    candidates$plx$sigma <- predict(candidates$models$sigma,se=T)
    candidates
}
