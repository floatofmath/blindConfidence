#' Simulation of estimation after blinded interim analysis
#'
#' This package provides all code used to perform and summarize the
#' simulation study shown in our paper "Estimation after blinded
#' interim analysis". The code for all simulations is provided in
#' \code{vignette{'simulations'}}, simulation results for coverage
#' probabilities and biases of estimates are saved in
#' \code{data{gridsim}}; the results from the case study simulation
#' are saved in \code{data{casesim}}; the results for maximizing the
#' bias of the effect size and variance estimate as well as the
#' inflation of the non-coverage probability are saved in
#' \code{data{maxsim}}, \code{data{maxmeansim}},
#' \code{data{maxvariancesim}} and
#' \code{data{maxcoveragesim}}. Rebuilding the vignette will use the
#' saved results to redraw all the graphics from the paper. In order
#' to rerun the simulations the corresponding code chunks need to be
#' configured for evaluation. NOTE: although this package should run
#' on Windows it has been tested and optimized under Linux.
#' 
#'
#' @name blindConfidence-package
#' @aliases blindConfidence
#' @docType package
#' @title Simulation of estimation after blinded interim analysis
#' @author Florian Klinglmueller
#' @import data.table 
NULL

#' Compute sample size of z-test
#'
#' \code{zss} computes per-group the sample size for the z-test
#' given standard deviation \code{s} mean difference \code{d} and
#' Type I and Type II errors \code{alpha} and \code{beta}
#' 
#' @param s standard deviation
#' @param d mean difference
#' @param alpha Type I error
#' @param beta Type II error
#' @examples
#' zss(1,.5,.025,.2)
#' @export
zss <- function(s,d,alpha,beta){
  2 * s^2 / d^2 * (qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2
}

#' Compute standard deviation from planned first stage sample size
#'
#' \code{zsd} computes the standard deviation that would result in
#' a first stage sample size \code{n1} given mean difference \code{d}
#' Type I and and Type II error- \code{alpha} and \code{beta} and
#' assuming that the first stage should include the fraction \code{v} of
#' the total sample size.
#'
#' @param n1 first stage sample size
#' @param d mean difference
#' @param v fraction of the total sample size
#' @param alpha Type I error probability
#' @param beta Type II error probability
#' @examples
#' zsd(15,.5,3/4,.025,.2)
#' @export
zsd <- function(n1,d,v = 1/2,alpha,beta){
    n <- n1/(v*2)
    sqrt(n * d^2 / ((qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2))
}


#' Simulate estimation after blinded interim analyses 
#'
#' \code{simVBIA} simulates \code{runs} adaptive trials with
#' blinded interim analysis and assesses the bias of the variance 
#' and mean estimates. Sample size review is performed either using
#' an interim estimate of the pooled variance that is adjusted
#' based on the assumed effect size or the uncorrected lumped variance.
#'
#' @param delta true effect size
#' @param sigma true standard deviation
#' @param d assumed effect size
#' @param n1 first stage sample size in treatment group
#' @param alpha = .025 pre-specified alpha level
#' @param beta pre-specified beta level
#' @param v timing of interim analysis
#' @param s assumed standard deviation, is more or less the same as n1
#' @param runs number of simulation runs
#' @param n2 first stage sample size in control group
#' @param n2min minimum secondstage sample size
#' @param fulldata return full data 
#' @param cf correction term added to the sample size rule
#' @template simvbia
#' @examples
#' simVBIA(.5,1,1,n1=10,n2min=2,runs=1000)
#' @export
simVBIA <- function(delta, # true effect size
                   sigma, # true standard deviation
                   d, # assumed effect size
                   n1,# ceiling(v * zss(s,d,alpha,beta)), first stage sample size in treatment group
                   alpha = .025, # pre-specified alpha level
                   beta = .2, # pre-specified beta level
                   v = 1/2, # timing of interim analysis
                   s = zsd(n1,d,v,alpha,beta), # assumed standard deviation, is more or less the same as n1
                   runs, # number of simulation runs
                   n2 =n1, # first stage sample size in control group
                   n2min = 0, # minimum secondstage sample size
                   fulldata = FALSE, #  
                   cf = 0
                   ){
  ## pre-compute t-distribution quantiles
  qttable=qt(1-alpha,1:1000)
  xi=(qnorm(alpha,lower.tail=FALSE)+qnorm(beta,lower.tail=FALSE))^2/d^2
  ## total first stage sample size
  n=n1+n2 
  lamb=n1*n2/n ##???

  ## simulated first stage means and sum of squares
  ## group 1
  sim <- data.frame(x1=rnorm(runs,sd=sigma/sqrt(n1)))

  ## here we also generate a second independent first stage
  ## this is used for code testing purposes as estimates using
  ## these values should be unbiased, and if not point to errors
  ## in the calculation
  ## x1r=rnorm(runs,sd=sigma/sqrt(n1))
  ## group 2
  sim$x2 <- with(sim,rnorm(runs,mean=delta,sd=sigma/sqrt(n2)))
  ## x2r=rnorm(runs,mean=delta,sd=sigma/sqrt(n2))

  ## sum of squares
  ## group 1
  sim$sq1 <- with(sim,sigma^2*rchisq(runs,n1-1))
  ## sq1r=sigma^2*rchisq(runs,n1-1)

  ## group 2
  sim$sq2 <- with(sim,sigma^2*rchisq(runs,n2-1))
  ## sq2r=sigma^2*rchisq(runs,n2-1)

  # overall standard deviation s_adj following Kieser &Friede 2003
  # may become <0, we set it to zero in that case (?) check
  sim$ps <- with(sim, pmax((sq1+sq2+lamb*{x1-x2}*{x1-x2})/(n-1)- d^2*n1/(2*n-2),0))

  ## without adjustment for mean guess
  sim$psu <- with(sim,(sq1+sq2+lamb*{x1-x2}*{x1-x2})/(n-1))

  ## round of errors!
  ## m1=pmax(round(2*xi*ps)-n1,0)
  sim$m1 <- with(sim,pmax(ceiling(2*xi*ps)-n1+cf,n2min))
  sim$m1u <- with(sim, pmax(ceiling(2*xi*psu)-n1+cf,n2min))
  sim$ps <- with(sim,NULL)
  sim$psu <- with(sim,NULL)
  ## replace by external values
  sim$m2 <- with(sim,m1)
  sim$m2u <- with(sim,m1u)

  ## secondstage means
  sim$y1 <- with(sim,rnorm(runs,sd=sigma/sqrt(m1)))
  sim$y2 <- with(sim,rnorm(runs,sd=sigma/sqrt(m2),mean=delta))
  sim$y1u <- with(sim,rnorm(runs,sd=sigma/sqrt(m1u)))
  sim$y2u <- with(sim,rnorm(runs,sd=sigma/sqrt(m2u),mean=delta))
  
  # set the  means to zero if second stage sample size is 0 
  sim$y1[sim$m1<1] <- 0
  sim$y2[sim$m2<1] <- 0
  sim$y1u[sim$m1u<1] <- 0
  sim$y2u[sim$m2u<1] <- 0
  
  # sum of squares for the second stage
  sim$tq1 <- with(sim,sigma^2*rchisq(runs,df=pmax(m1-1,1)) )
  sim$tq2 <- with(sim,sigma^2*rchisq(runs,df=pmax(m2-1,1)))
  sim$tq1u <- with(sim,sigma^2*rchisq(runs,df=pmax(m1u-1,1)) )
  sim$tq2u <- with(sim,sigma^2*rchisq(runs,df=pmax(m2u-1,1)))
  
  # set the  means to zero if second stage sample size is 0 
  sim$tq1[sim$m1<2] <-  0
  sim$tq2[sim$m2<2] <-  0
  sim$tq1u[sim$m1u<2] <- 0
  sim$tq2u[sim$m2u<2] <- 0

  sim$ll <- with(sim,{n1+m1}+{n2+m2}-2)
  sim$llu <- with(sim,{n1+m1u}+{n2+m2u}-2)


  sim$sp1u <- with(sim,sqrt((tq1u+sq1+(y1u-x1)^2*n1*m1u/{n1+m1u}+tq2u+sq2+(y2u-x2)^2*n2*m2u/{n2+m2u})/llu))
  sim$sp1 <- with(sim,sqrt((tq1+sq1+(y1-x1)^2*n1*m1/{n1+m1}+tq2+sq2+(y2-x2)^2*n2*m2/{n2+m2})/ll))
  ## t-test denominator (see kieser & friede 2003)
  sim$sp <- with(sim,sp1*sqrt((1/{n1+m1}+1/{n2+m2})))
  sim$spu <- with(sim,sp1u*sqrt((1/{n1+m1u}+1/{n2+m2u})))

  ## t-test numerator
  sim$diffmu <- with(sim,(n2*x2+m2*y2)/{n2+m2}-(n1*x1+m1*y1)/{n1+m1})
  sim$bar <- with(sim,sp*qttable[pmin(ll,1000)])# hier werden die approx. quantile verwendet. hat das einen impac)
  
  sim$diffmuu <- with(sim,(n2*x2+m2u*y2u)/{n2+m2u}-(n1*x1+m1u*y1u)/{n1+m1u})
  sim$baru <- with(sim,spu*qttable[pmin(llu,1000)])# hier werden die approx. quantile verwendet. hat das einen impac)

  err=c(
      lower.prob=with(sim,sum({diffmu-bar}>delta)/runs),
      upper.prob=with(sim,sum({diffmu+bar}<delta)/runs),
      total.prob=with(sim,sum(({diffmu-bar}>delta)+({diffmu+bar}<delta))/runs),
      mean.bias=with(sim,sum(diffmu-delta)/runs),
      variance.bias=with(sim,sum(sp1^2-sigma^2)/runs), #? sp^2 -1]
      vm =with(sim,var(diffmu)),
      ev =with(sim, mean(sp^2)), #? var(diffmu)-sem^2]
      exv =with(sim,mean(sigma^2*(1/{n1+m1}+1/{n2+m2}))),
      root.mse=with(sim,sqrt(sum((diffmu-delta)^2)/runs)),
      mean.m1=with(sim,median(m1)),
      min.m1=with(sim,unlist(quantile(m1,.01,na.rm=TRUE))),
      N10 = with(sim,sum(m1==0)/runs),
### uncorrected    
      uc.lower.prob=with(sim,sum({diffmuu-baru}>delta)/runs),
      uc.upper.prob=with(sim,sum({diffmuu+baru}<delta)/runs),
      uc.total.prob=with(sim,sum(({diffmuu-baru}>delta)+({diffmuu+baru}<delta))/runs),
      uc.mean.bias=with(sim,sum(diffmuu-delta)/runs),
      uc.variance.bias=with(sim,sum(sp1u^2-sigma^2)/runs), #? sp^2 -]1
      uc.vm=with(sim,var(diffmuu)), #]?
      uc.ev =with(sim,mean(spu^2)),
      uc.exv =with(sim, mean(sigma^2*(1/{n1+m1u}+1/{n2+m2u}))),
      uc.root.mse=with(sim,sqrt(sum((diffmuu-delta)^2)/runs)),
      uc.mean.m1=with(sim,median(m1u)),
      uc.min.m1=with(sim,unlist(quantile(m1u,.01,na.rm=TRUE))),
      uc.n10 = with(sim,sum(m1==0)/runs),
      tn1 = with(sim,mean(n1))
      )
  # list of lower prob, upper prob, total prob, mean bias, variance bias, standard error bias, sqrt(mse)
  #list(err=err,m1=m1)
  ## data <- list("x2" = x2, "y2" = y2,"sq2" = sq2,"tq2" = tq2, "x1" = x1, "y1" = y1,
  ##              "sq1" = sq1, "tq1" = tq1, "m1" = m1, "ps" = ps,
  ##              "y2u" = y2u ,"tq2u" = tq2u, "y1u" = y1u,
  ##              "tq1u" = tq1u, "m1u" = m1u, "psu" = psu)
  if(fulldata) return(list("sim"=sim,
                           "delta"=delta,
                           "runs"=runs,
                           "sigma"=sigma,
                           "n1"=n1,
                           "n2"=n2))
  return(err)
}

##' Combine the fulldata results of two simulation runs 
##'
##' @title Combine simulation results
##' @param ... any number of simulation results from \code{simVBIA} with \code{fulldata} set to \code{TRUE}
##' @return a simulation result object
##' @author float
combine_simVBIA <- function(...){
    param <- sapply(list(...),tail,5)
    test <- apply(param,1,function(th) all(th[-1] %in% th[1]))
    if(!all(test)){
        stop("Trying to combine simulation results from different parameter settings")
    }
    data <- lapply(list(...),`[[`,"sim")
    c(list("sim"=rbindlist(data)),tail(list(...)[[1]],-1))
}

##' Computes summary statistics for simulation results computed with \code{simVBIA}
##' 
##' @title Summarize simulation results
##' @param simresults result of \code{simVBIA} with option fulldata set to \code{TRUE}
##' @return summary statistics over simulation results
##' @template simvbia
##' @author float
summary_simVBIA <- function(simresults){
    c(lower.prob=simresults$sim[,sum({diffmu-bar}>simresults$delta)/simresults$runs],
      upper.prob=simresults$sim[,sum({diffmu+bar}<simresults$delta)/simresults$runs],
      total.prob=simresults$sim[,sum(({diffmu-bar}>simresults$delta)+({diffmu+bar}<simresults$delta))/simresults$runs],
      mean.bias=simresults$sim[,sum(diffmu-simresults$delta)/simresults$runs],
      variance.bias=simresults$sim[,sum(sp1^2-simresults$sigma^2)/simresults$runs], #? sp^2 -1]
      vm =sim[,var(diffmu)],
      ev =simresults$sim[, mean(sp^2)], #? var(diffmu)-sem^2]
      exv =simresults$sim[,mean(simresults$sigma^2*(1/{simresults$n1+m1}+1/{simresults$n2+m2}))],
      root.mse=simresults$sim[,sqrt(sum((diffmu-simresults$delta)^2)/simresults$runs)],
      mean.m1=simresults$sim[,median(m1)],
      min.m1=simresults$sim[,unlist(quantile(m1,.01,na.rm=TRUE))],
      N10 = simresults$sim[,sum(m1==0)/simresults$runs],
### uncorrected    
      uc.lower.prob=simresults$sim[,sum({diffmuu-baru}>simresults$delta)/simresults$runs],
      uc.upper.prob=simresults$sim[,sum({diffmuu+baru}<simresults$delta)/simresults$runs],
      uc.total.prob=simresults$sim[,sum(({diffmuu-baru}>simresults$delta)+({diffmuu+baru}<simresults$delta))/simresults$runs],
      uc.mean.bias=simresults$sim[,sum(diffmuu-simresults$delta)/simresults$runs],
      uc.variance.bias=simresults$sim[,sum(sp1u^2-simresults$sigma^2)/simresults$runs], #? sp^2 -]1
      uc.vm =sim[,var(diffmuu)],
      uc.ev =simresults$sim[,mean(spu^2)],
      uc.exv =simresults$sim[, mean(simresults$sigma^2*(1/{simresults$n1+m1u}+1/{simresults$n2+m2u}))],
      uc.root.mse=simresults$sim[,sqrt(sum((diffmuu-simresults$delta)^2)/simresults$runs)],
      uc.mean.m1=simresults$sim[,median(m1u)],
      uc.min.m1=simresults$sim[,unlist(quantile(m1u,.01,na.rm=TRUE))],
      uc.N10 = simresults$sim[,sum(m1==0)/simresults$runs],
      tn1 = simresults$sim[,mean(simresults$n1)]
      )
}
    
#' Compute a lower bound for the variance bias
#'
#' \code{lowerBound} computes a lower bound for the
#' negative bias of the variance estimate following blinded
#' sample size reestimation. 
#'
#' @param n1 per-group first stage sample size
#' @param d assumed effect size
#' @param alpha Type I error probability
#' @param beta Type II error probability
#' @examples
#' lowerBound(10,.5)
#' @export
lowerBound <- function(n1,d,alpha=.025,beta=.2){
  -(2*n1 - 1)/((2*n1 - 3)*2*(qnorm(1-alpha)+qnorm(1-beta))^2/d^2)
}


#' Plot the results of one simulation run
#'
#' \code{plotSim} plots the coverage probabilities
#' of the upper, lower and two-sided confidence intervals
#' for simulation \code{seq} simulation.
#'
#' @param seq a parameter sequence plotted on the x-axis
#' @param res the results a sequence of simulation runs along parameter \code{seq}
#' @param ... further arguments passed to plot
#' @examples
#' meandif <- seq(0.1,1,.1)
#' simulation <- sapply(meandif,function(m) simVBIA(delta=.25,sigma=2,m,n1=10,n2min=2,runs=1000))
#' plotSim(meandif,simulation)
#' @export
plotSim <- function(seq,res,...){
  lb <- min(c(1,1,1/2)*res[1:3,])
  ub <- max(c(1,1,1/2)*res[1:3,])
  plot(seq,unlist(res['upper.prob',]),type='b',pch='u',col='gray',ylim=c(lb,ub),...)
  lines(seq,unlist(res['lower.prob',]),type='b',pch='l',col='gray')
  lines(seq,unlist(res['total.prob',])/2,type='b')
  abline(h=.025)
}

#' Expected bias of the mean differrence conditional on the observed first stage variance
#' 
#' \code{cond.bias} computes the expected bias of the
#' mean difference conditional on the interim variance
#' estimate, as shown in Theorem 2 of \cite{posch2014estimation}. 
#'
#' @param S1os Interim variance estimate
#' @param n1 first stage per group sample size
#' @param delta true difference of means
#' @param sigma true standard deviation
#' @return expected minus true mean-difference
#' @export
cond.bias <- function(S1os,n1,delta,sigma){
    F1 <- function(x,S1os,n1,delta,sigma)
        dnorm(x,delta,sqrt({2*sigma^2}/n1)) * dchisq((S1os*(2*n1-1) - {x^2*n1/2})/(sigma^2),2*n1-2)
    b <- sqrt(2*S1os*(2*n1 - 1)/n1)
    K <- try(integrate(F1,-b,b,S1os=S1os,n1=n1,delta=delta,sigma=sigma,rel.tol=.Machine$double.eps,stop.on.error=FALSE)$value)
    F2 <- function(x,S1os,n1,delta,sigma)
        x*F1(x,S1os=S1os,n1=n1,delta=delta,sigma=sigma)
    Em <- try(integrate(F2,-b,b,S1os=S1os,n1=n1,delta=delta,sigma=sigma,rel.tol=.Machine$double.eps,stop.on.error=FALSE)$value)
    Em <- Em/K
    Em-delta
}

#' Maximum bias of the estimate mean difference following blinded interim analyses 
#'
#' \code{simMBIA} computes, via Monte-Carlo simulation,
#' the bias for estimates of the mean difference
#' if the blinded sample size reassessment rule is chosen
#' such that the bias is maximised.
#'
#' @param delta true mean difference
#' @param n1 first stage per group sample size
#' @param sigma true standard deviation
#' @param runs number of simulted trials
#' @return A vector following components:
#' \item{m.bias}{Maximum positiv mean bias that can be attained for a given parameter setting}
#' \item{m.bias.n}{Maximum negative mean bias that can be attained}
#' @export
simMBIA <- function(delta=0,n1=2,sigma=1,runs=100)
{
  ## means, variances group 1 and 2
  ma=rnorm(runs,0,sigma/sqrt(n1))
  mb=rnorm(runs,delta,sigma/sqrt(n1))
  sa=sigma^2*rchisq(runs,n1-1)
  sb=sigma^2*rchisq(runs,n1-1)
  ## mean differnce
  md = mb-ma
  ## interim variance estimate
  S1os=(sa+sb+(n1/2)*md^2)/(2*n1-1)
  ## expected bias
  biasv=simplify2array(lapply(S1os,cond.bias,n1=n1,delta=delta,sigma=sigma))
  
  ## if the expected bias is positive stop otherwise make second stage infinitely large
  est1=ifelse(biasv>0,md-delta,0)
  ## otherway around
  est2=ifelse(biasv<0,md-delta,0)
  ## stop  always
  est3=md-delta
  c(m.bias = mean(est1),
    m.bias.n = mean(est2))
}

#compMBIA <- function(delta=0,n1=2,sigma=1,runs=100)
    

