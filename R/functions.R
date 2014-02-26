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
zss <- function(s,d,alpha,beta){
  2 * s^2 / d^2 * (qnorm(alpha,lower=FALSE)+qnorm(beta,lower=FALSE))^2
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
zsd <- function(n1,d,v = 1/2,alpha,beta){
    n <- n1/(v*2)
    sqrt(n * d^2 / ((qnorm(alpha,lower=FALSE)+qnorm(beta,lower=FALSE))^2))
}


#' Simulate estimation after blinded interim analyses
#'
#' \code{simBIA} simulates \code{runs} adaptive trials with
#' blinded interim analysis and assesses the bias of the variance 
#' and mean estimates. Sample size review is performed either using
#' an interim estimate of the pooled variance that is adjusted
#' based on the assumed effect size or the uncorrected lumped variance.
#'
#' @param delta true effect size
#' @param sigma true standard deviation
#' @param d assumed effect size
#' @param alpha = .025 pre-specified alpha level
#' @param beta pre-specified beta level
#' @param v timing of interim analysis
#' @param s assumed standard deviation, is more or less the same as n1
#' @param n1 first stage sample size in treatment group
#' @param runs number of simulation runs
#' @param n2 first stage sample size in control group
#' @param n2min minimum secondstage sample size
#' @param testdata debugging option not for general use
#' @param fulldata return full data 
#' @param cf correction term added to the sample size rule
#' @param ...
#' @examples
#' simBIA(.5,1,1)
#' @export
simBIA <- function(delta, # true effect size
                   sigma, # true standard deviation
                   d, # assumed effect size
                   alpha = .025, # pre-specified alpha level
                   beta = .2, # pre-specified beta level
                   v = 1/2, # timing of interim analysis
                   s = zsd(n1,d,v,alpha,beta), # assumed standard deviation, is more or less the same as n1
                   n1 = ceiling(v * zss(s,d,alpha,beta)), # first stage sample size in treatment group
                   runs, # number of simulation runs
                       n2 =n1, # first stage sample size in control group
                   n2min = 0, # minimum secondstage sample size
                   testdata = FALSE, #  
                   fulldata = FALSE, #  
                   wurzel = FALSE, #  
                   cf = 0, #  
                   ...){
  ## pre-compute t-distribution quantiles
  qttable=qt(1-alpha,1:1000)
  xi=(qnorm(alpha,lower=FALSE)+qnorm(beta,lower=FALSE))^2/d^2
  n=n1+n2 # total first stage sample size
  lamb=n1*n2/n ##???
  # Simulation der first stage means und sum of squares
  x1=rnorm(runs,s=sigma/sqrt(n1))
  x1r=rnorm(runs,s=sigma/sqrt(n1))
  
  x2=rnorm(runs,m=delta,s=sigma/sqrt(n2))
  x2r=rnorm(runs,m=delta,s=sigma/sqrt(n2))

  sq1=sigma^2*rchisq(runs,n1-1)
  sq1r=sigma^2*rchisq(runs,n1-1)
  sq2=sigma^2*rchisq(runs,n2-1)
  sq2r=sigma^2*rchisq(runs,n2-1)

  ## replace by external values
  if(testdata & exists('florian')){
    x1 <- florian$x1
    x2 <- florian$x2
    sq1 <- florian$sq1
    sq2 <- florian$sq2
  }
  diff=x1-x2
  diffr=x1r-x2r
  # berechnung der overall standard deviation s_adj nach Kieser &Friede 2003
  # may become <0, we set it to zero in that case (?) check
  ps = pmax((sq1+sq2+lamb*diff*diff)/(n-1)- d^2*n1/(2*n-2),0) #check
  psu = (sq1+sq2+lamb*diff*diff)/(n-1)

  ## round of errors!
  #m1=pmax(round(2*xi*ps)-n1,0)
  m1 <- pmax(ceiling(2*xi*ps)-n1+cf,n2min)
  m1u <- pmax(ceiling(2*xi*psu)-n1+cf,n2min)

  ## replace by external values
  if(testdata & exists('florian')){
    other <- which(m1 != florian$m1)
    cat("Mismatching sample sizes: ",length(other),"\n")
    m1 <- florian$m1
  }
  m2=m1
  m2u=m1u
  #Berechnung der secondstage means
  y1=rnorm(runs,s=sigma/sqrt(m1))
  y2=rnorm(runs,s=sigma/sqrt(m2),m=delta)
  y1u=rnorm(runs,s=sigma/sqrt(m1u))
  y2u=rnorm(runs,s=sigma/sqrt(m2u),m=delta)
  
  # wir setzen die means auf null falls die second stage fallzahl 0 ist
  y1[m1<1]=0 
  y2[m2<1]=0
  y1u[m1u<1]=0 
  y2u[m2u<1]=0
  
  # Berechnung der sum of squares second stage
  tq1=sigma^2*rchisq(runs,df=pmax(m1-1,1)) 
  tq2=sigma^2*rchisq(runs,df=pmax(m2-1,1))
  tq1u=sigma^2*rchisq(runs,df=pmax(m1u-1,1)) 
  tq2u=sigma^2*rchisq(runs,df=pmax(m2u-1,1))

  # wir setzen diese  auf null falls die second stage fallzahl <2 ist
  tq1[m1<2]=0
  tq2[m2<2]=0
  tq1u[m1u<2]=0
  tq2u[m2u<2]=0

  ## replace by external values
  if(testdata & exists('florian')){
    y1 <- florian$y1
    y2 <- florian$y2
    tq1 <- florian$tq1
    tq2 <- florian$tq2
    y1u <- florian$y1
    y2u <- florian$y2
    tq1u <- florian$tq1
    tq2u <- florian$tq2    
  }
     
  
  # Berechnung der Overall means ond varianzen und konfidence bounds
  # hier muesste man checken, ob die formeln auch für kleine m1,m2, das heißt m1,m2=0,1 stimmen
  l1=n1+m1
  l2=n2+m2
  ll=l1+l2-2
  l1u=n1+m1u
  l2u=n2+m2u
  llu=l1u+l2u-2
  # siehe kieser & friede 2003
  ## t-test denominator
  sp1=sqrt((tq1+sq1+(y1-x1)^2*n1*m1/l1+tq2+sq2+(y2-x2)^2*n2*m2/l2)/ll)
  sp1r=sqrt((tq1+sq1r+(y1-x1r)^2*n1*m1/l1+tq2+sq2r+(y2-x2r)^2*n2*m2/l2)/ll)
  sp=sp1*sqrt((1/l1+1/l2))

  sp1u=sqrt((tq1u+sq1+(y1u-x1)^2*n1*m1u/l1u+tq2u+sq2+(y2u-x2)^2*n2*m2u/l2u)/llu)
  sp1ur=sqrt((tq1u+sq1r+(y1u-x1r)^2*n1*m1u/l1u+tq2u+sq2r+(y2u-x2r)^2*n2*m2u/l2u)/llu)
  spu=sp1u*sqrt((1/l1u+1/l2u))

  ## t-test numerator
  diffmu=(n2*x2+m2*y2)/l2-(n1*x1+m1*y1)/l1
  bar=sp*qttable[pmin(ll,1000)]# hier werden die approx. quantile verwendet. hat das einen impact
  cl=diffmu-bar
  cu=diffmu+bar
  
  diffmuu=(n2*x2+m2u*y2u)/l2u-(n1*x1+m1u*y1u)/l1u
  baru=spu*qttable[pmin(llu,1000)]# hier werden die approx. quantile verwendet. hat das einen impact
  clu=diffmuu-baru
  cuu=diffmuu+baru

  err=c(
      lower.prob=sum(cl>delta)/runs,
      upper.prob=sum(cu<delta)/runs,
      total.prob=sum((cl>delta)+(cu<delta))/runs,
      mean.bias=sum(diffmu-delta)/runs,
      variance.bias=sum(sp1^2-sigma^2)/runs, #? sp^2 -1
      variance.ubias=sum(sp1r^2-sigma^2)/runs, #? sp^2 -1
      se.bias=(sd(diffmu)-sqrt(mean(sp^2))), #? 
      root.mse=sqrt(sum((diffmu-delta)^2)/runs),
      mean.m1=median(m1),
      min.m1=unlist(quantile(m1,.01,na.rm=TRUE)),
      n10 = sum(m1==0)/runs,
### uncorrected    
      uc.lower.prob=sum(clu>delta)/runs,
      uc.upper.prob=sum(cuu<delta)/runs,
      uc.total.prob=sum((clu>delta)+(cuu<delta))/runs,
      uc.mean.bias=sum(diffmuu-delta)/runs,
      uc.variance.bias=sum(sp1u^2-sigma^2)/runs, #? sp^2 -1
      uc.variance.ubias=sum(sp1ur^2-sigma^2)/runs, #? sp^2 -1
      uc.se.bias=(sd(diffmuu)-sqrt(mean(spu^2))), #? 
      uc.root.mse=sqrt(sum((diffmuu-delta)^2)/runs),
      uc.mean.m1=median(m1u),
      uc.min.m1=unlist(quantile(m1u,.01,na.rm=TRUE)),
      uc.n10 = sum(m1==0)/runs,
      tn1 = mean(n1)
      )
  # list of lower prob, upper prob, total prob, mean bias, variance bias, standard error bias, sqrt(mse)
  if(testdata & exists("florian"))
    print(cbind(err,florian$err))
  #list(err=err,m1=m1)
  data <- list("x2" = x2, "y2" = y2,"sq2" = sq2,"tq2" = tq2, "x1" = x1, "y1" = y1,
               "sq1" = sq1, "tq1" = tq1, "m1" = m1, "ps" = ps,
               "y2u" = y2u ,"tq2u" = tq2u, "y1u" = y1u,
               "tq1u" = tq1u, "m1u" = m1u, "psu" = psu)
  data$err <- err
  if(fulldata) return(data)
  return(err)
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
#' @examples
#' meandif <- seq(0,1,.1)
#' simulation <- sapply(meandif,function(m) simBIA(delta=.25,sigma=2,m,runs=1000))
#' plotSim(meandif,simulation)
plotSim <- function(seq,res,...){
  lb <- min(c(1,1,1/2)*res[1:3,])
  ub <- max(c(1,1,1/2)*res[1:3,])
  plot(seq,unlist(res['upper.prob',]),type='b',pch='u',col='gray',ylim=c(lb,ub),...)
  lines(seq,unlist(res['lower.prob',]),type='b',pch='l',col='gray')
  lines(seq,unlist(res['total.prob',])/2,type='b')
  abline(h=.025)
}
