#' @return If fulldata is \code{FALSE} simVBIA returns a list with
#' the following items, summarizing the simulation results for
#' designs with sample size reassessment based on the adjusted
#' interim variance estimates:
#' \item{lower.prob}{Coverage probability of the lower confidence bound}
#' \item{upper.prob}{Coverage probability of the upper confidence bound}
#' \item{total.prob}{Coverage probability of the two sided confidence interval}
#' \item{mean.bias}{Bias of the mean estimate}
#' \item{variance.bias}{Bias of the variance estimate}
#' \item{ev}{Variance of the mean estimate as estimated from the adaptive trian}
#' \item{exv}{Variance of the mean estimate of a fixed sample trial with corresponding sample size}
#' \item{vm}{Variance of the mean estimate in the monte carlo sample}
#' \item{root.mse}{Root mean squared error of the effect estimate}
#' \item{mean.m1}{Mean second stage sample size}
#' \item{low.m1}{Lower 10 percen quantile of second stage sample sizes}
#' \item{n10}{Probability that the second stage sample size is zero}
#' as well as results for designs with sample size reassessment
#' based on unadjusted interim estimates:
#' \item{uc.mean.bias}{Bias of the mean estimate}
#' \item{uc.variance.bias}{Bias of the variance estimate}
#' \item{uc.ev}{Variance of the mean estimate as estimated from the adaptive trian}
#' \item{uc.exv}{Variance of the mean estimate of a fixed sample trial with corresponding sample size}
#' \item{uc.vm}{Variance of the mean estimate in the monte carlo sample}
#' \item{uc.root.mse}{Root mean squared error of the effect estimate}
#' \item{uc.mean.m1}{Mean second stage sample size}
#' \item{uc.low.m1}{Lower 10 percen quantile of second stage sample sizes}
#' \item{uc.N10}{Probability that the second stage sample size is zero}
#' \item{tn1}{mean(n1)}
