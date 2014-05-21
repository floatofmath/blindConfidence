\name{cond.bias}
\alias{cond.bias}
\title{Expected bias of the mean differrence conditional on the observed first stage variance}
\usage{
cond.bias(S1os, n1, delta, sigma)
}
\arguments{
  \item{S1os}{Interim variance estimate}

  \item{n1}{first stage per group sample size}

  \item{delta}{true difference of means}

  \item{sigma}{true standard deviation}
}
\value{
expected minus true mean-difference
}
\description{
\code{cond.bias} computes the expected bias of the mean
difference conditional on the interim variance estimate, as
shown in Theorem 2 of \cite{posch2014estimation}.
}

