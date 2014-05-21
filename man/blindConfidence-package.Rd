\docType{package}
\name{blindConfidence-package}
\alias{blindConfidence}
\alias{blindConfidence-package}
\title{Simulation of estimation after blinded interim analysis}
\description{
Simulation of estimation after blinded interim analysis
}
\details{
This package provides all code used to perform and
summarize the simulation study shown in our paper
"Estimation after blinded interim analysis". The code for
the general simulations is provided in
\code{vignette{'SimulateBias'}}, the corresponding
simulation results are saved in \code{data{gridsim}}; the
code for the case study simulation is in
\code{vignette{'CaseStudy'}} teh corresponding simulation
results are saved in \code{data{casesim}}. Rebuilding the
vignettes will use the saved results to redraw all the
graphics from the paper. In order to rerun the simulations
the corresponding code chunks need to be configured for
evaluation. Note that although this package will run on
Windows recomputing the simulation results requires
additional packages (e.g. parallel) that are not supported
on Windows.
}
\author{
Florian Klinglmueller
}

