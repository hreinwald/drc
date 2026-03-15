\name{C.dubia}

\alias{C.dubia}

\docType{data}

\title{Offsprings resulting from a toxicity test}

\description{
Results from a chronic reproduction toxicity test with seven different concentrations of waste water. The response was the number of offspring produced by the water flea \emph{Ceriodaphnia dubia}.
}

\usage{data(C.dubia)}

\format{
  A data frame with 50 observations of the following 2 variables.
  \describe{
    \item{\code{conc}}{a numeric vector giving waste water in percentage}
    \item{\code{number}}{a numeric vector}
  }
}

\source{
A. J. Bailer and J. T. Oris (1997). Estimating inhibition concentrations for different response scales
using generalized linear models. Environmental Toxicology and Chemistry, \bold{16}:1554--1559.
}


\examples{
library(drc)

## Displaying the data
head(C.dubia)

## Fitting a three-parameter log-logistic model
C.dubia.m1 <- drm(number ~ conc, data = C.dubia, fct = LL.3())
summary(C.dubia.m1)

## Plotting fitted curve together with the original data
plot(C.dubia.m1, xlab = "Concentration (\%)", ylab = "Number of offspring")
}
\keyword{datasets}
