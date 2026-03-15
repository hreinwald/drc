\name{arbovirus}

\alias{arbovirus}

\docType{data}

\title{arbovirus}

\description{Data from a dose-response experiment with an arbovirus involving two treatment groups (FP and SP). For each dose level, the total number of subjects and the numbers of dead and defective subjects were recorded.}

\usage{data(arbovirus)}

\format{
  A data frame with 9 observations on the following 5 variables.
  \describe{
    \item{\code{dose}}{a numeric vector}
    \item{\code{total}}{a numeric vector}
    \item{\code{dead}}{a numeric vector}
    \item{\code{def}}{a numeric vector}
    \item{\code{trt}}{a categorical vector}
  }
}

\examples{
library(drc)

## Displaying the data
head(arbovirus)

## Fitting a two-parameter log-logistic model for binomial response
arbovirus.m1 <- drm(dead/total ~ dose, trt, weights = total,
data = arbovirus, fct = LL.2(), type = "binomial")
summary(arbovirus.m1)

## Plotting the fitted curves
plot(arbovirus.m1, xlab = "Dose", ylab = "Proportion dead")
}
\keyword{datasets}
