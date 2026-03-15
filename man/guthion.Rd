\name{guthion}

\alias{guthion}

\docType{data}

\title{guthion}

\description{Data from an acute toxicity test with the insecticide guthion (azinphos-methyl). For each dose level in two treatment groups, the numbers of alive, moribund, and dead subjects were recorded.}

\usage{data(guthion)}

\format{
  A data frame with 6 observations on the following 6 variables.
  \describe{
    \item{\code{trt}}{a categorial vector}
    \item{\code{dose}}{a numeric vector}
    \item{\code{alive}}{a numeric vector}
    \item{\code{moribund}}{a numeric vector}
    \item{\code{dead}}{a numeric vector}
    \item{\code{total}}{a numeric vector}
  }
}

\examples{
library(drc)

## Displaying the data
head(guthion)

## Fitting a two-parameter log-logistic model for binomial response
guthion.m1 <- drm(dead/total ~ dose, trt, weights = total,
data = guthion, fct = LL.2(), type = "binomial")
summary(guthion.m1)

## Plotting the fitted curves
plot(guthion.m1, xlab = "Dose", ylab = "Proportion dead")
}
\keyword{datasets}
