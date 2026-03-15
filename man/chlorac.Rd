\name{chlorac}

\alias{chlorac}

\docType{data}

\title{chlorac}

\description{Data from an acute toxicity test where organisms were exposed to different concentrations of chloroacetaldehyde. The number of dead subjects out of a total were recorded for each concentration.}

\usage{data(chlorac)}

\format{
  A data frame with 6 observations on the following 3 variables.
  \describe{
    \item{\code{conc}}{a numeric vector}
    \item{\code{total}}{a numeric vector}
    \item{\code{num.dead}}{a numeric vector}
  }
}

\examples{
library(drc)

## Displaying the data
head(chlorac)

## Fitting a two-parameter log-logistic model for binomial response
chlorac.m1 <- drm(num.dead/total ~ conc, weights = total,
data = chlorac, fct = LL.2(), type = "binomial")
summary(chlorac.m1)

## Plotting the fitted curve
plot(chlorac.m1, xlab = "Concentration", ylab = "Proportion dead")
}

\keyword{datasets}
