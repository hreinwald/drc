\name{nicotine}

\alias{nicotine}

\docType{data}

\title{nicotine}

\description{Data from an acute toxicity test with nicotine. For each of several concentrations, the total number of subjects and the number of dead subjects were recorded.}

\usage{data(nicotine)}

\format{
  A data frame with 12 observations on the following 3 variables.
  \describe{
    \item{\code{conc}}{a numeric vector}
    \item{\code{total}}{a numeric vector}
    \item{\code{num.dead}}{a numeric vector}
  }
}

\examples{
library(drc)

## Displaying the data
head(nicotine)

## Fitting a two-parameter log-logistic model for binomial response
nicotine.m1 <- drm(num.dead/total ~ conc, weights = total,
data = nicotine, fct = LL.2(), type = "binomial")
summary(nicotine.m1)

## Plotting the fitted curve
plot(nicotine.m1, xlab = "Concentration", ylab = "Proportion dead", ylim = c(0, 1))
}
\keyword{datasets}
