\name{red.fescue}

\alias{red.fescue}

\docType{data}

\title{Red fescue}

\description{Data from a dose-response experiment with red fescue (\emph{Festuca rubra}). Biomass was measured at different dose levels and at two time points (day 0 and day 16).}

\usage{data(red.fescue)}

\format{
  A data frame with 26 observations on the following 3 variables.
  \describe{
    \item{\code{dose}}{a numeric vector}
    \item{\code{day}}{a numeric vector}
    \item{\code{biomass}}{a numeric vector}
  }
}

\examples{
\dontrun{
library(drc)

## Displaying the data
head(red.fescue)

## Fitting a four-parameter log-logistic model with separate curves per day
red.fescue.m1 <- drm(biomass ~ dose, day, data = red.fescue, fct = LL.4())
summary(red.fescue.m1)

## Plotting the fitted curves
plot(red.fescue.m1, xlab = "Dose", ylab = "Biomass")
}
}
\keyword{datasets}
