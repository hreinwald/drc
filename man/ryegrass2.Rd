\name{ryegrass2}

\alias{ryegrass2}

\docType{data}

\title{Ryegrass}

\description{Data from a dose-response experiment with ryegrass (\emph{Lolium} sp.). Biomass was measured at different dose levels and at two time points.}

\usage{data(ryegrass2)}

\format{
  A data frame with 27 observations on the following 3 variables.
  \describe{
    \item{\code{dose}}{a numeric vector}
    \item{\code{biomass}}{a numeric vector}
    \item{\code{day}}{a categorial vector}
  }
}

\examples{
library(drc)

## Displaying the data
head(ryegrass2)

## Fitting a four-parameter log-logistic model with separate curves per day
ryegrass2.m1 <- drm(biomass ~ dose, day, data = ryegrass2, fct = LL.4())
summary(ryegrass2.m1)

## Plotting the fitted curves
plot(ryegrass2.m1, xlab = "Dose", ylab = "Biomass")
}
\keyword{datasets}
