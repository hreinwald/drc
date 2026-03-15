\name{fluoranthene}

\alias{fluoranthene}

\docType{data}

\title{Death of fathead minnow larvae after exposure to fluoranthene}

\description{Fathead minnow larvae were exposed to fluoranthene, a polycyclic aromatic hydrocarbon, under two different algal densities resulting in different levels of ambient ultraviolet radiation. Number of dead larvaes were reported.
}

\usage{data(fluoranthene)}

\format{
  A data frame with 24 observations on the following 4 variables.
  \describe{
    \item{\code{algalconc}}{a numeric vector}
    \item{\code{conc}}{a numeric vector}
    \item{\code{totalnum}}{a numeric vector}
    \item{\code{mortality}}{a numeric vector}
  }
}

%\details{}

\source{
M. W. Wheeler, R. M. Park, and A. J. Bailer (2006). Comparing median lethal concentration values using confidence interval overlap or ratio tests. Environmental Toxicology and Chemistry, \bold{25}:1441--1444.
}

%\references{}

\examples{
library(drc)

## Displaying the data
head(fluoranthene)

## Fitting a two-parameter log-logistic model for binomial response
## with different curves per algal concentration
fluoranthene.m1 <- drm(mortality/totalnum ~ conc, algalconc, weights = totalnum,
data = fluoranthene, fct = LL.2(), type = "binomial")
summary(fluoranthene.m1)

## Plotting the fitted curves
plot(fluoranthene.m1, xlab = "Fluoranthene concentration",
ylab = "Proportion dead")
}


\keyword{datasets}
