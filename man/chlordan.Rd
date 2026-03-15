\name{chlordan}

\alias{chlordan}

\docType{data}

\title{Chlordan}

\description{Data from a chronic toxicity test measuring the reproduction of \emph{Daphnia} exposed to different concentrations of chlordane at two time points. The response measured was the number of offspring (repro) per replicate.}

\usage{data(chlordan)}

\format{
  A data frame with 60 observations on the following 5 variables.
  \describe{
    \item{\code{replicate}}{a numeric vector}
    \item{\code{conc}}{a numeric vector}
    \item{\code{repro}}{a numeric vector}
    \item{\code{time}}{a numeric vector}
  }
}

\examples{
\dontrun{
library(drc)

## Displaying the data
head(chlordan)

## Fitting a three-parameter log-logistic model
chlordan.m1 <- drm(repro ~ conc, data = chlordan, fct = LL.3())
summary(chlordan.m1)

## Plotting the fitted curve
plot(chlordan.m1, xlab = "Concentration", ylab = "Reproduction")
}
}
\keyword{datasets}
