\name{lemna}

\alias{lemna}

\docType{data}

\title{Lemna}

\description{Data from a dose-response experiment with the aquatic plant \emph{Lemna minor} (duckweed). The response measured was the frond number (count) at different concentrations of a test substance.}

\usage{data(lemna)}

\format{
  A data frame with 44 observations on the following 2 variables.
  \describe{
    \item{\code{conc}}{a numeric vector containing the concentration.}
    \item{\code{frond.num}}{a numeric vector containing the response (count).}
    }
}

\examples{
library(drc)

## Displaying the data
head(lemna)

## Fitting a four-parameter log-logistic model
lemna.m1 <- drm(frond.num ~ conc, data = lemna, fct = LL.4())
summary(lemna.m1)

## Plotting the fitted curve
plot(lemna.m1, xlab = "Concentration", ylab = "Frond number")
}
\keyword{datasets}
