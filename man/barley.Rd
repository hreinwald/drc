\name{barley}

\alias{barley}

\docType{data}

\title{Barley}

\description{Data from a dose-response experiment measuring the weight of barley (\emph{Hordeum vulgare}) at different dose levels of a substance.}

\usage{data(barley)}

\format{
  A data frame with 18 observations of the following 2 variables.
  \describe{
    \item{\code{Dose}}{a numeric vector}
    \item{\code{weight}}{a numeric vector}
    }
}

\examples{
\dontrun{
library(drc)

## Displaying the data
head(barley)

## Fitting a four-parameter log-logistic model
barley.m1 <- drm(weight ~ Dose, data = barley, fct = LL.4())
summary(barley.m1)

## Plotting the fitted curve
plot(barley.m1, xlab = "Dose", ylab = "Weight")
}
}

\keyword{datasets}
