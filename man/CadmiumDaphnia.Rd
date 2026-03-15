\name{CadmiumDaphnia}

\alias{CadmiumDaphnia}

\docType{data}

\title{Cadmium Daphnia Data}

\description{
  Data from an acute toxicity test exposing \emph{Daphnia} to cadmium over time. The endpoint measured was mortality (number of dead organisms) at each dose and time point.
}

\usage{data(CadmiumDaphnia)}

\format{
  A data frame with 58 observations on the following 6 variables.
  \describe{
    \item{\code{Dose}}{a numeric vector of dose values}
    \item{\code{Time}}{a numeric vector}
    \item{\code{Total}}{a numeric vector}
    \item{\code{Start}}{a numeric vector}
    \item{\code{End}}{a numeric vector}
    \item{\code{Dead}}{a numeric vector}
  }
}

\examples{
\dontrun{
library(drc)

## Displaying the data
head(CadmiumDaphnia)

## Fitting a two-parameter log-logistic model for binomial response at a single time point
CadmiumDaphnia.sub <- CadmiumDaphnia[CadmiumDaphnia$Time == 7, ]
CadmiumDaphnia.m1 <- drm(Dead/Total ~ as.numeric(as.character(Dose)), weights = Total,
data = CadmiumDaphnia.sub, fct = LL.2(), type = "binomial")
summary(CadmiumDaphnia.m1)

## Plotting the fitted curve
plot(CadmiumDaphnia.m1, xlab = "Cadmium dose", ylab = "Proportion dead")
}
}
\keyword{datasets}
