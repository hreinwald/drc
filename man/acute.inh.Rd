\name{acute.inh}

\alias{acute.inh}

\docType{data}

\title{Acute inhalation}

\description{Data from an acute inhalation toxicity test. For each of several dose levels, the total number of subjects and the number of dead subjects were recorded.}

\usage{data(acute.inh)}

\format{
  A data frame with 6 observations on the following 3 variables.
  \describe{
    \item{\code{dose}}{a numeric vector}
    \item{\code{total}}{a numeric vector}
    \item{\code{num.dead}}{a numeric vector}
  }
}

\examples{
\dontrun{
library(drc)

## Displaying the data
head(acute.inh)

## Fitting a two-parameter log-logistic model for binomial response
acute.inh.m1 <- drm(num.dead/total ~ dose, weights = total,
data = acute.inh, fct = LL.2(), type = "binomial")
summary(acute.inh.m1)

## Plotting the fitted curve
plot(acute.inh.m1, xlab = "Dose", ylab = "Proportion dead")
}
}


\keyword{datasets}
