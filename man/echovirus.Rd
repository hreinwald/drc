\name{echovirus}

\alias{echovirus}

\docType{data}

\title{Infections as response to exposure with \emph{Echovirus 12}}

\description{
  For each of four doses of a pathogen, \emph{Echovirus 12}, the number of exposed and infected human volunteers are reported.
}

\usage{data(echovirus)}

\format{
  A data frame with 4 observations on the following 3 variables.
  \describe{
    \item{\code{dose}}{a numeric vector reporting the dose in plague forming units (pfu)}
    \item{\code{total}}{a numeric vector}
    \item{\code{infected}}{a numeric vector}
  }
}

%\details{}

\source{
   H. Moon, S. B. Kim, J. J. Chen, N. I. George, and R. L. Kodell (2013). Model uncertainty
and model averaging in the estimation of infectious doses for microbial pathogens. Risk
Analysis, \bold{33(2)}:220-231.
}

%\references{}

\examples{
library(drc)

## Displaying the data
head(echovirus)

## Fitting a two-parameter log-logistic model for binomial response
echovirus.m1 <- drm(infected/total ~ dose, weights = total,
data = echovirus, fct = LL.2(), type = "binomial")
summary(echovirus.m1)

## Plotting the fitted curve
plot(echovirus.m1, xlab = "Dose (pfu)", ylab = "Proportion infected")
}
\keyword{datasets}
