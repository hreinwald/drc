\name{carbendazim}

\alias{carbendazim}

\docType{data}

\title{Damage of lymphocyte cells}

\description{
  For each of 13 dose levels the number of damaged lymphocyte cells were reported. Each dose level consisted of a total of 2000 cells.
}

\usage{data(carbendazim)}

\format{
  A data frame with 13 observations of the following 3 variables.
  \describe{
    \item{\code{dose}}{a numeric vector}
    \item{\code{total}}{a numeric vector}
    \item{\code{damage}}{a numeric vector}
  }
}


\source{
Bentley, K. S. and Kirkland, D. and Murphy, M. and Marshall, R. (2000). Evaluation of   thresholds for benomyl- and carbendazim-induced aneuploidy in cultured human lymphocytes using fluorescence in situ hybridization, \emph{Mutation Research/Genetic Toxicology and Environmental Mutagenesis}, \bold{464}, 41--51.}

\examples{
\dontrun{
library(drc)

## Displaying the data
head(carbendazim)

## Fitting a two-parameter log-logistic model for binomial response
carbendazim.m1 <- drm(damage/total ~ dose, weights = total,
data = carbendazim, fct = LL.2(), type = "binomial")
summary(carbendazim.m1)

## Plotting the fitted curve
plot(carbendazim.m1, xlab = "Dose", ylab = "Proportion damaged")
}
}

\keyword{datasets}
