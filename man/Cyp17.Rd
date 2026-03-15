\name{Cyp17}
\alias{Cyp17}
\docType{data}
\title{Cyp17 expression data}
\description{Observed Cyp17 gene expression measured at several dose levels across multiple experimental runs. CYP17 is a key enzyme in steroid hormone biosynthesis, and changes in its expression can indicate endocrine-disrupting effects.}
\usage{data(Cyp17)}
\format{
  A data frame with 63 observations on the following 3 variables.
  \describe{
    \item{\code{run}}{ID of 3 different runs}
    \item{\code{dose}}{5 dose levels (0, 0.1, 10, 100, 500)}
    \item{\code{expression}}{observed expression}
  }
}


\examples{
data(Cyp17)

## Display the structure of the data
head(Cyp17)

## Log-transform the expression values
Cyp17$logexpression <- log(Cyp17$expression) + 5

## Fit a four-parameter log-logistic model (ignoring run effects)
Cyp17.m1 <- drm(logexpression ~ dose, data = Cyp17, fct = LL.4())
summary(Cyp17.m1)
plot(Cyp17.m1, main = "Cyp17 dose-response")
}



\keyword{datasets}
