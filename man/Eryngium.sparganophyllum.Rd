\name{Eryngium.sparganophyllum}
\alias{Eryngium.sparganophyllum}
\alias{Eryngium.sparganophyllum0}

\docType{data}

\title{Germination of Eryngium sparganophyllum}

\description{Germination data from an experiments investigating the effect of different concentration of gibberellic acid on germination of Eryngium sparganophyllum seeds. Two datasets are provided: one resembling how data are entered in the first place ("Eryngium.sparganophyllum0") and one formatted and ready-to-use for the statistical analysis ("Eryngium.sparganophyllum")}

\usage{Eryngium.sparganophyllum}

\format{
  A data frame with 583 observations on the following variables.
  \describe{
    \item{\code{Treat}}{a factor with 15 levels denoting the concentration of gibberellic acid (in ppm)}
    \item{\code{Type}}{a factor with two levels denoting the type of treatment (gibberellic acid or temperature)}
    \item{\code{Day}}{a numeric vector recording time (in days) since the beginning of the experiment}
    \item{\code{Germ}}{a numeric vector of counts of germinated seeds}
    \item{\code{Start}}{a numeric vector of starting time points of monitoring intervals}
    \item{\code{End}}{a numeric vector of ending time points of monitoring intervals}
    \item{\code{Germinated}}{a numeric vector of counts of germinated seeds in a given interval}
    \item{\code{Rep}}{a numeric vector corresponding to the replicated sub-experiments; it is only a unique enumeration for the dataset "Eryngium.sparganophyllum"}
  }
}

\references{
Wolkis, D., Blackwell, S., Kaninaualiʻi Villanueva, S. (2020). Conservation seed physiology of the ciénega endemic, Eryngium sparganophyllum (Apiaceae). Conservation Physiology, 8, coaa017. https://doi.org/10.1093/conphys/coaa017
}

\keyword{datasets}
