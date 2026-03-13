\name{bees}

\alias{bees}

\docType{data}

\title{bees}

\description{
  Data are from a binary mixture experiment that involves multiple single-dose factorial designs where the insecticide imidacloprid is combined with each of 7 pesticides in turn. 
}

\usage{data(bees)}

\format{
  A data frame with 66 observations on the following 5 variables.
  \describe{
    \item{\code{mixture}}{Indicator of single-dose experiment (or control)}
    \item{\code{treat}}{Treatment or combination of treatments}
    \item{\code{rep}}{Replication number (there were 3 replicates per group)}
    \item{\code{dead0h}}{Number of dead bees initially}
    \item{\code{dead48h}}{Number of dead bees after 48 hours}
  }
}

\details{
Imidacloprid is a widely used insecticide. In a recent study potential synergistic effects on mortality of honey bees exposed to the insectide in binary mixtures with seven pesticides from different classes: acephate, λ-cyhalothrin, oxamyl, tetraconazole, sulfoxaflor, glyphosate, and clothianidin were investigated. Bees were reared in cages (25 insects per cage), with three cages per treatment group, and exposed to the mixture treatments for 48h. Mortality after 48h was the response. 
}
\source{
Data were retrieved from PLoS ONE repository.
}

\references{
Zhu YC, Yao J, Adamczyk J and Luttrell R, Synergistic toxicity and physiological impact of imidacloprid alone and binary mixtures with seven representative pesticides on honey bee (Apis mellifera). PLoS ONE 12: e0176837 (2017). https://doi.org/10.1371/journal.pone.0176837
  }

\examples{}
\keyword{datasets}
