\name{Daphnia}

\alias{Daphnia}

\docType{data}

\title{Daphnia}

\description{
Data are from a binary mixture experiment that was based on a fixed-ratio design involving 5 rays: the 2 rays for the pesticides prochloraz and alpha-cypermethrin and 3 mixture rays corresponding to virtual mixture proportions of 25:75, 50:50, and 75:25. 
}

\usage{data(Daphnia)}

\format{
  A data frame with 140 observations on the following 6 variables.
  \describe{
    \item{\code{dose.a}}{Dose of alpha-cypermethrin (mu g/L)}
    \item{\code{dose.p}}{Dose of prochloraz (mu g/L)}
    \item{\code{dose}}{Total dose in the mixture (mu g/L)}
    \item{\code{mix.frac}}{Mixture fraction}
    \item{\code{total}}{Total number of Daphnia}
    \item{\code{immob48}}{Number of immobile Daphnia after 48 hours}
  }
}

\details{
Synergistic and antagonistic effects of binary mixtures between a number of fungicides and the pyrethroid insecticide alpha-cypermethrin were investigated using a standard test system. Only data for the specific binary mixture of prochloraz and alpha-cypermethrin are provided. Data were obtained from a Daphnia acute immobilisation test where the test organisms were divided into groups of five, placed in containers, exposed to a dose (either a mixture dose or a dose from one of the two pesticides), and followed for 48h. 
}

\source{
Data were kindly provided by N. Cedergreen.
}

\references{
Noergaard KB and Cedergreen N, Pesticide cocktails can interact synergistically on aquatic crustaceans. Environ Sci Pollut Res 17: 957-967 (2010). https://doi.org/10.1007/s11356-009-0284-4
  }

\examples{}
\keyword{datasets}
