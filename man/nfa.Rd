\name{nfa}
\alias{nfa}
\docType{data}
\title{Network Formation Assay Data}
\description{Neurotoxicity test using a network formation assay studying
  the inhibition of network formation at acrylamide exposure.}
\usage{data(nfa)}
\format{
  A data frame with 45 observations on the following 4 variables.
  \describe{
    \item{\code{chip}}{chip ID}
    \item{\code{conc}}{7 concentrations of acrylamide, ranging from
      0-5mM}
    \item{\code{experiment}}{factor with levels 1 or 2 denoting two
      consecutive experiments}
    \item{\code{response}}{Number of connections [\%]}
  }
}
\references{
  Frimat, JP, Sisnaiske, J, Subbiah, S, Menne, H, Godoy, P, Lampen, P,
  Leist, M, Franzke, J, Hengstler, JG, van Thriel, C, West, J. The
  network formation assay: a spatially standardized neurite outgrowth
  analytical display for neurotoxicity screening. Lab Chip 2010; 10:701-709.
}

\examples{
data(nfa)

## Fit a four-parameter log-logistic model
nfa.m1 <- drm(response ~ conc, data = nfa, fct = LL.4())
summary(nfa.m1)
plot(nfa.m1, main = "NFA dose-response")
}

\keyword{datasets}
