\name{mdra}
\alias{mdra}
\docType{data}
\title{3T3 mouse fibroblasts and NRU assay}
\description{
The toxicity of sodium valproate was tested, using the 3T3 mouse fibroblasts and neutral red uptake (NRU) assay. 22 different experiments were performed independently in six laboratories, using eight concentration levels, each with six replicates on a 96-well plate. In addition, twelve measurements were taken for the solvent control.
}
\usage{data("mdra")}
\format{
  A data frame with 1320 observations on the following 4 variables.
  \describe{
    \item{\code{LabID}}{a factor with levels \code{A} \code{B} \code{C} \code{D} \code{E} \code{F}}
    \item{\code{ExperimentID}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10} \code{11} \code{12} \code{13} \code{14} \code{15} \code{16} \code{17} \code{18} \code{19} \code{20} \code{21} \code{22}}
    \item{\code{Concentration}}{a numeric vector}
    \item{\code{Response}}{a numeric vector}
  }
}
\source{
http://biostatistics.dkfz.de/download/mdra/MDRA_ExampleData.csv
}
\references{
Clothier, R., Gomez-Lechon, M. J., Kinsner-Ovaskainen, A., Kopp-Schneider, A., O'Connor, J. E., Prieto, P., and Stanzel, S. (2013). Comparative analysis of eight cytotoxicity assays evaluated within the ACuteTox Project. Toxicology in vitro, 27(4):1347--1356.
}
\examples{
data(mdra)

## Fit a three-parameter log-logistic model
mdra.m1 <- drm(Response ~ Concentration, data = mdra, fct = LL.3())
summary(mdra.m1)
plot(mdra.m1, main = "MDRA dose-response")
}
\keyword{datasets}
