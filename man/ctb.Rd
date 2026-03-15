\name{ctb}
\alias{ctb}
\docType{data}
\title{CellTiter-Blue Cell Viability Assay Data}
\description{Neurotoxicity test using the CellTiter-Blue Cell Viability
  Assay on SH-SY5Y cells for increasing concentrations of acrylamide.}
\usage{data(ctb)}
\format{
  A data frame with 647 observations on the following 5 variables.
  \describe{
    \item{\code{well}}{well ID of a 96 well plate}
    \item{\code{conc}}{12 concentrations of acrylamide, ranging from
      0-500mM}
    \item{\code{fluorescence}}{measured fluorescence after adding the
      resazurin reagent into the wells}
    \item{\code{day}}{integer denoting 3 different days}
    \item{\code{plate}}{factor with 7 levels representing the plate ID}
  }
}
\references{
  Frimat, JP, Sisnaiske, J, Subbiah, S, Menne, H, Godoy, P, Lampen, P,
  Leist, M, Franzke, J, Hengstler, JG, van Thriel, C, West, J. The
  network formation assay: a spatially standardized neurite outgrowth
  analytical display for neurotoxicity screening. Lab Chip 2010; 10:701-709.
}


\examples{
data(ctb)
ctb$day <- as.factor(ctb$day)

## Fit a four-parameter log-logistic model
ctb.m1 <- drm(fluorescence ~ conc, data = ctb, fct = LL.4())
summary(ctb.m1)
plot(ctb.m1, main = "CTB dose-response")
}



\keyword{datasets}
