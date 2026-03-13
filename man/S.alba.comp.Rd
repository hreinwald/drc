\name{S.alba.comp}

\alias{S.alba.comp}

\docType{data}

\title{Potency of two herbicides}

\description{
  Data are from an experiment, comparing the potency of the two herbicides glyphosate and bentazone in
  white mustard \emph{Sinapis alba}.
}

\usage{data(S.alba.comp)}

\format{
  A data frame with 141 observations on the following 8 variables.
  \describe{
  \item{\code{exp}}{a factor with levels ben1, ben2, gly1, gly2 indicating which experiment each observation belongs to.}
    \item{\code{herbicide}}{a factor with levels \code{Bentazone} \code{Glyphosate} (the two herbicides applied).}
    \item{\code{dose}}{a numeric vector containing the dose in g/ha.}
    \item{\code{drymatter}}{a numeric vector containing the response (dry matter in g/pot).}
    \item{\code{Tf}}{a numeric vector .}
    \item{\code{area}}{a numeric vector .}
    \item{\code{Fo}}{a numeric vector .}
    \item{\code{Fm}}{a numeric vector .}
    }
}

\details{
  The lower and upper limits for the two herbicides can be assumed identical, whereas slopes and ED50 values 
  are different (in the log-logistic model).
}

\source{
  Christensen, M. G. and Teicher, H. B., and Streibig, J. C. (2003) Linking fluorescence 
  induction curve and biomass in herbicide screening, \emph{Pest Management Science},
  \bold{59},  1303--1310.
}


\keyword{datasets}
