\name{auxins}

\alias{auxins}

\docType{data}

\title{Effect of technical grade and commercially formulated auxin herbicides}

\description{
  MCPA, 2,4-D, mecorprop and dichorlprop were applied either as technical grades
  materials or as commercial formulations. 
  Each experimental unit consisted of five 1-week old seedlings grown together
  in a pot of nutrient solution during 14 days. 
}


\usage{data(auxins)}

\format{
  A data frame with 150 observations on the following 5 variables.
  \describe{
    \item{\code{dryweight}}{a numeric vector}
    \item{\code{dose}}{a numeric vector}
    \item{\code{replicate}}{a factor with 3 levels}
    \item{\code{herbicide}}{a factor with 5 levels}
    \item{\code{formulation}}{a factor with 2 levels}
  }
}

\details{
  Data are parts of a larger joint action experiment with various herbicides.

  The eight herbicide preparations are naturally grouped into four pairs (herbicide:formulation) + control, and each pair of herbicides should have the same active ingredients but different formulation
  constituents, which were assumed to be biologically inert. The data consist
  of the 150 observations of dry weights, each observation being the weight
  of five plants grown in the same pot. All the eight herbicide preparations have
  essentially the same mode of action in the plant; they all act like the plant 
  auxins, which are plant regulators that affect cell enlongation an other 
  essential metabolic pathways. One of the objects of the experiment was to test
  if the response functions were identical except for a multiplicative factor in 
  the dose. This is a necessary, but not a sufficient, condition for a similar
  mode of action for the herbicides. 
}

\source{
  Streibig, J. C. (1987). Joint action of root-absorbed mixtures of auxin 
  herbicides in Sinapis alba L. and barley (Hordeum vulgare L.) 
  \emph{Weed Research}, \bold{27}, 337--347.
}

\references{
  Rudemo, M., Ruppert, D., and Streibig, J. C. (1989). Random-Effect Models 
  in Nonlinear Regression with Applications to Bioassay.
  \emph{Biometrics}, \bold{45}, 349--362.
}


\keyword{datasets}
