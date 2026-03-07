#' Searching through a range of initial parameter values to obtain convergence
#'
#' \code{searchdrc} provides a facility for searching through a range of parameter values
#' (one-dimensional) in order to obtain convergence of the estimation procedure.
#'
#' The function goes through the range with increments such that in total at most \code{len}
#' sets of parameter values are used as initial values for the estimation procedure. You would
#' need to identify the parameter which is most likely to cause problems for the estimation
#' procedure.
#'
#' @param object an object of class 'drc'. The object can be from a model that could not be fitted.
#' @param which a character string containing the parameter name.
#' @param range a numeric vector of length 2 specifying the interval endpoints for the range.
#' @param len numeric. The number of points in the interval.
#'
#' @return An object of class 'drc'.
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"searchdrc" <- function(object, which, range, len = 50)
{
   sv <- object$start
 
   parNames <- object$parNames[[2]]
   whichInd <- regexpr(paste("^", which, ":", sep = ""), parNames)
   whichInd <- ((1:length(parNames))[whichInd>0])[1]
   
   if (length(whichInd)<1) {stop(paste("No such parameter ", which, sep = ""))}

   found <- FALSE
   oldWarn <- getOption("warn")
   on.exit(options(warn = oldWarn), add = TRUE)
   for (i in seq(range[1], range[2], length.out = len))
   {
       sv[whichInd] <- i 

       options(warn = -1)
       modelFit <- try(update(object, start = sv, control = drmc(noMessage = TRUE)), silent=TRUE)
       options(warn = oldWarn)
       if (!inherits(modelFit, "try-error")) {found <- TRUE; break}
   }
   
   if (found) 
   {
       return(modelFit)
   } else {
       warning("Convergence failed.", call. = FALSE)
   }
}
