#' @title Sets control arguments
#'
#' @description Set control arguments in the control argument in the function \code{\link{drm}}.
#'
#' @param constr logical. If \code{TRUE} optimisation is constrained, only yielding non-negative
#'   parameters.
#' @param errorm logical specifying whether failed convergence in \code{\link{drm}} should
#'   result in an error or only a warning.
#' @param maxIt numeric. The maximum number of iterations in the optimisation procedure.
#' @param method character string. The method used in the optimisation procedure. See
#'   \code{\link{optim}} for available methods.
#' @param noMessage logical, specifying whether or not messages should be displayed.
#' @param relTol numeric. The relative tolerance in the optimisation procedure.
#' @param rmNA logical. Should \code{NA}s be removed from sum of squares used for estimation?
#'   Default is \code{FALSE} (not removed).
#' @param useD logical. If \code{TRUE} derivatives are used for estimation (if available).
#' @param trace logical. If \code{TRUE} the trace from \code{\link{optim}} is displayed.
#' @param otrace logical. If \code{TRUE} the output from \code{\link{optim}} is displayed.
#' @param warnVal numeric. If equal to 0 then the warnings are stored and displayed at the end.
#'   See under \sQuote{warn} in \code{\link{options}}. The default results in suppression of
#'   warnings.
#' @param dscaleThres numeric value specifying the threshold for dose scaling.
#' @param rscaleThres numeric value specifying the threshold for response scaling.
#' @param conCheck logical, switching on/off handling of control measurements.
#'
#' @return A list with components corresponding to each of the above arguments.
#'
#' @seealso \code{\link{drm}}, \code{\link{optim}}
#'
#' @examples
#' ## Displaying the default settings
#' drmc()
#'
#' ## Using the 'method' argument
#' model1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#' model2 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
#'   control = drmc(method = "Nelder-Mead"))
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"drmc" <- function(constr = FALSE, errorm = TRUE, maxIt = 500, method = "BFGS", 
noMessage = FALSE, relTol = 1e-7, rmNA = FALSE, useD = FALSE, trace = FALSE, 
otrace = FALSE, warnVal = -1, dscaleThres = 1e-15, rscaleThres = 1e-15, conCheck = TRUE)
{
    return(list(
                constr = constr,
                errorm = errorm,
                maxIt = maxIt, 
                method = method,
                noMessage = noMessage,
                relTol = relTol,
                rmNA = rmNA, 
                useD = useD,
                trace = trace,
                otrace = otrace,
                warnVal = warnVal,
                dscaleThres = dscaleThres,
                rscaleThres = rscaleThres,
                conCheck = conCheck
                ))
}
