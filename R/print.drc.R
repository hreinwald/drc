#' @title Printing key features
#'
#' @description
#' \code{print} displays brief information on an object of class 'drc'.
#'
#' @param x an object of class 'drc'.
#' @param ... additional arguments.
#' @param digits an integer giving the number of digits of the parameter coefficients. Default is 3.
#'
#' @return The object is returned invisibly.
#'
#' @author Christian Ritz
#'
#' @examples
#' ## Fitting a four-parameter log-logistic model
#' ryegrass.m1 <- drm(rootl ~conc, data = ryegrass, fct = LL.4())
#'
#' ## Displaying the model fit
#' print(ryegrass.m1)
#' ryegrass.m1  # gives the same output as the previous line
#'
#' @keywords models nonlinear
"print.drc" <- function(x, ..., digits = max(3, getOption("digits") - 3)) 
{
    object <- x

    classList <- class(object)
    cat(paste("\n", "A 'drc' model.", "\n", sep=""))
    
    ## Borrowing from print.lm
    cat("\nCall:\n", deparse(object$"call"), "\n\n", sep = "")
    if (length(coef(object))>0) 
    {
        cat("Coefficients:\n")
        print.default(format(coef(object), digits = digits), print.gap = 2, quote = FALSE)
    } else {
        cat("No coefficients\n")
    }
    cat("\n")
    
    invisible(object)
}
