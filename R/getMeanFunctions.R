#' Display available dose-response models
#'
#' Display information about available, built-in dose-response models.
#' The arguments \code{noParm} and \code{fname} can be combined.
#'
#' @param noParm numeric specifying the number of parameters of the models to be displayed.
#'   The default (NA) results in display of all models, regardless of number of parameters.
#' @param fname character string or vector of character strings specifying the short name(s)
#'   of the models to be displayed (need to match exactly).
#' @param flist list of built-in functions to be displayed.
#' @param display logical indicating whether or not the requested models should be displayed
#'   on the R console.
#'
#' @return An invisible list of functions or a list of strings with brief function descriptions.
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @examples
#' ## Listing all functions
#' getMeanFunctions()
#'
#' ## Listing all functions with 4 parameters
#' getMeanFunctions(4)
#'
#' ## Listing all (log-)logistic functions
#' getMeanFunctions(fname = "L")
#'
#' ## Listing all three-parameter (log-)logistic or Weibull functions
#' getMeanFunctions(3, fname = c("LL", "W"))
#'
#' @keywords models nonlinear
"getMeanFunctions" <- function(noParm = NA, fname = NULL, flist = NULL, display = TRUE)
{
    if (is.null(flist))
    {
        fctList <- list(
        LL.2(), LL.3(), LL.3u(), LL.4(), LL.5(),
        W1.2(), W1.3(), W1.3u(), W1.4(),
        W2.2(), W2.3(), W2.3u(), W2.4(),
        BC.4(), BC.5(), 
        LL2.2(), LL2.3(), LL2.3u(), LL2.4(), LL2.5(),                        
        AR.2(), AR.3(),
        EXD.2(), EXD.3(),
        G.2(), G.3(), G.3u(), G.4(),
        LN.2(), LN.3(), LN.3u(), LN.4(),
        L.3(), L.4(), L.5(),
        MM.2(), MM.3(),
        suppressWarnings(CRS.4a()), suppressWarnings(CRS.4b()), suppressWarnings(CRS.4c()),
        suppressWarnings(CRS.5a()), suppressWarnings(CRS.5b()), suppressWarnings(CRS.5c()),
        CRS.6(),
        UCRS.4a(fctName = "UCRS.4a"), UCRS.4b(fctName = "UCRS.4b"), UCRS.4c(fctName = "UCRS.4c"),
        UCRS.5a(fctName = "UCRS.5a"), UCRS.5b(fctName = "UCRS.5b"), UCRS.5c(fctName = "UCRS.5c"),
        NEC.2(), NEC.3(), NEC.4(),
        FPL.4(-1, 1),
        W2x.3(), W2x.4(),
        baro5(), gammadr()
        )
    } else {
        fctList <- flist
    }

    textVec <- NULL
    lapFct <- function(x) {c(x$"name", x$"text")} 
    if (!is.null(fname))
    {
        textVec <- fname
    }
    

    displayFunction <- function(object)
    {
        if ( ((is.na(noParm)) || (noParm == object$"noParm")) 
        && (is.null(textVec) || (any(textVec %in% object$"name"))) )
        {
            if (display)
            {
                cat(object$"text", "\n")
                cat(paste("(", object$"noParm", " parameters)", sep = ""), "\n")
                cat("In 'drc': ", object$"name", "\n\n")
            }
            return(object)
        } else {
            return(NULL)
        }
    }
    lapList1 <- lapply(fctList, displayFunction)
    
    ## Removing NULL elements in the list 'lapList1'
    lapList2 <- list()
    counter <- 1
    for (i in 1:length(lapList1))
    {
        if (!is.null(lapList1[[i]]))
        {
            lapList2[[counter]] <- lapList1[[i]]
            counter <- counter + 1
        }
    }
    
    if (!is.null(textVec) || !is.na(noParm))
    {
        invisible(lapList2)
    } else {
        invisible(lapply(lapList2, lapFct))
    }
}
