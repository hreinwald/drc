#' @title ANOVA for dose-response model fits
#'
#' @description
#' \code{anova} produces an analysis of variance table for one or two
#' non-linear model fits.
#'
#' @param object an object of class 'drc'.
#' @param ... additional arguments.
#' @param details logical indicating whether or not details on the models
#'   compared should be displayed. Default is TRUE (details are displayed).
#' @param test a character string specifying the test statistic to be applied.
#'   Use \code{"od"} to assess overdispersion for binomial data.
#'
#' @return An object of class 'anova'.
#'
#' @details
#' Specifying only a single object gives a test for lack-of-fit, comparing the
#' non-linear regression model to a more general one-way or two-way ANOVA
#' model.
#'
#' If two objects are specified a test for reduction from the larger to the
#' smaller model is given. (This only makes statistical sense if the models are
#' nested, that is: one model is a submodel of the other model.)
#'
#' @examples
#' ## Comparing Weibull three- and four-parameter models using an F test
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
#' ryegrass.m2 <- drm(rootl ~ conc, data = ryegrass, fct = W1.3())
#' anova(ryegrass.m2, ryegrass.m1)
#'
#' anova(ryegrass.m2, ryegrass.m1, details = FALSE)  # without details
#'
#' @author Christian Ritz
#'
#' @keywords models nonlinear
"anova.drc" <-
function(object, ..., details = TRUE, test = NULL)
{
    if (length(list(object, ...)) > 1) 
    {
        return(anova.drclist(object, ..., details = details, test = test))
    } else {
        stop("Use the function modelFit()")
    }
}
