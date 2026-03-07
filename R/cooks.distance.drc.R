#' @title Cook's distance for nonlinear dose-response models
#'
#' @description
#' Cook's distance values are provided for nonlinear dose-response model fits using the
#' same formulas as in linear regression but based on the corresponding approximate quantities
#' available for nonlinear models.
#'
#' @param model an object of class 'drc'.
#' @param ... additional arguments (not used).
#'
#' @return A vector of Cook's distance values, one value per observation.
#'
#' @author Christian Ritz
#'
#' @examples
#' ryegrass.LL.4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#'
#' cooks.distance(ryegrass.LL.4)
#'
#' @keywords models nonlinear
cooks.distance.drc <- function(model, ...)
{
    hatVal <- hatvalues(model)
    
    mse <- (rse(model))^2
    if (is.na(mse)) {mse <- 1}  
    # default for generalized linear models (assuming no overdispersion!)
    
    (residuals(model)^2) * (hatVal/((1-hatVal)^2)) / ((length(hatVal)- df.residual(model)) * mse)
}
