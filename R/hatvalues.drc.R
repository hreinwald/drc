#' @title Model diagnostics for nonlinear dose-response models
#'
#' @description
#' Hat values (leverage values) are provided for nonlinear dose-response model fits using the
#' same formulas as in linear regression but based on the corresponding approximate quantities
#' available for nonlinear models.
#'
#' @param model an object of class 'drc'.
#' @param ... additional arguments (not used).
#'
#' @details
#' Hat values are calculated using the formula given by Cook et al. (1986) and
#' McCullagh and Nelder (1989). The output values can be assessed in the same way as
#' in linear regression.
#'
#' @return A vector of leverage values (hat values), one value per observation.
#'
#' @references
#' Cook, R. D. and Tsai, C.-L. and Wei, B. C. (1986)
#' Bias in Nonlinear Regression,
#' \emph{Biometrika} \bold{73}, 615--623.
#'
#' McCullagh, P. and Nelder, J. A. (1989)
#' \emph{Generalized Linear Models},
#' Second edition, Chapman & Hall/CRC.
#'
#' @author Christian Ritz
#'
#' @examples
#' ryegrass.LL.4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#'
#' hatvalues(ryegrass.LL.4)
#'
#' @keywords models nonlinear
hatvalues.drc <- function(model, ...)
{
    xmat <- model$der
    diag(xmat %*% ginv(t(xmat) %*% xmat) %*% t(xmat))
#    names(hvector) <- as.character(1:length(hvector))
#    hvector
}

