#' @title Testing if there is a dose effect at all
#'
#' @description
#' A significance test is provided for the comparison of the dose-response model considered
#' and the simple linear regression model with slope 0 (a horizontal regression line
#' corresponding to no dose effect).
#'
#' @param object an object of class 'drc'.
#'
#' @details Perhaps useful for screening purposes.
#'
#' @return The likelihood ratio test statistic and the corresponding degrees of freedom
#'   and p-value are reported.
#'
#' @author Christian Ritz
#'
#' @examples
#' ryegrass.LL.4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#'
#' noEffect(ryegrass.LL.4)
#' # p-value < 0.0001: there is a highly significant dose effect!
#'
#' @keywords models nonlinear
noEffect <- function(object)
{
    if (identical(object$"type", "binomial"))
    {
        respVec <- object$"dataList"$resp
        weiVec <- object$"dataList"$weights
        llNull <- logLik(glm(cbind(respVec*weiVec, (1-respVec)*weiVec) ~ 1, family = binomial))
    }
    
    if (identical(object$"type", "Poisson"))
    {
        respVec <- object$"dataList"$resp
        llNull <- logLik(glm(respVec ~ 1, family = poisson))
    }

    if (identical(object$"type", "continuous"))
    {
        llNull <- logLik(lm(object$dataList$resp ~ 1))
    }
    lldrc <- logLik(object)
    lrt <- -2*(llNull - lldrc)
    dfDiff <- attr(lldrc, "df") - attr(llNull, "df")

    retVec <- c(lrt, dfDiff, 1 - pchisq(lrt, dfDiff))
    names(retVec) <- c("Chi-square test", "Df", "p-value")
    retVec
}