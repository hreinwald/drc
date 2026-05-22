#' @title Updating and re-fitting a model
#'
#' @description
#' \code{update} updates and re-fits a model on the basis of an object of class 'drc'.
#'
#' @param object an object of class 'drc'.
#' @param ... arguments to alter in object.
#' @param evaluate logical. If TRUE model is re-fit; otherwise an unevaluated call is returned.
#'
#' @return An object of class 'drc'.
#'
#' @author Christian Ritz
#'
#' @examples
#' ## Fitting a four-parameter Weibull model
#' model1 <- drm(ryegrass, fct = W1.4())
#'
#' ## Updating 'model1' by fitting a three-parameter Weibull model instead
#' model2 <- update(model1, fct = W1.3())
#' anova(model2, model1)
#'
#' @keywords models nonlinear
"update.drc" <- function (object, ..., evaluate = TRUE) 
{
    call <- object$call
    if (is.null(call)) 
        stop("need an object with call component")
        
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras) > 0) 
    {
        glsa <- names(as.list(args(drm)))
        names(extras) <- glsa[pmatch(names(extras), glsa[-length(glsa)])]
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) 
        {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate)
    {
        ## If call$data refers to a symbol that cannot be resolved in the
        ## calling frame (e.g. .x inside purrr::map()), fall back to the
        ## data stored in the fitted object.
        eval_env <- parent.frame()
        if (!is.null(call$data) && !is.null(object$origData))
        {
            data_resolvable <- tryCatch(
            {
                is.data.frame(eval(call$data, eval_env))
            },
            error = function(e) FALSE)
            if (!data_resolvable)
            {
                eval_env <- list2env(
                    list(.drc_stored_data__ = object$origData),
                    parent = eval_env)
                call$data <- quote(.drc_stored_data__)
            }
        }
        eval(call, eval_env)
    } else call
}
