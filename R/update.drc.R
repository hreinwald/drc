#' @title Updating and re-fitting a model
#'
#' @description
#' \code{update} updates and re-fits a model on the basis of an object of class 'drc'.
#'
#' @param object an object of class 'drc'.
#' @param ... arguments to alter in object.
#' @param evaluate logical. If TRUE model is re-fit; otherwise an unevaluated call is returned.
#'
#' @details
#' The arguments stored in the call of the fit (e.g. \code{formula = myFormula},
#' \code{data = myData}) are evaluated in the environment in which \code{drm()} was
#' called, so \code{update()} - and functions that refit models through it, such as
#' \code{\link{mselect}} and \code{\link{maED}} - also work for fits created inside
#' functions or in a separate environment, and never pick up an object of the same
#' name from the global environment by accident. New arguments supplied in \code{...}
#' are evaluated in the environment \code{update()} is called from.
#'
#' If \code{data} is not supplied and the original \code{data} expression can no
#' longer be found or now refers to different data, the data stored in the fit
#' (\code{object$origData}) is used (with a warning in the latter case).
#'
#' Fits made inside a function keep a reference to that function's environment
#' (as \code{lm()} does through its formula), which is also saved by \code{saveRDS()}.
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
    if (!evaluate) return(call)

    ## Arguments stored in the call are resolved in the environment the model was
    ## fitted in (attribute set by drm()); the new arguments in '...' keep being
    ## resolved in the caller's environment. Objects without the attribute (saved
    ## by older versions) are evaluated in the caller's environment as before.
    callerEnv <- parent.frame()
    fitEnv <- attr(object, "callEnv")
    if (!is.environment(fitEnv) || identical(fitEnv, callerEnv))
    {
        eval_env <- callerEnv
    } else {
        eval_env <- new.env(parent = fitEnv)
        notFound <- new.env()  # unique sentinel
        for (nm in unique(unlist(lapply(extras, all.names))))
        {
            val <- tryCatch(get(nm, envir = callerEnv), error = function(e) notFound)
            if (identical(val, notFound)) next
            fitVal <- tryCatch(get(nm, envir = fitEnv), error = function(e) notFound)
            if (!identical(fitVal, notFound) && identical(fitVal, val)) next
            assign(nm, val, envir = eval_env)
        }
    }

    ## Refits must use the data of the fit: fall back to the stored data if the
    ## 'data' expression cannot be resolved (e.g. .x inside purrr::map()) or if it
    ## now refers to different data than the fit was based on.
    if (!is.null(call$data) && !is.null(object$origData))
    {
        dataVal <- tryCatch(eval(call$data, eval_env), error = function(e) NULL)
        useStored <- !is.data.frame(dataVal)
        if (!useStored && !("data" %in% names(extras)) && !identical(dataVal, object$origData))
        {
            warning("'", paste(deparse(call$data), collapse = " "),
                    "' differs from the data the model was fitted to; ",
                    "refitting with the stored data (pass 'data' to update() to use other data)",
                    call. = FALSE)
            useStored <- TRUE
        }
        if (useStored)
        {
            eval_env <- list2env(
                list(.drc_stored_data__ = object$origData),
                parent = eval_env)
            call$data <- quote(.drc_stored_data__)
        }
    }
    eval(call, eval_env)
}
