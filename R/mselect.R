#' @title Dose-response model selection
#'
#' @description
#' Re-fits a dose-response model with alternative dose-response functions and
#' compares the resulting fits. The original model is always included.
#'
#' @param object A fitted object of class \code{"drc"}.
#' @param fctList An optional list of dose-response functions to fit to the same
#'   data. The original model is included even when it is not repeated in this
#'   list.
#' @param nested Logical; if \code{TRUE}, report p-values from tests comparing
#'   each model with the model immediately above it in the returned table.
#'   These tests are meaningful only when each pair is nested and is otherwise
#'   comparable. The first row has no comparison and its p-value is \code{NA}.
#' @param sorted Criterion used to order the result. \code{"IC"} (the default)
#'   orders by the information criterion in ascending order; \code{"Res var"}
#'   orders by residual variance in ascending order; \code{"Lack of fit"} orders
#'   by lack-of-fit p-value in descending order; and \code{"no"} preserves the
#'   input order. Missing values are placed last.
#' @param linreg Logical; if \code{TRUE}, also fit linear, quadratic, and cubic
#'   polynomial regressions. This option is available only for continuous data
#'   with a single dose-response curve.
#' @param icfct Function used to calculate the information criterion. Defaults
#'   to \code{\link{AIC}}; \code{\link{BIC}} can also be used. The function must
#'   return one numeric value for a fitted model.
#'
#' @details
#' The returned matrix reports the log likelihood, the selected information
#' criterion, the lack-of-fit p-value, and (for continuous data) the residual
#' variance. Smaller information-criterion and residual-variance values
#' indicate a better fit; a larger lack-of-fit p-value indicates less evidence
#' of lack of fit. A lack-of-fit p-value is \code{NA} when the fit type does not
#' provide this test. Log likelihoods should only be compared for models fit to
#' the same data under comparable likelihoods.
#'
#' Candidate models are re-fitted to the data stored in \code{object}, via
#' \code{\link{update.drc}}: the formula and other arguments of the original
#' call are evaluated in the environment the model was fitted in, so objects of
#' the same name in the global environment are never used. A
#' candidate that cannot be fitted is retained in the result with \code{NA}
#' statistics and a warning. \code{mselect()} returns a comparison table, not
#' the fitted candidate objects. With the default \code{sorted = "IC"}, the
#' first row is the model with the smallest information criterion.
#'
#' For nested comparisons, each p-value compares a row with the row immediately
#' above it after sorting. Thus, the selected \code{sorted} criterion determines
#' the comparison sequence. Failed fits have \code{NA} statistics; a nested
#' comparison is omitted if either model in that pair failed. A test error or
#' invalid p-value is reported as \code{NA} with a warning. These tests should
#' not be used to compare non-nested models.
#'
#' @return A numeric matrix with one row per model (including the original
#'   model) and columns for the log likelihood, information criterion,
#'   lack-of-fit p-value, and, for continuous data, residual variance. If
#'   \code{nested = TRUE}, a nested-test p-value column is also included unless
#'   \code{linreg = TRUE}.
#'
#' @author Christian Ritz
#'
#' @examples
#' ## Fit a starting model to the ryegrass data
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#'
#' ## Compare candidate dose-response functions; the smallest AIC is ranked first
#' candidate_fcts <- list(LL.3 = LL.3(), W1.4 = W1.4(), W2.4 = W2.4())
#' selection <- mselect(ryegrass.m1, fctList = candidate_fcts)
#' selection
#'
#' ## Identify and refit the AIC-selected model
#' best_name <- rownames(selection)[1]
#' best_fit <- if (best_name == ryegrass.m1$fct$name) {
#'   ryegrass.m1
#' } else {
#'   update(ryegrass.m1, fct = candidate_fcts[[best_name]])
#' }
#' best_name
#' best_fit$fct$name
#'
#' ## Add polynomial regressions to the comparison
#' mselect(ryegrass.m1, fctList = candidate_fcts, linreg = TRUE)
#'
#' ## Compare adjacent models in the returned order
#' mselect(ryegrass.m1, fctList = list(LL.5()), nested = TRUE)
#'
#' @keywords models nonlinear
"mselect" <-
function(object, fctList = NULL, nested = FALSE,
         sorted = c("IC", "Res var", "Lack of fit", "no"),
         linreg = FALSE, icfct = AIC)
{
    if (!inherits(object, "drc")) {
        stop("'object' must be a fitted object of class \"drc\".", call. = FALSE)
    }
    if (!is.logical(nested) || length(nested) != 1L || is.na(nested)) {
        stop("'nested' argument takes only the values: FALSE, TRUE (as a single non-missing logical value)",
             call. = FALSE)
    }
    if (!is.logical(linreg) || length(linreg) != 1L || is.na(linreg)) {
        stop("'linreg' must be a single non-missing logical value.", call. = FALSE)
    }
    isFct <- function(candidate) {
        is.function(candidate) ||
            (is.list(candidate) && is.function(candidate$fct))
    }
    if (!is.null(fctList) &&
        (!is.list(fctList) || isFct(fctList) ||
         any(!vapply(fctList, isFct, logical(1))))) {
        stop("'fctList' must be NULL or a list of dose-response functions; use list(LL.3(), ...).",
             call. = FALSE)
    }
    if (!is.function(icfct)) {
        stop("'icfct' must be a function such as AIC or BIC.", call. = FALSE)
    }

    sorted <- match.arg(sorted)
    contData <- identical(object$type, "continuous")
    if (sorted == "Res var" && !contData) {
        stop("'sorted = \"Res var\"' is only available for continuous data.",
             call. = FALSE)
    }
    if (linreg && !contData) {
        stop("'linreg = TRUE' is only available for continuous data.",
             call. = FALSE)
    }
    if (linreg && length(unique(object$dataList$curveid)) > 1L) {
        stop("'linreg = TRUE' is only available for a single dose-response curve.",
             call. = FALSE)
    }

    lenFL <- length(fctList)
    nCols <- 3L + as.integer(contData) + as.integer(nested)
    retMat <- matrix(NA_real_, lenFL + 1L, nCols)
    cnames <- c("logLik", "IC", "Lack of fit")
    if (contData) {
        cnames <- c(cnames, "Res var")
    }
    if (nested) {
        cnames <- c(cnames, "Nested F test")
    }
    colnames(retMat) <- cnames

    getModelName <- function(model, index = NULL) {
        name <- if (is.list(model)) model$name else NULL
        if (is.character(name) && length(name) == 1L &&
            !is.na(name) && nzchar(name)) {
            return(name)
        }
        if (!is.null(index) && !is.null(names(fctList))) {
            listName <- names(fctList)[index]
            if (!is.na(listName) && nzchar(listName)) {
                return(listName)
            }
        }
        if (is.null(index) && is.list(object$fct) &&
            is.character(object$fct$name) && length(object$fct$name) == 1L &&
            !is.na(object$fct$name) && nzchar(object$fct$name)) {
            return(object$fct$name)
        }
        if (is.null(index)) "Initial model" else paste("Model", index)
    }

    validatePValue <- function(value, description) {
        if (is.numeric(value) && length(value) == 1L &&
            is.na(value) && !is.nan(value)) {
            return(NA_real_)
        }
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
            value < 0 || value > 1) {
            warning(sprintf("%s is not a valid p-value; reporting NA.", description),
                    call. = FALSE)
            return(NA_real_)
        }
        as.numeric(value)
    }

    modelNames <- c(getModelName(object$fct), vapply(
        seq_len(lenFL),
        function(i) getModelName(fctList[[i]], i),
        character(1)
    ))
    rownames(retMat) <- make.unique(modelNames)
    fittedModels <- vector("list", lenFL + 1L)
    fittedModels[[1L]] <- object

    metricValue <- function(fun, fit, metric, modelName, required = FALSE) {
        value <- tryCatch(fun(fit), error = function(e) e)
        if (inherits(value, "error")) {
            message <- sprintf("Could not calculate %s for model '%s': %s",
                               metric, modelName, conditionMessage(value))
            if (required) {
                stop(message, call. = FALSE)
            }
            warning(paste0(message, "; reporting NA."), call. = FALSE)
            return(NA_real_)
        }
        if (!is.numeric(value) || length(value) != 1L || is.nan(value)) {
            message <- sprintf(
                "The calculation of %s for model '%s' must return one numeric value.",
                metric, modelName
            )
            if (required) {
                stop(message, call. = FALSE)
            }
            warning(paste0(message, " Reporting NA."), call. = FALSE)
            return(NA_real_)
        }
        as.numeric(value)
    }

    getLackOfFit <- function(fit, modelName) {
        result <- tryCatch(modelFit(fit), error = function(e) e)
        if (inherits(result, "error")) {
            warning(sprintf(
                "Could not calculate the lack-of-fit p-value for model '%s': %s; reporting NA.",
                modelName, conditionMessage(result)
            ), call. = FALSE)
            return(NA_real_)
        }
        if (is.null(result) || !is.matrix(result) && !is.data.frame(result) ||
            nrow(result) < 2L || !"p value" %in% colnames(result)) {
            return(NA_real_)
        }
        validatePValue(result[2L, "p value"],
                       sprintf("The lack-of-fit p-value for model '%s'", modelName))
    }

    getResidualVariance <- function(fit, modelName) {
        metricValue(function(x) summary(x)$resVar, fit, "residual variance",
                    modelName)
    }

    setModelMetrics <- function(row, fit, modelName, initial = FALSE) {
        retMat[row, "logLik"] <<- metricValue(logLik, fit, "log likelihood",
                                              modelName, required = initial)
        retMat[row, "IC"] <<- metricValue(icfct, fit, "information criterion",
                                          modelName, required = initial)
        retMat[row, "Lack of fit"] <<- getLackOfFit(fit, modelName)
        if (contData) {
            retMat[row, "Res var"] <<- getResidualVariance(fit, modelName)
        }
    }

    setModelMetrics(1L, object, rownames(retMat)[1L], initial = TRUE)

    if (lenFL > 0L) {
        for (i in seq_len(lenFL)) {
            row <- i + 1L
            modelName <- rownames(retMat)[row]
            tempObj <- tryCatch(
                update(object, fct = fctList[[i]], data = object[["origData"]]),
                error = function(e) e
            )

            if (inherits(tempObj, "error")) {
                warning(sprintf("Could not fit candidate model '%s': %s; reporting NA statistics.",
                                modelName, conditionMessage(tempObj)),
                        call. = FALSE)
                next
            }
            if (!inherits(tempObj, "drc")) {
                warning(sprintf(
                    "Candidate model '%s' did not return a fitted object of class \"drc\"; reporting NA statistics.",
                    modelName
                ), call. = FALSE)
                next
            }

            setModelMetrics(row, tempObj, modelName)
            fittedModels[[row]] <- tempObj
        }
    }

    if (linreg) {
        drcData <- as.data.frame(object$data[, c(2L, 1L), drop = FALSE])
        names(drcData) <- c("yVec", "xVec")

        linFitList <- tryCatch(
            list(
                lm(yVec ~ xVec, data = drcData),
                lm(yVec ~ xVec + I(xVec * xVec), data = drcData),
                lm(yVec ~ xVec + I(xVec * xVec) + I(xVec * xVec * xVec),
                   data = drcData)
            ),
            error = function(e) e
        )
        if (inherits(linFitList, "error")) {
            stop(sprintf("Could not fit the polynomial regression models: %s",
                         conditionMessage(linFitList)), call. = FALSE)
        }
        linModMat <- matrix(NA_real_, 3L, 4L)
        for (i in seq_along(linFitList)) {
            linModMat[i, ] <- c(
                metricValue(logLik, linFitList[[i]], "log likelihood",
                            c("Lin", "Quad", "Cubic")[i], required = TRUE),
                metricValue(icfct, linFitList[[i]], "information criterion",
                            c("Lin", "Quad", "Cubic")[i], required = TRUE),
                NA_real_,
                metricValue(function(x) summary(x)$sigma^2, linFitList[[i]],
                            "residual variance", c("Lin", "Quad", "Cubic")[i],
                            required = TRUE)
            )
        }
        rownames(linModMat) <- c("Lin", "Quad", "Cubic")
        linColNames <- cnames[seq_len(4L)]
        colnames(linModMat) <- linColNames

        if (nested) {
            retMat <- retMat[, seq_len(4L), drop = FALSE]
        }
        retMat <- rbind(retMat, linModMat)
    }

    if (sorted != "no") {
        decreasing <- identical(sorted, "Lack of fit")
        sortOrder <- order(retMat[, sorted], decreasing = decreasing,
                           na.last = TRUE)
        retMat <- retMat[sortOrder, , drop = FALSE]
        if (nested && !linreg) {
            fittedModels <- fittedModels[sortOrder]
        }
    }

    if (nested && !linreg) {
        retMat[, "Nested F test"] <- NA_real_
        if (length(fittedModels) > 1L) {
            for (i in 2:length(fittedModels)) {
                previousModel <- fittedModels[[i - 1L]]
                currentModel <- fittedModels[[i]]
                if (is.null(previousModel) || is.null(currentModel)) {
                    next
                }

                previousName <- rownames(retMat)[i - 1L]
                currentName <- rownames(retMat)[i]
                nestedResult <- tryCatch(
                    anova(previousModel, currentModel, details = FALSE),
                    error = function(e) e
                )
                if (inherits(nestedResult, "error")) {
                    warning(sprintf(
                        "Could not compare adjacent models '%s' and '%s': %s; reporting NA for the nested-test p-value.",
                        previousName, currentName, conditionMessage(nestedResult)
                    ), call. = FALSE)
                } else if (is.matrix(nestedResult) || is.data.frame(nestedResult)) {
                    pCol <- match("p value", colnames(nestedResult))
                    if (!is.na(pCol) && nrow(nestedResult) >= 2L) {
                        retMat[i, "Nested F test"] <- validatePValue(
                            nestedResult[2L, pCol],
                            sprintf("The nested-test p-value for models '%s' and '%s'",
                                    previousName, currentName)
                        )
                    }
                }
            }
        }
    }
    retMat
}
