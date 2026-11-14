#' Estimation of ED values using model-averaging
#'
#' Estimates and confidence intervals for ED values are estimated using
#' model-averaging.
#'
#' Model-averaging of individual estimates is carried out as described by
#' Buckland \emph{et al.} (1997) and Kang \emph{et al.} (2000) using
#' AIC-based weights. The two approaches differ w.r.t. the calculation of
#' confidence intervals: Buckland \emph{et al.} (1997) provide an approximate
#' variance formula under the assumption of perfectly correlated estimates
#' (so, confidence intervals will tend to be too wide). Kang \emph{et al.}
#' (2000) use the model weights to calculate confidence limits as weighted
#' means of the confidence limits for the individual fits.
#'
#' @param object an object of class \code{drc}.
#' @param fctList a list of non-linear functions to be compared. Models with a log(ED)
#'   parameterisation (\code{LL2.*}, \code{LN.*(loge = TRUE)}) are not allowed (an error with
#'   instructions is raised); use their dose-scale equivalents \code{LL.*} / \code{LN.*}.
#' @param respLev a numeric vector containing the response levels.
#' @param interval character string specifying the type of confidence intervals
#'   to be supplied. The default is \code{"none"}. The choices \code{"buckland"}
#'   and \code{"kang"} are explained in the Details section.
#' @param linreg logical indicating whether or not additionally a simple linear
#'   regression model should be fitted.
#' @param clevel character string specifying the curve id in case estimates for
#'   a specific curve or compound are requested. By default estimates are shown
#'   for all curves.
#' @param level numeric. The confidence level. Must be a single value strictly
#'   between 0 and 1. The default is \code{0.95}. It applies to both interval types; for
#'   \code{"kang"} it is also the level of the delta-method limits of the individual models
#'   (including the linear model if \code{linreg = TRUE}) that are averaged.
#' @param type character string. Whether the specified response levels are
#'   absolute or relative (default).
#' @param display logical. If \code{TRUE} results are displayed. Otherwise they
#'   are not (useful in simulations).
#' @param na.rm logical indicating whether or not \code{NA} values occurring
#'   during model fitting should be excluded from subsequent calculations.
#' @param extended logical specifying whether or not an extended output
#'   (including fit summaries) should be returned.
#' @param onExclude character string, what to do if models cannot contribute to the
#'   average because they could not be re-fitted or have non-finite ED values:
#'   \code{"warn"} (default) excludes them and issues one summary warning; \code{"stop"}
#'   raises an error instead (recommended for scripted pipelines, where a silently
#'   reduced model set would change the scientific result).
#'
#' @return If \code{extended = FALSE}, a matrix with two or more columns
#'   containing the model-averaged estimates and the corresponding estimated
#'   standard errors and, optionally, lower and upper confidence limits.
#'   If \code{extended = TRUE}, a list with components:
#'   \describe{
#'     \item{estimates}{Matrix of model-averaged ED estimates and intervals.}
#'     \item{fits}{Matrix of per-model ED estimates and AIC-based weights.}
#'     \item{models}{Data frame with one row per model (original model, candidates and,
#'       if requested, the linear model) and columns \code{model}, \code{AIC},
#'       \code{weight} and \code{status} (\code{"used"}, \code{"excluded: fit failed"}
#'       or \code{"excluded: non-finite ED"}). For several curves a column \code{curve}
#'       is added and the matrices are stacked.}
#'   }
#'   If at least one model was excluded, the returned matrix (\code{estimates} for
#'   \code{extended = TRUE}) additionally carries the attributes \code{"models"} (the data
#'   frame described above), \code{"n_used"} and \code{"n_requested"} (the number of models
#'   that contributed and that were requested; named by curve for several curves). The
#'   attributes are absent if all models were used.
#'
#' @details
#' The models in \code{fctList} are re-fitted with \code{\link{update.drc}}, i.e.
#' the formula and data of the original call are evaluated in the environment
#' \code{object} was fitted in (and \code{object$origData} is used if the data
#' was changed afterwards). A model whose \emph{fit} fails (e.g. no convergence) or
#' which yields non-finite ED values is excluded from the average (see
#' \code{onExclude}). If instead the call stored in \code{object} cannot be
#' re-evaluated at all (e.g. it refers to objects that no longer exist), no candidate
#' model can be fitted and \code{maED()} stops with an explanatory error rather than
#' returning the estimate of the original model as if it were a model average.
#'
#' @references
#'   Buckland, S. T. and Burnham, K. P. and Augustin, N. H. (1997)
#'   Model Selection: An Integral Part of Inference,
#'   \emph{Biometrics} \bold{53}, 603--618.
#'
#'   Kang, Seung-Ho and Kodell, Ralph L. and Chen, James J. (2000)
#'   Incorporating Model Uncertainties along with Data Uncertainties in
#'   Microbial Risk Assessment,
#'   \emph{Regulatory Toxicology and Pharmacology} \bold{32}, 68--72.
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso The function \code{\link{mselect}} provides a summary of fit
#'   statistics for several models fitted to the same data.
#'
#' @examples
#' ## Fitting an example dose-response model
#' ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#'
#' ## Model-averaging with default settings (no confidence intervals)
#' maED(
#'   ryegrass.m1,
#'   list(LL.5(), LN.4(), W1.4(), W2.4(), FPL.4(-1, 1), FPL.4(-2, 3), FPL.4(-0.5, 0.5)),
#'   c(10, 50, 90)
#' )
#'
#' ## Model-averaging with Buckland confidence intervals
#' maED(
#'   ryegrass.m1,
#'   list(LL.5(), LN.4(), W1.4(), W2.4()),
#'   c(10, 50, 90),
#'   interval = "buckland"
#' )
#'
#' ## Model-averaging with Kang confidence intervals
#' maED(
#'   ryegrass.m1,
#'   list(LL.5(), LN.4(), W1.4(), W2.4()),
#'   c(10, 50, 90),
#'   interval = "kang"
#' )
#'
#' @keywords models nonlinear
#' @export
maED <- function(
    object,
    fctList  = NULL,
    respLev  = c(10,20,50),
    interval = c("none", "buckland", "kang"),
    linreg   = FALSE,
    clevel   = NULL,
    level    = 0.95,
    type     = c("relative", "absolute"),
    display  = TRUE,
    na.rm    = FALSE,
    extended = FALSE,
    onExclude = c("warn", "stop")
) {
  
  ## --- Input validation -------------------------------------------------------
  
  if (!inherits(object, "drc")) {
    stop("'object' must be of class 'drc'")
  }
  if (!is.numeric(respLev) || length(respLev) == 0) {
    stop("'respLev' must be a non-empty numeric vector")
  }
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("'level' must be a single numeric value strictly between 0 and 1")
  }
  if (!is.logical(linreg) || length(linreg) != 1) {
    stop("'linreg' must be a single logical value")
  }
  if (!is.logical(display) || length(display) != 1) {
    stop("'display' must be a single logical value")
  }
  if (!is.logical(na.rm) || length(na.rm) != 1) {
    stop("'na.rm' must be a single logical value")
  }
  if (!is.logical(extended) || length(extended) != 1) {
    stop("'extended' must be a single logical value")
  }
  
  ## --- Resolve enumerated arguments -------------------------------------------
  
  interval <- match.arg(interval)
  type     <- match.arg(type)
  onExclude <- match.arg(onExclude)
  
  .checkNoLogEDModels(object, fctList)
  
  ## --- Handling multiple curves in a single dataset ---------------------------
  
  # When the parameter matrix has more than one column (i.e., multiple curves)
  # and no specific curve has been requested, recurse over each curve
  # individually and bind the results.
  ncolPM <- ncol(object[["parmMat"]])
  
  if (!identical(ncolPM, 1L) && is.null(clevel)) {
    curveIds   <- colnames(object[["parmMat"]])
    resultList <- vector("list", ncolPM)
    
    # Per-curve results are always requested in the extended form so that the
    # model status table can be stacked for all curves.
    for (i in seq_len(ncolPM)) {
      resultList[[i]] <- maED(
        object   = object,
        fctList  = fctList,
        respLev  = respLev,
        interval = interval,
        linreg   = linreg,
        clevel   = curveIds[i],
        level    = level,
        type     = type,
        display  = display,
        na.rm    = na.rm,
        extended = TRUE,
        onExclude = onExclude
      )
    }
    
    return(.maEDCombineCurves(resultList, curveIds, extended))
  }
  
  ## --- Model selection summary ------------------------------------------------
  
  msMat <- do.call(mselect, list(object = object, fctList = fctList, sorted = "no"))
  
  ## --- Pre-allocate ED estimate and SE matrices -------------------------------
  
  lenfl   <- length(fctList)
  lenrl   <- length(respLev)
  numRows <- lenfl + 1L
  
  edEst <- matrix(NA, numRows + linreg, lenrl)
  edSe  <- matrix(NA, numRows + linreg, lenrl)
  
  # Track which models failed to fit (try-error in ED computation)
  fitFailed <- logical(numRows + linreg)
  
  # Confidence limit matrices are always initialised to avoid undefined
  # variable errors in the 'kang' result-construction block. They have one row
  # per model including the optional linear model, like edEst and the weights.
  edCll <- matrix(NA, numRows + linreg, lenrl)
  edClu <- matrix(NA, numRows + linreg, lenrl)
  
  ## --- Set interval argument for individual ED calls --------------------------
  
  # Delta-method intervals are required for the Kang approach; otherwise no
  # per-model interval is needed because Buckland uses the SE directly.
  interval2 <- if (identical(interval, "kang")) "delta" else "none"
  
  ## --- ED estimates for the original model ------------------------------------
  
  edMat      <- ED(object, respLev, interval2, clevel, level = level, type = type, display = FALSE)
  edEst[1, ] <- edMat[, 1]
  edSe[1, ]  <- edMat[, 2]
  
  if (identical(interval2, "delta")) {
    edCll[1, ] <- edMat[, 3]
    edClu[1, ] <- edMat[, 4]
  }
  
  ## --- ED estimates for each model in fctList ---------------------------------
  
  # A failing re-fit is either a genuine fit failure of this candidate (e.g. no
  # convergence: exclude the model) or a sign that the call stored in 'object'
  # cannot be re-evaluated at all (affects every candidate: stop). To tell them
  # apart, the original model is re-fitted once after the first failure.
  callError   <- NULL
  callChecked <- FALSE
  
  for (i in seq_len(lenfl)) {
    refit <- tryCatch(
      update(object, fct = fctList[[i]]),
      error = function(e) e
    )
    
    if (inherits(refit, "error") && !callChecked) {
      callChecked <- TRUE
      callError   <- .maEDCallError(object)
    }
    if (!is.null(callError)) {
      stop(
        "maED(): cannot re-fit the candidate models because the call stored in ",
        "'object' cannot be re-evaluated: ", conditionMessage(callError), ". ",
        "Refit 'object' so that its formula and data are available (see ?update.drc).",
        call. = FALSE
      )
    }
    
    edMati <- if (inherits(refit, "error")) {
      refit
    } else {
      try(
        ED(
          refit,
          respLev,
          interval2,
          clevel,
          level   = level,
          type    = type,
          display = FALSE
        ),
        silent = TRUE
      )
    }
    
    if (inherits(edMati, "try-error") || inherits(edMati, "error")) {
      edMati <- matrix(NA, length(respLev), 4)
      fitFailed[i + 1L] <- TRUE
    }
    
    edEst[i + 1L, ] <- edMati[, 1]
    edSe[i + 1L, ]  <- edMati[, 2]
    
    if (identical(interval2, "delta")) {
      edCll[i + 1L, ] <- edMati[, 3]
      edClu[i + 1L, ] <- edMati[, 4]
    }
  }
  
  ## --- Optional linear regression fit ----------------------------------------
  
  if (linreg) {
    linFit1              <- lm(object[["data"]][, 2:1])
    edLin                <- ED.lin(linFit1, respLev)
    edEst[lenfl + 2L, ] <- unlist(edLin[, 1])
    edSe[lenfl + 2L, ]  <- unlist(edLin[, 2])
    
    # Delta-method limits of the linear model at the requested level (for Kang)
    edLinCi <- confint.basic(
      cbind(unlist(edLin[, 1]), unlist(edLin[, 2])),
      level, "continuous", df.residual(linFit1), FALSE
    )
    edCll[lenfl + 2L, ] <- edLinCi[, 1]
    edClu[lenfl + 2L, ] <- edLinCi[, 2]
    
    # Include linear model AIC in the weight calculation.
    aicVec <- c(msMat[, 2], AIC(linFit1))
    expVec <- as.vector(exp(-aicVec / 2))
  } else {
    aicVec <- msMat[, 2]
    expVec <- as.vector(exp(-msMat[, 2] / 2))
  }
  
  modelNames <- if (linreg) c(rownames(msMat), "Lin") else rownames(msMat)
  
  ## --- Filter out models with non-finite ED estimates or fitting failures ----
  
  # Save original ED estimates for display (so excluded models still show
  # their Inf/NaN/NA values in the fit summary, making the reason for
  # exclusion visible).
  edEstDisplay <- edEst
  
  # Identify models where any ED estimate is non-finite (Inf or NaN).
  nonFiniteMask <- apply(edEst, 1, function(x) any(is.infinite(x) | is.nan(x)))
  
  # Combined exclusion mask: non-finite ED values OR model fitting failures.
  excludeMask <- nonFiniteMask | fitFailed
  
  modelStatus <- rep("used", length(excludeMask))
  modelStatus[nonFiniteMask] <- "excluded: non-finite ED"
  modelStatus[fitFailed]     <- "excluded: fit failed"
  
  if (any(excludeMask)) {
    reasons <- vapply(which(excludeMask), function(k) {
      if (fitFailed[k]) {
        "model fitting or ED estimation failed"
      } else {
        badIdx <- which(is.infinite(edEst[k, ]) | is.nan(edEst[k, ]))
        paste0(
          "non-finite ED value(s) detected (",
          paste0("ED", respLev[badIdx], "=", edEst[k, badIdx], collapse = ", "), ")"
        )
      }
    }, character(1))
    
    exclMsg <- paste0(
      sum(excludeMask), " of ", length(excludeMask), " models excluded from model-averaging: ",
      paste0("'", modelNames[excludeMask], "' (", reasons, ")", collapse = "; ")
    )
    if (identical(onExclude, "stop")) {
      stop(exclMsg, call. = FALSE)
    }
    warning(exclMsg, call. = FALSE)
    
    edEst[excludeMask, ] <- NA
    edSe[excludeMask, ]  <- NA
    edCll[excludeMask, ] <- NA
    edClu[excludeMask, ] <- NA
  }
  
  ## --- AIC-based model weights ------------------------------------------------
  
  # Excluded models (non-finite ED or fitting failures) always get zero
  # weight, regardless of the na.rm parameter.
  expVec[excludeMask] <- 0
  
  # When models were excluded, na.rm must be TRUE in downstream sums so
  # that the NA placeholders left above do not propagate.  For the
  # remaining (non-excluded) models, the user-supplied na.rm still governs
  # how NA values from fitting failures are handled.
  effectiveNaRm <- na.rm || any(excludeMask)
  
  totalWeight <- sum(expVec, na.rm = effectiveNaRm)
  
  if (totalWeight == 0) {
    warning("No valid models remaining for model-averaging", call. = FALSE)
    wVec  <- rep(0, length(expVec))
    edVec <- rep(NA_real_, lenrl)
  } else {
    wVec  <- expVec / totalWeight
    edVec <- apply(edEst * wVec, MARGIN = 2, FUN = sum, na.rm = effectiveNaRm)
  }
  
  ## --- Construct result matrix ------------------------------------------------
  
  if (identical(interval, "none")) {
    retMat <- as.matrix(cbind(edVec))
    colnames(retMat) <- colnames(edMat)[1]
    
  } else if (identical(interval, "buckland")) {
    if (totalWeight == 0) {
      seVec <- rep(NA_real_, lenrl)
    } else {
      seVec <- apply(
        sqrt(edSe^2 + (t(t(edEst) - edVec))^2) * wVec,
        MARGIN = 2,
        FUN    = sum,
        na.rm  = effectiveNaRm
      )
    }
    quantVal <- qnorm(1 - (1 - level) / 2) * seVec
    retMat   <- as.matrix(cbind(edVec, seVec, edVec - quantVal, edVec + quantVal))
    colnames(retMat) <- c(colnames(edMat)[c(1, 2)], "Lower", "Upper")
    
  } else {
    if (totalWeight == 0) {
      seVec  <- rep(NA_real_, lenrl)
      retMat <- matrix(NA_real_, lenrl, 4)
    } else {
      seVec <- apply(
        sqrt(edSe^2 + (t(t(edEst) - edVec))^2) * wVec,
        MARGIN = 2,
        FUN    = sum,
        na.rm  = effectiveNaRm
      )
      retMat <- as.matrix(cbind(
        apply(edEst * wVec, MARGIN = 2, FUN = sum, na.rm = effectiveNaRm),
        seVec,
        apply(edCll * wVec, MARGIN = 2, FUN = sum, na.rm = effectiveNaRm),
        apply(edClu * wVec, MARGIN = 2, FUN = sum, na.rm = effectiveNaRm)
      ))
    }
    colnames(retMat) <- c(colnames(edMat)[c(1, 2)], "Lower", "Upper")
  }
  
  rownames(retMat) <- rownames(edMat)
  
  ## --- Model status table -----------------------------------------------------
  
  modelsDf <- data.frame(
    model  = modelNames,
    AIC    = as.numeric(aicVec),
    weight = as.numeric(wVec),
    status = modelStatus,
    stringsAsFactors = FALSE
  )
  
  # Attributes are only attached if models were excluded, to keep the printed
  # output of regular results unchanged.
  if (any(excludeMask)) {
    attr(retMat, "models")      <- modelsDf
    attr(retMat, "n_used")      <- sum(modelStatus == "used")
    attr(retMat, "n_requested") <- length(modelStatus)
  }  
  ## --- Construct fit summary matrix -------------------------------------------
  
  # Use original (unfiltered) ED estimates for display so that excluded
  # models show their Inf/NaN values alongside their zero weight.
  disMat           <- as.matrix(cbind(edEstDisplay, wVec))
  colnames(disMat) <- c(paste0("ED", respLev), "Weight")
  rownames(disMat) <- if (linreg) c(rownames(msMat), "Lin") else rownames(msMat)
  
  ## --- Optional display -------------------------------------------------------
  
  if (display) {
    print(disMat)
    cat("\n")
  }
  
  ## --- Return -----------------------------------------------------------------
  
  if (extended) {
    return(list(estimates = retMat, fits = disMat, models = modelsDf))
  }
  
  return(retMat)
}

## Re-fits the original model from its stored call. Returns NULL if that works
## and the error condition if the call cannot be re-evaluated.
.maEDCallError <- function(object) {
  res <- tryCatch(suppressWarnings(update(object)), error = function(e) e)
  if (inherits(res, "error")) res else NULL
}

## Assembles the per-curve results (each from maED(extended = TRUE)) of a
## multi-curve fit into the usual return value.
.maEDCombineCurves <- function(resultList, curveIds, extended) {
  estimates <- do.call(rbind, lapply(resultList, function(r) {
    est <- r$estimates
    attributes(est) <- attributes(est)[c("dim", "dimnames")]
    est
  }))
  
  models <- do.call(rbind, Map(function(r, id) {
    cbind(curve = id, r$models, stringsAsFactors = FALSE)
  }, resultList, curveIds))
  
  if (any(models$status != "used")) {
    nUsed <- vapply(curveIds, function(id) sum(models$status[models$curve == id] == "used"), numeric(1))
    nReq  <- vapply(curveIds, function(id) sum(models$curve == id), numeric(1))
    attr(estimates, "models")      <- models
    attr(estimates, "n_used")      <- nUsed
    attr(estimates, "n_requested") <- nReq
  }
  
  if (!extended) return(estimates)
  
  fits <- do.call(rbind, Map(function(r, id) {
    f <- r$fits
    rownames(f) <- paste(id, rownames(f), sep = ":")
    f
  }, resultList, curveIds))
  
  list(estimates = estimates, fits = fits, models = models)
}
