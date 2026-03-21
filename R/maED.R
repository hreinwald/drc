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
#' @param fctList a list of non-linear functions to be compared.
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
#'   between 0 and 1. The default is \code{0.95}.
#' @param type character string. Whether the specified response levels are
#'   absolute or relative (default).
#' @param display logical. If \code{TRUE} results are displayed. Otherwise they
#'   are not (useful in simulations).
#' @param na.rm logical indicating whether or not \code{NA} values occurring
#'   during model fitting should be excluded from subsequent calculations.
#' @param extended logical specifying whether or not an extended output
#'   (including fit summaries) should be returned.
#'
#' @return If \code{extended = FALSE}, a matrix with two or more columns
#'   containing the model-averaged estimates and the corresponding estimated
#'   standard errors and, optionally, lower and upper confidence limits.
#'   If \code{extended = TRUE}, a list with components:
#'   \describe{
#'     \item{estimates}{Matrix of model-averaged ED estimates and intervals.}
#'     \item{fits}{Matrix of per-model ED estimates and AIC-based weights.}
#'   }
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
    extended = FALSE
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
  
  ## --- Handling multiple curves in a single dataset ---------------------------
  
  # When the parameter matrix has more than one column (i.e., multiple curves)
  # and no specific curve has been requested, recurse over each curve
  # individually and bind the results.
  ncolPM <- ncol(object[["parmMat"]])
  
  if (!identical(ncolPM, 1L) && is.null(clevel)) {
    curveIds   <- colnames(object[["parmMat"]])
    resultList <- vector("list", ncolPM)
    
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
        extended = extended
      )
    }
    
    return(do.call(rbind, resultList))
  }
  
  ## --- Model selection summary ------------------------------------------------
  
  msMat <- do.call(mselect, list(object = object, fctList = fctList, sorted = "no"))
  
  ## --- Pre-allocate ED estimate and SE matrices -------------------------------
  
  lenfl   <- length(fctList)
  lenrl   <- length(respLev)
  numRows <- lenfl + 1L
  
  edEst <- matrix(NA, numRows + linreg, lenrl)
  edSe  <- matrix(NA, numRows + linreg, lenrl)
  
  # Confidence limit matrices are always initialised to avoid undefined
  # variable errors in the 'kang' result-construction block.
  edCll <- matrix(NA, numRows, lenrl)
  edClu <- matrix(NA, numRows, lenrl)
  
  ## --- Set interval argument for individual ED calls --------------------------
  
  # Delta-method intervals are required for the Kang approach; otherwise no
  # per-model interval is needed because Buckland uses the SE directly.
  interval2 <- if (identical(interval, "kang")) "delta" else "none"
  
  ## --- ED estimates for the original model ------------------------------------
  
  edMat      <- ED(object, respLev, interval2, clevel, type = type, display = FALSE)
  edEst[1, ] <- edMat[, 1]
  edSe[1, ]  <- edMat[, 2]
  
  if (identical(interval2, "delta")) {
    edCll[1, ] <- edMat[, 3]
    edClu[1, ] <- edMat[, 4]
  }
  
  ## --- ED estimates for each model in fctList ---------------------------------
  
  for (i in seq_len(lenfl)) {
    edMati <- try(
      ED(
        update(object, fct = fctList[[i]]),
        respLev,
        interval2,
        clevel,
        type    = type,
        display = FALSE
      ),
      silent = TRUE
    )
    
    if (inherits(edMati, "try-error")) {
      edMati <- matrix(NA, length(respLev), 4)
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
    
    # Include linear model AIC in the weight calculation.
    expVec <- as.vector(exp(-c(msMat[, 2], AIC(linFit1)) / 2))
  } else {
    expVec <- as.vector(exp(-msMat[, 2] / 2))
  }
  
  ## --- Filter out models with non-finite ED estimates -----------------------
  
  # Save original ED estimates for display (so excluded models still show
  # their Inf/NaN values in the fit summary, making the reason for exclusion
  # visible).
  edEstDisplay <- edEst
  
  # Identify models where any ED estimate is non-finite (Inf or NaN).
  # NA values from model fitting failures are NOT flagged here; those are
  # governed by the 'na.rm' parameter instead.
  excludeMask <- apply(edEst, 1, function(x) any(is.infinite(x) | is.nan(x)))
  
  if (any(excludeMask)) {
    modelNames <- if (linreg) c(rownames(msMat), "Lin") else rownames(msMat)
    for (k in which(excludeMask)) {
      badIdx <- which(is.infinite(edEst[k, ]) | is.nan(edEst[k, ]))
      warning(
        "Model '", modelNames[k], "' excluded from model-averaging: ",
        "non-finite ED value(s) detected (",
        paste0("ED", respLev[badIdx], "=", edEst[k, badIdx], collapse = ", "), ")",
        call. = FALSE
      )
    }
    edEst[excludeMask, ] <- NA
    edSe[excludeMask, ]  <- NA
    excludeCI <- excludeMask[seq_len(numRows)]
    edCll[excludeCI, ] <- NA
    edClu[excludeCI, ] <- NA
  }
  
  ## --- AIC-based model weights ------------------------------------------------
  
  # Excluded models (non-finite ED) always get zero weight, regardless of
  # the na.rm parameter.
  expVec[excludeMask] <- 0
  
  # When models were excluded, na.rm must be TRUE in downstream sums so
  # that the NA placeholders left above do not propagate.  For the
  # remaining (non-excluded) models, the user-supplied na.rm still governs
  # how NA values from fitting failures are handled.
  effectiveNaRm <- na.rm || any(excludeMask)
  
  wVec  <- expVec / sum(expVec, na.rm = effectiveNaRm)
  edVec <- apply(edEst * wVec, MARGIN = 2, FUN = sum, na.rm = effectiveNaRm)
  
  ## --- Construct result matrix ------------------------------------------------
  
  if (identical(interval, "none")) {
    retMat <- as.matrix(cbind(edVec))
    colnames(retMat) <- colnames(edMat)[1]
    
  } else if (identical(interval, "buckland")) {
    seVec <- apply(
      sqrt(edSe^2 + (t(t(edEst) - apply(edEst, MARGIN = 2, FUN = mean, na.rm = effectiveNaRm)))^2) * wVec,
      MARGIN = 2,
      FUN    = sum,
      na.rm  = effectiveNaRm
    )
    quantVal <- qnorm(1 - (1 - level) / 2) * seVec
    retMat   <- as.matrix(cbind(edVec, seVec, edVec - quantVal, edVec + quantVal))
    colnames(retMat) <- c(colnames(edMat)[c(1, 2)], "Lower", "Upper")
    
  } else {
    retMat <- as.matrix(cbind(
      apply(edEst * wVec, MARGIN = 2, FUN = sum, na.rm = effectiveNaRm),
      apply(edCll * wVec, MARGIN = 2, FUN = sum, na.rm = effectiveNaRm),
      apply(edClu * wVec, MARGIN = 2, FUN = sum, na.rm = effectiveNaRm)
    ))
    colnames(retMat) <- colnames(edMat)[c(1, 3, 4)]
  }
  
  rownames(retMat) <- rownames(edMat)
  
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
    return(list(estimates = retMat, fits = disMat))
  }
  
  return(retMat)
}
