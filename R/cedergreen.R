# HELPER -------------------------------------------------------------------

# Define the edfct and maxfct helper functions here. Needed for the cedergreen function to work properly. 
# These functions are used for calculating effective doses and finding the maximum hormesis, respectively.

## -- Defining the edfct (Refactored) -- ##
#' Calculate Effective Dose for the Cedergreen-Ritz Hormesis Model
#'
#' @description
#' An internal helper function to calculate the effective dose (ED) and its
#' derivatives for the Cedergreen-Ritz five-parameter hormesis model. It uses
#' `uniroot` to find the dose for a given response level.
#'
#' @param parm A numeric vector of the non-fixed model parameters.
#' @param all_params A numeric vector template for all model parameters (b,c,d,e,f).
#' @param not_fixed A logical or integer vector indicating the non-fixed parameters.
#' @param alpha A numeric value for the hormesis model's alpha shape parameter.
#' @param respl The response level to calculate the dose for (e.g., 50 for ED50).
#' @param reference A character string ("control" or "absolute") for calculating the response.
#' @param type A character string specifying the type of ED calculation.
#' @param lower The lower bound of the dose interval for the root-finding search.
#' @param upper The upper bound of the dose interval for the root-finding search.
#'
#' @return A list containing the calculated effective dose and a vector of its
#'   partial derivatives with respect to the non-fixed parameters.
#'   
#' @author Hannes Reinwald
#' 
#' @keywords internal
#'
cedergreen_edfct <- function( 
    parm,         # Vector of non-fixed parameters
    all_params,   # Full parameter vector (template)
    not_fixed,    # Index/logical of non-fixed params
    alpha,        # Hormesis shape parameter
    respl,        # Response level (e.g., 50)
    reference,    # Reference for EDhelper
    type,         # Type for EDhelper
    lower = 1e-4,
    upper = 10000
){
  # 1. Self-contained: Reconstruct the full parameter vector
  all_params[not_fixed] <- parm
  
  # 2. Readability: Use named parameters
  p_named <- list(b = all_params[1], c = all_params[2], d = all_params[3], 
                  e = all_params[4], f = all_params[5])
  
  # Calculate the target response proportion
  p_percent <- EDhelper(all_params, respl, reference, type, TRUE)
  target_prop <- (100 - p_percent) / 100
  
  # Define the dose-response model with clear names
  response_model <- function(dose, p, alpha) {
    p$c + (p$d - p$c + p$f * exp(-1 / (dose^alpha))) / (1 + exp(p$b * (log(dose) - log(p$e))))
  }
  
  # Find the dose at which the maximum response occurs to handle non-monotonicity
  dose_grid <- exp(seq(log(lower), log(upper), length.out = 1000))
  response_values <- response_model(dose_grid, p_named, alpha)
  dose_at_max_resp <- dose_grid[which.max(response_values)]
  
  # Define the equation to solve for the effective dose (ED)
  root_eqn <- function(dose) {
    # Rearranged model: F(dose, params) = 0
    target_prop * (1 + exp(p_named$b * (log(dose) - log(p_named$e)))) - 
      (1 + p_named$f * exp(-1 / (dose^alpha)) / (p_named$d - p_named$c))
  }
  
  # 3. Robustness: Solve for ED with error handling
  effective_dose <- tryCatch({
    uniroot(root_eqn, lower = dose_at_max_resp, upper = upper)$root
  }, error = function(e) {
    warning(paste("Root finding failed for ED", respl, ". Returning NA.", sep=""))
    return(NA)
  })
  
  if (is.na(effective_dose)) {
    return(list(NA, rep(NA, sum(not_fixed))))
  }
  
  # 4. Clarity in Derivatives: Calculate derivatives for the Delta Method
  # (This part remains complex, but breaking it down would be the next step)  
  # Note: The derivative calculation is kept brief here for demonstration.
  # In a real refactoring, each term of derParm and derDose would be calculated separately.
  tempVal1 <- exp(p_named$b * (log(effective_dose) - log(p_named$e)))
  tempVal2 <- p_named$d - p_named$c
  
  derParm <- c(target_prop*tempVal1*(log(effective_dose)-log(p_named$e)), 
               -p_named$f*exp(-1/(effective_dose^alpha))/((tempVal2)^2),
               p_named$f*exp(-1/(effective_dose^alpha))/((tempVal2)^2), 
               -target_prop*tempVal1*p_named$b/p_named$e,
               -exp(-1/(effective_dose^alpha))/tempVal2)
  
  derDose <- target_prop*tempVal1*p_named$b/effective_dose - p_named$f/tempVal2*exp(-1/(effective_dose^alpha))/(effective_dose^(1+alpha))*alpha 
  
  # Correct application of Implicit Function Theorem
  ed_derivatives <- -derParm / derDose
  return(list(effective_dose, ed_derivatives[not_fixed]))
}


## -- Defining the maxfct function (Refactored) -- ##
#' Find the Dose and Response at Maximum Hormesis
#'
#' @description
#' This function finds the dose that elicits the maximum hormetic (stimulatory)
#' response for the Cedergreen-Ritz model and the response value at that dose.
#'
#' @param all_params A named list of all model parameters (b, c, d, e, f).
#' @param alpha The hormesis alpha shape parameter.
#' @param lower The lower bound of the dose interval to search for the maximum.
#' @param upper The upper bound of the dose interval to search for the maximum.
#'
#' @return A numeric vector containing two values: the dose at the maximum
#'   response, and the maximum response value itself. Returns `c(NA, NA)` on failure.
#'   
#' @author Hannes Reinwald
#'
#' @keywords internal
#'
#' @importFrom stats optimize
cedergreen_maxfct <- function(all_params, alpha, lower = 1e-6, upper = 1000, .optimize_fn = stats::optimize)
{
  # Define the dose-response model using named parameters for clarity
  response_model <- function(dose, p, alpha) {
    p$c + (p$d - p$c + p$f * exp(-1 / (dose^alpha))) / (1 + exp(p$b * (log(dose) - log(p$e))))
  }

  # Use optimize() to directly find the dose that maximizes the response.
  # It is more robust than finding the root of the derivative.
  # We search for the maximum by telling optimize to maximize=TRUE.
  opt_result <- tryCatch({
    .optimize_fn(
      f = response_model,
      interval = c(lower, upper),
      p = all_params,
      alpha = alpha,
      maximum = TRUE
    )
  }, error = function(e) {
    warning("Optimization failed to find a maximum hormesis dose.")
    return(NULL)
  })

  if (is.null(opt_result)) {
    return(c(maxDose = NA, maxResponse = NA))
  }

  # Return the dose at the maximum and the value of the function at that maximum
  return(c(maxDose = opt_result$maximum, maxResponse = opt_result$objective))
}


# MAIN ---------------------------------------------------------------------

#' @title Cedergreen-Ritz-Streibig Model
#' @description Provides the Cedergreen-Ritz-Streibig function, a five-parameter model 
#'   for describing dose-response curves that exhibit hormesis (a stimulatory or 
#'   beneficial effect at low doses). This function generates a model object suitable 
#'   for use with non-linear regression functions like \code{\link[drc]{drm}}.
#'
#' @details 
#' The Cedergreen-Ritz-Streibig model is defined by the following equation:
#' \deqn{f(x) = c + \frac{d - c + f \exp(-1/x^{\alpha})}{1 + \exp(b(\log(x) - \log(e)))}}
#' The parameter \eqn{f} determines the size of the hormetic effect (stimulation). 
#' If \eqn{f=0}, the model simplifies to the standard four-parameter log-logistic model.
#' The parameter \eqn{\alpha} is a shape parameter that must be specified by the user.
#'
#' @param fixed A numeric vector of length 5 specifying any parameters to be held fixed 
#'   during the estimation. The order is \code{c(b, c, d, e, f)}. Use \code{NA} for 
#'   parameters that should be estimated. The default is to estimate all parameters.
#' @param names A character vector of length 5 providing names for the parameters. 
#'   The default is \code{c("b", "c", "d", "e", "f")}.
#' @param method A character string specifying the method for the self-starter function 
#'   to use for finding initial parameter values. Options are \code{"loglinear"}, 
#'   \code{"anke"}, \code{"method3"}, and \code{"normolle"}. This is only used if \code{ssfct} is \code{NULL}.
#' @param ssfct A custom self-starter function. If \code{NULL} (the default), a 
#'   self-starter is automatically generated by calling \code{\link{cedergreen.ssf}} 
#'   with the specified \code{method}, \code{fixed}, and \code{alpha} arguments.
#' @param alpha A mandatory numeric value specifying the fixed shape parameter \eqn{\alpha}. 
#'   The function will stop if this is not provided.
#' @param fctName An optional character string to name the function object.
#' @param fctText An optional character string providing a descriptive text for the model.
#'
#' @return A list of class \code{mllogistic}, containing the model function (\code{fct}), 
#'   the self-starter function (\code{ssfct}), parameter names (\code{names}), and other 
#'   components required for use with modeling functions like \code{\link[drc]{drm}}.
#'
#' @seealso \code{\link[drc]{drm}} for model fitting, and \code{\link{cedergreen.ssf}} for the 
#'   underlying self-starter function.
#'   
#' @author Christian Ritz, Hannes Reinwald
#' 
#' @keywords models nonlinear
#'
#' @examples
#' dose <- c(0, 0.1, 0.5, 1, 5, 10, 20)
#' response <- c(100, 102, 95, 80, 40, 25, 20)
#' my_data <- data.frame(dose = dose, response = response)
#' model_fit <- drm(response ~ dose, data = my_data,
#'                  fct = cedergreen(alpha = 0.5))
#' summary(model_fit)
#'
#' @export
"cedergreen" <- function(
    fixed  = c(NA, NA, NA, NA, NA), 
    names  = c("b", "c", "d", "e", "f"), 
    method = c("loglinear", "anke", "method3", "normolle"), 
    ssfct  = NULL, 
    alpha, 
    fctName, 
    fctText 
){
  ## Checking arguments and setting up fixed parameter logic
  numParm <- 5
  if (!is.character(names) || !(length(names) == numParm)) {stop("Not correct 'names' argument")}
  if (!(length(fixed) == numParm)) {stop("Not correct 'fixed' argument")}    
  if (missing(alpha)) {stop("'alpha' argument must be specified")}
  
  # Determine if fixed parameters are being used. This will be passed to ssfct.
  useFixed <- !all(is.na(fixed))
  
  # Match the method argument
  method   = match.arg(method)
  notFixed = is.na(fixed)
  parmVec  = rep(0, numParm)
  parmVec[!notFixed] <- fixed[!notFixed]
  
  ## Defining the non-linear function
  fct <- function(dose, parm){
    parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
    parmMat[, notFixed] <- parm
    # The Cedergreen-Ritz-Streibig model equation
    parmMat[,2] + (parmMat[,3] - parmMat[,2] + parmMat[,5] * exp(-1 / (dose^alpha))) / 
      (1 + exp(parmMat[,1] * (log(dose) - log(parmMat[,4]))))
  }
  
  ## --- Correctly defining the self-starter ---
  # If a custom self-starter is not provided, use our robust cedergreen.ssf
  if (is.null(ssfct)) {
    # This is the correct way to call the external self-starter. 
    # It passes the method, the fixed vector, alpha, and the useFixed flag.
    ssfct <- cedergreen.ssf(method = method, fixed = fixed, alpha = alpha, useFixed = useFixed)
  }
  
  ## Defining names for the parameters to be estimated
  names <- names[notFixed]
  
  ## Specifying the derivatives (derivatives code is kept as is)
  deriv1 <- function(dose, parm)
  {
    parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
    parmMat[, notFixed] <- parm
    
    t0 = exp(-1/(dose^alpha))
    t1 = parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*t0
    t2 = exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
    t3 = 1 + t2                          
    t4 = (1 + t2)^(-2)
    
    # A helper function 'xlogx' would need to be defined for this to work
    cbind(
      -t1*xlogx(dose/parmMat[, 4], parmMat[, 1])*t4, 
      1 - 1/t3, 
      1/t3, 
      1*t2*(parmMat[, 1]/parmMat[, 4])*t4, 
      t0/t3
    )[, notFixed]
  }
  
  
  ## Defining the ED function: wrapper closure that delegates to cedergreen_edfct
  ## The framework calls edfct(parm, respl, reference, type, ...) so we bind
  ## parmVec, notFixed, and alpha from the enclosing scope.
  edfct <- function(parm, respl, reference, type, lower = 1e-4, upper = 10000, ...)
  {
    cedergreen_edfct(parm, parmVec, notFixed, alpha, respl, reference, type, lower, upper)
  }
  
  ## Finding the maximal hormesis: wrapper closure that delegates to cedergreen_maxfct
  ## The framework calls maxfct(parm, lower, upper) where parm is the non-fixed
  ## parameter vector. We reconstruct the full named-list and bind alpha.
  maxfct <- function(parm, lower = 1e-3, upper = 1000, .optimize_fn = stats::optimize)
  {
    parmVec[notFixed] <- parm
    all_params <- list(b = parmVec[1], c = parmVec[2], d = parmVec[3],
                       e = parmVec[4], f = parmVec[5])
    cedergreen_maxfct(all_params, alpha, lower, upper, .optimize_fn)
  }
  
  # Return results
  returnList <- list(
    fct    = fct, 
    ssfct  = ssfct, 
    names  = names, 
    deriv1 = deriv1, # Note: deriv1 is incomplete without xlogx
    deriv2 = NULL,
    edfct  = edfct,
    maxfct = maxfct,
    name   = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
    text   = ifelse(missing(fctText), "Cedergreen-Ritz-Streibig", fctText),     
    noParm = sum(is.na(fixed))
  )
  class(returnList) <- "mllogistic"
  invisible(returnList)
}


# WRAPPER -------------------------------------------------------------------

#' Wrapper for 5-parameter Cedergreen-Ritz-Streibig Model
#' 
#' @author Hannes Reinwald
#'
#' @description
#' A convenience wrapper for the \code{drc::cedergreen} function, preset for a 
#' 5-parameter model. It provides flexible handling for the alpha parameter.
#'
#' @details
#' This function simplifies the creation of a 5-parameter Cedergreen-Ritz-Streibig 
#' model by setting sensible defaults for the parameter names. It allows the 
#' alpha parameter to be specified either by a predefined character shortcut 
#' ('a', 'b', 'c') or by a direct numeric value.
#' 
#' By default the function runs with `alpha=1`, which corresponds to the `CRS.4a` model. 
#' Setting `alpha=0.5` corresponds to the `CRS.4b` model, and `alpha=0.25` corresponds to the `CRS.4c` model.
#' 
#' By default, all parameters are set to be estimated (i.e., \code{fixed} is all \code{NA}), 
#' but users can specify any parameters to be held constant during estimation. 
#' The self-starter function is automatically generated based on the specified method and 
#' fixed parameters, ensuring that initial values are appropriately calculated for the model fitting process.
#'
#' The function automatically generates a model name (`fctName`) and description 
#' (`fctText`) unless they are explicitly provided by the user.
#'
#' @param names A character vector of length 5 specifying the names of the model 
#'   parameters. Default is \code{c("b", "c", "d", "e", "f")}.
#' @param fixed A numeric vector of length 5. Use \code{NA} for parameters to be 
#'   estimated and a numeric value for parameters to be fixed. Default is all 
#'   \code{NA}.
#' @param alpha_type A character or a numeric value. Can be one of 'a' (alpha=1), 
#'   'b' (alpha=0.5), 'c' (alpha=0.25), or a specific numeric value for alpha.
#' @param fctName An optional character string to name the model function. If 
#'   \code{NULL} (the default), a name is generated automatically.
#' @param fctText An optional character string describing the model. If 
#'   \code{NULL} (the default), a description is generated automatically.
#' @param ... Additional arguments to be passed to \code{drc::cedergreen}, such 
#'   as \code{data}.
#'
#' @return A \code{drc} model object of class \code{cedergreen}. If the underlying 
#'   \code{drc::cedergreen} call fails, it issues a warning and returns \code{NULL}.
#'
#' @examples
#' # Create a CRS.5 model specification
#' crs_model_a <- CRS.5()
#'
#' # Fix the lower limit to 0 and use a custom numeric alpha
#' crs_model_custom <- CRS.5(
#'   fixed = c(NA, 0, NA, NA, NA), alpha_type = 0.75
#' )
#' 
#' @export
CRS.5 = function(names = c("b", "c", "d", "e", "f"), 
                 fixed = c(NA, NA, NA, NA, NA), 
                 alpha_type = "a", # one of 'a', 'b' or 'c' or a numeric value specifiying alpha
                 fctName = NULL,
                 fctText = NULL,
                 ... ){
  
  # Input sanity check
  if (!is.character(names) | !(length(names) == 5)) {
    stop("Not correct 'names' argument")
  }
  
  # Specify respective alpha values 
  alpha_value = list(a = 1, b = 0.5, c = 0.25)
  
  # Check alpha_type and calculate the numeric alpha value 'a'
  if(!is.numeric(alpha_type)){
    # --- STABILITY IMPROVEMENT ---
    # Check if the provided character is a valid key before trying to access it
    if (!(alpha_type %in% names(alpha_value))) {
      stop("Invalid 'alpha_type'. Must be one of 'a', 'b', 'c', or a numeric value.")
    }
    a = alpha_value[[alpha_type]]
  } else {
    a = alpha_type
    # Modify alpha_type for unique function naming if a numeric is given
    alpha_type = paste0(".alpha:", alpha_type)
  }
  stopifnot(is.numeric(a)) # This check is now safe
  
  # Propper variable handling for optional arguments
  if(is.null(fctName)) fctName = paste0(as.character(match.call()[[1]]), alpha_type)
  if(is.null(fctText)) fctText = paste0("Cedergreen-Ritz-Streibig (alpha=", a, ")")
  
  # Run cedergreen() within tryCatch for proper error handling
  res = tryCatch({
    drc::cedergreen(fixed = fixed, names = names, alpha = a, 
                    fctName = fctName, fctText = fctText, ...)
  }, error = function(e) {
    warning(paste("The cedergreen() model call failed with an error:", conditionMessage(e)))
    return(NULL)
  })
  
  return(res)
}


# DEPRECATED -----------------------------------------------------------

## 4 Parametric --------------------------------------------------------

#' @title Cedergreen-Ritz-Streibig Model with Lower Limit Fixed at 0 and Alpha = 1
#' (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated as of version 3.3.0. Please use [CRS.5()] instead,
#' which provides a more general and flexible interface.
#'
#' A four-parameter Cedergreen-Ritz-Streibig (CRS) hormesis model where the lower
#' asymptote (`c`) is fixed at 0 and the alpha parameter controlling the steepness
#' of the hormetic component is fixed at 1. The four free parameters are `b`, `d`,
#' `e`, and `f`.
#'
#' @param names A character vector of length 5 specifying the names of the model
#'   parameters in the following order:
#'   \describe{
#'     \item{`b`}{Hill slope (steepness of the dose-response curve).}
#'     \item{`c`}{Lower asymptote (fixed at 0 via the `fixed` argument).}
#'     \item{`d`}{Upper asymptote.}
#'     \item{`e`}{Effective dose producing a response midway between `c` and `d`
#'       (ED50).}
#'     \item{`f`}{Hormesis parameter controlling the magnitude of the stimulatory
#'       effect at low doses.}
#'   }
#'   Defaults to `c("b", "c", "d", "e", "f")`.
#'
#' @param fixed A numeric vector of length 5 specifying fixed (non-estimated)
#'   parameter values. Use `NA` for parameters that should be estimated freely.
#'   Defaults to `c(NA, 0, NA, NA, NA)`, which fixes the lower asymptote `c`
#'   at 0.
#'
#' @param ... Additional arguments passed to [cedergreen()].
#'
#' @return A list of class `"drcMean"` as returned by [cedergreen()], containing
#'   the model definition including the mean function, its gradient, parameter
#'   names, and fixed values. This object is intended for use as the `fct`
#'   argument in [drm()].
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [CRS.5()] — the recommended replacement for this deprecated function.
#' * [cedergreen()] — the underlying model constructor.
#' * [CRS.5a()] — the five-parameter CRS model with alpha = 1.
#' * [UCRS.4a()] — the unconstrained four-parameter CRS model with alpha = 1.
#'
#' @examples
#' # NOTE: CRS.4a() is deprecated. Use CRS.5() instead.
#' # The example below is retained for backward compatibility illustration only.
#'
#' lettuce.crsm1 <- drm( lettuce[, c(2, 1)], fct = CRS.4a() )
#' summary(lettuce.crsm1)
#' ED(lettuce.crsm1, c(50))
#'
#' # Recommended replacement:
#' fct_spec <- CRS.5(alpha_type = "a", fixed = c(NA, 0, NA, NA, NA))
#' lettuce.crs5 <- drm(lettuce[, c(2, 1)], fct = fct_spec)
#' summary(lettuce.crs5)
#' ED(lettuce.crs5, c(50))
#'
#' @keywords models nonlinear
#' @export
"CRS.4a" <- function(
    names = c("b", "c", "d", "e", "f"), 
    fixed = c(NA, 0, NA, NA, NA), # fixed c = 0
    ... ){
  
  # Deprecated warning
  lifecycle::deprecate_warn(
    when    = "3.3.0",
    what    = "CRS.4a()",
    with    = "CRS.5()"
  )
  
  return(
    cedergreen(
      fixed = fixed, names = names, alpha = 1, 
      fctName = as.character(match.call()[[1]]), 
      fctText = "Cedergreen-Ritz-Streibig with lower limit 0 (alpha=1)", 
      ...
    )
  )
}


#' @title Alias for CRS.4a (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is a deprecated alias for [CRS.4a()], itself deprecated as of
#' version 3.3.0. Please use [CRS.5()] instead, which provides a more general
#' and flexible interface.
#'
#' @inherit CRS.4a params return
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [CRS.5()] — the recommended replacement for this deprecated function.
#' * [CRS.4a()] — the function this alias points to.
#' * [cedergreen()] — the underlying model constructor.
#'
#' @examples
#' # NOTE: ml3a() is a deprecated alias for CRS.4a(). Use CRS.5() instead.
#' # The example below is retained for backward compatibility illustration only.
#'
#' lettuce.crsm1 <- drm( lettuce[, c(2, 1)], fct = ml3a() )
#' summary(lettuce.crsm1)
#' ED(lettuce.crsm1, c(50))
#'
#' # Recommended replacement:
#' fct_spec <- CRS.5(alpha_type = "a", fixed = c(NA, 0, NA, NA, NA))
#' lettuce.crs5 <- drm(lettuce[, c(2, 1)], fct = fct_spec)
#' summary(lettuce.crs5)
#' ED(lettuce.crs5, c(50))
#'
#' @keywords models nonlinear
#' @export
ml3a <- CRS.4a


#' @title Cedergreen-Ritz-Streibig Model with Lower Limit Fixed at 0 and Alpha = 0.5
#' (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated as of version 3.3.0. Please use [CRS.5()] instead,
#' which provides a more general and flexible interface.
#'
#' A four-parameter Cedergreen-Ritz-Streibig (CRS) hormesis model where the lower
#' asymptote (`c`) is fixed at 0 and the alpha parameter controlling the steepness
#' of the hormetic component is fixed at 0.5. The four free parameters are `b`, `d`,
#' `e`, and `f`.
#'
#' @param names A character vector of length 5 specifying the names of the model
#'   parameters in the following order:
#'   \describe{
#'     \item{`b`}{Hill slope (steepness of the dose-response curve).}
#'     \item{`c`}{Lower asymptote (fixed at 0 via the `fixed` argument).}
#'     \item{`d`}{Upper asymptote.}
#'     \item{`e`}{Effective dose producing a response midway between `c` and `d`
#'       (ED50).}
#'     \item{`f`}{Hormesis parameter controlling the magnitude of the stimulatory
#'       effect at low doses.}
#'   }
#'   Defaults to `c("b", "c", "d", "e", "f")`.
#'
#' @param fixed A numeric vector of length 5 specifying fixed (non-estimated)
#'   parameter values. Use `NA` for parameters that should be estimated freely.
#'   Defaults to `c(NA, 0, NA, NA, NA)`, which fixes the lower asymptote `c`
#'   at 0.
#'
#' @param ... Additional arguments passed to [cedergreen()].
#'
#' @return A list of class `"drcMean"` as returned by [cedergreen()], containing
#'   the model definition including the mean function, its gradient, parameter
#'   names, and fixed values. This object is intended for use as the `fct`
#'   argument in [drm()].
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [CRS.5()] — the recommended replacement for this deprecated function.
#' * [cedergreen()] — the underlying model constructor.
#' * [CRS.4a()] — the four-parameter CRS model with lower limit fixed at 0 and alpha = 1.
#' * [CRS.5b()] — the five-parameter CRS model with alpha = 0.5.
#'
#' @examples
#' # NOTE: CRS.4b() is deprecated. Use CRS.5() instead.
#' # The example below is retained for backward compatibility illustration only.
#'
#' lettuce.crsm2 <- drm( lettuce[, c(2, 1)], fct = CRS.4b() )
#' summary(lettuce.crsm2)
#' ED(lettuce.crsm2, c(50))
#'
#' # Recommended replacement:
#' fct_spec <- CRS.5(alpha_type = "b", fixed = c(NA, 0, NA, NA, NA))
#' lettuce.crs5 <- drm(lettuce[, c(2, 1)], fct = fct_spec)
#' summary(lettuce.crs5)
#' ED(lettuce.crs5, c(50))
#'
#' @keywords models nonlinear
#' @export
"CRS.4b" <- function(
    names = c("b", "c", "d", "e", "f"), 
    fixed = c(NA, 0, NA, NA, NA), # fixed c = 0
    ... ){
  
  # Deprecated warning
  lifecycle::deprecate_warn(
    when    = "3.3.0",
    what    = "CRS.4b()",
    with    = "CRS.5()"
  )
  
  return(
    cedergreen(
      fixed = fixed, names = names, alpha = 0.5, 
      fctName = as.character(match.call()[[1]]), 
      fctText = "Cedergreen-Ritz-Streibig with lower limit 0 (alpha=0.5)", 
      ...
    )
  )
}

#' @title Alias for CRS.4b (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is a deprecated alias for [CRS.4b()], itself deprecated as of
#' version 3.3.0. Please use [CRS.5()] instead, which provides a more general
#' and flexible interface.
#'
#' @inherit CRS.4b params return
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [CRS.5()] — the recommended replacement for this deprecated function.
#' * [CRS.4b()] — the function this alias points to.
#' * [cedergreen()] — the underlying model constructor.
#'
#' @examples
#' # NOTE: ml3b() is a deprecated alias for CRS.4b(). Use CRS.5() instead.
#' # The example below is retained for backward compatibility illustration only.
#'
#' lettuce.crsm2 <- drm( lettuce[, c(2, 1)], fct = ml3b() )
#' summary(lettuce.crsm2)
#' ED(lettuce.crsm2, c(50))
#'
#' # Recommended replacement:
#' fct_spec <- CRS.5(alpha_type = "b", fixed = c(NA, 0, NA, NA, NA))
#' lettuce.crs5 <- drm(lettuce[, c(2, 1)], fct = fct_spec)
#' summary(lettuce.crs5)
#' ED(lettuce.crs5, c(50))
#'
#' @keywords models nonlinear
#' @export
ml3b <- CRS.4b


#' @title Cedergreen-Ritz-Streibig Model with Lower Limit Fixed at 0 and Alpha = 0.25
#' (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated as of version 3.3.0. Please use [CRS.5()] instead,
#' which provides a more general and flexible interface.
#'
#' A four-parameter Cedergreen-Ritz-Streibig (CRS) hormesis model where the lower
#' asymptote (`c`) is fixed at 0 and the alpha parameter controlling the steepness
#' of the hormetic component is fixed at 0.25. The four free parameters are `b`, `d`,
#' `e`, and `f`.
#'
#' @param names A character vector of length 5 specifying the names of the model
#'   parameters in the following order:
#'   \describe{
#'     \item{`b`}{Hill slope (steepness of the dose-response curve).}
#'     \item{`c`}{Lower asymptote (fixed at 0 via the `fixed` argument).}
#'     \item{`d`}{Upper asymptote.}
#'     \item{`e`}{Effective dose producing a response midway between `c` and `d`
#'       (ED50).}
#'     \item{`f`}{Hormesis parameter controlling the magnitude of the stimulatory
#'       effect at low doses.}
#'   }
#'   Defaults to `c("b", "c", "d", "e", "f")`.
#'
#' @param fixed A numeric vector of length 5 specifying fixed (non-estimated)
#'   parameter values. Use `NA` for parameters that should be estimated freely.
#'   Defaults to `c(NA, 0, NA, NA, NA)`, which fixes the lower asymptote `c`
#'   at 0.
#'
#' @param ... Additional arguments passed to [cedergreen()].
#'
#' @return A list of class `"drcMean"` as returned by [cedergreen()], containing
#'   the model definition including the mean function, its gradient, parameter
#'   names, and fixed values. This object is intended for use as the `fct`
#'   argument in [drm()].
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [CRS.5()] — the recommended replacement for this deprecated function.
#' * [cedergreen()] — the underlying model constructor.
#' * [CRS.4a()] — the four-parameter CRS model with lower limit fixed at 0 and alpha = 1.
#' * [CRS.5c()] — the five-parameter CRS model with alpha = 0.25.
#'
#' @examples
#' # NOTE: CRS.4c() is deprecated. Use CRS.5() instead.
#' # The example below is retained for backward compatibility illustration only.
#'
#' lettuce.crsm3 <- drm( lettuce[, c(2, 1)], fct = CRS.4c() )
#' summary(lettuce.crsm3)
#' ED(lettuce.crsm3, c(50))
#'
#' # Recommended replacement:
#' fct_spec <- CRS.5(alpha_type = "c", fixed = c(NA, 0, NA, NA, NA))
#' lettuce.crs5 <- drm(lettuce[, c(2, 1)], fct = fct_spec)
#' summary(lettuce.crs5)
#' ED(lettuce.crs5, c(50))
#'
#' @keywords models nonlinear
#' @export
"CRS.4c" <- function(
    names = c("b", "c", "d", "e", "f"), 
    fixed = c(NA, 0, NA, NA, NA), # fixed c = 0
    ... ){
  
  # Deprecated warning
  lifecycle::deprecate_warn(
    when    = "3.3.0",
    what    = "CRS.4c()",
    with    = "CRS.5()"
  )
  
  return(
    cedergreen(
      fixed = fixed, names = names, alpha = 0.25, 
      fctName = as.character(match.call()[[1]]), 
      fctText = "Cedergreen-Ritz-Streibig with lower limit 0 (alpha=0.25)", 
      ...
    )
  )
}

#' @title Alias for CRS.4c (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is a deprecated alias for [CRS.4c()], itself deprecated as of
#' version 3.3.0. Please use [CRS.5()] instead, which provides a more general
#' and flexible interface.
#'
#' @inherit CRS.4c params return
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [CRS.5()] — the recommended replacement for this deprecated function.
#' * [CRS.4c()] — the function this alias points to.
#' * [cedergreen()] — the underlying model constructor.
#'
#' @examples
#' # NOTE: ml3c() is a deprecated alias for CRS.4c(). Use CRS.5() instead.
#' # The example below is retained for backward compatibility illustration only.
#'
#' lettuce.crsm3 <- drm( lettuce[, c(2, 1)], fct = ml3c() )
#' summary(lettuce.crsm3)
#' ED(lettuce.crsm3, c(50))
#'
#' # Recommended replacement:
#' fct_spec <- CRS.5(alpha_type = "c", fixed = c(NA, 0, NA, NA, NA))
#' lettuce.crs5 <- drm(lettuce[, c(2, 1)], fct = fct_spec)
#' summary(lettuce.crs5)
#' ED(lettuce.crs5, c(50))
#'
#' @keywords models nonlinear
#' @export
ml3c <- CRS.4c


## 5 Parametric --------------------------------------------------------

#' @title Cedergreen-Ritz-Streibig Five-Parameter Model with Alpha = 1
#' (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated as of version 3.3.0. Please use [CRS.5()] instead,
#' which provides a more general and flexible interface.
#'
#' A five-parameter Cedergreen-Ritz-Streibig (CRS) hormesis model where the alpha
#' parameter controlling the steepness of the hormetic component is fixed at 1.
#' All five parameters `b`, `c`, `d`, `e`, and `f` are freely estimated.
#'
#' @param names A character vector of length 5 specifying the names of the model
#'   parameters in the following order:
#'   \describe{
#'     \item{`b`}{Hill slope (steepness of the dose-response curve).}
#'     \item{`c`}{Lower asymptote (freely estimated).}
#'     \item{`d`}{Upper asymptote.}
#'     \item{`e`}{Effective dose producing a response midway between `c` and `d`
#'       (ED50).}
#'     \item{`f`}{Hormesis parameter controlling the magnitude of the stimulatory
#'       effect at low doses.}
#'   }
#'   Defaults to `c("b", "c", "d", "e", "f")`.
#'
#' @param fixed A numeric vector of length 5 specifying fixed (non-estimated)
#'   parameter values. Use `NA` for parameters that should be estimated freely.
#'   Defaults to `c(NA, NA, NA, NA, NA)`, meaning all five parameters are
#'   freely estimated.
#'
#' @param ... Additional arguments passed to [cedergreen()].
#'
#' @return A list of class `"drcMean"` as returned by [cedergreen()], containing
#'   the model definition including the mean function, its gradient, parameter
#'   names, and fixed values. This object is intended for use as the `fct`
#'   argument in [drm()].
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [CRS.5()] — the recommended replacement for this deprecated function.
#' * [cedergreen()] — the underlying model constructor.
#' * [CRS.4a()] — the four-parameter CRS model with lower limit fixed at 0 and alpha = 1.
#' * [UCRS.5a()] — the unconstrained five-parameter CRS model with alpha = 1.
#'
#' @examples
#' # NOTE: CRS.5a() is deprecated. Use CRS.5() instead.
#' # The example below is retained for backward compatibility illustration only.
#'
#' lettuce.m1 <- drm( lettuce[, c(2, 1)], fct = CRS.5a() )
#' summary(lettuce.m1)
#' ED(lettuce.m1, c(50))
#'
#' # Recommended replacement:
#' lettuce.crs5 <- drm( lettuce[, c(2, 1)], fct = CRS.5(alpha_type = "a") )
#' summary(lettuce.crs5)
#' ED(lettuce.crs5, c(50))
#'
#' @keywords models nonlinear
#' @export
"CRS.5a" <- function(
    names = c("b", "c", "d", "e", "f"), 
    fixed = c(NA, NA, NA, NA, NA),
    ... ){
  
  # Deprecated warning
  lifecycle::deprecate_warn(
    when    = "3.3.0",
    what    = "CRS.5a()",
    with    = "CRS.5()"
  )
  
  return(
    cedergreen(
      fixed = fixed, names = names, alpha = 1, 
      fctName = as.character(match.call()[[1]]), 
      fctText = "Cedergreen-Ritz-Streibig (alpha=1)", 
      ...
    )
  )
}


#' @title Alias for CRS.5a (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is a deprecated alias for [CRS.5a()], itself deprecated as of
#' version 3.3.0. Please use [CRS.5()] instead, which provides a more general
#' and flexible interface.
#'
#' @inherit CRS.5a params return
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [CRS.5()] — the recommended replacement for this deprecated function.
#' * [CRS.5a()] — the function this alias points to.
#' * [cedergreen()] — the underlying model constructor.
#'
#' @examples
#' # NOTE: ml4a() is a deprecated alias for CRS.5a(). Use CRS.5() instead.
#' # The example below is retained for backward compatibility illustration only.
#'
#' lettuce.m1 <- drm( lettuce[, c(2, 1)], fct = ml4a() )
#' summary(lettuce.m1)
#' ED(lettuce.m1, c(50))
#'
#' # Recommended replacement:
#' lettuce.crs5 <- drm( lettuce[, c(2, 1)], fct = CRS.5(alpha_type = "a") )
#' summary(lettuce.crs5)
#' ED(lettuce.crs5, c(50))
#'
#' @keywords models nonlinear
#' @export
ml4a <- CRS.5a


#' @title Cedergreen-Ritz-Streibig Five-Parameter Model with Alpha = 0.5
#' (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated as of version 3.3.0. Please use [CRS.5()] instead,
#' which provides a more general and flexible interface.
#'
#' A five-parameter Cedergreen-Ritz-Streibig (CRS) hormesis model where the alpha
#' parameter controlling the steepness of the hormetic component is fixed at 0.5.
#' All five parameters `b`, `c`, `d`, `e`, and `f` are freely estimated.
#'
#' @param names A character vector of length 5 specifying the names of the model
#'   parameters in the following order:
#'   \describe{
#'     \item{`b`}{Hill slope (steepness of the dose-response curve).}
#'     \item{`c`}{Lower asymptote (freely estimated).}
#'     \item{`d`}{Upper asymptote.}
#'     \item{`e`}{Effective dose producing a response midway between `c` and `d`
#'       (ED50).}
#'     \item{`f`}{Hormesis parameter controlling the magnitude of the stimulatory
#'       effect at low doses.}
#'   }
#'   Defaults to `c("b", "c", "d", "e", "f")`.
#'
#' @param fixed A numeric vector of length 5 specifying fixed (non-estimated)
#'   parameter values. Use `NA` for parameters that should be estimated freely.
#'   Defaults to `c(NA, NA, NA, NA, NA)`, meaning all five parameters are
#'   freely estimated.
#'
#' @param ... Additional arguments passed to [cedergreen()].
#'
#' @return A list of class `"drcMean"` as returned by [cedergreen()], containing
#'   the model definition including the mean function, its gradient, parameter
#'   names, and fixed values. This object is intended for use as the `fct`
#'   argument in [drm()].
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [CRS.5()] — the recommended replacement for this deprecated function.
#' * [cedergreen()] — the underlying model constructor.
#' * [CRS.4b()] — the four-parameter CRS model with lower limit fixed at 0 and alpha = 0.5.
#' * [CRS.5a()] — the five-parameter CRS model with alpha = 1.
#'
#' @examples
#' # NOTE: CRS.5b() is deprecated. Use CRS.5() instead.
#' # The example below is retained for backward compatibility illustration only.
#'
#' lettuce.m2 <- drm( lettuce[, c(2, 1)], fct = CRS.5b() )
#' summary(lettuce.m2)
#' ED(lettuce.m2, c(50))
#'
#' # Recommended replacement:
#' lettuce.crs5 <- drm( lettuce[, c(2, 1)], fct = CRS.5(alpha_type = "b") )
#' summary(lettuce.crs5)
#' ED(lettuce.crs5, c(50))
#'
#' @keywords models nonlinear
#' @export
"CRS.5b" <- function(
    names = c("b", "c", "d", "e", "f"), 
    fixed = c(NA, NA, NA, NA, NA),
    ... ){
  
  # Deprecated warning
  lifecycle::deprecate_warn(
    when    = "3.3.0",
    what    = "CRS.5b()",
    with    = "CRS.5()"
  )
  
  return(
    cedergreen(
      fixed = fixed, names = names, alpha = 0.5, 
      fctName = as.character(match.call()[[1]]), 
      fctText = "Cedergreen-Ritz-Streibig (alpha=0.5)", 
      ...
    )
  )
}


#' @title Alias for CRS.5b (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is a deprecated alias for [CRS.5b()], itself deprecated as of
#' version 3.3.0. Please use [CRS.5()] instead, which provides a more general
#' and flexible interface.
#'
#' @inherit CRS.5b params return
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [CRS.5()] — the recommended replacement for this deprecated function.
#' * [CRS.5b()] — the function this alias points to.
#' * [cedergreen()] — the underlying model constructor.
#'
#' @examples
#' # NOTE: ml4b() is a deprecated alias for CRS.5b(). Use CRS.5() instead.
#' # The example below is retained for backward compatibility illustration only.
#'
#' lettuce.m2 <- drm( lettuce[, c(2, 1)], fct = ml4b() )
#' summary(lettuce.m2)
#' ED(lettuce.m2, c(50))
#'
#' # Recommended replacement:
#' lettuce.crs5 <- drm( lettuce[, c(2, 1)], fct = CRS.5(alpha_type = "b") )
#' summary(lettuce.crs5)
#' ED(lettuce.crs5, c(50))
#'
#' @keywords models nonlinear
#' @export
ml4b <- CRS.5b


#' @title Cedergreen-Ritz-Streibig Five-Parameter Model with Alpha = 0.25
#' (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated as of version 3.3.0. Please use [CRS.5()] instead,
#' which provides a more general and flexible interface.
#'
#' A five-parameter Cedergreen-Ritz-Streibig (CRS) hormesis model where the alpha
#' parameter controlling the steepness of the hormetic component is fixed at 0.25.
#' All five parameters `b`, `c`, `d`, `e`, and `f` are freely estimated.
#'
#' @param names A character vector of length 5 specifying the names of the model
#'   parameters in the following order:
#'   \describe{
#'     \item{`b`}{Hill slope (steepness of the dose-response curve).}
#'     \item{`c`}{Lower asymptote (freely estimated).}
#'     \item{`d`}{Upper asymptote.}
#'     \item{`e`}{Effective dose producing a response midway between `c` and `d`
#'       (ED50).}
#'     \item{`f`}{Hormesis parameter controlling the magnitude of the stimulatory
#'       effect at low doses.}
#'   }
#'   Defaults to `c("b", "c", "d", "e", "f")`.
#'
#' @param fixed A numeric vector of length 5 specifying fixed (non-estimated)
#'   parameter values. Use `NA` for parameters that should be estimated freely.
#'   Defaults to `c(NA, NA, NA, NA, NA)`, meaning all five parameters are
#'   freely estimated.
#'
#' @param ... Additional arguments passed to [cedergreen()].
#'
#' @return A list of class `"drcMean"` as returned by [cedergreen()], containing
#'   the model definition including the mean function, its gradient, parameter
#'   names, and fixed values. This object is intended for use as the `fct`
#'   argument in [drm()].
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [CRS.5()] — the recommended replacement for this deprecated function.
#' * [cedergreen()] — the underlying model constructor.
#' * [CRS.4c()] — the four-parameter CRS model with lower limit fixed at 0 and alpha = 0.25.
#' * [CRS.5b()] — the five-parameter CRS model with alpha = 0.5.
#'
#' @examples
#' # NOTE: CRS.5c() is deprecated. Use CRS.5() instead.
#' # The example below is retained for backward compatibility illustration only.
#'
#' lettuce.m3 <- drm( lettuce[, c(2, 1)], fct = CRS.5c() )
#' summary(lettuce.m3)
#' ED(lettuce.m3, c(50))
#'
#' # Recommended replacement:
#' lettuce.crs5 <- drm( lettuce[, c(2, 1)], fct = CRS.5(alpha_type = "c") )
#' summary(lettuce.crs5)
#' ED(lettuce.crs5, c(50))
#'
#' @keywords models nonlinear
#' @export
"CRS.5c" <- function(
    names = c("b", "c", "d", "e", "f"), 
    fixed = c(NA, NA, NA, NA, NA),
    ... ){
  
  # Deprecated warning
  lifecycle::deprecate_warn(
    when    = "3.3.0",
    what    = "CRS.5c()",
    with    = "CRS.5()"
  )
  
  return(
    cedergreen(
      fixed = fixed, names = names, alpha = 0.25, 
      fctName = as.character(match.call()[[1]]), 
      fctText = "Cedergreen-Ritz-Streibig (alpha=0.25)", 
      ...
    )
  )
}

#' @title Alias for CRS.5c (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is a deprecated alias for [CRS.5c()], itself deprecated as of
#' version 3.3.0. Please use [CRS.5()] instead, which provides a more general
#' and flexible interface.
#'
#' @inherit CRS.5c params return
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @seealso
#' * [CRS.5()] — the recommended replacement for this deprecated function.
#' * [CRS.5c()] — the function this alias points to.
#' * [cedergreen()] — the underlying model constructor.
#'
#' @examples
#' # NOTE: ml4c() is a deprecated alias for CRS.5c(). Use CRS.5() instead.
#' # The example below is retained for backward compatibility illustration only.
#'
#' lettuce.m3 <- drm( lettuce[, c(2, 1)], fct = ml4c() )
#' summary(lettuce.m3)
#' ED(lettuce.m3, c(50))
#'
#' # Recommended replacement:
#' lettuce.crs5 <- drm( lettuce[, c(2, 1)], fct = CRS.5(alpha_type = "c") )
#' summary(lettuce.crs5)
#' ED(lettuce.crs5, c(50))
#'
#' @keywords models nonlinear
#' @export
ml4c <- CRS.5c