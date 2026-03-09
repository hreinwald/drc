# Some test code for the cedergreen2 function

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
#' @keywords internal
#'
cedergreen_maxfct <- function(all_params, alpha, lower = 1e-6, upper = 1000)
{
  # Define the dose-response model using named parameters for clarity
  response_model <- function(dose, p, alpha) {
    p$c + (p$d - p$c + p$f * exp(-1 / (dose^alpha))) / (1 + exp(p$b * (log(dose) - log(p$e))))
  }

  # Use optimize() to directly find the dose that maximizes the response.
  # It is more robust than finding the root of the derivative.
  # We search for the maximum by telling optimize to maximize=TRUE.
  opt_result <- tryCatch({
    optimize(
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
#' @export
#' @examples
#' \donttest{
#' # This example requires the 'drc' package to be installed.
#' if (requireNamespace("drc", quietly = TRUE)) {
#' 
#'   library(drc)
#' 
#'   # Sample data exhibiting a hormetic effect (initial stimulation)
#'   dose <- c(0, 0.1, 0.5, 1, 5, 10, 20)
#'   response <- c(100, 102, 95, 80, 40, 25, 20)
#'   my_data <- data.frame(dose = dose, response = response)
#' 
#'   # 1. Basic model fit, specifying the mandatory 'alpha' parameter
#'   # The self-starter will find initial values for all 5 parameters.
#'   model_fit <- drm(response ~ dose, data = my_data,
#'                    fct = cedergreen(alpha = 0.5))
#'   
#'   summary(model_fit)
#'   plot(model_fit, main = "Cedergreen Model Fit")
#' 
#'   # 2. Model fit with a fixed lower limit (c = 20)
#'   # The self-starter will now only find initial values for b, d, e, and f.
#'   model_fit_fixed <- drm(response ~ dose, data = my_data,
#'                          fct = cedergreen(alpha = 0.5, fixed = c(NA, 20, NA, NA, NA)))
#'         
#'   summary(model_fit_fixed)
#' }
#'}
#' 
"cedergreen2" <- function(
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
    

    # Add edfct and maxfct functions here ... 
    # >>> REVIEW LATER <<<
    # My goal would be to define cedergreen_edfct and cedergreen_maxfct as internal helper functions, 
    # and then assign them to edfct and maxfct respectively. This keeps the main function clean and focused on its primary role.
    # Defined the helper functions above already but they need to be integrated here properly. For now, I'll just assign them directly.

    ## Defining the ED function <<< replace by cedergreen_edfct >>>
    # edfct <- function(parm, respl, reference, type, lower = 1e-4, upper = 10000, ...)
    # {
    # #        if (is.null(upper)) {upper <- 1000}
    # #        if (missing(upper2)) {upper2 <- 1000}
    # interval <- c(lower, upper) 
    # parmVec[notFixed] <- parm
    # p <- EDhelper(parmVec, respl, reference, type, TRUE)  # FALSE)  Changed 2010-06-02 after e-mail from Claire       
    # tempVal <- (100-p)/100
    
    # helpFct <- function(dose) {parmVec[2]+(parmVec[3]-parmVec[2]+parmVec[5]*exp(-1/(dose^alpha)))/(1+exp(parmVec[1]*(log(dose)-log(parmVec[4]))))}
    # #        doseVec <- exp(seq(-upper2, upper2, length=1000))
    # doseVec <- exp(seq(log(interval[1]), log(interval[2]), length=1000))
    # maxAt <- doseVec[which.max(helpFct(doseVec))]
    # #        print(maxAt)
    # #        print(upper)
    
    # eqn <- function(dose) {tempVal*(1+exp(parmVec[1]*(log(dose)-log(parmVec[4]))))-(1+parmVec[5]*exp(-1/(dose^alpha))/(parmVec[3]-parmVec[2]))}
    # EDp <- uniroot(eqn, lower=maxAt, upper=upper)$root
    
    # EDdose <- EDp
    # tempVal1 <- exp(parmVec[1]*(log(EDdose)-log(parmVec[4])))
    # tempVal2 <- parmVec[3]-parmVec[2]
    # derParm <- c(tempVal*tempVal1*(log(EDdose)-log(parmVec[4])), -parmVec[5]*exp(-1/(EDdose^alpha))/((tempVal2)^2),
    #              parmVec[5]*exp(-1/(EDdose^alpha))/((tempVal2)^2), -tempVal*tempVal1*parmVec[1]/parmVec[4],
    #              -exp(-1/(EDdose^alpha))/tempVal2)
    # derDose <- tempVal*tempVal1*parmVec[1]/EDdose-parmVec[5]/tempVal2*exp(-1/(EDdose^alpha))/(EDdose^(1+alpha))*alpha 
    
    # EDder <- derParm/derDose
    
    # return(list(EDp, EDder[notFixed]))
    # }

    ## Finding the maximal hormesis <<< replace by cedergreen_maxfct >>>
    # maxfct <- function(parm, lower = 1e-3, upper = 1000)
    # {
    # #        if (is.null(upper)) {upper <- 1000}
    # #        if (is.null(interval)) {interval <- c(1e-3, 1000)}            
    # #        alpha <- 0.5
    # parmVec[notFixed] <- parm
    
    # optfct <- function(t)
    # {
    #   expTerm1 <- parmVec[5]*exp(-1/(t^alpha))
    #   expTerm2 <- exp(parmVec[1]*(log(t)-log(parmVec[4])))
      
    #   return(expTerm1*alpha/(t^(alpha+1))*(1+expTerm2)-(parmVec[3]-parmVec[2]+expTerm1)*expTerm2*parmVec[1]/t)
    # }
    
    # ED1 <- edfct(parm, 1, lower, upper)[[1]]
    
    # doseVec <- exp(seq(log(1e-6), log(ED1), length = 100))
    # #        print((doseVec[optfct(doseVec)>0])[1])
    
    # maxDose <- uniroot(optfct, c((doseVec[optfct(doseVec)>0])[1], ED1))$root
    # return(c(maxDose, fct(maxDose, matrix(parm, 1, length(names)))))
    # }
    
    # Return results
    returnList <- list(
        fct    = fct, 
        ssfct  = ssfct, 
        names  = names, 
        deriv1 = deriv1, # Note: deriv1 is incomplete without xlogx
        deriv2 = NULL,
        edfct  = cedergreen_edfct,  # Assigning the robust ED function; replaced edfct with cedergreen_edfct
        maxfct = cedergreen_maxfct, # Assigning the robust maxfct function; replaced maxfct with cedergreen_maxfct
        name   = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
        text   = ifelse(missing(fctText), "Cedergreen-Ritz-Streibig", fctText),     
        noParm = sum(is.na(fixed))
    )
    class(returnList) <- "mllogistic"
    invisible(returnList)
}