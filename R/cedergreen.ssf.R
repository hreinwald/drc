#' @title Self-starter for the Cedergreen-Ritz-Streibig Dose-Response Model
#' @description A self-starting function for the Cedergreen-Ritz-Streibig model, 
#'   used to find initial parameter estimates for non-linear regression (e.g., with `nls` or `drc`).
#'
#' @details This function is a closure that returns another function. The returned 
#'   function takes a data frame and calculates initial values for the model parameters 
#'   (b, c, d, e, f). This self-starter relies on several helper functions 
#'   (e.g., `findcd`, `findbe1`, `findbe2`, `findbe3`) which must be available in the 
#'   calling environment.
#'
#' @param method A character string specifying the method for estimating initial 
#'   'b' and 'e' parameters. Using descriptive names is preferred.
#' @param fixed A numeric vector of fixed parameter values, with `NA` for 
#'   parameters that need to be estimated. The required order is `c(b, c, d, e, f)`.
#' @param alpha A numeric value for the alpha parameter, which is treated as a known
#'   constant during the estimation of the other initial parameters.
#' @param useFixed A logical value. If `TRUE`, the function will use the non-NA 
#'   values provided in the `fixed` argument as fixed parameters and only estimate the others.
#'
#' @return A numeric vector of initial parameter estimates for the model parameters 
#'   that were not specified as `fixed`.
#' @keywords internal
#' @export
"cedergreen.ssf" <- function(method = c("loglinear", "anke", "method3", "normolle"), fixed, alpha, useFixed = FALSE) {
    # Note: The dot (.) has a special meaning in R. It is used to separate a generic function from its method in the S3 object-oriented system.
    # Therefore it is not recomended to use the "." string in function naming. The make sure that R does not falsly interpret the "." symbol I put
    # the function asignment into quotes here. Hoping this makes it more robut. 
    method <- match.arg(method)
    
    ## Helper functions for transformations and calculations
    # Transformation of the response variable
    y_transform <- function(y, c_val, d_val) {
        log((d_val - y) / (y - c_val))
    }
    # Function to calculate the 'b' parameter
    b_function <- function(x, y, c_val, d_val, e_val) {
        y_transform(y, c_val, d_val) / log(x / e_val)
    }
    # Function to calculate the 'e' parameter (ED50)
    e_function <- function(x, y, b_val, c_val, d_val) {
        x * exp(-y_transform(y, c_val, d_val) / b_val)
    }
    
    ## Assign the chosen method for finding initial 'b' and 'e' parameter values.
    ## This relies on external helper functions: findbe1, findbe2, findbe3.
    find_be_method <- switch(method,
        "loglinear" = findbe1(function(x) {
            # Safely calculate log of dose, returning NA for non-positive values
            log_dose <- rep(NA, length(x))
            log_dose[x > 0] <- log(x[x > 0])
            log_dose
        }, y_transform),
        "anke" = findbe2(b_function, e_function, "Anke"),
        "method3" = findbe3(),
        "normolle" = findbe2(b_function, e_function, "Normolle")
    )
    
    # This is the actual self-starter function returned by the closure
    function(dframe) {
        dose <- dframe[, 1]
        response <- dframe[, 2]

        # --- Parameter Initialization ---
        # The 'fixed' vector order is c(b, c, d, e, f)
        
        # Initial values for c (lower limit) and d (upper limit)
        if (useFixed && !is.na(fixed[2]) && !is.na(fixed[3])) {
            c_init <- fixed[2]
            d_init <- fixed[3]
        } else {
            # Calculate if not fixed. Relies on external findcd().
            initial_cd <- findcd(dose, response)
            c_init <- if (useFixed && !is.na(fixed[2])) fixed[2] else initial_cd[1]
            d_init <- if (useFixed && !is.na(fixed[3])) fixed[3] else initial_cd[2]
        }

        # --- Robustness Check ---
        # The y_transform requires response values to be strictly between c and d.
        if (any(response <= c_init) || any(response >= d_init)) {
            warning("Response values detected outside the initial (c, d) asymptotes. Adjusting asymptotes slightly to prevent math errors.")
            # Adjust c and d to encompass all data points, with a small buffer
            c_init <- min(c_init, min(response) - 0.01 * abs(min(response)))
            d_init <- max(d_init, max(response) + 0.01 * abs(max(response)))
        }

        # Initial values for b and e
        if (useFixed && !is.na(fixed[1]) && !is.na(fixed[4])) {
            b_init <- fixed[1]
            e_init <- fixed[4]
        } else {
            initial_be <- find_be_method(dose, response, c_init, d_init)
            b_init <- if (useFixed && !is.na(fixed[1])) fixed[1] else initial_be[1]
            e_init <- if (useFixed && !is.na(fixed[4])) fixed[4] else initial_be[2]
        }
        
        # Initial value for f
        if (useFixed && !is.na(fixed[5])) {
            f_init <- fixed[5]
        } else {
            # This calculation is based on the model's properties at a specific point.
            # It ensures the curve shape is reasonable based on the median response.
            f_init <- (2 * (median(response) - c_init) - (d_init - c_init)) * exp(1 / (e_init^alpha))
        }
        
        # Assemble a named vector of all initial parameter estimates
        initial_estimates <- c(b = b_init, c = c_init, d = d_init, e = e_init, f = f_init)
        
        # Return only the estimates for parameters that are NOT fixed (i.e., are NA in the 'fixed' vector)
        return(initial_estimates[is.na(fixed)])
    }
}

## -- OLD FUNCTION CALL BELOW -- ##
# #' @title Self-starter for Cedergreen model
# #' @keywords internal
# "cedergreen.ssf" <- function(method = c("1", "2", "3", "4"), fixed, alpha, useFixed = FALSE)
# {
#     method <- match.arg(method)
    
#     ## Defining helper functions (used below)
#     ytrans <- function(y, cVal, dVal) {log((dVal - y)/(y - cVal))}
#     bfct <- function(x, y, cVal, dVal, eVal) {ytrans(y, cVal, dVal) / log(x / eVal)}
#     efct <- function(x, y, bVal, cVal, dVal) {x * exp(-ytrans(y, cVal, dVal)/bVal)}
    
#     ## Assigning function for finding initial b and e parameter values    
#     findbe <- switch(method,
#     "1" = findbe1(function(x) {rVec <- log(x); rVec[!x>0] <- NA; rVec}, ytrans),
#     "2" = findbe2(bfct, efct, "Anke"),
#     "3" = findbe3(),
#     "4" = findbe2(bfct, efct, "Normolle"))
    
#     function(dframe)
#     {
#         x <- dframe[, 1]
#         y <- dframe[, 2]

#         ## Finding initial values for c and d parameters
#         cdVal <- findcd(x, y)
#         if (useFixed) {}  # not implemented at the moment
    
#         ## Finding initial values for b and e parameters    
#         beVal <- findbe(x, y, cdVal[1], cdVal[2])       
    
#         ## Finding initial value for f parameter
#         fVal <- (2*(median(y) - cdVal[1]) - (cdVal[2] - cdVal[1])) * exp(1/(beVal[2]^alpha))
        
#         return(c(beVal[1], cdVal, beVal[2], fVal)[is.na(fixed)])
#     }
# }
