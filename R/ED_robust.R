#' Select Appropriate Confidence Interval Method for a drc Model
#'
#' This function determines the recommended confidence interval calculation method
#' ('type' argument in drc::ED) based on the model family of a 'drc' object.
#'
#' @param model A drc model object or a character string specifying the model name (e.g., "LL.4").
#' @param small_n A logical value. If TRUE, the t-distribution-based Fieller's method ("tfls")
#'   is used for small samples for applicable models. If FALSE, the normal-distribution-based
#'   method ("fls") is used. Defaults to TRUE.
#' @param fls_pattern A regular expression character string. This pattern is used to identify
#'   model families for which the "fls" or "tfls" method is appropriate. The default
#'   covers standard log-logistic, log-normal, Brain-Cousens, and Cedergreen-Ritz-Streibig models.
#' @param verbose A logical value. If TRUE, a message is printed when the function
#'   resorts to its default choice because the model type was not explicitly matched.
#'   Defaults to TRUE.
#'
#' @return A character string: "tfls", "fls", or "delta", representing the
#'   recommended interval type for use in `drc::ED()`.
#'
#' @author Hannes Reinwald
#'
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' # Assuming 'ryegrass_model' is a drc object created with LL.4()
#' # library(drc)
#' # ryegrass_model <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
#'
#' # 1. Pass the model object directly
#' get_ed_interval(ryegrass_model)
#' #> [1] "tfls"
#'
#' # 2. Pass the model name as a string
#' get_ed_interval("LL.4")
#' #> [1] "tfls"
#'
#' # 3. Example with a Weibull model
#' get_ed_interval("W1.4")
#' #> [1] "delta"
#'
#' # 4. Example with a large sample size assumption
#' get_ed_interval(ryegrass_model, small_n = FALSE)
#' #> [1] "fls"
#'
#' # 5. Example of a model that falls through to the default
#' #    (e.g., a hypothetical linear model 'LIN.1')
#' get_ed_interval("LIN.1")
#' #> Defaulting to 'tfls' for model type: LIN.1
#' #> [1] "tfls"
#' }
#'
get_ed_interval <- function(
    model,
    small_n = TRUE,
    fls_pattern = "^LL|^LN|^BC|^CRS",
    verbose = FALSE
) {
  # --- Input Validation and Name Extraction ---
  if (inherits(model, "drc")) {
    # If 'model' is a drc object, extract the function name
    model_name <- as.character(model$call$fct)[1]
  } else if (is.character(model) && length(model) == 1) {
    # If 'model' is a character string
    model_name <- model
  } else {
    stop("Input 'model' must be a 'drc' object or a single character string.")
  }
  
  # --- Core Logic ---
  if (grepl(fls_pattern, model_name, ignore.case = TRUE)) {
    # Log-logistic, Log-normal, and other user-defined families
    return(ifelse(small_n, "tfls", "fls"))
  } else if (grepl("^W", model_name, ignore.case = TRUE)) {
    # Weibull models
    return("delta")
  } else {
    # Default for other models (e.g., linear, quadratic)
    if (verbose) {
      message(paste("Defaulting to 'tfls' for model type:", model_name))
    }
    return("tfls")
  }
}


#' Robust Calculation of Effective Doses (ED)
#'
#' @description
#' This function serves as a robust wrapper for `drc::ED`. It calculates 
#' effective doses (EDs) for multiple specified response levels. Its primary 
#' feature is the ability to gracefully handle cases where an ED value is not 
#' mathematically estimable from the model (e.g., the requested response is 
#' outside the model's asymptotes). Instead of throwing an error, it returns a 
#' row of `NA` values for that specific response level, ensuring the overall 
#' analysis can proceed.
#'
#' @param mod An object of class 'drc', representing the fitted dose-response model.
#' @param respLev A numeric vector specifying the response levels for which to 
#'   calculate ED values (e.g., `c(10, 50)` for ED10 and ED50).
#' @param interval A character string specifying the method for calculating 
#'   confidence intervals. Defaults to the output of `get_ed_interval()`. 
#'   Common options include "delta", "tfls", or "buckland".
#' @param CI_level A numeric value between 0 and 1 indicating the confidence 
#'   level for the intervals (e.g., 0.95 for a 95% CI).
#' @param ... Additional arguments to be passed directly to `drc::ED`.
#'
#' @return 
#' A `data.table` where each row corresponds to a requested response level. 
#' The table includes the ED estimate, standard error, confidence interval 
#' (Lower, Upper), and metadata about the calculation (confidence level, method, 
#' model name, and EC level). Rows for non-estimable EDs are populated with `NA`.
#'
#' @author Hannes Reinwald
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Load necessary packages for the example
#' library(drc)
#' library(dplyr)
#' library(data.table)
#'
#' # Use a dataset with hormesis where some ED levels are not reachable
#' data(lettuce)
#'
#' # Run a hormesis model (4-parameter log-logistic with a non-zero lower asymptote)
#' m = drm(weight~conc, data = lettuce, fct = BC.4())
#'
#' # create plot header
#' p = modelFit(m)$`p value`[2] %>% round(.,4)
#' aic = AIC(m) %>% round(.,2)
#' plot_title = paste0(m$fct$name, " (LoF p = ", p, ", AIC = ", aic, ")")
#'
#' # plot data
#' plot(m, type = "all", main = plot_title)
#'
#' # Get the EC values robustly, including levels that may not be estimable
#' ED_robust(m, respLev = c(10, 50, 98), CI_level = 0.95)
#'
#' # Expected output will show calculations for ED10 and ED50, and a row of NAs for ED98.
#' }
#'
ED_robust <- function(mod, respLev = c(10, 20, 50), 
                      interval = get_ed_interval(mod$fct$name, small_n = TRUE), 
                      CI_level = 0.95, verbose = FALSE, ...) {
  # Use lapply to iterate over each response level.
  # This will return a list of data frames (one for each level).
  results_list <- lapply(respLev, function(ec) {
    #message("Calculating ED for response level: ", ec, "%")
    
    # This is the structure of a row for a failed calculation
    na_row <- data.frame(
      Estimate = NA_real_,
      stderr = NA_real_,
      Lower = NA_real_,
      Upper = NA_real_,
      confint_level = CI_level,
      confint_method = interval,
      model = drm_name(mod),
      EC = ec
    )
    
    # Use tryCatch to run ED() and capture any errors
    ed_result = tryCatch({
      
      # Attempt to calculate the ED value
      res = drc::ED(mod, respLev = ec, interval = interval, 
                    level = CI_level, display = FALSE, ... )
      
      # Additional check: Is the estimate positive or NA?
      if (is.na(res[1, "Estimate"]) || res[1, "Estimate"] <= 0){ return(NULL) }
      
      # If successful and positive, return the result
      if(verbose) message("Successfully calculated ED for response level: ", ec, "%")
      res # Return the result to be processed into a data frame
    }, error = function(e) {
      # If an error occurs (like the uniroot error), return NULL
      if(verbose) message("Error calculating ED for response level: ", ec, "% - ", e$message)
      return(NULL)
    }) # end of tryCatch
    
    
    # If ed_result is NULL (due to error or non-positive estimate), return the NA row
    if ( is.null(ed_result) ) {
      return(na_row)
    } else {
      # If successful, process the result into a clean data frame
      if(verbose) message("Appending info ...")
      as.data.frame(ed_result) %>%
        rename(stderr = "Std. Error") %>%
        mutate(
          confint_level = CI_level,
          confint_method = interval,
          model = drm_name(mod),
          EC = as.numeric(sub("^e.*[:]", "", rownames(ed_result)))
        )
    }
  })
  
  # Combine the list of single-row data frames into one final data frame
  return( data.table::rbindlist(results_list) )
}



#' Robust Calculation of Model-Averaged Effective Doses
#'
#' @description
#' This function serves as a robust wrapper for `drc::maED`. It calculates 
#' model-averaged effective doses (EDs) for specified response levels. The key 
#' feature is its resilience to errors; it iterates through each response level 
#' individually and handles failures gracefully by returning `NA` values for that 
#' level, rather than terminating the entire operation.
#'
#' @details
#' The function enhances `drc::maED` by introducing a robust calculation loop. 
#' It iterates over each element of `respLev` and calls `drc::maED` within a 
#' `tryCatch` block. This approach isolates failures, preventing an error at one 
#' response level (e.g., an EC99 that cannot be estimated) from halting the 
#' calculation of others.
#'
#' Furthermore, after a successful calculation, the function checks if the 
#' resulting 'Estimate' is positive. If the estimate is `NA`, non-positive, or 
#' if the `tryCatch` block catches an error, the function returns a structured 
#' row of `NA`s for that response level, ensuring a consistent output format.
#'
#' @param mod A model object of class 'drc', which serves as the base model for
#'   the averaging.
#' @param fct_ls A list of alternative dose-response functions (e.g., `LL.3()`, 
#'   `W1.4()`) to be used in the model averaging process. The list should be 
#'   named.
#' @param respLev A numeric vector specifying the response levels (in 
#'   percentages) for which to calculate the EDs (e.g., `c(10, 50)` for EC10 
#'   and EC50).
#' @param interval A character string specifying the type of confidence interval 
#'   to be supplied. The default is "buckland". See `drc::maED` for other options.
#' @param CI_level A numeric value between 0 and 1 specifying the confidence 
#'   level for the confidence intervals. Default is 0.95.
#' @param verbose A logical value. If `TRUE`, the function will print status 
#'   messages about the calculation progress and any errors encountered for each 
#'   response level. Default is `FALSE`.
#' @param ... Additional arguments to be passed to the underlying `drc::maED` 
#'   function.
#'
#' @return A `data.frame` with one row for each response level specified in 
#'   `respLev`. The columns are:
#'   \item{Estimate}{The estimated model-averaged effective dose.}
#'   \item{stderr}{The standard error of the estimate.}
#'   \item{Lower}{The lower bound of the confidence interval.}
#'   \item{Upper}{The upper bound of the confidence interval.}
#'   \item{confint_level}{The confidence level used for the interval.}
#'   \item{confint_method}{The method used for the confidence interval calculation.}
#'   \item{model}{A character string listing the models used for averaging.}
#'   \item{EC}{The response level (as a percentage).}
#'   If the calculation for a specific response level fails or results in a 
#'   non-positive estimate, the corresponding row will contain `NA` values for 
#'   `Estimate`, `stderr`, `Lower`, and `Upper`.
#'
#' @seealso \code{\link[drc]{maED}}
#' 
#' @author Hannes Reinwald
#'
#' @export
#' @importFrom dplyr %>% rename mutate
#' @importFrom data.table rbindlist
#'
#' @examples
#' \dontrun{
#' # Load necessary packages
#' library(drc)
#' library(dplyr)
#' library(data.table)
#' 
#' # Use a sample dataset from the drc package
#' data(spinach)
#' base_model <- drm(SLOPE ~ DOSE, data = spinach, fct = LL.4())
#' 
#' # 1. Fit a base model (e.g., four-parameter log-logistic)
#' # Use a dataset with hormesis where some ED levels might not be reachable
#' data(lettuce)
#' #' Run a hormesis model (4-parameter log-logistic with a non-zero lower asymptote)
#' base_model = drm(weight~conc, data = lettuce, fct = BC.5())
#' plot(base_model, type = "all", main = base_model$fct$name)
#' 
#' # 2. Define a named list of alternative models for averaging
#' model_list <- list(LL.4(), W1.4(), W2.4(), CRS.4c())
#' names(model_list) <- c("LL.4", "W1.4", "W2.4", "CRS.5c")
#' 
#' # 3. Inspect the model comparison
#' model_comparison = mselect(base_model, model_list, nested = TRUE)
#' head(model_comparison)
#' 
#' # 3. Calculate model-averaged ED values for multiple response levels
#' # This includes a level (EC99) that might be difficult to estimate.
#' model_names = c("W2.4", "CRS.5c")
#' ma_eds <- maED_robust(base_model, 
#'                       fct_ls = model_list[model_names], 
#'                       respLev = c(10, 50, 99),
#'                       verbose = TRUE)
#' # Print the results
#' # Note how EC99 results in NAs without stopping the calculation for EC10/EC50.
#' print(ma_eds)
#' 
#' # Plot model average EC50 value 
#' ec50 = ma_eds[ma_eds$EC == 50,"Estimate"]
#' plot(base_model, type = "all", main = base_model$fct$name)
#' abline(v = ec50, col = "darkred", lty = 2, lwd = 2)
#' 
#' # Example with a different confidence level and more models
#' # Note how the previous error disappears when using differen EC intervals.
#' ma_eds_90ci <- maED_robust(base_model, 
#'                            fct_ls = model_list, 
#'                            respLev = c(10, 20, 50),
#'                            verbose = TRUE,
#'                            CI_level = 0.90)
#' print(ma_eds_90ci)
#' }
#' 
maED_robust <- function(mod, fct_ls = NULL, respLev = c(10, 20, 50), 
                        interval = "buckland", 
                        CI_level = 0.95, verbose = FALSE, ...) {
  
  # Use lapply to iterate over each response level.
  # This will return a list of data frames (one for each level).
  results_list <- lapply(respLev, function(ec) {
    
    # Pre-calculate the model name string, as it's needed for the NA row.
    # This logic is taken directly from your original my_maED function.
    model_name <- paste0(sub("[:].*$", "", c(mod$fct$name, names(fct_ls))), collapse = "/")
    
    # This is the structure of a row for a failed calculation.
    na_row <- data.frame(
      Estimate = NA_real_,
      stderr = NA_real_,
      Lower = NA_real_,
      Upper = NA_real_,
      confint_level = CI_level,
      confint_method = interval,
      model = model_name,
      EC = ec
    )
    
    # Use tryCatch to run maED() and capture any errors.
    ma_ed_result <- tryCatch({
      
      # Attempt to calculate the model-averaged ED value for the single response level.
      res <- drc::maED(mod, fctList = fct_ls, respLev = ec, interval = interval, 
                       level = CI_level, display = FALSE, na.rm = TRUE, ...)
      
      # Additional check: Is the estimate positive and not NA?
      if (is.na(res[1, "Estimate"]) || res[1, "Estimate"] <= 0) {
        return(NULL)
      }
      
      # If successful and positive, return the result.
      if (verbose) message("Successfully calculated maED for response level: ", ec, "%")
      res # Return the result to be processed into a data frame.
      
    }, error = function(e) {
      # If an error occurs, return NULL.
      if (verbose) message("Error calculating maED for response level: ", ec, "% - ", e$message)
      return(NULL)
    }) # end of tryCatch
    
    # If ma_ed_result is NULL (due to error or non-positive estimate), return the NA row.
    if (is.null(ma_ed_result)) {
      return(na_row)
    } else {
      # If successful, process the result into a clean data frame.
      if (verbose) message("Appending info ...")
      as.data.frame(ma_ed_result) %>%
        rename(stderr = "Std. Error") %>%
        mutate(
          confint_level = CI_level,
          confint_method = interval,
          model = model_name,
          EC = as.numeric(sub("^e.*[:]", "", rownames(ma_ed_result)))
        )
    }
  })
  
  # Combine the list of single-row data frames into one final data frame.
  return(data.table::rbindlist(results_list))
}