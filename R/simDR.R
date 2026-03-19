#' Simulating ED values under various scenarios
#'
#' Simulating ED values for a given model and given dose values.
#'
#' The arguments \code{mpar} and \code{sigma} are typically obtained from a
#' previous model fit. Only dose-response models assuming normally distributed
#' errors can be used.
#'
#' @param mpar numeric vector of model parameters.
#' @param sigma numeric specifying the residual standard deviation.
#' @param fct list supplying the chosen dose-response mean function (e.g., \code{LL.4()}).
#' @param noSim numeric giving the number of simulations. Defaults to \code{1000}.
#' @param conc numeric vector of concentration/dose values. Must contain at least 5 values.
#' @param edVec numeric vector of ED levels to estimate in each simulation. Defaults to
#'   \code{c(10, 50)}.
#' @param seedVal numeric giving the seed used to initialise the random number generator.
#'   Defaults to \code{20070723}.
#'
#' @return Invisibly returns a list with one element:
#'   \describe{
#'     \item{\code{se}}{A 3D array of dimensions
#'       \code{(length(conc) - 4) x 6 x length(edVec)} containing empirical
#'       standard deviations of the estimated ED values. Rows correspond to the
#'       number of concentration levels used (starting from 5). Columns correspond
#'       to the number of replicates per concentration (1 to 6). The third dimension
#'       corresponds to each ED level in \code{edVec}.}
#'   }
#'   The array values are also printed to the console during execution.
#'
#' @author Christian Ritz, Hannes Reinwald
#'
#' @examples
#' ryegrass.m1 <- drm(ryegrass, fct = LL.4())
#'
#' simDR(
#'   mpar     = coef(ryegrass.m1),
#'   sigma    = sqrt(summary(ryegrass.m1)$resVar),
#'   fct      = LL.4(),
#'   noSim    = 2,
#'   conc     = c(1.88, 3.75, 7.50, 0.94, 15, 0.47, 30, 0.23, 60),
#'   seedVal  = 20070723
#' )
#'
#' @keywords models nonlinear
simDR <- function(mpar, sigma, fct, noSim = 1000, conc, edVec = c(10, 50), seedVal = 20070723)
{
  set.seed(seedVal)
  
  ## Calculating the true ED values
  n_ed    <- length(edVec)
  ed_true <- rep(0, n_ed)
  for (i in 1:n_ed)
  {
    ed_true[i] <- fct$edfct(mpar, edVec[i], type = "relative")[[1]]
  }
  
  ## Run simulations
  ed_array  <- array(NA, c(length(conc) - 4, 6, n_ed))
  sim_biases <- matrix(NA, noSim, n_ed)
  
  for (i in 5:length(conc))
  {
    conc_sorted <- sort(conc[1:i])
    for (j in 1:6)
    {
      conc_rep <- rep(conc_sorted, rep(j, i))
      sim_data <- rdrm(noSim, fct, mpar, conc_rep, ypar = sigma)
      
      for (k in 1:noSim)
      {
        fit <- try(drm(sim_data$y[k, ] ~ sim_data$x[k, ], fct = fct), silent = TRUE)
        if (!inherits(fit, "try-error"))
        {
          ed_estimates     <- ED(fit, edVec, display = FALSE)
          sim_biases[k, ] <- ed_estimates[, 1] - ed_true
        }
      }
      ed_array[i - 4, j, ] <- apply(sim_biases, 2, sd, na.rm = TRUE)
    }
  }
  
  ## Display results
  cat("Concentrations used:", conc, "\n\n")
  for (i in 1:n_ed)
  {
    result_matrix             <- ed_array[, , i]
    colnames(result_matrix)   <- 1:6
    rownames(result_matrix)   <- 5:length(conc)
    
    cat("ED value considered:", edVec[i], "\n")
    cat("Conc. no.\\Replicates:", "\n")
    print(result_matrix)
    cat("\n\n")
  }
  
  invisible(list(se = ed_array))
}