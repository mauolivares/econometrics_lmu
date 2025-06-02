#' @title Generate Data for Monte Carlo
#' 
#' @description This function generates data for the simulation experiments.
#' 
#' @param n   Integer. Sample size. The default value is 100.
#' @param dgp String. Different models that generate iid data. The default option corresponds with the simple linear regression model with homoskedastic errors. The other two options generate data for the linear regression model with heteroskedastic errors and the linear regression model with errors from a distribution with heavy tails.
#' 
#' @return List. A list with the dependent variable Y and the regressors X for an iid sample of size n.
#' 
#' @author Mauricio Olivares
#' @import MASS Matrix


gen_data <- function(n = 100, dgp = c("model 1", "model 2", "model 3")) {

  # Match the arguments 
  dgp <- match.arg(dgp)

  # Generate the parameter of interest, beta
  beta <- c(0.5, 2, -1)
  # Generate Regressors
  x_1 <- rep(1, n)
  x_2 <- rnorm(n)
  x_3 <- runif(n)
  x <- cbind(x_1, x_2, x_3)

  # --------------- #
  if (dgp == "model 1") {
    # Model 1 generares data for the simple linear
    # regression model with homoskedastic errors.
    sigma <- 1
    # Generate the dependent variable
    y <- x %*% beta + rnorm(n, 0, sigma)
  } else if (dgp == "model 2") {
    # Model 2 generares data for the linear regression model
    # with heteroskedastic errors
    sigma <- exp(0.2 * x_2)
    # Generate the dependent variable
    y <- x %*% beta + rnorm(n, 0, sigma)
  } else if (dgp == "model 3") {
    # Model 3 generares data for the linear regression model
    # with errors from a distribution with heavy tails
    # Generate the dependent variable
    y <- x %*% beta + rt(n, 2+10E-2)
  } else {
    stop(paste("Unrecognized option '", dgp, "'", sep = ""))
  }

  # Label Data and Export it as a list
  output <- list()
  output$Y <- y
  colnames(output$Y) <- "Y"
  output$X <- cbind(x_1, x_2, x_3)
  colnames(output$X) <- c("Intercept", "X1", "X2")
  return(output)
}