# -------------------------------------------------------#
# Large Sample Behavior in the Linear Regression Model:  #
# a Monte Carlo Illustration                             #
# Econometric Theory, LMU Munich                         #
# Mauricio Olivares                                      #
# -------------------------------------------------------#


# ---------------------- #
#  Clean the work space

rm(list = ls())
cat("\014") 
local({r <- getOption("repos"); r["CRAN"] <- "http://cran.r-project.org"; options(repos=r)}) #set repo

# ---------------------------------- #
# Set working directory and set seed
# for random generating numbers

setwd("/Users/lehmann/Google Drive/My Drive/Teaching/EcTheory-LMU/econometrics_lmu/")
set.seed(5, kind = "L'Ecuyer-CMRG")


# ---------------- #
# Load Packages

pkg <- list("Matrix", "ggplot2", "dplyr", "sandwich", "lmtest", "car", "nlme")
lapply(pkg, require, character.only = TRUE)

# ------------------------#
# Load Internal Functions
# These functions live inside the
# working directory defined earlier.

functions <- list("scripts/gen_data.R", "scripts/plots_estimates.R")
lapply(functions, source)
rm(functions)

# ----------------------------- #
# Parameters of the simulations

sample_sizes <- c(50, 100, 1000)  # Different sample sizes
num_simulations <- 5000  # Number of simulations
simulate_ols <- matrix(NA, nrow = length(sample_sizes), ncol = num_simulations)  # Matrix to store the OLS estimators

# ----------------------------- #
for (t in seq_along(sample_sizes)){
  for (s in 1:num_simulations){
    # Generate data
    data <- gen_data(n = sample_sizes[t], dgp = "model 2")
    Y <- data$Y
    X <- data$X
    # Estimate the OLS estimator 
    betaHat <- solve(t(X) %*% X) %*% t(X) %*% Y
    # Store the OLS estimator for beta_1
    simulate_ols[t, s] <- betaHat[2]
  }
}

# ----------------------------- #
#        LLN Illustration       #
# ----------------------------- #

# Histograms of the estimated slopes to showcase the LLN
my_plots(sample_sizes, simulate_ols)
plot_asympt(sample_sizes, simulate_ols)


# ----------------------------- #
#        CLT Illustration       #
# ----------------------------- #

# The idea behind the CLT is that if we
# recenter and properly scale by the square
# root of the sample size, the distribution
# of the OLS estimator should converge to 
# a normal distribution.

simulate_ols_stand <- matrix(NA, nrow = length(sample_sizes), ncol = num_simulations)  # Matrix to store the centered OLS estimators
for(t in seq_along(sample_sizes)){
  simulate_ols_stand[t, ] <-  sqrt(sample_sizes[t]) * (simulate_ols[t, ] - 2)
}

# Histograms of the estimated slopes to showcase the CLT
plot_asympt(sample_sizes, simulate_ols_stand)

# rmarkdown::render("tutorials/HAC_inference.Rmd", output_dir = "html")

