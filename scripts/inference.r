# -------------------------------------------------------#
# Asymptotic Inference in the Linear Regression Model    #
# Econometric Theory, Summer Semester 2025               #
# Mauricio Olivares, June 20th, 2025                     #
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

pkg <- list("Matrix", "ggplot2", "sandwich", "lmtest", "car")
lapply(pkg, require, character.only = TRUE)

# ------------------------#
# Load Internal Functions
# These functions live inside the
# working directory defined earlier.

functions <- list("scripts/gen_data.R", "scripts/plots_estimates.R")
lapply(functions, source)
rm(functions)

# ----------------------------------- #
#    OLS Estimation and Inference     #
# ----------------------------------- #

 # Generate data
data <- gen_data(n = 500, dgp = "model 2")
Y <- data$Y
X <- data$X

# Base R includes an intercept by default. 
# We need to remove it from our model because
# gen_data.R returns a matrix with an intercept.
fitModel <- lm(Y ~ X - 1)
summary(fitModel)

# Estimate of the covariance matrix
Vhat <- vcovHC(fitModel, type = "HC0")

# -------------------------------- #
#    Testing Linear Hypotheses     #
# -------------------------------- #

# Recall the true values of the slope parameters
# beta0 = 0.5
# beta1 = 2
# beta2 = -1

# Test the hypothesis using linearHypothesis function
# Hypothesis: H0: beta1 = 2
test_result <- linearHypothesis(fitModel, "XX1 = 2", vcov = vcovHC(fitModel, type = "HC0"), test = "Chisq")
# Display the test result
print(test_result)
# Joint hypothesis: H0: beta1 = 1 and beta2 = 0
joint_test_result <- linearHypothesis(fitModel, c("XX1 = 1", "XX2 = 0"), vcov = vcovHC(fitModel, type = "HC0"), test = "Chisq")
# Display the test result
print(joint_test_result)


# --------------------------------- #
#    Testing Non-Linear Hypotheses  #
# --------------------------------- #

# Our goal is to test the following hypotheses:
# H_0: beta0*beta1+beta2=0
#
# And then the following hypothesis:
# H_0: beta0*beta1+beta2=0 & beta1^2=3


# Case I: The nonlinear hypothesis
# H_0: beta0*beta1+beta2=0
# is a scalar, so one option is to use the delta method.
# We can achieve this using the deltaMethod() function
# from the car package to calculate the delta method

deltaMethod(fitModel, "b1*b2+b3", parameterNames = c("b1", "b2", "b3"))

# Alternatively, we could use the Wald statistic.

# Before we proceed, let's define the functions that
# calculate theta hat matrix A hat (see Assumption A5)
# for these particular hypotheses (single and multiple)

A_jacobian <- function(a, type = c("single", "multiple")) {
type <- match.arg(type)
  if (type == "single") {
    return(c(a[2], a[1], 1))
}
  if (type == "multiple") {
    A <- c(a[2], a[1], 1, 0, 2 * a[2], 0 )
    return(matrix(A, nrow = 3, ncol = 2, byrow = FALSE))
  }
  stop("Unrecognized type. Use 'single' or 'multiple'.")
}

theta <- function(a, type = c("single", "multiple")) {
  type <- match.arg(type)
  if (type == "single") {
    return(c(a[1] * a[2] + a[3]))
  }
  if (type == "multiple") {
    return(c(a[1]*a[2] + a[3], a[2]^2 - 3))
  }
  stop("Unrecognized type. Use 'single' or 'multiple'.")
}

# Calculate theta hat in the single hypothesis case
theta_hat_1 <- theta(coef(fitModel), type = "single")
# Matrix A is the Jacobian matrix of the transformation
Ahat_1 <- A_jacobian(coef(fitModel), type = "single")

# Calculate the Wald statistic
Wn_1 <- t(theta_hat_1) %*% solve(t(Ahat_1) %*% Vhat %*% Ahat_1) %*% theta_hat_1

# Rejection Rule: reject H0 if p-value<0.05
pval <- pchisq(Wn, df = 1, lower.tail = FALSE)
pval


# Let's now test the multiple null hypothesis
# H_0: beta0*beta1+beta2=0 and beta1^2=3

theta_hat_2 <- theta(coef(fitModel), type = "multiple")
Ahat_2 <- A_jacobian(coef(fitModel), type = "multiple")
Wn <- t(theta_hat_2) %*% solve(t(Ahat_2) %*% Vhat %*% Ahat_2) %*% theta_hat_2
Wn
# Rejection Rule: reject H0 if p-value<0.05
pval <- pchisq(Wn, df = 2, lower.tail = FALSE)
pval


rmarkdown::render("tutorials/Asymptotic_Inference.Rmd", output_dir = "html")
