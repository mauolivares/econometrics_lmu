# -------------------------------------------------------#
# Large Sample Behavior in the Linear Regression Model:  #
# a Monte Carlo Illustration                             #
# Econometric Theory, Summer Semester 2024               #
# Mauricio Olivares, June 14th, 2024                     #
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
# Alternatively
fitModel_b <- lm(data$Y ~ data$X[,2:3])
summary(fitModel_b)

# Standard errors of the OLS estimator
V0 <- summary(fitModel)$coefficients[, "Std. Error"]

# To adjust for heteroskedasticity, you can use the
# sandwich package to obtain robust standard errors.

# Heteroskedasticity-Consistent (HC) covariance matrix estimators
hc0_se <- sqrt(diag(vcovHC(fitModel, type = "HC0")))
hc1_se <- sqrt(diag(vcovHC(fitModel, type = "HC1")))
hc2_se <- sqrt(diag(vcovHC(fitModel, type = "HC2")))
hc3_se <- sqrt(diag(vcovHC(fitModel, type = "HC3")))

# Display the robust standard errors
robust_se <- data.frame(V0=V0, HC0 = hc0_se, HC1 = hc1_se, HC2 = hc2_se, HC3 = hc3_se)
print(robust_se)

# These formulas provide the estimates using the
# estimators V0, HC0, HC1, HC2, and HC3, we covered in class.
# For example, if want to derive the HC0 estimator by hand:

# Step 1. Calculate the OLS residuals
betaHat <- solve(t(X) %*% X) %*% t(X) %*% Y
# Step 2. Calculate the variance of the OLS estimator
residuals <- Y - X %*% betaHat
residuals <- as.vector(residuals)
vHat0 <- solve(t(X) %*% X) %*% t(X) %*% diag(residuals^2) %*% X %*% solve(t(X) %*% X)

# Step 3. Calculate the standard errors
robust_se$HC0hand <- sqrt(diag(vHat0))
# Display the robust standard errors
print(robust_se)

# To incorporate heteroskedasticity-consistent (HC) estimators
# when reporting standard errors in the summary() output, you can
# use the coeftest() function from the lmtest package along with
# the vcovHC() function from the sandwich package. This approach 
# will allow you to display the model coefficients with robust standard errors.

# Summary of the model with standard errors adjusted for heteroskedasticity
# HC0 type
coeftest(fitModel, vcov = vcovHC(fitModel, type = "HC0"))
# HC1 type
coeftest(fitModel, vcov = vcovHC(fitModel, type = "HC1"))
# HC2 type
coeftest(fitModel, vcov = vcovHC(fitModel, type = "HC2"))
# HC3 type
coeftest(fitModel, vcov = vcovHC(fitModel, type = "HC3"))

# Summarize the model with different HC estimators using
# the summarize_with_robust_se() function we created

summarize_with_robust_se(fitModel)


# -------------------------------- #
#    Testing Linear Hypotheses     #
# -------------------------------- #

# Test the hypothesis using linearHypothesis function
# Hypothesis: H0: beta1 = 2
test_result <- linearHypothesis(fitModel, "XX1 = 2", vcov = vcovHC(fitModel, type = "HC0"), test = "Chisq")
# Display the test result
print(test_result)
# Joint hypothesis: H0: beta1 = 1 and beta2 = 0
joint_test_result <- linearHypothesis(fitModel, c("XX1 = 1", "XX2 = 0"), vcov = vcovHC(fitModel, type = "HC0"), test = "Chisq")
# Display the test result
print(joint_test_result)
