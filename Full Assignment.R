#1a. 
#px+1 = 0.99, px+2 = 0.985, 3px+2 = 0.95, qx+4 = 0.02
#npx+nqx = 1, px+4 = 0.98
#npx = lx+n/lx
#3px+1 = lx+4/lx+1
#3px+1 = (1/px+4)*3px+2*px+1
#3px+1 = (1/0.98)*0.95*0.99
#Let a be 3px+1
a = (1/0.98)*0.95*0.99
a
#let b be 2px+1
#2px+1 = lx+3/lx+1
#2px+1 = px+2*px+1
b = 0.99*0.985
b

#1b. 
EY1 = 1
EY2 = 2
EY3 = 3
VarY1 = 1
VarY2 = 2
VarY3 = 3
CovY1Y2 = 0
CovY1Y3 = 0
CovY2Y3 = 1
#W = Y3/(Y1+Y2)
#Given Y1 has 0 covariance with Y3 and Y1 has 0 covariance with Y3, it is not necessarily implied 
# as independence, but we can apporximate it as independence. 
#Independece, implies E(Y1+Y2) = E(Y1)+E(Y2)
EW = EY3/(EY1+EY2)
EW
#Similarly, VarW = E(W^2) - E(W)^2
#VarW =  (Var[Y3] + E[Y3]^2)/((Var[Y1] + E[Y1]^2) + (Var[Y2] + E[Y2]^2) + 2Cov(Y1,Y2)) - 1
VarW = 1 - 5/9
VarW

#2a. 
# Set seed and simulate data
set.seed(1)
xdata <- rnorm(1000, mean = 1, sd = 1)

# Function to calculate the first order partial derivatives of the 
#summed log-likelihoods of a Normal distribution
first_derivatives <- function(mu, sigma, xdata) {
  n <- length(xdata)
  dmu <- sum(xdata - mu) / sigma^2
  dsigma <- sum((xdata - mu)^2 / sigma^3) - n / sigma
  return(c(dmu, dsigma))
}

# Function to calculate the second order partial derivatives of the 
#summed log-likelihoods of a Normal distribution
second_derivatives <- function(mu, sigma, xdata) {
  n <- length(xdata)
  d2mu_dmu2 <- -n / sigma^2
  d2mu_dsigma2 <- sum(-1 / sigma^2 * (xdata - mu))
  d2sigma_dmu2 <- d2mu_dsigma2
  d2sigma_dsigma2 <- sum(3 / sigma^4 * (xdata - mu)^2) - n / sigma^2
  return(matrix(c(d2mu_dmu2, d2mu_dsigma2, d2sigma_dmu2, d2sigma_dsigma2), 
                nrow = 2, ncol = 2))
}

# Calculate the first and second order partial derivatives at population true values
mu_true <- 1
sigma_true <- 1
first_derivatives(mu_true, sigma_true, xdata)
second_derivatives(mu_true, sigma_true, xdata)

#2b. 
#Install packages and set data.
install.packages("maxLik", lib = '/private/var/folders/jc/kbh7y6jj2_q0dbkbq75wnld40000gn/T/RtmpUuq79r/downloaded_packages')
set.seed(1)
xdata <- rnorm(1000, mean = 1, sd = 1)
#First, create the log-likelihood function of both the mu and sigma functions
logL <- function(params, data) {
  mu <- params[1]
  sigma <- params[2]
  -sum(dnorm(data, mean = mu, sd = sigma, log = TRUE))}
#Fit data using maxLik function and use estimates
fit <- (maxLik(logL, start = c(mu = 0, sigma = 1), data = xdata))
mu_hat <- fit$estimate[1]
sigma_hat <- fit$estimate[2]
#Provide the confidence intervals
ci_mu <- confint(fit, "mu", level = 0.95)
confint(mu_hat, level = 0.95)
ci_sigma <- confint(fit, "sigma", level = 0.95)
print("mu_hat =", mu_hat, "\n")
cat("95% CI for mu: [", ci_mu[1], ", ", ci_mu[2], "]\n")
cat("sigma_hat =", sigma_hat, "\n")
cat("95% CI for sigma: [", ci_sigma[1], ", ", ci_sigma[2], "]")

3a.
#Download necessary packages for the data
library("KMsurv")
library("survival")
data(burn)
summary(burn)
str(burn)

#Time = t1, status = z11, gender = z2, race = z3, percentage burn = z4
fit <- coxph(Surv(T1, Z11) ~ Z2 + Z3 + Z4, data = burn)
summary(fit)
#The R^2 value of 0.245 indicates that 20% of status can be explained by time. 

#3b. 
# Fit the Cox PH model
fitcoxph <- coxph(Surv(T1, Z11) ~ Z2 + Z3 + Z4, data = burn)

# Check the goodness of fit of the model
summary(fitcoxph)
plot(fit, conf.int = TRUE, xlab = "Time (days)", ylab = "Probability of survival")

# Create a new dataset with the covariate values of interest
subdata <- subset(burn, Z3 = "1" & Z2 == "0" & tbsa > 50)
fit_sub <- survfit(Surv(T1, Z11) ~ 1, data = subdata)

# Find the estimated S(30)
S30 <- summary(fit_sub, T1 = '30')$surv[1]
S30

#3c. 
summary(fit)
#Based on the output we can see that the coefficient for gender was -0.49, 
#for race -0.03, and 0.04 for total surface area burned.
#This is suggestive of the fact that females have a higher hazard to males, there is
#no difference in race and that the total surface area increases as a percentage. 

#4a. 
library("KMsurv")
library("survival")
data(bmt)

# Subset data for ALL disease group
all_data <- subset(bmt, disease ="ALL")

# Perform survival analysis using KME
kme <- survfit(Surv(t1, group) ~ 1, data = all_data)

# Estimated S(500)
est_S_500 <- summary(kme, times = 500)$surv[1]

est_S_500

#4b. 
# Subset data for NAE disease group
nae_data <- subset(bmt, disease ="NAE")

# Perform survival analysis using KME
kme2 <- survfit(Surv(t1, group) ~ 1, data = nae_data)

# Estimated S(500)
est_S_500_2 <- summary(kme2, times = 500)$surv[1]

est_S_500_2

#4c. 
# Subset data for AML low and high risk groups
aml_low_data <- subset(bmt, disease ="AML" & tx=="Low")
aml_high_data <- subset(bmt, disease ="AML" & tx=="High")

# Perform survival analysis using KME for AML low risk group
kme_aml_low <- survfit(Surv(t1, group) ~ 1, data = aml_low_data)
kme_aml_low_ci <- survfit(Surv(t1, group) ~ 1, data = aml_low_data, conf.int = 0.95)
kme_aml_low_ci

# Perform survival analysis using KME for AML high risk group
kme_aml_high <- survfit(Surv(t1, group) ~ 1, data = aml_high_data)
kme_aml_high_ci <- survfit(Surv(t1, group) ~ 1, data = aml_high_data, conf.int = 0.95)
kme_aml_high_ci

# Plot estimated survival probabilities for both groups
plot(kme_aml_low, main = "Kaplan-Meier Survival Plot for AML Low and High Risk Groups", 
     xlab = "Time (days)", ylab = "Survival Probability", col = "black")
lines(kme_aml_low_ci, col = "red", lty = 6)
lines(kme_aml_high, col = "blue")
lines(kme_aml_high_ci, col = "blue", lty = 2)
legend("topright", legend = c("AML Low Risk", "AML High Risk"), 
       col = c("red", "blue"), lty = 1)
#The low CI and the high CI are the same, so the lines overlap.