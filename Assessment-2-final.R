############################################################################
# Bayesian Inference and Computation for Data Scientists (ZZSC5960)        #
# Assessment 2                                                             #
# Monte Carlo algorithms project: DNA methylation of cytosines             #
# Submitted by - Rushmila Islam                                            #
############################################################################

##############################
# Q 1, 2, 3
# Accept/Reject Sampling
#############################

# Define the density function
f <- function(x) {
  1/3 * dbeta(x, 1, 5) + 1/3 * dbeta(x, 3, 5) + 1/3 * dbeta(x, 10, 5)
}

# Q1: Plot the density function f(x)
par(mfrow=c(1,1))
th=seq(0,1,length=10000)
plot(th, f(th), type = "l", lwd = 2, 
     xlab = expression(theta), ylab = "Density", 
     main = "Density of mixture beta distributions",
     col='blue')

set.seed(1234)
N=10000
theta = runif(N,0,1) #generate samples from proposal distribution 
u <- runif(N,0,1) #generate u; 10,000 values from a uniform distribution
theta_max=0.0 #value of theta where density point is maximum

K <- f(theta_max) #maximum value of the density
K


#############################
# A/R Algorithm 1 
#############################
theta2=c() #theta accepted
for(i in 1:N){
  #theta_prop <- runif(1,0,1)
  acc_prob <- f(theta[i])/K 
  #u <- runif(1,0,1)
  if(u[i]<acc_prob){
    theta2 <- c(theta2, theta[i]) 
  }
}



# Plotting accept-reject sampling from algorithm 1 data against f(x) 
ymax <-2
par(mfrow=c(1,1))
plot(theta, u, xlab=expression(theta), ylab = "", 
     main = "Sample", xlim=c(0,1), ylim=c(0,ymax))
lines(th, f(th), col="blue", lwd=2)


###
par(mfrow=c(1,1))
hist(theta2, probability=T, xlab=expression(theta), ylab = "Density", 
     main = "Accept-reject sampling (Algorithm 1)", xlim=c(0,1), ylim=c(0,ymax), col='lightyellow')
lines(th, f(th), col="blue", lwd=2)

# Acceptance rate for algorithm 1
acceptance_rate_1 <- 1/K #length(theta2)/N
cat("Acceptance rate (Algorithm 1) = ", acceptance_rate_1, "\n") #0.6

# Plotting accept-reject sampling from algorithm 1
par(mfrow=c(1,1))
hist(theta2, probability=T, xlab=expression(theta), ylab = "Density", 
     main = "A/R sampling (Algorithm 1)", xlim=c(0,1), ylim=c(0,ymax), col='lightyellow')
lines(th, f(th), col="blue", lwd=2)

##############################
# Q 4, 5
# Importance Sampling
#############################

w=f(theta)
W=w/sum(w)
d=density(theta, weights=W, from=0, to=1)

# Plotting IS 
par(mfrow=c(1,1))
plot(th,f(th), type="l", lwd=2, col="blue",
     xlab=expression(theta), ylab="Density", 
     main="Importance sampling")
lines(d, col="magenta", lwd=2, lty=2)


par(mfrow=c(1,1), mar=c(4,4,2,2))
hist(theta2, probability=T, 
     xlab=expression(theta), ylab="Density", main="Accept-reject sampling vs Importance sampling",
     xlim=c(0,1),ylim=c(0,ymax), col="lightyellow") #AR
lines(d, col="magenta", lwd=2, lty=3) #IS
lines(th,f(th),col="blue",lwd=2) #Target density
legend("topright", legend = c("Target density", "Approximation IS", "Approximation A/R"), 
       col = c("blue", "magenta", "lightyellow"), lty = 1, lwd = 2)



# Mean calculation
sample_true_mean_ <- sum(theta)/length(theta)
sample_true_mean_
ar_approximation <- mean(theta2) #A/R approximation: 0.4017751
ar_approximation
is_approximation <- sum(w * theta)/sum(w) #IS approximation: 0.4048396
is_approximation

# Quantile
par(mfrow=c(1,1), mar=c(4,4,2,2))
plot(th,f(th),col="blue",type="l", lwd=2,
     xlab=expression(theta), ylab="Density",
     main="Approximation of IS vs A/R") #Target density

theta_IS <- sample(theta, size=N, prob=W, replace=T)
lines(density(theta_IS), col="magenta")
quantile(theta_IS, probs=c(0.25, 0.50, 0.75))

theta_AR <- sample(theta2, size=N, replace=T)
lines(density(theta_AR), col="black")
quantile(theta_AR, probs=c(0.25, 0.50, 0.75))


legend("topright", legend = c("Approximation IS", "Approximation A/R"), 
       col = c("magenta", "black"), lty = 1, lwd = 2)


