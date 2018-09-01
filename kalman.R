# Author: Eric Kalosa-Kenyon
# Fiter a constant from observations with random noise all in 1d
#   i.e. use the simplest possible model
#
# Process model:
#   x_k = x_(k-1) + w_k         i.e. AR(1), state transition = A = 1
# Observation model:
#   z_k = x_k + v_k             i.e. H = 1

rm(list=ls())
set.seed(33)

# Set parameters
const = pi
Q = 10^(-5)     # process noise
R = 1/2         # measurement noise
x0 = rnorm(1, mean=const, sd=Q)

nk = 100 # 100 points indexed by k \in 1:100

# Allocate simulation objects
xs = matrix(NA, nk, 1)  # Process
xs[1] = x0
zs = matrix(NA, nk, 1)  # Observations
zs[1] = rnorm(1, mean=xs[1], sd=R)

# Perform simulation
for(i in 2:nk){
    xs[i] = rnorm(1, mean=xs[i-1], sd=Q)    # propagate process
    zs[i] = rnorm(1, mean=xs[i], sd=R)      # observation of propagated process
}

# Filter the process observations
xspr = matrix(NA, nk, 1)    # prior guess for observation
xspo = matrix(NA, nk, 1)    # posterior guess for observation
erpr = matrix(NA, nk, 1)    # prior error estimate
erpo = matrix(NA, nk, 1)    # posterior error estimate
ks = matrix(NA, nk, 1)      # kalman gain filters

# Source: Welch & Bishop, An Introduction to the Kalman Filter (2006)

xspr[1] = xs[1]/2   # initial guesses (in reality xs are completely unavailable)
erpr[1] = Q + R
ks[1] = erpr[1] * (erpr[1] + R)^(-1)
xspo[1] = xspr[1] + ks[1]*(zs[1] - xspr[1])
erpo[1] = (1-ks[1])*erpr[1]

for(i in 2:nk){
    xspr[i] = xspo[i-1]                     # prior guess based on proc. model
    erpr[i] = erpo[i-1] + Q                 # error prior
    ks[i] = erpr[i] * (erpr[i] + R)^(-1)    # kalman gain
    innov = zs[i] - xspr[i]                 # difference between prior and obsv
    xspo[i] = xspr[i] + ks[i]*innov         # kalman filter
    erpo[i] = (1-ks[i])*erpr[i]             # posterior estimate of obsv error
}

# TODO: plot everything
# TODO: write blog post about this
# TODO: track a stock price from this
