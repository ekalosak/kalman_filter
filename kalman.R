# Author: Eric Kalosa-Kenyon
# Fiter a constant from observations with random noise all in 1d
#   i.e. use the simplest possible model
#
# Process model:
#   x_k = x_(k-1) + w_k         i.e. AR(1), state transition = A = 1
# Observation model:
#   z_k = x_k + v_k             i.e. H = 1

## @knitr setup
set.seed(33)

library(ggplot2)
library(latex2exp)
library(mvtnorm)
# library(gganimate)

## @knitr param_simple

# Set parameters
const = pi
Q = 10^(-5)     # process noise
R = 1/2         # measurement noise
x0 = rnorm(1, mean=const, sd=Q)

nk = 100 # 100 points indexed by k \in 1:100

## @knitr simulate_simple

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

## @knitr plot_simple
pltdf = data.frame(
                   time=1:nk,
                   truth=xs,
                   observations=zs,
                   prior=xspr,
                   posterior=xspo
                   )

sz = 1.5
plt_simple = ggplot(pltdf, aes(x=time)) +
    geom_line(aes(y=truth, color="blue"), size=sz) + # process
    geom_point(aes(y=observations, color="red"), alpha=0.6) + # observations
    geom_line(aes(y=posterior, color="green"), alpha=0.7, size=sz) + # filtered
    labs(
         title=paste("Kalman filter of a noisy constant,", const),
         x="Time", y="Value", color=""
         ) +
    scale_color_manual(
                       labels=c("Truth", "Posterior", "Observation"),
                       values=c("steelblue", "maroon", "orange")
                       )

plt_simple

## @knitr plot_error_profile_1d
pltdf$err = with(pltdf, sqrt((truth-posterior)^2))

plt_simple_error = ggplot(pltdf, aes(x=time, y=err)) +
    geom_col(fill="goldenrod3") +
    labs(title=paste("Error of posterior"), x="Time", y="RSE")
plt_simple_error

## @knitr simulate_movement

# Parameterize
dimx = 6 # 2 pos, 2 vel, 2 accel
dimz = 2 # observe only 2 position components
x0 = c(0,0, 0.3,2, 0,-1/5)
stopifnot(length(x0) == dimx)
proc_noise_pos = 0.15
proc_noise_vel = 0.1
obsv_noise_pos = 0.1

# Create noise and process objects
proc_sig = matrix(0, dimx, dimx)
proc_sig[1,1] = proc_sig[2,2] = proc_noise_pos
proc_sig[3,3] = proc_sig[4,4] = proc_noise_vel

obsv_sig = diag(obsv_noise_pos, 2)

A = matrix(0, dimx, dimx) # A is dynamics matrix
A[1,1] = A[1,3] = 1 # x_(t+1) = x_(t) + vel(x)_(t)
A[2,2] = A[2,4] = 1 # same for y
A[3,3] = A[3,5] = 1 # velocity_next = velocity_current + acel_current
A[4,4] = A[4,6] = 1 # velocity_next = velocity_current + acel_current
A[5,5] = 1 # acceleration is constant
A[6,6] = 1

B = matrix(0, dimz, dimx) # B is observation matrix
B[1,1] = B[2,2] = 1

# Allocate and initialize simulation objects
xs = matrix(NA, nk, dimx)  # Process
xs[1,] = x0
zs = matrix(NA, nk, dimz)  # Observations
z0 = B %*% x0 + t(rmvnorm(1, mean=c(0,0), sigma=obsv_sig))
zs[1,] = z0

# Perform simulation
nstopped = nk
for(i in 2:nk){
    proc_noise = t(rmvnorm(1, mean=rep(0, dimx), sigma=proc_sig))
    xnext = A %*% xs[i-1,] + proc_noise

    if(xnext[2] < 0){ # if projectile hit ground
        nk = i
        break
    }
    xs[i,] = xnext

    obsv_noise = t(rmvnorm(1, mean=rep(0, dimz), sigma=obsv_sig))
    znext = B %*% xnext + obsv_noise
    zs[i,] = znext
}

xs = na.omit(xs)
zs = na.omit(zs)

# TODO: plot the animated projectile observations

# Filter the process observations
xspr = matrix(NA, nk, dimz)    # prior guess for observation
xspo = matrix(NA, nk, dimz)    # posterior guess for observation
erpr = matrix(NA, nk, dimz)    # prior error estimate
erpo = matrix(NA, nk, dimz)    # posterior error estimate
ks = matrix(NA, nk, dimz)      # kalman gain filters

# Source: Welch & Bishop, An Introduction to the Kalman Filter (2006)

# TODO:
xspr[1,] = c(0,0)
# erpr[1] = Q + R
# ks[1] = erpr[1] * (erpr[1] + R)^(-1)
# xspo[1] = xspr[1] + ks[1]*(zs[1] - xspr[1])
# erpo[1] = (1-ks[1])*erpr[1]

# for(i in 2:nk){
#     xspr[i] = xspo[i-1]                     # prior guess based on proc. model
#     erpr[i] = erpo[i-1] + Q                 # error prior
#     ks[i] = erpr[i] * (erpr[i] + R)^(-1)    # kalman gain
#     innov = zs[i] - xspr[i]                 # difference between prior and obsv
#     xspo[i] = xspr[i] + ks[i]*innov         # kalman filter
#     erpo[i] = (1-ks[i])*erpr[i]             # posterior estimate of obsv error
# }
#
#
# ## @knitr the_rest
#
# # TODO: grid with different observation noises and initial guess errors,
# #       test convergence rates and plot histograms for credible bands on those
# #       rates.
#
# # TODO: write blog post about this
# # TODO: track a stock price from this
