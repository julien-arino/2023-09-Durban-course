library(deSolve)
library(latex2exp)

# Load useful_functions.R, which contains the crop function
source("useful_functions.R")

# Functions for the incidence function
h = function(E, params) {
  # Use a Michaelis Menten type growth
  OUT = params$g_max * E / (params$g_half+E)
  return(OUT)
}
p_t = function(t, params) {
  angle = 2*pi/params$p_period
  OUT = cos(angle*t) # Make the base cos wave
  OUT = OUT/2*(params$p_max-params$p_min) # Scale
  OUT = OUT-min(OUT)+params$p_min # Shift up
  return(OUT)
} 
g = function(E, params, t) {
  OUT = params$N * params$beta * p_t(t, params) * h(E,params)
  return(OUT)
}

# Right hand side of the ODE
rhs_Capasso_ODE = function(t, x, params) {
  with(as.list(c(x, params)), {
    dE = c_H*H-d_E*E
    dH = g(E, params, t)-gamma_H*H
    list(c(dE, dH))
  })
}

R0 = function(params) {
  with(as.list(params), {
    R0 = list()
    R0$min = N*beta*p_min*g_max*c_H / (g_half*d_E*gamma_H)
    R0$max = N*beta*p_max*g_max*c_H / (g_half*d_E*gamma_H)
    return(R0)
  })
}

# Put parameters in a list
params = list()
params$N = 1000       # Total population
params$gamma_H = 1/10 # Infectious period
params$d_E = 1/5      # Lifetime agent
params$c_H = 0.1      # Flow from humans
# Human characteristics and behaviour
params$beta = 0.2 # Fraction susceptible
params$p = 0.1    # Probability of having a "snack"
# Growth function
params$g_max = 10
params$g_half = 100
# Periodic forcing parameters (on proba of having a snack)
params$p_period = 30
params$p_min = 0.001
params$p_max = 0.005
# Final time
params$t_f = 365

# Initial conditions. It is useful to have named variables here (S=.., I=..),
# so that the ODE algorithm knows what variables are passed to it (otherwise, 
# the RHS would need S=x[1], etc.)
IC <- c(E = 10, H = 0)

# Set time span. Can be either end points or specified time points. We do the latter
# (when fitting, which we are not doing here, this is useful)
tspan = seq(from = 0, to = params$t_f, by = 0.1)

sol_ODE = ode(y = IC,
              func = rhs_Capasso_ODE,
              times = tspan,
              parms = params)


# Are we plotting for a dark background
plot_blackBG = TRUE
# Plot (and save) the result
png("../FIGS/Capasso_ETP_2.png", width = 1000, height = 600)
if (plot_blackBG) {
  par(bg = 'black', fg = 'white') # set background to black, foreground white
  colour = "white"
} else {
  colour = "black"
}

y_max = max(max(sol_ODE[,"H"]), max(sol_ODE[,"E"]))
plot(sol_ODE[,"time"], sol_ODE[,"H"],
     type = "l", lwd = 2,
     xlab = "Time (days)", ylab = "Value",
     ylim = c(0, y_max),
     col.axis = colour, cex.axis = 2,
     col.lab = colour, cex.lab = 2)
lines(sol_ODE[,"time"], sol_ODE[,"E"], 
      lwd = 2, lty = 3)
legend("bottomright", legend = c("H(t)", "E(t)"),
       lwd = c(2,2), lty = c(1,3), inset = 0.01)
dev.off()
#crop_figure("../FIGS/Capasso_ETP_1.png")

tmp = R0(params)
print(tmp)