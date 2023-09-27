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
g = function(E, params) {
  OUT = params$N * params$beta * params$p * h(E,params)
  return(OUT)
}

# Right hand side of the ODE
rhs_Capasso_ODE = function(t, x, params) {
  with(as.list(c(x, params)), {
    dE = c_H*H-d_E*E
    dH = g(E, params)-gamma_H*H
    list(c(dE, dH))
  })
}

R0 = function(params) {
  with(as.list(params), {
    R0 = N*beta*p*g_max*c_H / (g_half*d_E*gamma_H)
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
params$p = 0.1    # Probability of having "snack"
# Growth function
params$g_max = 10
params$g_half = 100
# If we fix the R0, what do we fix it to?
params$R_0 = 1.5
# Final time
params$t_f = 150
# R0 is (beta/gamma)*S0, so beta=R0*gamma/S0
# params$beta = params$R_0*params$gamma/(params$Pop-params$I_0)

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
png("../FIGS/Capasso_ETP_1.png", width = 1000, height = 600)
if (plot_blackBG) {
  par(bg = 'black', fg = 'white') # set background to black, foreground white
  colour = "white"
} else {
  colour = "black"
}

plot(sol_ODE[,"time"], sol_ODE[,"H"],
     type = "l", lwd = 2,
     xlab = "Time (days)", ylab = "Value",
     col.axis = colour, cex.axis = 2,
     col.lab = colour, cex.lab = 2)
lines(sol_ODE[,"time"], sol_ODE[,"E"], 
      lwd = 2, lty = 3)
legend("bottomright", legend = c("H(t)", "E(t)"),
       lwd = c(2,2), lty = c(1,3), inset = 0.01)
dev.off()
#crop_figure("../FIGS/Capasso_ETP_1.png")

R0(params)