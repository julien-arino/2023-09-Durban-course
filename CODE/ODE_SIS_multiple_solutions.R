library(deSolve)
library(latex2exp)

# Load a few useful functions
source("useful_functions.R")

# The right hand side of the ODE
rhs_SIS <- function(t, x, p) {
  with(as.list(c(x, p)), {
    N = S + I
    dS = d*(N-S) + gamma * I - d * S - beta * S * I / N
    dI = beta * S * I / N - (d + gamma) * I
    return(list(c(dS, dI)))
  })
}

# Compute the EPs
valeur_PE = function(params) {
  with(as.list(c(params)), {
    OUT = list()
    if (R_0<1) {
      OUT$S_EP = Pop
      OUT$E_EP = 0
      OUT$col = "dodgerblue4"
    } else {
      OUT$S_EP = 1/R_0*Pop
      OUT$E_EP = (1-1/R_0)*Pop
      OUT$col = "darkorange4"
    }
    return(OUT)
  })
}

# Put the parameters in a list
params = list()
params$Pop = 100
params$d = 1/(80 * 365.25)
params$gamma = 1/14
params$t_f = 150
params$I_0 = 2
# Note that we did not set R_0 or beta. This is done in a loop

# Initial conditions. It is useful to have named variables ((S=.., 
# I=..), here, so that the algorithm knows which variables it is being passed 
# (otherwise, you will need to write something like S=x[1], etc.)
IC <- c(S = (params$Pop-params$I_0), I = params$I_0)

# Times at which the solution will be returned.
tspan = seq(from = 0, to = params$t_f, by = 0.1)

# Now simulate the ODE. Loop on several values of R_0
R_0 = c(0.8, 1.5, 2.5)
# Save results in a list together with EP values
sol_ODE = list()
EP = list()
# Now loop on R_0
for (r0 in R_0) {
  # Name for list entry
  entry_name = sprintf("$R_0$=%1.1f",r0)
  # Keep the current value of R_0 to compute EPs
  params$R_0 = r0
  # R0=(beta/(d+gamma)) => beta=R0*(d+gamma)
  params$beta = r0*(params$d+params$gamma)
  # Call numerical integrator
  sol_ODE[[entry_name]] = ode(y = IC,
                              func = rhs_SIS,
                              times = tspan,
                              parms = params)
  EP[[entry_name]] = valeur_PE(params)
  EP[[entry_name]]$lty = which(r0 == R_0)
}

# Get maximum value of I across all simulations for plot. Note the use of lapply.
max_I = max(unlist(lapply(sol_ODE, function(x) max(x[,"I"]))))

# Plot
png(file = "../SLIDES/FIGS/sol_endemic_SIS_several.png",
    width = 1200, height = 800, res = 200)
y_axis = plot_hr_yaxis(sol_ODE[[1]][,"time"], sol_ODE[[1]][,"I"],
                       y_range = c(0, max_I),
                       type = "l", lwd = 5, col = EP[[1]]$col, lty = EP[[1]]$lty,
                       xlab = "Temps (jours)", ylab = "Prévalence")
points(x = params$t_f, y = EP[[1]]$E_EP*y_axis$factor, 
       col = EP[[1]]$col, pch = 19, cex = 2)
for (i in 2:length(sol_ODE)) {
  lines(sol_ODE[[i]][,"time"], sol_ODE[[i]][,"I"]*y_axis$factor,
        type = "l", lwd = 5, col = EP[[i]]$col, lty = EP[[i]]$lty)
  points(x = params$t_f, y = EP[[i]]$E_EP*y_axis$factor, 
         col = EP[[i]]$col, pch = 19, cex = 2)
}
legend("topleft", legend = TeX(names(EP)), cex = 0.8,
       col = unlist(lapply(EP, function(x) x$col)),
       lty = unlist(lapply(EP, function(x) x$lty)),
       lwd = c(3,3,3))
dev.off()
crop_figure(file = "../SLIDES/FIGS/sol_endemic_SIS_several.png")

