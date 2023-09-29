# Example simulation of a simple SIS model
library(GillespieSSA2)
library(adaptivetau)
library(ABM)

# Source a file with a few helpful functions for plotting (nice axes labels, crop figure)
source("useful_functions.R")

# Total population
Pop = 1000
# Initial number of infectious
I_0 = 7
# Parameters
gamma = 1/5
# R0 would be (beta/gamma)*S0, so beta=R0*gamma/S0
beta = 1.5*gamma/(Pop-I_0)
b = 0.1 
d = 0.09
pi_I = 1/(10*365.25)
pi_T = 1/(15*365.25)
delta_A = 1/180
delta_T = 1/(5*365.25)
alpha = 0.1
detect = 1/(180)
# Final time
t_f = 100



IC <- c(S = (Pop-I_0), I = I_0, A = 0, T = 0)
params <- c(gamma = gamma, beta = beta, b = b, d = d, pi_I = pi_I,
            pi_T = pi_T, delta_A = delta_A, delta_T = delta_T,
            alpha = alpha, detect = detect)
reactions <- list(
  # propensity function effects name for reaction
  reaction("b*(S+I+A+T)", c(S=+1), "birth"),
  reaction("d*S", c(S=-1), "death_S"),
  reaction("d*I", c(I=-1), "death_I"),
  reaction("d*T", c(T=-1), "death_T"),
  reaction("d*A", c(A=-1), "death_A"),
  reaction("beta*S*I", c(S=-1,I=+1), "new_infection_by_I"),
  reaction("beta*S*alpha*T", c(S=-1,I=+1), "new_infection_by_T"),
  reaction("pi_I*I", c(A=+1,I=-1), "I_to_A"),
  reaction("detect*I", c(I=-1,T=+1), "I_to_T"),
  reaction("pi_T*T", c(I=-1,A=+1), "T_to_A"),
  reaction("delta_A*A", c(A=-1), "delta_A"),
  reaction("delta_T*T", c(T=-1), "delta_T")
)

transitions = list(c(S=+1),
                   c(S=-1),
                   c(I=-1),
                   c(T=-1),
                   c(A=-1),
                   c(S=-1,I=+1),
                   c(S=-1,I=+1),
                   c(A=+1,I=-1),
                   c(I=-1,T=+1),
                   c(I=-1,A=+1),
                   c(A=-1),
                   c(T=-1)) # predator dies
# Function to calculate transition rates, given variables and parameters
lvrates <- function(x, params, t) {
  with(as.list(c(x, params)), {
    return(c(b*(S+I+A+T),
             d*S,d*I,
             d*T,
             d*A,
             beta*S*I,
             beta*S*alpha*T,
             pi_I*I,
             detect*I,
             pi_T*T,
             delta_A*A,
             delta_T*T)
    )
  })
}

set.seed(NULL)
sol = ssa.adaptivetau(IC,
                      transitions, 
                      lvrates, 
                      params, 
                      tf = 400)

plot(sol[,"time"], sol[,"I"], type ="l", col = "red")
lines(sol[,"time"], sol[,"T"], type ="l", col = "blue")


# sol <- ssa(
#     initial_state = IC,
#     reactions = reactions,
#     params = params,
#     method = ssa_exact(),
#     final_time = t_f,
#     #log_firings = TRUE
# )
# 
# # Prepare y-axis for human readable form
# y_axis = make_y_axis(c(0, max(sol$state[,"I"])))
# 
# # Are we plotting for a dark background
# plot_blackBG = TRUE
# if (plot_blackBG) {
#   colour = "white"
# } else {
#   colour = "black"
# }
# 
# png(file = "../FIGS/one_CTMC_sim.png",
#     width = 1200, height = 800, res = 200)
# if (plot_blackBG) {
#   par(bg = 'black', fg = 'white') # set background to black, foreground white
# }
# plot(sol$time, sol$state[,"I"],
#      type = "l",
#      col.axis = colour, cex.axis = 1.25,
#      col.lab = colour, cex.lab = 1.1,
#      yaxt = "n",
#      xlab = "Time (days)", ylab = "Prevalence")
# axis(2, at = y_axis$ticks, labels = y_axis$labels, 
#      las = 1,
#      col.axis = colour,
#      cex.axis = 0.75)
# dev.off()
# crop_figure(file = "../FIGS/one_CTMC_sim.png")
# 
# 
# nb_sims = 50
# sol = list()
# tictoc::tic()
# for (i in 1:nb_sims) {
#   writeLines(paste("Start simulation", i))
#   set.seed(NULL)
#   sol[[i]] <-
#     ssa(
#       initial_state = IC,
#       reactions = reactions,
#       params = params,
#       method = ssa_exact(),
#       final_time = t_f,
#     )
# }
# tictoc::toc()
# 
# # Determine maximum value of I for plot
# I_max = max(unlist(lapply(sol, function(x) max(x$state[,"I"], na.rm = TRUE))))
# # Prepare y-axis for human readable form
# y_axis = make_y_axis(c(0, I_max))
# 
# # We want to show trajectories that go to zero differently from those that go endemic,
# # so we do a bit of preprocessing, adding a colour field each solution
# for (i in 1:nb_sims) {
#   idx_last_I = dim(sol[[i]]$state)[1]
#   val_last_I = as.numeric(sol[[i]]$state[idx_last_I,"I"])
#   if (val_last_I == 0) {
#     sol[[i]]$colour = "dodgerblue4"
#     sol[[i]]$lwd = 1
#   } else {
#     sol[[i]]$colour = ifelse(plot_blackBG, "white", "black")
#     sol[[i]]$lwd = 0.5
#   }
# }
# 
# # Now do the plot
# png(file = "../FIGS/several_CTMC_sims.png",
#     width = 1200, height = 800, res = 200)
# if (plot_blackBG) {
#   par(bg = 'black', fg = 'white') # set background to black, foreground white
# }
# plot(sol[[1]]$time, sol[[1]]$state[,"I"]*y_axis$factor,
#      xlab = "Time (days)", ylab = "Prevalence",
#      type = "l",
#      xlim = c(0, t_f), ylim = c(0, I_max), 
#      col.axis = colour, cex.axis = 1.25,
#      col.lab = colour, cex.lab = 1.1,
#      yaxt = "n",
#      col = sol[[1]]$colour, lwd = sol[[1]]$lwd)
# for (i in 2:nb_sims) {
#   lines(sol[[i]]$time, sol[[i]]$state[,"I"]*y_axis$factor,
#         type = "l",
#         col = sol[[i]]$colour, lwd = sol[[i]]$lwd)
# }
# axis(2, at = y_axis$ticks, labels = y_axis$labels, 
#      las = 1,
#      col.axis = colour,
#      cex.axis = 0.75)
# dev.off()
# crop_figure("../FIGS/several_CTMC_sims.png")
# 
# 


gamma = newExpWaitingTime(0.2) # the recovery rate
beta = 0.4 # the transmission rate
sigma=0.5 # the rate leaving the latent stage and becoming infectious
n = 10000 # the population size
sim = Simulation$new(n, function(i) if (i <= 5) "I" else "S")
sim$addLogger(newCounter("S", "S"))
sim$addLogger(newCounter("I", "I"))
sim$addLogger(newCounter("R", "R"))
m = newRandomMixing()
sim$addContact(m)
sim$addTransition("E"->"I", sigma)
sim$addTransition("I"->"R", gamma)
sim$addTransition("I" + "S" -> "I" + "E" ~ m, beta) #, changed_callback=changed)
T = 15 # time when the control measure is implemented
p = 0.5 # the reduction in transmission rate due to control measure
beta.t = function(time) {
  if (time>=T) return(rexp(1, beta * p))
  t = rexp(1, beta)
  if (time + t <= T) t else T - time + rexp(1, beta * p)
}
sim$addTransition("I" + "S" -> "I" + "E" ~ m, beta.t) #, changed_callback=changed)
result = sim$run(0:200)
print(result)