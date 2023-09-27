# Tetanus model of Cvjetanovic

library(deSolve)

# Load useful_functions.R, which contains the crop function
source("useful_functions.R")

tetanus_Cvjetanovic = function(t, y, params) {
  with(as.list(c(y, params)), {
    T = S+L_b+L+I_b+I+R+R_b
    dD = pi_IbD*gamma_b*I_b+pi_ID*gamma*I
    delta_T = dD/T
    dS_b = b*T
    dS = b*(1-lambda_b)*(T-R)+nu*R+nu_b*I+pi_IbS*gamma_b*I_b +
      pi_IS*gamma*I-(lambda+d-delta_T)*S
    dL_b = lambda_b*b*(T-R)-(epsilon_b+d-delta_T)*L_b
    dL = lambda*S-(epsilon+d-delta_T)*L
    dI_b = epsilon_b*L_b-(gamma_b+d-delta_T)*I
    dI = epsilon*L-(gamma+d-delta_T)*I
    dR = pi_IbR*gamma_b*I_b+pi_IR*gamma*I-(nu+d-delta_T)*R
    dR_b = b*R-(nu_b+d-delta_T)*R_b
    list(c(S_b+dS_b,S+dS,L_b+dL_b,L+dL,I_b+dI_b,I+dI,R+dR,R_b+dR_b,D+dD))
  })
}

params = list()
params$epsilon_b = 0.1667
params$epsilon = 0.125
params$gamma_b = 1/3
params$gamma = 0.0714
params$nu = 0.000274
params$nu_b = 0.005479
params$b = 0.00009889
params$d = 0.0000411

params$pi_IbS = 0.05
params$pi_IS = 0.3
params$pi_IbR = 0.05
params$pi_IR = 0.3
params$pi_IbD = 0.9
params$pi_ID = 0.4

params$lambda_b = 0.01
params$lambda = 0.01

IC = c(S_b = 0,
       S = 100000,
       L_b = 0,
       L = 0,
       I_b = 0,
       I = 0,
       R = 0,
       R_b = 0,
       D = 0)
tspan = 0:30 
sol <- ode(func = tetanus_Cvjetanovic, y = IC, times = tspan, 
           parms = params, method = "iteration")

T = sol[,"S"]+sol[,"L_b"]+sol[,"L"]+sol[,"I_b"]+sol[,"I"]+sol[,"R"]+sol[,"R_b"]

plot(sol[,"time"], T,
     type = "b")
