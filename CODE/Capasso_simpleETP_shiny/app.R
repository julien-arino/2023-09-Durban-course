#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(deSolve)

# Define UI
ui <- fluidPage(
    # Application title
    titlePanel("Simple ETP model of Capasso"),
    # Sidebar with slider inputs for some parameters
    sidebarLayout(
        sidebarPanel(
          sliderInput("inv_gamma_H",
                      "Average infectious period (days):",
                      min = 1,
                      max = 30,
                      value = 10),
          sliderInput("inv_d_E",
                      "Average lifetime agent (days):",
                      min = 1,
                      max = 10,
                      value = 5),
          sliderInput("c_H",
                      "Flow from humans:",
                      min = 0.01,
                      max = 2,
                      value = 0.1),
          sliderInput("beta",
                      "Fraction susceptible:",
                      min = 0,
                      max = 1,
                      value = 0.2),
          sliderInput("p",
                      "Probability of having a snack:",
                      min = 0,
                      max = 1,
                      value = 0.1),
          sliderInput("g_max",
                      "g_max:",
                      min = 0,
                      max = 100,
                      value = 10),
          sliderInput("g_half",
                      "g_half:",
                      min = 1,
                      max = 10,
                      value = 5),
          sliderInput("tf",
                      "Final time:",
                      min = 1,
                      max = 200,
                      value = 100,
                      step = 1),
        ),
        # Show a plot of the generated solution
        mainPanel(
            plotOutput("a_odePlot",width="750px")
        )
    )
)

# The functions needed for the computation (described in the simple code)
h = function(E, params) {
  OUT = params$g_max * E / (params$g_half+E)
  return(OUT)
}
g = function(E, params) {
  with(as.list(params), {
  OUT = N*beta*p*h(E,params)
  return(OUT)
  })
}
rhs_Capasso_ODE = function(t, x, params) {
  with(as.list(c(x, params)), {
    dE = c_H*H-d_E*E
    dH = g(E, params)-gamma_H*H
    tmp = g(E, params)
    list(c(dE, dH))
  })
}
R0 = function(params) {
  with(as.list(params), {
    R0 = N*beta*p*g_max*c_H / (g_half*d_E*gamma_H)
    return(R0)
  })
}

# Define server logic required to draw the result
server <- function(input, output) {
  ##
  ## Expression that generates the plot
  ##
  output$a_odePlot <- renderPlot({
    params <- list()
    params$N = 1000 # We could let this vary, we don't here..
    for (param_name in names(input)) {
      if (grepl("inv_", param_name)) {
        new_param_name = gsub("inv_", "", param_name)
        params[[new_param_name]] = 1/input[[param_name]]
      } else {
        params[[param_name]] = input[[param_name]]
      }
    }
    # Initial conditions and time span
    IC <- c(E = 10, H = 0)
    tspan <- seq(from = 0, to = params$tf, by = 0.1)
    # Compute solution
    sol_ODE = ode(y = IC,
                  func = rhs_Capasso_ODE,
                  times = tspan,
                  parms = params)
    # Make the plot
    y_max = max(max(sol_ODE[,"H"]),sol_ODE[,"E"])
    plot(sol_ODE[,"time"], sol_ODE[,"H"],
         type = "l", lwd = 2,
         xlab = "Time (days)", ylab = "Value",
         ylim = c(0, y_max),
         main = sprintf("R_0=%1.2f", round(R0(params),2)))
    lines(sol_ODE[,"time"], sol_ODE[,"E"], 
          lwd = 2, lty = 3)
    legend("topleft", legend = c("H(t)", "E(t)"),
           lwd = c(2,2), lty = c(1,3), inset = 0.01)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
