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
library(latex2exp)

rhs_SIR_KMK <- function(t, x, p) {
  with(as.list(c(x, p)), {
    dS = - beta * S * I
    dI = beta * S * I - gamma * I
    dR = gamma * I
    return(list(c(dS, dI, dR)))
  })
}

final_size_eq = function(S_inf, S0 = 999, I0 = 1, R_0 = 2.5) {
  OUT = S0*(log(S0)-log(S_inf)) - (S0+I0-S_inf)*R_0
  return(OUT)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Two sims of the Kermack-McKendrick SIR epidemic model"),
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          sliderInput("R_01",
                      "R_0^1:",
                      min = 0.5,
                      max = 5,
                      value = 2.5),
          sliderInput("R_02",
                      "R_0^2:",
                      min = 0.5,
                      max = 5,
                      value = 2.5),
          sliderInput("inv_gamma",
                      "Average infectious period:",
                      min = 0.5,
                      max = 14,
                      value = 4),
          sliderInput("t_f",
                      "Final time:",
                      min = 100,
                      max = 365.25,
                      value = 200)
        ),
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$distPlot <- renderPlot({
      # Initial condition for S (to compute R_0)
      N = 1000
      I0 = 20
      S0 = N-I0
      R0 = 0
      # Get R_0 and gamma from the sliders
      R_01 = input$R_01
      R_02 = input$R_02
      gamma = 1/input$inv_gamma
      # Set beta1 so that R_01 = the wanted value
      beta1 = R_01 * gamma / (N-I0)
      params1 = list(gamma = gamma, beta = beta1)
      # Set beta2 so that R_02 = the wanted value
      beta2 = R_02 * gamma / (N-I0)
      params2 = list(gamma = gamma, beta = beta2)
      # Initial conditions
      IC = c(S = S0, I = I0, R = R0)
      times = seq(0, input$t_f, 1)
      # Results for the first simulation
      sol1_KMK <- ode(IC, times, rhs_SIR_KMK, params1)
      S_infty1 = uniroot(
        f = function(x) final_size_eq(S_inf = x,
                                      S0 = S0,
                                      I0 = I0,
                                      R_0 = R_01),
        interval = c(0.05, S0))[1]
      attack_rate1 = 
        (as.numeric(S0)-as.numeric(S_infty1)) / as.numeric(S0) * 100
      # Results for the second simulation
      sol2_KMK <- ode(IC, times, rhs_SIR_KMK, params2)
      S_infty2 = uniroot(
        f = function(x) final_size_eq(S_inf = x,
                                      S0 = S0,
                                      I0 = I0,
                                      R_0 = R_02),
        interval = c(0.05, S0))[1]
      attack_rate2 = 
        (as.numeric(S0)-as.numeric(S_infty2)) / as.numeric(S0) * 100
      # Find y limits for plots
      y_max = max(max(sol1_KMK[,"I"]), max(sol2_KMK[,"I"]))
      ## Plot results
      plot(sol1_KMK[, "time"], sol1_KMK[, "I"], 
           type = "l", lwd = 2, ylim = c(0, y_max),
           col = "red",
           main = TeX(sprintf("KMK SIR, attack rate 1=%3.2f%%, attack rate 2=%3.2f%%", attack_rate1, attack_rate2)),
           xlab = "Time (days)", ylab = "Prevalence")
      lines(sol2_KMK[, "time"], sol2_KMK[, "I"], 
            type = "l", lwd = 2, col = "blue")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
