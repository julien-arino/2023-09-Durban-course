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
    titlePanel("One sim of the Kermack-McKendrick SIR epidemic model"),
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          sliderInput("R_0",
                      "R_0:",
                      min = 0.5,
                      max = 5,
                      value = 0.5),
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
      R_0 = input$R_0
      gamma = 1/input$inv_gamma
      # Set beta so that R_0 = 1.5
      beta = R_0 * gamma / (N-I0)
      params = list(gamma = gamma, beta = beta)
      IC = c(S = S0, I = I0, R = R0)
      times = seq(0, input$t_f, 1)
      sol_KMK <- ode(IC, times, rhs_SIR_KMK, params)
      S_infty = uniroot(
        f = function(x) final_size_eq(S_inf = x,
                                      S0 = S0,
                                      I0 = I0,
                                      R_0 = R_0),
        interval = c(0.05, S0))[1]
      attack_rate = 
        (as.numeric(S0)-as.numeric(S_infty)) / as.numeric(S0) * 100
      ## Plot results
      plot(sol_KMK[, "time"], sol_KMK[, "I"], 
           type = "l", lwd = 2,
           main = TeX(sprintf("KMK SIR, $R_0=%1.2f$, $S_\\infty=%3.1f$, attack rate=%3.2f", R_0, S_infty, attack_rate)),
           xlab = "Time (days)", ylab = "Prevalence")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
