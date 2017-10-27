# load required libraries
library(shiny)
library(plotly)
library(RCurl)

# Load a list w/ function names
fun.url <- getURL("https://raw.githubusercontent.com/sylwesterf/optim/master/fun_names.csv")
fun.names <- read.csv(text = fun.url)
alg.url <- getURL("https://raw.githubusercontent.com/sylwesterf/optim/master/alg_names.csv")
alg.names <- read.csv(text = alg.url)

# Load scripts w/ definitions of algorithms and functions
script.obj <- getURL("https://raw.githubusercontent.com/sylwesterf/optim/master/obj_fun.R")
script.alg <- getURL("https://raw.githubusercontent.com/sylwesterf/optim/master/opt_alg.R")
eval(parse(text = script.obj))
eval(parse(text = script.alg))

# List of test functions
fun.list <- list(rastr,schaffer2,ackley,spheref,rosen,beale,
                 goldpr,booth,bukin6,matya,levy13,hump,easom,
                 crossit,egg,holder,mccorm,schaffer4,stybtang)

# List of algorithms
alg.list <- list(steepestDescent, simulatedAnnealing, 
                 geneticAlgorithm, particleSwarm)
# ui
ui <- fluidPage(
  
  headerPanel("Metaheuristic optimization"),
  
  sidebarPanel(
    
    # tagsform for conditional dropdown selector
    tags$form(
      
    # select an algorithm for analysis
    selectInput('input.alg', 'Choose an algorithm:', 
                as.vector(alg.names[,1])
    ),
    br(),
    
    # select a test function
    selectInput('input.fun', 'Choose a test function:', 
                as.vector(fun.names[,1])
    ),
    verbatimTextOutput('optimal'),
    br(),
    
    # set parameters
    uiOutput("ui"),
    
    # submit
    actionButton("button1", "Submit")
    )
  ),
  
  #output plot and summary data
  mainPanel(
    h4("Summary data:"),
    verbatimTextOutput('summary'),
    plotOutput('plot.c'),
    plotOutput('plot.m'),
    plotOutput('plot.h')
  )
  
)

# server 
server <- function(input, output) {
  
  # print an optimal solution for selected function
  output$optimal <- renderPrint({ fun.names[which(fun.names==input$input.fun),c("f.opt", "x.opt", "y.opt")] })
  
  # selection control
  output$ui <- renderUI({
    
    if (is.null(input$input.alg))
      return()
    
    # Depending on input$input.alg, we'll generate a different
    # UI component and send it to the client.
    switch(input$input.alg,
           "Steepest Descent" =  fluidRow(column(12, numericInput("sd.x", "Initial x:", value = 8),
                                                 numericInput("sd.y", "Initial y:", value = -9),
                                                 numericInput("sd.step", "Step size (10^):", value = -2),
                                                 numericInput("sd.crit", "Termination criterion (10^):", value = -13),
                                                 numericInput("sd.max", "Maximum number of iterations:", value = 100)
                                                 )
                                          ),
           
           "Simulated Annealing" = fluidRow(column(12, numericInput("sa.x", "Initial x:", value = -5),
                                                   numericInput("sa.y", "Initial y:", value = 7),
                                                   numericInput("sa.alpha", "Annealing schedule parameter:", value = 0.99),
                                                   numericInput("sa.temp", "Inital temperature:", value = 100),
                                                   numericInput("sa.delta", "Neighbourhood radius:", value = 0.2),
                                                   numericInput("sa.max", "Maximum number of iterations:", value = 100)
                                                   )
                                            ),
           
           "Genetic Algorithm" = fluidRow(column(12, numericInput("ga.minx", "Minimum value of x coordinate:", value = -20),
                                                 numericInput("ga.miny", "Minimum value of y coordinate:", value = -20),
                                                 numericInput("ga.maxx", "Maximum value of x coordinate:", value = 20),
                                                 numericInput("ga.maxy", "Maximum value of y coordinate:", value = 20),
                                                 numericInput("ga.col", "Coordinate encryption length:", value = 50),
                                                 numericInput("ga.pop", "Size of the population:", value = 50),
                                                 numericInput("ga.prob", "Probability of single genome mutation:", value = 0.05),
                                                 numericInput("ga.max", "Maximum number of iterations:", value = 100)
                                                 )
                                          ),
           
           "Particle Swarm Optimization"	= fluidRow(column(12, numericInput("pso.pop", "Population size:", value = 50),
                                                           numericInput("pso.lox", "Lower boundary X for initial particles:", value = -10),
                                                           numericInput("pso.loy", "Lower boundary Y for initial particles:", value = -10),
                                                           numericInput("pso.upx", "Upper boundary X for initial particles:", value = 10),
                                                           numericInput("pso.upy", "Upper boundary Y for initial particles:", value = 10),
                                                           numericInput("pso.wmin", "Intertia weight min:", value = 0.4),
                                                           numericInput("pso.wmax", "Intertia weight max:", value = 0.8),
                                                           numericInput("pso.c1", "c1 learning factor (individual experience):", value = 1.5),
                                                           numericInput("pso.c2", "c2 learning factor (social communication):", value = 0.9),
                                                           #numericInput("pso.crit", "MaxDistQuick stopping criterion threshold (10^):", value = -13),
                                                           numericInput("pso.max", "Maximum number of iterations:", value = 100)
                                                           )
                                                    )
    )
  })
  
  # calculate
  values <- reactiveValues(variable = NA)
  
  observe({
    
    if(input$button1 > 0) {
      
      values$alg.num <- which(alg.names==isolate(input$input.alg))
      values$fun.num <- which(fun.names==isolate(input$input.fun))
      
      if (values$alg.num == 1) {
        values$result <- isolate(alg.list[[values$alg.num]](fun.list[[values$fun.num]],
                                                            c(input$sd.x, input$sd.y),
                                                            10^(input$sd.step),
                                                            10^(input$sd.crit),
                                                            input$sd.max
                                                            )
                                 )

      } else if (values$alg.num == 2) {
        values$result <- isolate(alg.list[[values$alg.num]](fun.list[[values$fun.num]],
                                                            c(input$sa.x, input$sa.y),
                                                            input$sa.alpha,
                                                            input$sa.temp,
                                                            input$sa.delta,
                                                            input$sa.max
                                                            )
                                 )
        
      } else if (values$alg.num == 3) {
        values$result <- isolate(geneticAlgorithm(fun.list[[values$fun.num]],
                                                            c(input$ga.minx, input$ga.miny),
                                                            c(input$ga.maxx, input$ga.maxy),
                                                            input$ga.col,
                                                            input$ga.pop,
                                                            input$ga.prob,
                                                            input$ga.max
                                                            )
                                 )
        
      } else if (values$alg.num == 4) {
        values$result <- isolate(alg.list[[values$alg.num]](fun.list[[values$fun.num]],
                                                            input$pso.pop,
                                                            2,
                                                            c(input$pso.lox,input$pso.loy),
                                                            c(input$pso.upx,input$pso.upy), 
                                                            c(input$pso.wmin, input$pso.wmax),
                                                            input$pso.c1,
                                                            input$pso.c2,
                                                            input$pso.max
                                                            )
                                 )
        
      }
      
      # summary data
      output$summary <- renderPrint({ (values$result[c("f.opt", "x.opt", "epoch")]) })
      
      # convergence plot
      output$plot.c <- renderPlot({ 
        plot(values$result$f.hist, 
             col = 'blue', type = 'l', lwd=1, xlab = "Iteration number", ylab = "Function value", 
             main = paste("Convergence plot")
          )
        })
      
      # itertime histogram
      if (is.null(values$result$iter.time)){
        return(NULL)}
      else {
        output$plot.h <- renderPlot({
          hist(values$result$iter.time, 
              xlab = "Iteration time", main = paste("Iteration time histogram")
              )
        })
      }
      
      # plot optimization paths
      x_seq <- seq(min(values$result$x.hist)-1, max(values$result$x.hist)+1, length = 100)
      matrVal <- matrix(0, nrow = 100, ncol = 100)
      for(iRow in 1 : 100){
        for(iCol in 1 : 100){
          matrVal[iRow, iCol] <- fun.list[[values$fun.num]](c(x_seq[iRow], x_seq[iCol]))
        }
      }
      
      output$plot.m <- renderPlot({
        
        contour(x_seq, x_seq, matrVal)
        
        lines(values$result$x.hist, col = 'red', type = 'p', lwd=1)
      })
    }
  })
  
}

# run an app
shinyApp(ui, server)
