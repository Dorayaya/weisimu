

library(shinyWidgets)
library(weisimu)

params = c("n","shape","p")

ui <- fluidPage(
  fluidRow(
    column(4,
      numericRangeInput("n", "sample size", value = c(20, 100), min = 21),
      numericRangeInput("shape", "shape param", value = c(0.5,0.5), min = 0.1),
      numericRangeInput("p", "trimming prop.", value = c(0.05, 0.05), min = 0.01, max = 0.5),
      numericInput("S", "num. simulations", value = 100, min = 1),
      selectInput("xaxis", "Which param on x-axis?", params)
    ),
    column(4,
      plotOutput("density"))
  ),
  actionButton("simulate", "Simulate!"),
  plotOutput("plot")
)

server <- function(input, output) {
  reactplot = eventReactive(input$simulate, {
    if (input$xaxis=="n") simtrim_by(n=seq(min(input$n),max(input$n),by=1),
                                     shape=min(input$shape),scale=1,p=min(input$p),S=input$S,
                                     main=paste0("Simulation: MSE for trimmed vs. untrimmed mean (shape=",min(input$shape),"; p=",min(input$p),")"))
    if (input$xaxis=="shape") simtrim_by(n=min(input$n),
                                         shape=seq(min(input$shape),max(input$shape),length.out=100),
                                         scale=1,p=min(input$p),S=input$S,
                                         main=paste0("Simulation: MSE for trimmed vs. untrimmed mean (n=",min(input$n),"; p=",min(input$p),")"))
    if (input$xaxis=="p") simtrim_by(n=min(input$n),
                                     shape=min(input$shape),scale=1,
                                     p=seq(min(input$p),max(input$p),length.out=100),S=input$S,
                                     main=paste0("Simulation: MSE for trimmed vs. untrimmed mean (shape=",min(input$shape),"; n=",min(input$n),")"))    
  }
  )
  
  output$plot <- renderPlot({
    reactplot()
  }, height=450, width=700)
  
  output$density = renderPlot({
    curve(dweibull(x, shape=min(input$shape)),from=0,to=4,col="Blue",ylab="density",lwd=2,
          main="pdf of weibull with given shape parameters")
    curve(dweibull(x, shape=max(input$shape)),from=0,to=4,col="Red",ylab="density",lwd=2,add=T)
    legend("topright",legend=paste0("shape=",c(min(input$shape),max(input$shape))),
           col=c("Blue","Red"),lty=c(1,1),lwd=c(2,2))
  }, height=350, width=450)
}

shinyApp(ui, server)