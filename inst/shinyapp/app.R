library(shiny)
library(R6)

source("VariantProcessor.R") 
source("RustPipeline.R")  
source("DifferentialExpression.R") 
source("DataHandling.R")

# Create an instance of Classes
rustPipeline <- RustPipeline$new()

# Initialize shared helper classes
dataHandling <- DataHandling$new()

# Initialize and pass helpers in the constructor 
differentialExpression <- DifferentialExpression$new(
  dataHandling = dataHandling
)

variantProcessor <- VariantProcessor$new(
  dataHandling = dataHandling
)

# If i want to use upload like others
#options(shiny.maxRequestSize = 40960 * 1024^2)  # ~40 GB


# Define UI ----
ui <- fluidPage(
  titlePanel("TCR-Antigen screen minimal workflow"),
  tabsetPanel(
    variantProcessor$ui(),  
    rustPipeline$ui(),
    differentialExpression$ui()
  )
)

# Define server logic ----
server <- function(input, output, session) {
  variantProcessor$server(input, output)
  rustPipeline$server(input, output)
  differentialExpression$server(input, output)
  dataHandling$server(input, output)
}

# Run the app ----
shinyApp(ui = ui, server = server)
