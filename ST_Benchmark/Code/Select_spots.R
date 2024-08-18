# Install packages if not already installed
if (!require("shiny")) install.packages("shiny")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")

# Load the packages
library(shiny)
library(ggplot2)
library(dplyr)

ui <- fluidPage(
  titlePanel("Tissue Spot Selector"),
  sidebarLayout(
    sidebarPanel(
      helpText("Select points on the plot to see their rownames."),
      textOutput("selectedRows")
    ),
    mainPanel(
      plotOutput("plot", click = "plot_click", brush = "plot_brush")
    )
  )
)

server <- function(input, output, session) {
  # Read your dataset
  data <- layer_manual_MOB  # Correct variable name for consistency
  
  # Render the plot
  output$plot <- renderPlot({
    ggplot(data, aes(x = x, y = y)) +
      geom_point() +
      theme_minimal()
  })
  
  # Subset data based on brush selection
  selectedData <- reactive({
    brushedPoints(data, input$plot_brush)
  })
  
  # Output selected rownames as a string with each name in quotes and separated by commas
  output$selectedRows <- renderText({
    req(selectedData())  # Ensure there is data selected
    selected_rows <- rownames(selectedData())
    if (length(selected_rows) > 0) {
      # Formatting each row name in quotes and joining them with commas
      formatted_rows <- sapply(selected_rows, function(x) paste0('"', x, '"'))
      paste(formatted_rows, collapse = ", ")
    } else {
      "No points selected."
    }
  })
}
shinyApp(ui = ui, server = server)
