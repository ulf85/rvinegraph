
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(networkD3)


shinyUI(fluidPage(

  # Application title
  titlePanel("R-vine Tree Visualisation"),

  # Sidebar 
  sidebarLayout(
    sidebarPanel(
      HTML("Upload an .RData file with your RVM-object.<br>
           The App will plot the corresponding R-vine tree(s).<br>
           In the text field below you can choose which tree should be plotted.<br><br>"),
      # Upload data:
      fileInput("file", "Upload Rdata-file:", accept = ".RData"),
      textInput("RVMname", "Name of the RVM-object (A list of loaded objects is printed to the right.)"),
      numericInput("tree", "Which tree should be plotted?", min=1, value = 1),
      downloadButton('downloadPlot', 'Download Plot'),
      br(),
      HTML("TODO: change font size, color, etc.<br>"),
      HTML("TODO: export figure")
    ),

    # network plot
    mainPanel(
      simpleNetworkOutput("RVineTree"),
      tableOutput("table")#,
      #htmlOutput("text")
    )
  )
))
