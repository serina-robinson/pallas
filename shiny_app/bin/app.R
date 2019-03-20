---
output: html_document
runtime: shiny
---
library(shiny)
library(DECIPHER)
library(Biostrings)

# Define UI for data upload app ----
ui <- fluidPage(

  # App title ----
  titlePanel("PALLAS: Prediction tool for Adenylation, Ligase, and Luciferase Activity and Substrate specificity"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Select a file ----
      fileInput("file1", "Choose a FASTA File",
                multiple = FALSE),

      # Horizontal line ----
      tags$hr(),

      # Input: Select separator ----
      radioButtons("seqtype", "Sequence type",
                   choices = c(Nucleotide = "nuc",
                               Protein = "aa"),
                   selected = "aa")
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Data file ----
      tableOutput("contents")

    )

  )
)

# Define server logic to read selected file ----
server <- function(input, output) {

  output$contents <- renderTable({

    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    req(input$file1)

    # when reading semicolon separated files,
    # having a comma separator causes `read.csv` to error
    tryCatch(
      {
        if(input$seqtype == "aa") {
          return(readAAStringSet(input$file1$datapath))
        }
        else {
          return(readDNAStringSet(input$file1$datapath))
        }
        
      })
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )



}

# Create Shiny app ----
shinyApp(ui, server)
