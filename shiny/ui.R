library(shiny)
require(Seurat)

#source("global.R")
#strIDs <- gsub("\\.rds","",sort(list.files(wd,"rds$")))
wd <-'/gpfs/data01/glasslab/home/zhl022/shiny/CH_shiny/db'
strIDs <- sort(list.files(wd,"rds$"))  # list of rds files

shinyUI(fluidPage(
  
  # App title ----
  titlePanel("Hello scRNA-seq"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    
    sidebarPanel(
      # Input: Slider for the number of bins ----
      selectInput(inputId = "dataset",
                  label = "Choose a dataset:",
                  choices = c('select',strIDs),
                  selected= 'select'),
      actionButton('go_rds','Confirm dataset'),
      textOutput('status_load_rds'),      
      
      textInput("genelist", label = h3("enter genes names to visulize, separate by comma"), value = ""),
      actionButton('status_gene_names', "confirm gene names"),
      textOutput('status_gene_names')
      ),
      

    
    
    
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Histogram ----
      tabsetPanel(
        tabPanel("Cluster", plotOutput(outputId = "distPlot1",width = "90%",height = "800px")),
        tabPanel("Genes-tsne", plotOutput(outputId = "distPlot2",width = "90%",height = "800px")),
        tabPanel("Genes-violin", plotOutput(outputId = "distPlot3",width = "90%",height = "800px"))
        )
    )
  )
  )
)

