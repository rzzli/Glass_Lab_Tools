library(shiny)
require(dplyr)
require(Seurat)

#wd<-'/data/singleCell/RNA/'
#setwd(wd)
#strRDS <- sort(list.files(wd,"rds$"))
#strIDs <- gsub("\\.rds","",strRDS)
#names(strRDS) <- strIDs
# Define server logic required to draw a histogram
wd<-'/gpfs/data01/glasslab/home/zhl022/shiny/CH_shiny/db/'


shinyServer(function(input, output) {
  
  observe({
    if (input$go_rds !=0 & input$dataset != 'select') {
      #datasetInput <- reactive({
        #strRDS[input$dataset]
      #})
        X<<- readRDS(paste(wd,input$dataset,sep=""))
        allgenes<<-X@assays$RNA@counts@Dimnames[1]
        output$status_load_rds=renderText('rds loaded')
        output$distPlot1 <- renderPlot({
          g <- DimPlot(object = X, label = T, pt.size = 1)
          print(g)
        })
        observeEvent(input$status_gene_names, {
          names <-  unlist(strsplit(input$genelist, ","))
          overlap <-setdiff(names, allgenes)

          if(length(overlap)==0){
            output$status_gene_names=renderText("successfully upload gene names")

          } else if (length(names)>9) {
            output$status_gene_names=renderText("gene names must be <=9")
          }
          else {
            output$status_gene_names=renderText(toString(overlap))
          }
          
          output$distPlot2 <-renderPlot({
            print(FeaturePlot(X,features = unlist(strsplit(input$genelist, ","))))
          })
          output$distPlot3 <-renderPlot({
            print(VlnPlot(X, features = unlist(strsplit(input$genelist, ","))))
          })
        })
    }
  })
})

