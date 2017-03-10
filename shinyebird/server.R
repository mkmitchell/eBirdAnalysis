library(ggplot2)
library(scales)
library(diptest)
library(shiny)


shinyServer(function(input, output) {
  workspace = "D:/ebird"
  bcr = "BCR.csv"
  
  observeEvent(input$do, {
     print(as.numeric(input$do))
  })
  
  ebird = reactive({
    # Input ArcGIS Model csv file
    infile = input$species
    print(infile)
    inbird = paste(infile,".txt",sep="")
    print(inbird)
    ############################################################################
    # Read in ebird data
    tempin = read.delim(paste(workspace,inbird,sep="/"), sep="\t", header=TRUE, quote = "", stringsAsFactors = FALSE, na.strings=c(""))
    bcrData = read.csv(paste(workspace, bcr, sep="/"), header=TRUE)
    merge(tempin, bcrData, by.x = "BCR.CODE", by.y = "BCR")
  })
  
  output$selectedSpecies = renderUI({
    df = ebird()
    items = unique(df$BCRNAME)
    selectInput("bcr", "BCR:", items)
    })
  
})