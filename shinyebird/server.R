library(ggplot2)
library(scales)
library(diptest)
library(shiny)
library(TTR)


shinyServer(function(input, output) {
  workspace = "/data"
  
  ebird = reactiveValues((data = Null))
  
  observeEvent(input$do, {
     print(as.numeric(input$do))
     withProgress(message = 'Loading dataset', value = 0, {
      # Input ArcGIS Model csv file
      infile = input$species
      inbird = paste(infile,".csv",sep="")
      print(inbird)
      ############################################################################
      # Read in ebird data
      temp = read.csv(paste(workspace,inbird,sep="/"), sep=",", header=TRUE, quote = "", stringsAsFactors = FALSE, na.strings=c(""))
      incProgress(0.8, detail = "Finished pulling in eBird.  Formatting")
      # Drop extra columns and data
      
      # Reorder months
      temp$Month = factor(temp$Month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))
      temp$Week = factor(temp$Week, levels=c(35:53,1:17))
      #Set X as na
      temp$OBSERVATION.COUNT = ifelse(temp$OBSERVATION.COUNT == "X", 1, temp$OBSERVATION.COUNT)
      temp$OBSERVATION.COUNT = as.numeric(temp$OBSERVATION.COUNT)
      temp$BCRNUMNAME = paste(temp$BCR.CODE, temp$BCRNAME, sep="_")
    })
    ebird() = temp
  })
  
  output$selectedSpecies = renderUI({
    df = ebird()
    items = unique(df$BCRNUMNAME)
    selectInput("bcr", "BCR:", items)
    })
  
  output$whichSpecies = renderText({
    input$species
  })
  
  computeSummary = reactive({
    df = ebird()
    df = subset(df, df$BCRNUMNAME == input$bcr)
    aggMean = aggregate(df$OBSERVATION.COUNT, list(Week=df$Week, BCR=df$BCR.CODE), mean)
    #plot(aggMean$Week, aggMean$OBSERVATION.COUNT)
    ggplot(aggMean, aes(x=aggMean$Week, y=x)) +
      stat_summary(aes(y = x,group=1), fun.y=mean, colour="red", geom="line",group=1) + 
      labs(y="Mean Observation count", x="Week number from the first week in September until the last week in April") +
      ggtitle(paste("Figure 2. Observation count mean by BCR plotted over wintering period for ", input$species, sep=""))
  })
  computeSmooth = reactive({
    df = ebird()
    df = subset(df, df$BCRNUMNAME == input$bcr)
    aggMean = aggregate(df$OBSERVATION.COUNT, list(Week=df$Week, BCR=df$BCR.CODE), mean)
    aggMean$smooth = SMA(aggMean[, "x"],3)
    #plot(aggMean$Week, aggMean$OBSERVATION.COUNT)
    ggplot(aggMean, aes(x=aggMean$Week, y=aggMean$smooth)) +
      stat_summary(aes(y = aggMean$smooth,group=1), fun.y=mean, colour="red", geom="line",group=1) + 
      labs(y="Mean Observation count", x="Week number from the first week in September until the last week in April") +
      scale_x_discrete(labels=c("Sept", 36:43, "November", 45:53, "Jan", 2:17)) +
      ggtitle(paste("Figure 2. Smoothed Observation count mean by BCR plotted over wintering period for ", input$species, sep=""))
  })
  
  computePVal = reactive({
    df = ebird()
    df = subset(df, df$BCRNUMNAME == input$bcr)
    testsetup = aggregate(df$OBSERVATION.COUNT, list(Week=df$Week, BCR=df$BCR.CODE, BCRNUMNAME = df$BCRNUMNAME), mean)
    testsmooth = SMA(testsetup[, "x"], 3)
    test = dip.test(testsmooth)
    bcr_name = unique(testsetup$BCRNUMNAME)
    paste("P-value:", test$p.value[[1]]," / BCR:", bcr_name, sep=" ")
  })
  
  output$statsTable = renderPlot({
    if(input$do == 0) return(NULL)
    computeSummary()
  })
  output$smoothTable = renderPlot({
    if(input$do == 0) return(NULL)
    computeSmooth()
  })
  
  output$pVal = renderText({
    if(input$do == 0) return(NULL)
    computePVal()
  })
  
})