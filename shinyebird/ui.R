shinyUI(fluidPage(
  
  titlePanel("eBird analysis"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("species", "Species:", c("American green winged teal" = "agwt","Black Duck" = "abdu", "American Widg" = "amwi", "Blue winged teal" = "bwte", "Canvasback" = "canv", "Cin teal" = "citi", "Gadwall"= "gadw", "Great Scaup" = "grsc", "Lesser" = "lesc", "Mallard" = "mall", "Pintail" = "nopi", "Shoveler" = "nsho", "Readhead"= "redh", "Ringed neck" = "rndu", "Ruddy Duck" = "rudu", "Wood duck" = "wodu")),
      uiOutput("selectedSpecies"),
      actionButton("do", "Run Analysis")
    ),
      
    mainPanel(
      h3(textOutput("whichSpecies")),
      plotOutput(("statsTable")),
      textOutput(("pVal"))
    )
  )
))