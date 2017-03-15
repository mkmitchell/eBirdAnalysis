shinyUI(fluidPage(
  
  titlePanel("Waterfowl eBird analysis"),
  
  sidebarLayout(
    sidebarPanel(
      "Please choose a species and wait for the BCR list and species code to update.  Sometimes it takes awhile to load the datasets.  When ready, the fields will not be greyed out.  You only need to hit run Analysis once.  The rest will be dynamic", br(),br(),
      selectInput("species", "Species:", c("American green winged teal" = "agwt","Black Duck" = "abdu", "American Widg" = "amwi", "Blue winged teal" = "bwte", "Canvasback" = "canv", "Cin teal" = "citi", "Gadwall"= "gadw", "Great Scaup" = "grsc", "Lesser" = "lesc", "Mallard" = "mall", "Pintail" = "nopi", "Shoveler" = "nsho", "Readhead"= "redh", "Ringed neck" = "rndu", "Ruddy Duck" = "rudu", "Wood duck" = "wodu")),
      uiOutput("selectedSpecies"),
      actionButton("do", "Run Analysis")
    ),
      
    mainPanel(
      strong(
      div(style="display:inline-block","Species Code: "),  
      div(style="display:inline-block; text-transform:uppercase", textOutput("whichSpecies"))),
      plotOutput(("statsTable")),
      plotOutput(("smoothTable")),br(),br(),
      "Test statistic increases with multimodality so a lower value means you have a unimodal distribution.",br(),
      strong(textOutput(("pVal"))),br(),br(),
      "The test statistic comes from dip.test: https://cran.r-project.org/web/packages/diptest/diptest.pdf"
      
    )
  )
))