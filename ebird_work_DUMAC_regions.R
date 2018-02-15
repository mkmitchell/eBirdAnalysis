library(ggplot2)
library(scales)
library(diptest)
library(shiny)
library(TTR)

############################################################################
# Setup environment and variables
############################################################################
# Variable designiation
# Workspace directory
workspace = "D:/ebird/datadownload/post/"
#birdlist = list("abdu", "agwt", "amwi", "bago", "bbwd", "bcte", "blsc", "buff", "bwte", "canv", "cite", "coei", "cogo", "come", "eidr", "gadw", "home", "kiei", "lesc", "ltdu", "mall", "merg", "nopi", "nsho", "rbme", "redh", "rndu", "rudu", "scau", "scot", "susc", "wodu", "wwsc")
birdlist = list("abdu", "agwt", "amwi", "bago", "bbwd", "bcte", "blsc", "buff", "bwte", "canv", "cite", "coei", "cogo", "come", "eidr", "gadw", "home", "kiei", "lesc", "ltdu", "mall", "merg", "nopi", "nsho", "rbme", "redh", "rndu", "rudu", "scau", "scot", "susc", "wodu", "wwsc")
bcr = "D:/ebird/datadownload/BCR.csv"
bcr = read.csv(bcr, header=TRUE)
runextra = 0
for(sp in 1:length(birdlist)) {
  inbird = birdlist[[sp]]
  print(inbird)
  ############################################################################
  # Read in ebird data (16 species)
  ebird <- read.csv(paste(workspace,inbird,".csv",sep=""), header=TRUE)
  # Join arcData to coverData based on cover_type
  #ebird = merge(ebird, bcr, by.x = "BCR.CODE", by.y = "BCR")
  setwd(workspace)
  
  # Drop extra columns and data
  ebird = subset(ebird, ebird$COUNTRY_CODE == "MX")
  if (nrow(ebird) < 4) { next}
  #ebird = subset(ebird, ebird$APPROVED == "1")
  ebird$PROJECT.CODE = NULL
  ebird$PROTOCOL.TYPE = NULL
  ebird$SAMPLING.EVENT.IDENTIFIER = NULL
  ebird$FIRST.NAME = NULL
  ebird$LAST.NAME = NULL
  ebird$OBSERVER.ID = NULL
  ebird$BREEDING.BIRD.ATLAS.CODE = NULL
  ebird$GLOBAL.UNIQUE.IDENTIFIER = NULL
  ebird$SUBSPECIES.COMMON.NAME = NULL
  ebird$SUBSPECIES.SCIENTIFIC.NAME = NULL
  ebird$AGE.SEX = NULL
  ebird$COUNTRY = NULL
  ebird$IBA.CODE = NULL
  ebird$SUBNATIONAL1_CODE = NULL
  ebird$SUBNATIONAL2_CODE = NULL
  ebird$ATLAS.BLOCK = NULL
  ebird$LOCALITY = NULL
  ebird$LOCALITY.ID = NULL
  ebird$LOCALITY.TYPE = NULL
  ebird$REVIEWED = NULL
  ebird$REASON = NULL
  ebird$TRIP.COMMENTS = NULL
  ebird$EFFORT.AREA.HA = NULL
  ebird$EFFORT.DISTANCE.KM = NULL
  ebird$SPECIES.COMMENTS = NULL
  ebird$X = NULL
  ebird$APPROVED = NULL
  ebird$TIME.OBSERVATIONS.STARTED = NULL
  ebird$DURATION.MINUTES = NULL
  ebird$NUMBER.OBSERVERS = NULL
  ebird$ALL.SPECIES.REPORTED = NULL
  ebird$GROUP.IDENTIFIER = NULL
  ebird$TAXONOMIC.ORDER = NULL
  ebird$CATEGORY = NULL
  
  # Set date field and divide up by year, month, week for plotting
  ebird$Date = as.Date(ebird$OBSERVATION.DATE, "%Y-%m-%d")
  ebird$Year = strtoi(format(ebird$Date, "%Y"))
  to.month = function(x) as.integer(format(x, "%m"))
  to.day <- function(x) as.integer(format(x, "%d"))
  ebird$Month = to.month(ebird$Date)
  ebird$Day = to.day(ebird$Date)
  ebird$MonthDay = paste(ebird$Month, ebird$Day, sep="/")
  # Remove years < 2005 and months 6,7,8
  ebird = subset(ebird, ebird$Year > 2005 & ebird$Year <= 2016)
  ebird = subset(ebird, ebird$Month >= 9 | ebird$Month <= 4)
  
  ebird$Winter = ifelse(ebird$Month <= 4, paste(ebird$Year - 1,ebird$Year, sep = "/"),paste(ebird$Year, ebird$Year + 1, sep = "/"))
  ebird$Week = as.numeric(format(ebird$Date, "%U"))
  ebird = subset(ebird, ebird$Winter != "2005/2006" & ebird$Winter != "2016/2017")
  
  # Reorder months
  ebird$Month = factor(ebird$Month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))
  ebird$Week = factor(ebird$Week, levels=c(35:53,0:17))
  test = c()
  getFactor = function(x) {
    for (m in 9:12){
      if (m %in% c(9,11)) {
        maxVal = 30
      } else {
        maxVal = 31
      }
      for (d in 1:maxVal){
        test = append(test, paste(m,d,sep="/"))
      }
    }
    for (m in 1:4){
      if (m == 2) {
        maxVal = 28
      } else if (m == 4) {
        maxVal = 30
      } else {
        maxVal = 31
      }
      for (d in 1:maxVal){
        test = append(test, paste(m,d,sep="/"))
      }
    }
    return(test)
  }
  ebird$MonthDay = factor(ebird$MonthDay,levels=(getFactor()))
  
  
  #Set X as na
  ebird$OBSERVATION.COUNT = ifelse(ebird$OBSERVATION.COUNT == "X", 1, ebird$OBSERVATION.COUNT)
  ebird$OBSERVATION.COUNT = as.numeric(ebird$OBSERVATION.COUNT)
  
  #Agg all MX
  print("All BCR")
    aggMean = aggregate(ebird$OBSERVATION.COUNT, list(Week=ebird$MonthDay), mean)
    if (nrow(aggMean) > 3) {
      ss=smooth.spline(x=aggMean$Week, y=aggMean$x, spar=0.7, keep.data = TRUE)
      ss$x = aggMean$Week
      plot = plot(x=ss$x, y=ss$y, type="l", 
                  ylab="Average count",
                  xlab="Date",
                  cex.lab=1.5
      )
      title(paste("Mexico migration curve - ", inbird, sep=""))
      lines(x=ss$x,y=ss$y, col="red")
      dev.copy(jpeg,paste("D:/ebird/Mexico/Migration_", inbird, "_Mean.jpg",sep=""))
      dev.off()
    }
    else 
    {
      print(paste("not enough ", inbird, sep=""))
  }
  
  if (runextra == 1){
    # Graph mean
    gulfcoast = subset(ebird, ebird$BCR.CODE %in% c(36,49,52,55,65,66))
    print("Gulf Coast")
    if (nrow(gulfcoast) > 3){
      aggMean = aggregate(gulfcoast$OBSERVATION.COUNT, list(Week=gulfcoast$MonthDay), mean)
      if (nrow(aggMean > 3)) {
        ss=smooth.spline(x=aggMean$Week, y=aggMean$x, spar=0.7, keep.data = TRUE)
        ss$x = aggMean$Week
        plot = plot(x=ss$x, y=ss$y, type="l", 
             ylab="Average count",
             xlab="Date",
             cex.lab=1.5
        )
        title(paste("Gulf Coast of Mexico - ", inbird, sep=""))
        lines(x=ss$x,y=ss$y, col="red")
        dev.copy(jpeg,paste("D:/ebird/Mexico/Migration_", inbird, "_GulfCoast_Mean.jpg",sep=""))
        dev.off()
      }
    } else {
      print(paste("not enough ", inbird, sep=""))
    }
    
    pacmx = subset(ebird, ebird$BCR.CODE %in% c(33,40,43,44))
    print("Pac")
    if (nrow(pacmx) > 3){
      aggMean = aggregate(pacmx$OBSERVATION.COUNT, list(Week=pacmx$MonthDay), mean)
      if (nrow(aggMean) > 3)) {
        ss=smooth.spline(x=aggMean$Week, y=aggMean$x, spar=0.7, keep.data = TRUE)
        ss$x = aggMean$Week
        plot = plot(x=ss$x, y=ss$y, type="l", 
                    ylab="Average count",
                    xlab="Date",
                    cex.lab=1.5
        )
        title(paste("Pacific Coast of Mexico - ", inbird, sep=""))
        lines(x=ss$x,y=ss$y, col="red")
        dev.copy(jpeg,paste("D:/ebird/Mexico/Migration_", inbird, "_PacMX_Mean.jpg",sep=""))
        dev.off()
      }
    }
    nchigh = subset(ebird, ebird$BCR.CODE %in% c(34,35,46,47))
    print("Highlands")
    if (nrow(nchigh) > 3){
      aggMean = aggregate(nchigh$OBSERVATION.COUNT, list(Week=nchigh$MonthDay), mean)
      if (nrow(aggMean > 3)) {
        ss=smooth.spline(x=aggMean$Week, y=aggMean$x, spar=0.7, keep.data = TRUE)
        ss$x = aggMean$Week
        plot = plot(x=ss$x, y=ss$y, type="l", 
                    ylab="Average count",
                    xlab="Date",
                    cex.lab=1.5
        )
        title(paste("Northern and Central Highlands - ", inbird, sep=""))
        lines(x=ss$x,y=ss$y, col="red")
        dev.copy(jpeg,paste("D:/ebird/Mexico/Migration_", inbird, "_Highlands_Mean.jpg",sep=""))
        dev.off()
      }
    }
  }
}


