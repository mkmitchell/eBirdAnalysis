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
workspace = "D:/ebird/data"
species = "mall"
# Input ArcGIS Model csv file
inbird = "mall.csv"
############################################################################
# Read in ebird data (16 species)
ebird <- read.csv(paste(workspace,inbird,sep="/"), sep=",", header=TRUE, quote = "", stringsAsFactors = FALSE, na.strings=c(""))
setwd(workspace)

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
ebird$BCRNUMNAME = paste(ebird$BCR.CODE, ebird$BCRNAME, sep="_")
ebird = subset(ebird, ebird$STATE_PROVINCE == 'Missouri')
ebird$huntzone = ifelse(ebird$LATITUDE >= 39.2, 'N',ifelse(ebird$LATITUDE <=37.5, 'S', 'M'))

# Graph mean
aggMean = aggregate(ebird$OBSERVATION.COUNT, list(Week=ebird$MonthDay, Huntzone = ebird$huntzone), mean)
for(i in unique(aggMean$Huntzone)){
  print(i)
  temp = subset(aggMean, aggMean$Huntzone == i)
  ss=smooth.spline(x=temp$Week, y=temp$x, spar=0.7, keep.data = TRUE)
  ss$x = temp$Week
  plot(x=ss$x, y=ss$y, type="l", 
       ylab="Average count",
       xlab="Date",
       cex.lab=1.5
  )
  title(i)
  lines(x=ss$x,y=ss$y, col="red")
  df = data.frame(ss$x, ss$y)
  write.csv(df, file=paste("D:/MOMALLSmooth_", i, ".csv", sep=""))
}

