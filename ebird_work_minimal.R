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
popobj <- read.csv(paste(workspace,"PopObj.csv",sep="/"), sep=",", header=TRUE, quote = "", stringsAsFactors = FALSE, na.strings=c(""))
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

# Graph mean
sub = subset(ebird, ebird$BCR.CODE == 26)
aggMean = aggregate(sub$OBSERVATION.COUNT, list(Week=sub$MonthDay, BCR=sub$BCR.CODE, BCRName = sub$BCRNAME, BCRNUMNAME = sub$BCRNUMNAME), mean)
ss =smooth.spline(x=aggMean$Week, y=aggMean$x, spar=0.7, keep.data = TRUE)
ss$x = aggMean$Week
plot(x=ss$x, y=ss$y, type="l", 
     ylab="Average count",
     xlab="Date",
     cex.lab=1.5
)
lines(x=ss$x,y=ss$y, col="red")
test = dip.test(ss$yin)
print(test$p.value[[1]])

#ggsave(file=paste(species, "_Mean.jpg",sep=""), plot=plot, width=8, height=4)

print("P values < 0.05 indicate significant bimodality and values greater than 0.05 but less than 0.10 suggest bimodality with marginal significance")
for(i in unique(aggMean$BCR)) {
  testsetup = aggregate(ebird$OBSERVATION.COUNT, list(Week=ebird$Week, BCR=ebird$BCR.CODE, BCRName = ebird$BCRNAME), mean)
  testsetup = testsetup[which(testsetup$BCR == i),]
  test = dip.test(testsetup$x)
  bcr_name = unique(testsetup$BCRName)
  print(paste("P-value:", test$p.value[[1]]," / BCR:", bcr_name, sep=" "))
}
max(ss$y)
ss$ypct = ss$y/max(ss$y)*100
capspecies = toupper(species)
newsub = subset(popobj, popobj$species == capspecies)
dt = dip.test(ss$yin)
newOutput = list(species, ifelse(dt$p.value[[1]]<0.05, "4b", "4d"), list(ss$x), list(ss$y), list(ss$ypct))
noDF = data.frame(newOutput)
names(noDF)[1] = "species"
names(noDF)[2] = "code"
names(noDF)[3] = "ssX"
names(noDF)[4] = "ssY"
names(noDF)[5] = "poppct"
ebird = merge(newOutput, popobj, by.x = "species", by.y = "species")
