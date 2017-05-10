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
# list("abdu", "agwt", "amwi", "bwte", "canv", "cite", "gadw", "mall", "nopi", "nsho", "redh", "rndu", "rudu", "scau", "wodu")
birdlist = list("abdu", "agwt", "amwi", "bwte", "canv", "cite", "gadw", "mall", "nopi", "nsho", "redh", "rndu", "rudu", "scau", "wodu")
for (sp in 1:length(birdlist)) {
  species = birdlist[[sp]]
  print(species)
  # Input ArcGIS Model csv file
  inbird = paste(species, ".csv", sep="")
  ############################################################################
  # Read in ebird data (16 species)
  ebird <- read.csv(paste(workspace,inbird,sep="/"), sep=",", header=TRUE, quote = "", stringsAsFactors = FALSE, na.strings=c(""))
  popobj <- read.csv(paste(workspace,"PopObj2.csv",sep="/"), sep=",", header=TRUE, quote = "", stringsAsFactors = FALSE, na.strings=c(""))
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
  
  outTotal = data.frame()
  outCurve = data.frame()
  # Graph mean
  aggMean = aggregate(ebird$OBSERVATION.COUNT, list(Week=ebird$MonthDay, BCR=ebird$BCR.CODE, BCRName = ebird$BCRNAME, BCRNUMNAME = ebird$BCRNUMNAME), mean)
  
  for(i in unique(aggMean$BCR)) {
    sub = subset(aggMean, aggMean$BCR == i)
    if (nrow(sub) < 4) {
      next
    }
    ss =smooth.spline(x=sub$Week, y=sub$x, spar=0.7, keep.data = TRUE)
    ss$x = aggMean$Week
    
    #### NEW STUFF ####
    ss$ypct = ss$y/max(ss$y)*100
    capspecies = toupper(species)
    newsub = subset(popobj, popobj$species == capspecies)
    newsub = aggregate(cbind(LTAPopObj, X80percPopObj)~fips+species+state+countyname+BCR+CODE, data=newsub, sum, na.rm=TRUE)
    
    dt = dip.test(ss$yin)
    pval =  dt$p.value[[1]]
    newOutput = list(species, pval, list(ss$y), list(ss$ypct), i)
    noDF = data.frame(newOutput)
    #test = data.frame(species, ifelse(dt$p.value[[1]]<0.05, "4b", "4d"), ss$y, ss$ypct)
    names(noDF)[1] = "species"
    names(noDF)[2] = "pval"
    #names(noDF)[3] = "ssX"
    names(noDF)[3] = "ssY"
    names(noDF)[4] = "poppct"
    names(noDF)[5] = "BCR"
    noDF$species = toupper(noDF$species)
    # noDF$unique = paste(noDF$species, noDF$CODE, noDF$BCR)
    # newsub$unique = paste(newsub$species, newsub$CODE, newsub$BCR)
    spauc = merge(noDF, newsub)
    spauc$poppct = ifelse(spauc$poppct < 0, 0, spauc$poppct)
    spauc$LTAPopTot = spauc$poppct * .01 * spauc$LTAPopObj
    spauc$X80PopTot = spauc$poppct * .01 * spauc$X80percPopObj
    if (nrow(spauc) > 0) {
      aggItAll = aggregate(cbind(LTAPopTot, X80PopTot)~fips+species+state+countyname+BCR+CODE+pval, data=spauc, sum, na.rm=TRUE)
      outTotal = rbind(outTotal, aggItAll)
      outCurve = rbind(outCurve, data.frame(ss$y, ss$ypct, BCR = i))
    }
  }
  outTest = outTotal[!(outTotal$pval < 0.05 & outTotal$CODE == "4d"),]
  outTest = outTest[!(outTest$pval >= 0.05 & outTest$CODE == "4b"),]
  outTotal = outTest
  outTotal = subset(outTotal, (outTotal$LTAPopTot != 0) & (outTotal$X80PopTot != 0))
  outCurve$species = species
  write.csv(outCurve, file=paste(workspace, "/output/",species, "_Curve.csv", sep=""), row.names = F)
  write.csv(outTotal, file=paste(workspace, "/output/",species, "_Output.csv", sep=""), row.names = F)
}
