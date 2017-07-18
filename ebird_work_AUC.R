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
workspace = "D:/ebird/newdata"
# list("abdu", "agwt", "amwi", "buff", "bwte", "canv", "cite", "coei", "gadw", "kiei", "ltdu", "mall", "nopi", "nsho", "redh", "rndu", "rudu", "scau", "wodu")
birdlist = list("abdu", "agwt", "amwi", "buff", "bwte", "canv", "cite", "coei", "gadw", "kiei", "ltdu", "mall", "nopi", "nsho", "redh", "rndu", "rudu", "scau", "wodu")
data <- read.csv(paste("D:/ebird/Fleming/", "correctedDailyharv19992014.csv", sep=""), na.strings = "")
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
  # Mean of ebird data by Week
  aggMean = aggregate(ebird$OBSERVATION.COUNT, list(Week=ebird$MonthDay, BCR=ebird$BCR.CODE, BCRName = ebird$BCRNAME, BCRNUMNAME = ebird$BCRNUMNAME), mean)
  # Iterate through ebird data by BCR
  for(i in unique(aggMean$BCR)) {
    sub = subset(aggMean, aggMean$BCR == i)
    if (nrow(sub) < 4) {
      next
    }
    ss =smooth.spline(x=sub$Week, y=sub$x, spar=0.7, keep.data = TRUE)
    ss$x = aggMean$Week
    
    # Create ebird curve and subset county level stepdown data
    ss$ypct = ss$y/max(ss$y)*100
    capspecies = toupper(species)
    newsub = subset(popobj, popobj$species == capspecies)
    newsub = aggregate(cbind(LTAPopObj, X80percPopObj)~fips+species+state+countyname+BCR+CODE, data=newsub, sum, na.rm=TRUE)
    
    # Run dip test
    dt = dip.test(ss$yin)
    pval =  dt$p.value[[1]]
    
    # TEST AUC Comparison for b/d classification
    autumn = subset(sub, (as.Date(sub$Week, format="%m/%d") >= as.Date("9/1", format="%m/%d") & (as.Date(sub$Week, format="%m/%d")) <= as.Date("11/30", format="%m/%d")))
    if (nrow(autumn) < 1) {
      autumn = 0
    } else{
      autumn = aggregate(autumn$x, list(BCR=autumn$BCR), sum)
      autumn = autumn$x[1]
    }
    winter = subset(sub, (as.Date(sub$Week, format="%m/%d") >= as.Date("12/1", format="%m/%d") | (as.Date(sub$Week, format="%m/%d")) <= as.Date("1/31", format="%m/%d")))
    if (nrow(winter) < 1) {
      winter = 0
    } else {
      winter = aggregate(winter$x, list(BCR=winter$BCR), sum)
      winter = winter$x[1]
    }
    
    classification = ''
    if (autumn > winter) {
      classification = '4b'
    } else {
      classification = '4d'
    }
    
    ## Create classification based off of stepdowns
    sumharvest <- aggregate(data$Harvest, by=list(data$lumpedAOU, data$ST, data$fips,data$Season), sum)  # sum extracted harvest by county
    sumharvest = subset(sumharvest, sumharvest$Group.1 == capspecies)
    names(sumharvest) <- c("spp","ST", "fips","season", "harvest")
    if (nrow(sumharvest) > 0) {
      harvestfall = sumharvest[sumharvest$season=="fall",]
      harvestfall$fallharvest = harvestfall$harvest
      harvestwinter = sumharvest[sumharvest$season=="winter",]
      harvestwinter$winterharvest = harvestwinter$harvest
      mergeharvest = merge(harvestfall, harvestwinter, by=c("spp", "fips"))
      mergeharvest$HarvestCode = ifelse(mergeharvest$fallharvest > mergeharvest$winterharvest, "4b", "4d")
      newsub$HarvestCode = mergeharvest[match(newsub$fips, mergeharvest$fips),11]
      newsub$HarvestCode = ifelse(is.na(newsub$HarvestCode), 0, newsub$HarvestCode)  
    } else {
      newsub$HarvestCode = 0
    }
    
    newOutput = list(species, pval, list(ss$y), list(ss$ypct), i, classification)
    
    #Create new dataframe from the list of species curves by bcr
    noDF = data.frame(newOutput)
    #test = data.frame(species, ifelse(dt$p.value[[1]]<0.05, "4b", "4d"), ss$y, ss$ypct)
    names(noDF)[1] = "species"
    names(noDF)[2] = "pval"
    #names(noDF)[3] = "ssX"
    names(noDF)[3] = "ssY"
    names(noDF)[4] = "poppct"
    names(noDF)[5] = "BCR"
    names(noDF)[6] = "eBirdClass"
    noDF$species = toupper(noDF$species)
    # noDF$unique = paste(noDF$species, noDF$CODE, noDF$BCR)
    # newsub$unique = paste(newsub$species, newsub$CODE, newsub$BCR)
    # Merge curve data with the stepdown data
    spauc = merge(noDF, newsub, by=c("species", "BCR"))
    spauc$poppct = ifelse(spauc$poppct < 0, 0, spauc$poppct)
    spauc$LTAPopTot = spauc$poppct * .01 * spauc$LTAPopObj
    spauc$X80PopTot = spauc$poppct * .01 * spauc$X80percPopObj
    
    #Merge current speices/BCR to aggregate DF
    if (nrow(spauc) > 0) {
      aggItAll = aggregate(cbind(LTAPopTot, X80PopTot)~fips+species+state+countyname+BCR+CODE+eBirdClass+pval+HarvestCode, data=spauc, sum, na.rm=TRUE)
      outTotal = rbind(outTotal, aggItAll)
      outCurve = rbind(outCurve, data.frame(ss$y, ss$ypct, BCR = i))
    }
  }
  outTest = outTotal[(outTotal$CODE == outTotal$eBirdClass),]
  #outTest = outTest[!(outTest$pval >= 0.05 & outTest$CODE == "4b"),]
  outTotal = outTest
  outTotal = subset(outTotal, (outTotal$LTAPopTot != 0) & (outTotal$X80PopTot != 0))
  outCurve$species = species
  write.csv(outCurve, file=paste(workspace, "/output/",species, "_Curve.csv", sep=""), row.names = F)
  write.csv(outTotal, file=paste(workspace, "/output/",species, "_Output.csv", sep=""), row.names = F)
}