library(ggplot2)
library(scales)
library(shiny)
library(TTR)

############################################################################
# Setup environment and variables
############################################################################
# Variable designiation
# Workspace directory
workspace = "D:/ebird/Mexico/Analysis/"
# list("agwt", "amwi","bcte","gadw", "mall", "nopi", "nsho", "redh", "rndu", "rudu", "scau")
birdlist = list("agwt", "amwi","bcte","gadw", "mall", "nopi", "nsho", "redh", "rndu", "rudu", "scau")
# Fleming et al. R code input data.  I'm not sure the source of this file.  Canada data was added to the original file format schema.
harvestdata <- read.csv(paste("D:/ebird/Mexico/Analysis/", "MX_harvest.csv", sep=""), na.strings = "")
sumharvestin <- aggregate(harvestdata$Harvest, by=list(harvestdata$lumpedAOU, harvestdata$fips), sum)  # sum extracted harvest by county
CombineItAll = data.frame()

#loop that does work for each bird in the bird list specified above.  All Species.csv files were derived from Fleming et al. R code output.
for (sp in 1:length(birdlist)) {
  species = birdlist[[sp]]
  print(species)
  capspecies = toupper(species)
  # Input ArcGIS Model csv file.  BWTE and CITE are both treated as BCTE.  Work is done slightly differently for these two.
  inbird = paste(species, ".csv", sep="")
  ############################################################################
  # Read in ebird data (16 species)
  ebird <- read.csv(paste(workspace,"ebird",inbird,sep="/"), sep=",", header=TRUE, quote = "", stringsAsFactors = FALSE, na.strings=c(""))
  popobj <- read.csv(paste(workspace,"Output","Mexico_1_2_18.csv",sep="/"), sep=",", header=TRUE, quote = "", stringsAsFactors = FALSE, na.strings=c(""))
  setwd(workspace)
  #Setup harvest and bcr input
  sumharvest = subset(sumharvestin, sumharvestin$Group.1 == capspecies)
  names(sumharvest) <- c("spp","fips","harvest")
  popobj = subset(popobj, popobj$species == capspecies)
  ebird = subset(ebird, ebird$COUNTRY_CODE == "MX")
  harvestbcr = merge(sumharvest, popobj, by.x = "fips", by.y = "fips")
  
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
  print(unique(aggMean$BCR))
  for(i in unique(aggMean$BCR)) {
    # #if (i != 23) {
    # #  next
    # #}
    # if (i ==15 ){
    #   sub = subset(aggMean, aggMean$BCR == 32)
    # } else {
    #   sub = subset(aggMean, aggMean$BCR == i)
    # }
    # if (nrow(sub) < 4) {
    #   next
    # }
    ssorg =smooth.spline(x=aggMean$Week, y=aggMean$x, spar=0.7, keep.data = TRUE)
    ss = stats:::predict.smooth.spline(ssorg, unique(sort(c(seq(1,242, by = 1)))))
    ss$x = c(getFactor())[0:242]
    ss$y[ss$y < 0] = 0
    print(i)
    #Classify as 4b or 4d
    ## Create classification based off of stepdowns
    max = 0
    harvestbcrsp = aggregate(harvestbcr$harvest, list(Species=harvestbcr$spp, BCR=harvestbcr$BCR), sum)
    harvestclass = ''
    if (nrow(harvestbcrsp) > 0) {
      harvestclass = "4d"
      seasonsub = subset(ss$y, (as.Date(ss$x, format="%m/%d") >= as.Date("12/1", format="%m/%d") | (as.Date(ss$x, format="%m/%d")) <= as.Date("3/30", format="%m/%d")))
      #seasonspline = smooth.spline(x=seasonsub$Week, y=seasonsub$x, spar=0.7, keep.data = TRUE)
      max = max(seasonsub)
    } else {
      print("not enough harvestbcr")
      next
    }
   
    # Create ebird curve and subset county level stepdown data
    ss$ypct = ss$y/max*100
    newsub = subset(popobj, popobj$species == capspecies)
    newsub = aggregate(cbind(LTAPopObj, X80percPopObj)~fips+species+BCR+CODE, data=newsub, sum, na.rm=TRUE)
    
    newOutput = list(species, list(ss$y), list(ss$ypct), i)
    
    #Create new dataframe from the list of species curves by bcr
    noDF = data.frame(newOutput)
    #test = data.frame(species, ifelse(dt$p.value[[1]]<0.05, "4b", "4d"), ss$y, ss$ypct)
    names(noDF)[1] = "species"
    names(noDF)[2] = "ssY"
    names(noDF)[3] = "poppct"
    names(noDF)[4] = "BCR"
    noDF$species = toupper(noDF$species)
    # noDF$unique = paste(noDF$species, noDF$CODE, noDF$BCR)
    # newsub$unique = paste(newsub$species, newsub$CODE, newsub$BCR)
    # Merge curve data with the stepdown data
    spauc = merge(noDF, newsub, by=c("species", "BCR"))
    spauc$poppct = ifelse(spauc$poppct < 0, 0, spauc$poppct)
    spauc$LTAPopTot = spauc$poppct * .01 * spauc$LTAPopObj
    spauc$X80PopTot = spauc$poppct * .01 * spauc$X80percPopObj
    
    #Merge current species/BCR to aggregate DF
    if (nrow(spauc) > 0) {
      aggItAll = aggregate(cbind(LTAPopTot, X80PopTot)~fips+species+BCR+CODE, data=spauc, sum, na.rm=TRUE)
      outTotal = rbind(outTotal, aggItAll)
      outCurve = rbind(outCurve, data.frame(ss$y, ss$ypct, BCR = i))
    }
  } #End BCR Loop

  #Classify as 4b or 4d
  ## Create classification based off of stepdowns
  outTotal = subset(outTotal, (outTotal$LTAPopTot != 0) & (outTotal$X80PopTot != 0))
  outTotal =  merge(outTotal, popobj[,c("species", "fips", "LTAPopObj", "X80percPopObj")], by=c("species", "fips"))
  CombineItAll = rbind(CombineItAll, outTotal)
  outCurve$species = species
  write.csv(outCurve, file=paste(workspace, "/output/",species, "_Curve.csv", sep=""), row.names = F)
  write.csv(outTotal, file=paste(workspace, "/output/",species, "_Output.csv", sep=""), row.names = F)
} #End species loop
write.csv(CombineItAll, file=paste(workspace, "/output/All_species", "_Output.csv", sep=""), row.names = F)
