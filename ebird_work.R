library(ggplot2)
library(scales)
library(diptest)
############################################################################
# Setup environment and variables
############################################################################
# Variable designiation
# Workspace directory
workspace = "D:/ebird"
species = "AGWT"
# Input ArcGIS Model csv file
inbird = "ebd_US_agwtea1_relNov-2016.txt"
bcr = "BCR.csv"
############################################################################
# Read in ebird data (16 species)
ebird <- read.delim(paste(workspace,inbird,sep="/"), sep="\t", header=TRUE, quote = "", stringsAsFactors = FALSE, na.strings=c(""))
bcr = read.csv(paste(workspace, bcr, sep="/"), header=TRUE)
# Join arcData to coverData based on cover_type
ebird = merge(ebird, bcr, by.x = "BCR.CODE", by.y = "BCR")
setwd(workspace)

# Drop extra columns and data
ebird = subset(ebird, ebird$COUNTRY_CODE == "US")
ebird = subset(ebird, ebird$APPROVED == "1")
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
ebird$Month = strtoi(format(ebird$Date, "%m"))
ebird$Day = strtoi(format(ebird$Date, "%d"))
ebird$Month = ifelse(format(ebird$Date, "%m") == "09", 9, ebird$Month)
# Remove years < 2005 and months 6,7,8
ebird = subset(ebird, ebird$Year > 2005 & ebird$Year <= 2016)
ebird = subset(ebird, ebird$Month >= 9 | ebird$Month <= 4)

ebird$Winter = ifelse(ebird$Month <= 4, paste(ebird$Year - 1,ebird$Year, sep = "/"),paste(ebird$Year, ebird$Year + 1, sep = "/"))
ebird$Week = as.numeric(format(ebird$Date, "%U"))
#ebird$Month = ifelse(ebird$Month <= 5, substring(ebird$Month,2),ebird$Month)
ebird = subset(ebird, ebird$Winter != "2005/2006" & ebird$Winter != "2016/2017")

# Reorder months
ebird$Month = factor(ebird$Month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))
ebird$Week = factor(ebird$Week, levels=c(35:53,0:17))
#Set X as na
ebird$OBSERVATION.COUNT = ifelse(ebird$OBSERVATION.COUNT == "X", 1, ebird$OBSERVATION.COUNT)
ebird$OBSERVATION.COUNT = as.numeric(ebird$OBSERVATION.COUNT)

aggCount = aggregate(ebird$OBSERVATION.COUNT, list(Week=ebird$Week, BCR=ebird$BCR.CODE, Winter=ebird$Winter), sum)
# Graph one winter
plot = ggplot(data=aggCount, aes(x=Week, y=x, group=Winter, color=Winter)) +
  stat_smooth(se=FALSE, na.rm = TRUE) + 
  ggtitle(paste("Figure 1. Observation count sum by BCR plotted over wintering period for ", species, sep="")) +
  facet_wrap(~ BCR)
# + scale_y_continuous(labels = comma) # optional code to remove scientific notation
print(plot)
ggsave(file=paste(species, "_Pop.jpg",sep=""), plot=plot, width=7, height=4)


# Graph mean
aggMean = aggregate(ebird$OBSERVATION.COUNT, list(Week=ebird$Week, BCR=ebird$BCR.CODE), mean)
plot = ggplot(aggMean, aes(x=Week, y=x)) + geom_point() +
  ylim(5,500) +
  ggtitle(paste("Figure 2. Observation count mean by BCR plotted over wintering period for ", species, sep="")) +
  stat_summary(aes(y = x,group=1), fun.y=mean, colour="red", geom="line",group=1) +
  facet_wrap(~ BCR)
print(plot)
ggsave(file=paste(species, "_Mean.jpg",sep=""), plot=plot, width=8, height=4)

print("P values < 0.05 indicate significant bimodality and values greater than 0.05 but less than 0.10 suggest bimodality with marginal significance")
for(i in unique(aggMean$BCR)) {
  testsetup = aggregate(ebird$OBSERVATION.COUNT, list(Week=ebird$Week, BCR=ebird$BCR.CODE, BCRName = ebird$BCRNAME), mean)
  testsetup = testsetup[which(testsetup$BCR == i),]
  test = dip.test(testsetup$x)
  bcr_name = unique(testsetup$BCRName)
  print(paste("P-value:", test$statistic[[1]]," / BCR:", bcr_name, sep=" "))
}

