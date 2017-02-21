library(ggplot2)
library(scales)
# Read in ebird data (16 species)
ebird <- read.delim("D:/ebird/ebd_US_agwtea1_relNov-2016.txt", sep="\t", header=TRUE, quote = "", stringsAsFactors = FALSE, na.strings=c(""))
species = "AGWT"
setwd("D:/ebird")

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
ebird$Week <- as.Date(cut(ebird$Date, breaks = "week", start.on.monday = FALSE))
#ebird$Month = ifelse(ebird$Month <= 5, substring(ebird$Month,2),ebird$Month)
ebird = subset(ebird, ebird$Winter != "2005/2006" & ebird$Winter != "2016/2017")

# Reorder months
ebird$Month = factor(ebird$Month, levels=c(9, 10, 11, 12, 1, 2, 3, 4))
#Set X as na
ebird$OBSERVATION.COUNT = ifelse(ebird$OBSERVATION.COUNT == "X", 1, ebird$OBSERVATION.COUNT)
ebird$OBSERVATION.COUNT = as.numeric(ebird$OBSERVATION.COUNT)

aggCount = aggregate(ebird$OBSERVATION.COUNT, list(Month=ebird$Month, ebird$BCR.CODE, Winter=ebird$Winter), sum)
# Graph one winter
plot = ggplot(data=aggCount, aes(x=Month, y=x, group=Winter, color=Winter)) +
  stat_smooth(se=FALSE, na.rm = TRUE) + 
  ggtitle(paste("Figure 1. Observation count sum mapped over wintering period for ", species, sep="")) +  
  facet_wrap(~ Group.2)
# + scale_y_continuous(labels = comma) # optional code to remove scientific notation
print(plot)
ggsave(file=paste(species, "_Pop.jpg",sep=""), plot=plot, width=7, height=4)