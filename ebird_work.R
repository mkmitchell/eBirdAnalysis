library(ggplot2)
library(scales)
# Read in ebird data (16 species)
ebird <- read.delim("D:/ebird/ebd_US_agwtea1_relNov-2016.txt", sep="\t", header=TRUE, quote = "", stringsAsFactors = FALSE, na.strings=c(""))

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
ebird$Year = as.Date(cut(ebird$Date, breaks = "year"))
ebird$Month = as.Date(cut(ebird$Date, breaks = "month"))
ebird$Week <- as.Date(cut(ebird$Date, breaks = "week", start.on.monday = FALSE))

#Set X as na
ebird$OBSERVATION.COUNT[is.na(ebird$OBSERVATION.COUNT)] <- 1

# Graph by month
ggplot(data = ebird,
       aes(Month, OBSERVATION.COUNT, group=Year)) +
  geom_line(size=1) +
  stat_summary(fun.y = sum, # adds up all observations for the month
               geom = "bar") + # or "line"
  scale_x_date(
    labels=date_format("%m"),
    date_breaks = "1 month") # custom x-axis labels