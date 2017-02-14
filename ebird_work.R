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
ebird$Year = format(ebird$Date, "%Y")
ebird$Month = format(ebird$Date, "%m")
# Remove years < 2005 and months 5,6,7,8
ebird = subset(ebird, ebird$Year >= "2005")
ebird = subset(ebird, ebird$Month >= "09" | ebird$Month <= "04")

ebird$Month = ifelse(ebird$Month == "09", "9", ebird$Month)
ebird$Winter = ifelse(strtoi(ebird$Month) <= 4, paste(strtoi(ebird$Year) - 1,ebird$Year, sep = "/"),paste(ebird$Year, strtoi(ebird$Year) + 1, sep = "/"))
ebird$Week <- as.Date(cut(ebird$Date, breaks = "week", start.on.monday = FALSE))
ebird$Month = ifelse(strtoi(ebird$Month) <= 4, substring(ebird$Month,2),ebird$Month)


# Reorder months
ebird$Month = factor(ebird$Month, levels=c("9", "10", "11", "12", "1", "2", "3", "4"))
#Set X as na
ebird$OBSERVATION.COUNT = ifelse(ebird$OBSERVATION.COUNT == "X", 1, ebird$OBSERVATION.COUNT)
ebird$OBSERVATION.COUNT = as.numeric(ebird$OBSERVATION.COUNT)

aggCount = aggregate(ebird$OBSERVATION.COUNT, list(Month=ebird$Month, Winter=ebird$Winter), sum)
# Graph one winter
ggplot(data=aggCount, aes(x=Month, y=x, group=Winter, color=Winter)) +
  geom_line()

# Graph by month
ggplot(data=ebird, aes(x=Month, y=OBSERVATION.COUNT, group=Winter, color=Winter)) +
  geom_line()
