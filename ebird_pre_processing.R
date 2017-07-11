library(scales)

############################################################################
# Setup environment and variables
############################################################################
# Variable designiation
# Workspace directory
workspace = "D:/ebird/"
# Input ArcGIS Model csv file
inbird = "ebd_US_lessca_relNov-2016.txt"
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

ebird$OBSERVATION.COUNT = ifelse(ebird$OBSERVATION.COUNT == "X", 1, ebird$OBSERVATION.COUNT)
ebird$OBSERVATION.COUNT = as.numeric(ebird$OBSERVATION.COUNT)

write.csv(ebird, file=paste(workspace, paste(inbird, "post.csv", sep="_"), sep=""), quote=FALSE, row.names=FALSE)


