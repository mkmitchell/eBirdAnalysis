library(ggplot2)
# Read in ebird data (16 species)
ebird <- read.delim("D:/ebird/ebird_data.txt", sep="\t", header=TRUE, quote = "", stringsAsFactors = FALSE, na.strings=c(""))

# Set date field and divide up by year, month, week for plotting
ebird$Date = as.Date(ebird$OBSERVATION.DATE, "%Y-%m-%d")
ebird$Year = as.Date(cut(ebird$Date, breaks = "year"))
ebird$Month = as.Date(cut(ebird$Date, breaks = "month"))
ebird$Week <- as.Date(cut(ebird$Date, breaks = "week", start.on.monday = FALSE))

# Split by species
mall = ebird[ebird$COMMON.NAME] == 'Mallard'

# Graph by month
ggplot(data = mall, aes(Date, OBSERVATION.COUNT)) + geom_line()