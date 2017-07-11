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
workspace = "D:/ebird"
species = "MALL"
# Input ArcGIS Model csv file
inbird = "mall.txt"
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



###################################################################################
# Silvermans test
#
# http://www-bcf.usc.edu/~gourab/code-bmt/tables/table-2/
# silverman.test <-function(x,k,M=999,adjust=FALSE,digits=6)
###################################################################################
silverman.test <-function(x,k,M=999,adjust=FALSE,digits=6){
  # x: data
  # k: number of modes to be tested
  # M: number of bootstrap replications
  
  #check if seed is available (as done in boot package)
  #if so save it
  seedAvailable = exists(x=".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(seedAvailable)
    saved_seed = .Random.seed 
  else{
    rnorm(1)
    saved_seed = .Random.seed
  }
  
  #temp function for bootstrapping
  y.obs <- function(x,h,sig=sd(x)){
    mean(x) + (x-mean(x)+h*rnorm(length(x),0,1))/((1+h^2/sig^2)^(1/2))
    #(x+h*rnorm(length(x),0,1))/((1+h^2/sig^2)^(1/2))
  }
  
  #temp function for density calculation
  nor.kernel <- function(x,h){
    density(x,bw=h,kernel ="gaussian")$y
  }
  
  #start of the test
  h0 <- h.crit(x, k)
  n <- 0
  
  for (i in 1:M) {
    x.boot <- sort(y.obs(sample(x, replace=TRUE),h0))
    mod.temp <- nr.modes(nor.kernel(x.boot,h0))
    if (mod.temp > k){
      n <- n+1
    }
  }
  
  p <- n/M
  ptemp=p
  
  if(adjust==TRUE){
    if(k==1){
      #asymptotic levels of silvermantest by Hall/York
      x=c(0,0.005,0.010,0.020,0.030,0.040,0.050,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.25,0.30,0.35,0.40,0.50)
      y=c(0,0,0,0.002,0.004,0.006,0.010,0.012,0.016,0.021,0.025,0.032,0.038,0.043,0.050,0.057,0.062,0.07,0.079,0.088,0.094,0.102,0.149,0.202,0.252,0.308,0.423)
      sp = interpSpline(x,y)
      #adjusting the p-value
      if(p<0.005)
        p=0
      else{
        p = predict(sp,p)$y
        p = round(p,digits)
      }
      
    }
    else{
      print("The option to adjust the p-value is valid only for k=1")
    } 
  }
  
  #return(list(saved_seed=saved_seed,p_value=p))
  #test_obj = new("Silvermantest", data=x, p_value = p,saved_seed=saved_seed,k=k)
  return(p)
}

h.crit <-
  function(x,k,prec=6){
    
    #temp function
    nor.kernel <- function(x,h){
      density(x,bw=h,kernel ="gaussian")$y
    }
    
    
    digits=prec
    prec=10^(-prec)
    x <- sort(x)
    minh <- min(diff(x))		#minimal possible h
    maxh <- diff(range(x))/2	#maximal possible h
    a <- maxh
    b <- minh
    zaehler=0
    
    while (abs(b-a)>prec){
      m <- nr.modes(nor.kernel(x,a))
      
      b <- a
      if (m > k){
        minh <- a
        a <- (a + maxh)/2
      } 
      else {
        maxh <- a
        a <- (a - minh)/2
      }
    }
    
    a=round(a,digits)
    
    
    if(nr.modes( nor.kernel(x,a) ) <= k){
      #subtract until more than k modes
      while(nr.modes( nor.kernel(x,a) ) <= k){
        a = a - prec
      }
      a=a+prec
    }
    
    if(nr.modes( nor.kernel(x,a) ) > k){
      #add until nr. of moodes correct
      while(nr.modes( nor.kernel(x,a) ) > k){
        a = a + prec
      }
    }
    
    a
  }


nr.modes <-
  function(y){
    
    d1 <- diff(y)
    signs <- diff(d1/abs(d1))
    length(signs[signs==-2])
    
  }
###################################################################################
# End Silvermans test http://www-bcf.usc.edu/~gourab/code-bmt/tables/table-2/




# Graph mean
aggMean = aggregate(ebird$OBSERVATION.COUNT, list(Week=ebird$MonthDay, BCR=ebird$BCR.CODE, BCRName = ebird$BCRNAME), mean)
for(i in unique(aggMean$BCR)){
  temp = aggMean[which(aggMean$BCR == 26),]
  ss=smooth.spline(x=temp$Week, y=temp$x, spar=0.7, keep.data = TRUE)
  ss$x = temp$Week
  plot(x=ss$x, y=ss$y, type="l", 
       ylab="Average count",
       xlab="Date",
       cex.lab=1.5
  )
  title(i)
  lines(x=ss$x,y=ss$y, col="red")
  test = dip.test(ss$yin)
  print(i)
  print(paste("Diptest: ", test$p.value[[1]], sep=""))
  print(paste("Silverman: ", silverman.test(ss$y, 2, M=999, adjust=FALSE), sep=""))
}
#ggsave(file=paste(species, "_Mean.jpg",sep=""), plot=plot, width=8, height=4)

print("P values < 0.05 indicate significant bimodality and values greater than 0.05 but less than 0.10 suggest bimodality with marginal significance")
for(i in unique(aggMean$BCR)) {
  testsetup = aggregate(ebird$OBSERVATION.COUNT, list(Week=ebird$Week, BCR=ebird$BCR.CODE, BCRName = ebird$BCRNAME), mean)
  testsetup = testsetup[which(testsetup$BCR == i),]
  test = dip.test(testsetup$x)
  bcr_name = unique(testsetup$BCRName)
  print(paste("P-value:", test$p.value[[1]]," / BCR:", bcr_name, sep=" "))
}

