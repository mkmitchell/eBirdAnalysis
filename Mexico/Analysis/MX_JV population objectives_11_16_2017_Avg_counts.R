######################################################################################
###   Estimating midwinter population objectives using harvest and MWI             ###
###   Kathy Fleming, PHAB, last updated April 2017.                                ###
###   Edited by Mike Mitchell November 2017                                        ###
###                                                                                ###
###   This code calculates winter population objectives (for both LTA and          ###
###   80th percentile) for all JVs using 2 methods. A .csv format table is         ###
###   created for both the LTA and the 80th percentile over the LTA for each       ###
###   method.  Pdf files with Bar charts for each JV are created for methods       ###
###   4B and 4d.  Bar charts for each method can also be created in R for a        ###
###   single species by designating their 4-letter code below.                     ###
###                                                                                ###
###                                                                                ###
###   Four data files are required to run this code:  correcteddailyharv19992014.csv,# 
###   jvcounty.csv, MWandMO_objectives.csv, and JVabbrev.csv                       ###
###                                                                                ###
###   Adding Mexico population objectives to the end of this code requires         ###
###   additional files: 
######################################################################################
  
  
  library(reshape)
  
  ## specify path to access data tables and store results as well as the data input file names that reside in the path
  path <- "D:/ebird/Mexico/Analysis/2_19/"
  
  jvabbvinput <- "JVabbrev.csv"
  
  # Mexico data input
  mxharvestinput = "Mexico_counts.csv"
  mxmwinput = "Mexico_objectives.csv"
  mxJVcountyinput = "MX_jvcounty.csv"
  
  ## species to plot as bar graphs
  species <-  "MALL"
  
  
  
  ######################################### MEIXCO METHOD 4d #####################################################
  # harvest data:  1978-2014, winter period only but to end of Jan (Dec 1 - Jan 31)
  # MWI data: not used
  # Cont obj: LTA, 80 perc
  # winter survival rate:  
  # mxharvestinput = "Mexico_subset.csv"
  # mxmwinput = "Mexico_objectives.csv"
  #########################################################################################################
  
  ### Calculate proportion of harvest by county, using only harvest data for overwinter period:
  ### Dec 1 -- Jan 22
  
  data <- read.csv(paste(path, mxharvestinput, sep=""), na.strings = "")
  Hmonth <- substr(data$HDate, 1,2)
  Hmonth <- as.numeric(sub("/","",Hmonth))
  harvdata <- cbind(data, Hmonth)
  decharv <- subset(harvdata, Hmonth=="12")  # extract december harvest
  harvday <- ifelse (substr(decharv$HDate, 5,5)=="/",substr(decharv$HDate, 4,4),substr(decharv$HDate, 4,5))
  decharv <- cbind(decharv, harvday)
  decharv$harvday <- as.integer(decharv$harvday)
  sumharvest <- aggregate(decharv$Harvest, by=list(decharv$lumpedAOU, 
                                                   decharv$ST,
                                                   decharv$fips, 
                                                   decharv$Year), sum)  # sum extracted harvest by county
  names(sumharvest) <- c("spp","ST", "fips","Year", "harvest")
  
  ## Average count
  avgharvest <- aggregate(sumharvest$harvest, by=list(sumharvest$spp, sumharvest$fips), mean) # tot harv by spp in all states
  names(avgharvest)= c("spp", "fips", "harvest")
  totavg = aggregate(avgharvest$harvest, by=list(avgharvest$spp), sum)
  names(totavg) = c("spp", "totalharvest")
  # Sanity check.  aggregate sum of avgharvest$harvest by spp is the same as totavg
  aggregate(avgharvest$harvest, by=list(avgharvest$spp), sum)
  avgharv = merge(avgharvest, totavg, by="spp")
  avgharv$propharv <- avgharv$harvest/avgharv$totalharvest
  # Sanity check.  Aggregate sum of avgharvest$propharv by spp = 1
  aggregate(avgharv$propharv, by=list(avgharv$spp), sum)
  # End average count
  propharvest = avgharv
  ########## testing average count first
  propharvest[propharvest$propharv=="NaN",]$propharv <- rep(0,length(propharvest[propharvest$propharv=="NaN",]$propharv))
  ### Calculate proportion of LTA and 80 percentile in each county: 
  
  contobj <- read.csv(paste(path, mxmwinput, sep="")) # winter LTA and 80 perc LTA objectives by spp, corrected for Mx MWI
  contobj$spp = toupper(contobj$spp)
  alldata <- merge(propharvest,contobj, by="spp")
  LTApopobj <- alldata$propharv*alldata$LTA   # multiply each county prop of total harv by winter LTA obj
  perc80popobj <- alldata$propharv*alldata$X80p # multiply each county prop of total harv by 80 perc LTA obj
  popobj <- cbind(alldata, LTApopobj,perc80popobj)
  
  popobjtot <- aggregate(popobj[,8:9], by=list(popobj$spp),sum)  # should total to winter objectives by species and year
  #Sanity check
  aggregate(popobj[,8:9], by=list(popobj$spp),sum)  # should total to winter objectives by species and year
  #popobjstate <- aggregate(popobj[,8:9], by=list(popobj$spp), sum)  # total state winter objectives by species and year
  
  ###  Calculate population objectives by JV
  
  JVcounty <- read.csv(paste(path, mxJVcountyinput, sep=""))
  popobjJV <- merge(popobj, JVcounty, by=c("fips"))
  propLTApopobj <- popobjJV$propsqkm*popobjJV$LTApopobj
  propperc80popobj <- popobjJV$propsqkm*popobjJV$perc80popobj
  popobjJVtots <- cbind(popobjJV,propLTApopobj, propperc80popobj)
  #sanity check
  aggregate(popobjJVtots$propsqkm, by=list(popobjJVtots$spp, popobjJVtots$fips), sum) #should all add to 1
  aggregate(popobjJVtots$propLTApopobj, by=list(popobjJVtots$spp), sum)
  
  ## Create tables containing JV pop objectives
  
  JVtableLTA <- cast(popobjJVtots, JV~spp, value="propLTApopobj",sum)
  write.table(JVtableLTA, paste(path, "output/", "MX LTA popobj by JV Method 4d.csv", sep=""), sep=",", row.names=FALSE)
  JVtable80perc <- cast(popobjJVtots, JV~spp, value="propperc80popobj",sum)
  write.table(JVtable80perc, paste(path, "output/", "MX 80th perc popobj by JV Method 4d.csv", sep=""), sep=",", row.names=FALSE)
  
  
  ## Create table of county objectives for mapping in ArcGIS
  
  countypopLTA <- cast(popobjJVtots, fips~spp, value="propLTApopobj", sum)
  write.csv(countypopLTA, paste(path, "output/", "MX Meth4DcountypopLTA.csv", sep=""))
  countypop80 <- cast(popobjJVtots, fips~spp, value="propperc80popobj", sum)
  write.csv(countypop80, paste(path, "output/", "MX Meth4Dcountypop80.csv", sep=""))
  
  ## Create Bar Plots of LTA by JVs
  
  sppLTA <- subset(contobj, spp==species)
  LTA <- as.integer(sppLTA$MW_LTA)
  JVtot <- aggregate(popobjJVtots[,20:21], by=list(popobjJV$JV, popobjJV$spp), sum)
  names(JVtot) <- c("JV", "spp","propLTApopobj","propperc80popobj")
  JVabbrev <- read.csv(paste(path,jvabbvinput, sep=""))
  names(JVabbrev) <- c("JV", "JVabbrev")
  JVobj <- subset(JVtot, JVtot$spp==species)
  JVobjabbrev <- merge(JVobj, JVabbrev, by="JV")
  barplot(height=JVobjabbrev$propLTApopobj, col="gray", xlab="JV", ylab="Mexico Winter Population Objective", 
          names.arg=JVobjabbrev$JVabbrev, space=0, cex.names=0.8,
          main=unique(paste("Method 4d: ", JVobj$spp, ", using Mexico harvest data (Dec 1-Jan 31), no MWI, \nfall/wint surv rate 0.85, Cont Pop Obj =",
                            prettyNum(LTA,big.mark=",", scientific=F), sep=""),
                      ylim= c(0, max(JVobjabbrev$propLTApopobj))))
  
  # Plot bar chart for each JV with species totals -- LTA
  
  pdf(paste(path, "output/", "MX Method 4D JV bar charts LTA.pdf", sep=""), paper="a4r", width = 9)
  jvs <- unique(JVtot$JV)#[c(1:11,13:19)]
  
  for (i in 1:length(jvs))  {
    jvspptots <- JVtot[JVtot$JV==jvs[i],]
    jvmax <- signif(max(jvspptots$propLTApopobj),2)
    marks <- c(0, signif(jvmax/5), signif(jvmax*2/5),signif(jvmax*3/5),signif(jvmax*4/5),signif(jvmax)) 
    jvplot <-barplot(height=jvspptots$propLTApopobj, col="blue", xlab="Species", ylab="Mexico Winter Population Objective, LTA", 
                     names.arg=FALSE,main=unique(paste(jvs[i], " JV,\nPopulation Objectives by Species, LTA \nMethod 4D -- Dec 1-Jan 31 harvest data", sep="")),
                     axes=FALSE, ylim=c(0, (max(jvspptots$propLTApopobj)+(max(jvspptots$propLTApopobj)/5))))
    axis(2,at=marks,labels=format(marks,scientific=FALSE, big.mark=","), las=2, line=-1.5, hadj=0.9)
    jvplot
    text(x=jvplot[,1], y = 0, adj=c(1,2), jvspptots$spp, cex=0.8, srt=45, xpd=TRUE)
    
  }
  
  dev.off()
  
  # Plot bar chart for each JV with species totals -- 80th perc
  
  pdf(paste(path, "output/", "MX Method 4D JV bar charts 80th perc.pdf", sep=""), paper="a4r", width = 9)
  jvs <- unique(JVtot$JV)
  
  for (i in 1:length(jvs))  {
    jvspptots <- JVtot[JVtot$JV==jvs[i],]
    jvmax <- signif(max(jvspptots$propperc80popobj),2)
    marks <- c(0, signif(jvmax/5), signif(jvmax*2/5),signif(jvmax*3/5),signif(jvmax*4/5),signif(jvmax)) 
    jvplot <-barplot(height=jvspptots$propperc80popobj, col="blue", xlab="Species", ylab="Mexico Winter Population Objective, 80th Percentile", 
                     names.arg=FALSE,main=unique(paste(jvs[i], " JV,\nPopulation Objectives by Species, 80th percentile \nMethod 4D -- Dec 1-Jan 31 harvest data", sep="")),
                     axes=FALSE, ylim=c(0, (max(jvspptots$propperc80popobj)+(max(jvspptots$propperc80popobj)/5))))
    axis(2,at=marks, labels=format(marks,scientific=FALSE, big.mark=","), las=2, line=-1.5, hadj=0.9)
    jvplot
    text(x=jvplot[,1], y = 0, adj=c(1,2), jvspptots$spp, cex=0.8, srt=45, xpd=TRUE)
    
  }
  
  dev.off()