#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Mike Burrows
# Contact: michael.burrows@sams.ac.uk
#
######### Analyze Marine Heat Waves on MWHS sites for the eDNA project #########


####### SST 2x2 aggregated values daily ######

setwd("W:/sa01mb/SSTdaily")
# .libPaths("C:/Users/sa01mb/Documents/R/R-4x-library")

setwd("D:/Documents/Presentations/2023/MarineHWMallorca2023/DataAnalysis")

# unload()
# 
# install.packages("rlang",lib="C:/Users/sa01mb/Documents/R/R-4x-library")
# 
# install.packages("heatwaveR",lib="C:/Users/sa01mb/Documents/R/R-4x-library")
# install.packages("terra",lib="C:/Users/sa01mb/Documents/R/R-4x-library")
# install.packages("sp",lib="C:/Users/sa01mb/Documents/R/R-4x-library")
# install.packages("rgdal",lib="C:/Users/sa01mb/Documents/R/R-4x-library")
# install.packages("ncdf4",lib="C:/Users/sa01mb/Documents/R/R-4x-library")
# install.packages("data.table",lib="C:/Users/sa01mb/Documents/R/R-4x-library")
# install.packages("plyr",lib="C:/Users/sa01mb/Documents/R/R-4x-library")
# install.packages("pbapply",lib="C:/Users/sa01mb/Documents/R/R-4x-library")
# install.packages("RColorBrewer",lib="C:/Users/sa01mb/Documents/R/R-4x-library")
# install.packages("doParallel",lib="C:/Users/sa01mb/Documents/R/R-4x-library")


library(rlang)

library(terra)

#library(raster)
library(sp)
library(rgdal)
library(ncdf4)
library(data.table)
library(plyr)
library(pbapply)
library(parallel)
library(doParallel)
library(mgcv)

library(RColorBrewer)

library(heatwaveR)

col1<-colorRamp(c("red", "blue"),space="rgb")

# # http://data.ceda.ac.uk/badc/cmip5/data/cmip5/output1/NOAA-GFDL/GFDL-ESM2M/historicalNat/mon/ocean/Omon/r1i1p1/latest/thetao/thetao_Omon_GFDL-ESM2M_historicalNat_r1i1p1_198601-199012.nc.html
# 
# # 
# 
# Winfo<-nc_open("C:/Users/sa01mb/Downloads/thetao_Omon_GFDL-ESM2M_rcp26_r1i1p1_200601-201012.nc.nc")
#Winfo<-nc_open("sst.day.mean.2000.v2.nc")

#Winfo<-nc_open("D:/Downloads/sst.day.anom.2023.nc")

# D:/Downloads/sst.day.anom.2023.nc

#print(Winfo)

#test<-ncvar_get(Winfo, "sst")#,start=c(1,1,1),count=c(1440,720,))  # X-Y-time


# Test functions

agg5x5<-function(x) {
  rasx<-terra::rast(x)
  aggx<-terra::aggregate(rasx,fact=20,FUN=mean,na.rm=T)
  return<-aggx
}

agg2x2<-function(x) {
  rasx<-terra::rast(x)
  aggx<-terra::aggregate(rasx,fact=8,FUN=mean,na.rm=T)
  return<-aggx
}

#daily5x5<-pbapply(test,MARGIN=3,FUN=agg5x5) # works

#image(testx) # works

## S4 method for signature 'SpatRaster'

#writeCDF(r_c, paste0("sstdaily5x5_",year,".nc"), varname="SSTdaily", longname="", unit="Centigrade")


####### Read all daily files by year ########## 

for (year in 1982:2022) {
  print(year)
  if (year<2019){
    filetxt<-paste0("sst.day.mean.",year,".v2.nc") 
  } else {
    filetxt<-paste0("sst.day.mean.",year,".nc") 
  }
  Winfo<-nc_open(filetxt)
  test<-ncvar_get(Winfo, "sst")#=,start=c(1,1,1),count=c(1440,720,1))  # X-Y-time
  
  ## Aggregated daily data in a 5x5 grid
  
  daily5x5<-pbapply(test,MARGIN=3,FUN=agg5x5)
  daily5x5rast<-rast(daily5x5)
  writeCDF(daily5x5rast, paste0("sstdaily5x5_",year,".nc"), varname="SSTdaily", longname="", unit="Centigrade",
           overwrite=T)
  
}

####### Same using 2 x 2 lat/long ########## WORKS!!! #####

library(parallel)

years1<-1982:2018; years2=2019:2022

sstfiles<-c(paste0("sst.day.mean.",years1,".v2.nc"),paste0("sst.day.mean.",years2,".nc"))
lsstfiles<-as.list(sstfiles)

agg2x2nc<-function(filetxt) {
  Winfo<-ncdf4::nc_open(filetxt)
  yearx<-substr(filetxt,14,17)
  x<-ncdf4::ncvar_get(Winfo, "sst")#=,start=c(1,1,1),count=c(1440,720,1))  # X-Y-time
  rasx<-terra::rast(x)
  aggx<-terra::aggregate(rasx,fact=8,FUN=mean,na.rm=T)
  terra::writeCDF(aggx, paste0("sstdaily2x2_",yearx,".nc"), varname="SSTdaily", longname="", unit="Centigrade", overwrite=T)
  return<-aggx
}

cl <- makeCluster (12)

result <- parLapply (cl, lsstfiles, fun = agg2x2nc)

stopCluster (cl)

## Read and plot 2 x 2 files

yearx=2000
sstfile1<-paste0("sstdaily2x2_",yearx,".nc")
rasx<-terra::rast(sstfile1)

plot(rasx[[1]])
matrx<-as.matrix

# The structure of the OISST data after aggregation 2x2

rasx[[1]][]

latvals<-rep(seq(-89,89,2),180)
longvals<-rep(seq(359,1,-2),each=90)
longvals1<-ifelse(longvals>180,longvals-360,longvals)
idvals<-1:length(latvals)

gridinfo2x2<-data.frame(latvals,longvals,longvals1,idvals)

findval<-idvals[latvals==55 & longvals==351] #433 North Atlantic

sstrast<-rast(cbind(longvals1,latvals,rasx[[1]][]),type="xyz")

plot(sstrast)


############ 
# 
# # 
# 
# 
# # quants<-c(0.05,0.1,0.9,0.95)
# # 
# # 
# # for (year in 1982:2018) {
# #   print(year)
# #   filetxt<-paste0("sst.day.mean.",year,".v2.nc")
# #   Winfo<-nc_open(filetxt)
# #   test<-ncvar_get(Winfo, "sst")#=,start=c(1,1,1),count=c(1440,720,1))  # X-Y-time
# #   yeart<-pbapply(test,1:2, FUN=quantile, probs = c(0.05,0.1,0.9,0.95),  na.rm = TRUE)
# #   for (j in 1:4) {
# #     yearras<-raster(yeart[j,,])
# #     writeRaster(yearras,paste0("sstquant",quants[j],year,".nc"),overwrite=T)
# #   }
# # }
# 
# # #
# # 
# # q95exmap<-nctoras1(tlon,tlat,raster("sstquant0.951982.nc"))
# # plot(q95exmap)
# # writeRaster(trendmap,"sstquant0.951982.tif",overwrite=T)
# 
# 
# 
# 
# 
# # Combine yearly 5x5 
# 
# dayvalfiles<-list.files(pattern="sstdaily5x5_")
# 
# dayvalmaps<-pblapply(dayvalfiles,rast)
# dayvalrast<-rast(dayvalmaps)
# 
# writeCDF(dayvalrast, "allsstdaily5x5_1982to2022.nc", varname="SSTdaily", longname="", unit="Centigrade",
#          overwrite=F)

######

#dayvalsmatrix<-as.matrix(dayvalrast,wide=T)
years1<-1982:2018; years2=2019:2022


dayvalsdates<-seq(as.Date("1982-01-01"), as.Date("2022-12-31"), by="days")

sst2x2files<-c(paste0("sstdaily2x2_",c(years1,years2),".nc"))

lsst2x2files<-as.list(sst2x2files)

setwd("E:/SA01MB/DailySST") # faster read?? Yes!

daily2x2rast<-terra::rast(sst2x2files)

#### Calculate heatwave metrics ####

### Get heatwaves

climperiod<-c("1983-01-01","2012-12-31")

gridcells<-as.data.frame(daily2x2rast[[1]], xy=TRUE)

datesx=dayvalsdates
x=unlist(testts)

gethws<-function(x,datesx,climperiod) {
  
  mhwdata<-data.frame(t=datesx,temp=unlist(x))
  nadata<-sum(is.na(mhwdata))
  
  mhwclimatology<-heatwaveR::ts2clm(data=mhwdata,climatologyPeriod=c("1983-01-01","2012-12-31"))
  mhw1<-heatwaveR::detect_event(data=mhwclimatology)
  
  return(mhw1$event)
}

# Make data into a single matrix for passing to the parallel code - takes a little while

daily2x2vals<-as.matrix(daily2x2rast,wide=F)

# Set up parallel cluster

cl <- makeCluster(20) # only on Samphire

# Make sure library calls are accessible

clusterEvalQ(cl, .libPaths(c("C:/programs/rlib","C:/Users/sa01mb/Documents/R/R-4x-library")))

result <- parApply(cl=cl,daily2x2vals, MARGIN=1, FUN = gethws, dayvalsdates, climperiod)  # Works

stopCluster (cl)

saveRDS(result, file = "heatwaves2x2.Rds")

#heatwaveR::lolli_plot(result[[100]], metric = "intensity_max") # needs climatology not events



###### Get SST means and trends per cell ##########

# Make data into a single matrix for passing to the parallel code - takes a little while

daily2x2vals<-as.matrix(daily2x2rast,wide=F)

testsst<-daily2x2vals[433,]
# 
# sstdata<-data.frame(t=dayvalsdates,temp=testsst)


getssmeantrendquants<-function(x,datesx,climperiod) {
  
  sstdata<-data.frame(t=datesx,temp=unlist(x))
  nadata<-sum(!is.na(sstdata$temp) & !is.nan(sstdata$temp)) 
  if (nadata>0) {
    climperiodpsx<-as.Date(climperiod)
    sstdata<-subset(sstdata,t>climperiodpsx[1] & t<climperiodpsx[2])
    
    sstmean<-mean(sstdata$temp,na.rm=T)
    sstquantile<-quantile(sstdata$temp,c(0,0.1,0.5,0.9,1),na.rm = T)
    lmsstrend<-try(lm(temp ~ t, data = sstdata))
    trendb<-summary(lmsstrend)$coef
    rsquare<-summary(lmsstrend)$r.squared
    return(list(sstmean,sstquantile,trendb,rsquare)) 
  } else {
    return(NA)
  }
}

getsstyearts<-function(x,datesx,climperiod) {
  
  sstdata<-data.table(t=datesx,temp=unlist(x))
  nadata<-sum(!is.na(sstdata$temp) & !is.nan(sstdata$temp)) 
  if (nadata>0) {
    climperiodpsx<-as.Date(climperiod)
    sstdata$sstyear<-year(sstdata$t)
    sstyearmean<-sstdata[,list(yearmean=mean(temp,na.rm=T)),by="sstyear"]
    sstyearmean1<-subset(sstyearmean,sstyear>=year(climperiodpsx[1]) & sstyear<=year(climperiodpsx[2]))
    sstyearmeanclim<-mean(sstyearmean1$yearmean)
    sstyearmean$yanom<-sstyearmean$yearmean-sstyearmeanclim
    # sstmean<-mean(sstdata$temp,na.rm=T)
    # sstquantile<-quantile(sstdata$temp,c(0,0.1,0.5,0.9,1),na.rm = T)
    # lmsstrend<-try(lm(temp ~ t, data = sstdata))
    # trendb<-summary(lmsstrend)$coef
    # rsquare<-summary(lmsstrend)$r.squared
    return(sstyearmean) 
  } else {
    return(NA)
  }
  
}


# testres<-getssmeantrendquants(testsst, dayvalsdates, climperiod)

# Set up parallel cluster

detectCores(logical = FALSE)
cl <- makeCluster(20) # >12 only on Samphire

# Make sure library calls are accessible

clusterEvalQ(cl, .libPaths(c("C:/programs/rlib","C:/Users/sa01mb/Documents/R/R-4x-library")))

resultclm <- parApply(cl=cl,daily2x2vals, MARGIN=1, FUN = getssmeantrendquants, dayvalsdates, climperiod)  # Works

stopCluster (cl)

saveRDS(resultclm, file = "climate2x2.Rds")

testfunsst<-getsstyearts(testsst, dayvalsdates, climperiod)

sstyearlytimeseries<-pbapply(daily2x2vals, MARGIN=1, FUN = getsstyearts, dayvalsdates, climperiod)

saveRDS(sstyearlytimeseries, file = "sstyears2x2.Rds")

# Read processed data

resultclm<-readRDS("climate2x2a.Rds")

sstyearlytimeseries<-readRDS("sstyears2x2.Rds")

# Convert lists to data frames


# Climate trends

climdata <- sapply(resultclm, function(x) length(x)>1) # a vector with true for sst cells and false otherwise

resultclm1<-resultclm[climdata]
resultclm2<-lapply(resultclm1,unlist)
resultclm3<-as.data.frame(do.call(rbind,resultclm2)) # do.call(cbind, my_list)

names(resultclm3)[c(1,7:15)]<-c("meansst","trendint","trendintse","trend","trendse",
                                "trendintt","trendt","trendintp","trendp","trendrsq")

latvals<-rep(seq(-89,89,2),180)
longvals<-rep(seq(359,1,-2),each=90)
longvals1<-ifelse(longvals>180,longvals-360,longvals)
idvals<-1:length(latvals)

resultclm3$lat<-latvals[climdata]
resultclm3$lon<-longvals1[climdata]

# Yearly SST values: means and anomalies

sstyearlytimeseries1<-sstyearlytimeseries[climdata]
sstyearly2x2<-rbindlist(sstyearlytimeseries1,idcol=T)
sstyearly2x2$lon<-rep(longvals1[climdata],each=41)
sstyearly2x2$lat<-rep(latvals[climdata],each=41)

# Check outcomes

meansstmap<-with(resultclm3,rast(cbind(lon,lat,meansst),type="xyz"))
plot(meansstmap) # Correct


testmapdata<-subset(sstyearly2x2,sstyear==2000)
yearsstmap<-with(sstyearly2x2,rast(cbind(lon,lat,yearmean),type="xyz"))
plot(yearsstmap)  #   *** CORRECT*****


##### Process MHW 2x2 data #####

result<-readRDS("heatwaves2x2.Rds")

getmhwfreqyearts<-function(mhweventb,yearspan) {
  
  mhwevents<-data.table(mhweventb)
  mhwevents$monthpeak<-month(mhwevents$date_peak)
  
  mhwevents$year<-year(mhwevents$date_peak)
  mhwevents$demidec<-cut(mhwevents$year,breaks=seq(yearspan[1],yearspan[2],5))
  mhwevents$sumwint<-ifelse(mhwevents$month %in% 6:11,2,1)
  
  mhwyears<-mhwevents[,list(mhwfrq=.N,
                            avintmn=mean(intensity_mean),
                            avintmx=mean(intensity_max),
                            avdur=mean(duration),
                            sumdur=sum(duration),
                            sumcum=sum(intensity_cumulative)),by=c("year","sumwint")]
  
  vars<-names(mhwyears)[3:8]
  
  mhwyearsc<-dcast(mhwyears,year~sumwint,value.var =vars)
  
  mhwyearsonly<-data.frame(year=seq(yearspan[1],yearspan[2]))
  mhwyearsc1<-merge(mhwyearsonly,mhwyearsc,all.x=T)
  mhwyearsc1[is.na(mhwyearsc1)]<-0
  
  return(mhwyearsc1)
  
}

testdt<-getmhwfreqyearts(result[[433]],c(1983,2022))

allmhws<-pblapply(result,getmhwfreqyearts,c(1983,2022))

# saveRDS(allmhws, file = "allmhws2x2.Rds")

allmhws<-readRDS("allmhws2x2.Rds")

# test output

allmhws[[433]]

# Merge all together

allmhwsdt<-rbindlist(allmhws,fill=T,idcol=T)

latvals<-rep(seq(-89,89,2),180)
longvals<-rep(seq(359,1,-2),each=90)
idvals<-1:length(latvals)

mhwlatlon<-data.table(lon=longvals1,lat=latvals,id=idvals)

allmhwsdt$totalhwfq<-rowSums(allmhwsdt[,c("mhwfrq_1","mhwfrq_2")],na.rm=T)
allmhwsdt$totalcum<-rowSums(allmhwsdt[,c("sumcum_1","sumcum_2")],na.rm=T)
allmhwsdt$totaldur<-rowSums(allmhwsdt[,c("sumdur_1","sumdur_2")],na.rm=T)
allmhwsdt$totalmaxint<-mapply(max,allmhwsdt$avintmx_1,allmhwsdt$avintmx_2,na.rm=T)

allmhwsdt1<-merge(mhwlatlon,allmhwsdt,by.x="id",by.y=".id")

allmhwsdt1$intlat02<-2*floor(allmhwsdt1$lat/2)
allmhwsdt1$intlong02<-2*floor(allmhwsdt1$lon/2)




# trend summaries


# Get trends

library(plyr)

allmhwsdt[is.na(allmhwsdt)]<-0



# All HWs

hwallfqtrends <- dlply(allmhwsdt, .variables=c(".id"), function(df) 
  try(lm(totalhwfq ~ year, data = df)))
cellhwtrends<-ldply(hwallfqtrends, coef)

hwallmaxinttrends <- dlply(allmhwsdt, .variables=c(".id"), function(df) 
  try(lm(totalmaxint ~ year, data = df)))
cellhwmxinttrends<-ldply(hwallmaxinttrends, coef)


# Winter HWs

hwwinfqtrends <- dlply(allmhwsdt, .variables=c(".id"), function(df) 
  try(lm(mhwfrq_1 ~ year, data = df)))
celltrendswin<-ldply(hwwinfqtrends, coef)

hwwinmxtrends <- dlply(allmhwsdt, .variables=c(".id"), function(df) 
  try(lm(avintmx_1 ~ year, data = df)))
celltrendsmxwin<-ldply(hwwinmxtrends, coef)

# Summer HWs

hwsumfqtrends <- dlply(allmhwsdt, .variables=c(".id"), function(df) 
  try(lm(mhwfrq_2 ~ year, data = df)))
celltrendssumr<-ldply(hwsumfqtrends, coef)

hwsummxtrends <- dlply(allmhwsdt, .variables=c(".id"), function(df) 
  try(lm(avintmx_2 ~ year, data = df)))
celltrendsmxsum<-ldply(hwsummxtrends, coef)


##### Maps of MHW and climate data, including trends ###########

# MHW data #

mhwtrendmap<-rast(cbind(longvals1,latvals,cellhwtrends$year),type="xyz")

mhwmaxtrendmap<-rast(cbind(longvals1,latvals,cellhwmxinttrends$year),type="xyz")


mhwwintrendmap<-rast(cbind(longvals1,latvals,celltrendswin$year),type="xyz")
mhwsumtrendmap<-rast(cbind(longvals1,latvals,celltrendssumr$year),type="xyz")


trendbreaks=c(-Inf,seq(-0.2,0.2,0.02),Inf)
colbreaks= colorRampPalette(c("blue","white","red"),space="rgb",interpolate = c("linear"))(length(trendbreaks)-1)

plot(mhwtrendmap,col=colbreaks,breaks=trendbreaks,main="All HWs")

plot(mhwmaxtrendmap,col=colbreaks,breaks=trendbreaks,main="All HWs: Maximum intensity")

plot(mhwwintrendmap,col=colbreaks,breaks=trendbreaks,main="N Hemisphere Winter HWs (Dec-May)")
plot(mhwsumtrendmap,col=colbreaks,breaks=trendbreaks,main="N Hemisphere Summer HWs (Jun-Nov)")

writeRaster(mhwtrendmap,"D:/Documents/Presentations/2023/MarineHWMallorca2023/mhwtrendmap.tif",overwrite=T)

# SST data #

ssttrendmap<-with(resultclm3,rast(cbind(lon,lat,trend),type="xyz"))
ssttrendsemap<-with(resultclm3,rast(cbind(lon,lat,trendse),type="xyz"))

sstmeanmap<-with(resultclm3,rast(cbind(lon,lat,meansst),type="xyz"))
sstp90map<-with(resultclm3,rast(cbind(lon,lat,`90%`),type="xyz"))
sstp10map<-with(resultclm3,rast(cbind(lon,lat,`10%`),type="xyz"))

sstmeantrendvals<-data.frame(with(resultclm3,cbind(lon,lat,meansst,trend)))

plot(ssttrendmap,main=paste0("SST trend ",climperiod[1], " ",climperiod[2]))

plot(ssttrendsemap,main=paste0("SST trend SE ",climperiod[1], " ",climperiod[2]))
plot(sstp90map)
plot(sstp10map)
plot(sstmeanmap)

writeRaster(ssttrendmap,"D:/Documents/Presentations/2023/MarineHWMallorca2023/DataAnalysis/ssttrendmap.tif",overwrite=T)