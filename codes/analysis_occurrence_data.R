#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Mike Burrows
# Contact: michael.burrows@sams.ac.uk
#
########## Read and analyse eDNA data from UNESCO World Heritage Sites #########

# Load packages
library(dplyr)
library(stringr)
library(purrr)
library(arrow)
library(data.table)
fs::dir_create("figures")


####

speciestherm<-read_parquet("results/species_tsummaries.parquet")

speciesthermdt<-data.table(speciestherm)

speciesthermsst<-speciesthermdt[depth=="depthsurf" & variant=="mean",,]
speciesthermsstc<-dcast(speciesthermdt,species~metric,value.var = "value",fun=mean)

## Site by species occurrence table

speciesthermsite<-read_parquet("results/tsummaries_aggregated.parquet")

speciesthermsitedt<-data.table(speciesthermsite)

### CTI by higher geography areas #####

names(speciesthermsitedt)
# [1] "species"            "AphiaID"            "source_obis"        "source_gbif"        "source_dna"        
# [6] "higherGeography"    "q_0.5"              "q_0.75"             "q_0.9"              "q_0.95"            
# [11] "q_0.99"             "q_1"                "mean"               "sd"                 "q_0"               
# [16] "q_0.01"             "q_0.05"             "q_0.1"              "q_0.25"             "site_current"      
# [21] "site_ssp126_dec100" "site_ssp126_dec50"  "site_ssp245_dec100" "site_ssp245_dec50"  "site_ssp370_dec100"
# [26] "site_ssp370_dec50"  "site_ssp460_dec100" "site_ssp460_dec50"  "site_ssp585_dec100" "site_ssp585_dec50" 
# [31] "records"            "where"

ctihighergeog<-speciesthermsitedt[,list(ctiavg=mean(q_0.5),
                                        sdcti=sd(q_0.5),
                                        str=mean(q_0.9-q_0.1),
                                        sstavg=mean(site_current),
                                        nspp=.N),by=c("higherGeography","where")]

# Save for later use
write_parquet(ctihighergeog, "results/cti_highergeo.parquet")

# Save figures
pdf("figures/AllSpeciesThermalMetrics.pdf",height=3.5,width=8)
par(mfrow=c(1,3))
plot(data=ctihighergeog,ctiavg~sstavg,pch=16,col=unclass(as.factor(where)),cex=1.5,main="All species",
     xlab="Sea surface temperature (°C)", ylab="Community Temperature Index (°C)")
abline(a=0,b=1,lty=3)
legend("topleft",legend=levels(as.factor(ctihighergeog$where)),inset=0.02,col=1:3,pt.cex=1.5,pch=16)

plot(data=ctihighergeog,sdcti~sstavg,pch=16,col=unclass(as.factor(where)),cex=1.5,main="All species",
     xlab="Sea surface temperature (°C)", ylab="SD CTI (°C)")
#abline(a=0,b=1,lty=3)
legend("bottomleft",legend=levels(as.factor(ctihighergeog$where)),inset=0.02,col=1:3,pt.cex=1.5,pch=16)

plot(data=ctihighergeog,str~sstavg,pch=16,col=unclass(as.factor(where)),cex=1.5,main="All species",
     xlab="Sea surface temperature (°C)", ylab="Species thermal range (average, °C)")
#abline(a=0,b=1,lty=3)
legend("bottomleft",legend=levels(as.factor(ctihighergeog$where)),inset=0.02,col=1:3,pt.cex=1.5,pch=16)

### Add taxon information

specieslist<-read_parquet("results/species_list.parquet")

speciesthermsitedt1<-merge(speciesthermsitedt,specieslist,by="species")

fishthermsitedt<-speciesthermsitedt1[group=="fish",,]

### CTI by fish only #####


fishctihighergeog<-fishthermsitedt[,list(ctiavg=mean(q_0.5),
                                         sdcti=sd(q_0.5),
                                         str=mean(q_0.9-q_0.1),
                                         sstavg=mean(site_current),
                                         pgt90now=sum(site_current>q_0.9)/.N,
                                         pgt90ssp2=sum(site_ssp245_dec50>q_0.9)/.N,
                                         nspp=.N),by=c("higherGeography","where")]

par(mfrow=c(1,3))
plot(data=fishctihighergeog,ctiavg~sstavg,pch=16,col=unclass(as.factor(where)),cex=1.5,
     xlab="Sea surface temperature (°C)", ylab="Community Temperature Index (°C)",main="Fish")
abline(a=0,b=1,lty=3)
legend("topleft",legend=levels(as.factor(fishctihighergeog$where)),inset=0.02,col=1:3,pt.cex=1.5,pch=16)

plot(data=fishctihighergeog,sdcti~sstavg,pch=16,col=unclass(as.factor(where)),cex=1.5,
     xlab="Sea surface temperature (°C)", ylab="SD CTI (°C)",main="Fish")
#abline(a=0,b=1,lty=3)
legend("bottomleft",legend=levels(as.factor(fishctihighergeog$where)),inset=0.02,col=1:3,pt.cex=1.5,pch=16)

plot(data=fishctihighergeog,str~sstavg,pch=16,col=unclass(as.factor(where)),cex=1.5,
     xlab="Sea surface temperature (°C)", ylab="Species thermal range (average, °C)",main="Fish")
#abline(a=0,b=1,lty=3)
legend("bottomleft",legend=levels(as.factor(fishctihighergeog$where)),inset=0.02,col=1:3,pt.cex=1.5,pch=16)
dev.off()

### Proportion of species living above their T90 temperatures

par(mfrow=c(1,2))
plot(data=fishctihighergeog,pgt90now~sstavg,pch=16,col=unclass(as.factor(where)),cex=1.5,ylim=c(0,1),
     xlab="Sea surface temperature (°C)", ylab="Proportion of species with SST>T90",main="Fish")
abline(a=0,b=1,lty=3)
legend("topleft",legend=levels(as.factor(fishctihighergeog$where)),inset=0.02,col=1:3,pt.cex=1.5,pch=16)

plot(data=fishctihighergeog,pgt90ssp2~sstavg,pch=16,col=unclass(as.factor(where)),cex=1.5,ylim=c(0,1),
     xlab="Sea surface temperature (°C)", ylab="Proportion of species with SST>T90 under SSP2",main="Fish")
abline(a=0,b=1,lty=3)
legend("topleft",legend=levels(as.factor(fishctihighergeog$where)),inset=0.02,col=1:3,pt.cex=1.5,pch=16)




