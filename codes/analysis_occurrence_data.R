#################### eDNA expeditions - scientific analysis ####################
########################## Environmental data download #########################
# January of 2024
# Author: Silas Principe, Mike Burrows
# Contact: s.principe@unesco.org; michael.burrows@sams.ac.uk
#
############## Plot eDNA data from UNESCO World Heritage Sites #################

# Load packages
library(dplyr)
library(stringr)
library(purrr)
library(arrow)
library(data.table)
library(ggplot2)
fs::dir_create("figures")


speciestherm<-read_parquet("results/species_tsummaries.parquet")

speciesthermdt<-data.table(speciestherm)

speciesthermsst<-speciesthermdt[depth=="depthsurf" & variant=="mean",,]
speciesthermsstc<-dcast(speciesthermdt,species~metric,value.var = "value",fun=mean)

## Site by species occurrence table

speciesthermsite<-read_parquet("results/tsummaries_aggregated.parquet")

speciesthermsitedt<-data.table(speciesthermsite)

speciesthermsitedt$where <- ifelse(speciesthermsitedt$where == "Both" | 
                                     speciesthermsitedt$where == "OBIS/GBIF", "Databases", "eDNA")

### CTI by higher geography areas #####
ctihighergeog<-speciesthermsitedt[,list(ctiavg=mean(q_0.5),
                                        sdcti=sd(q_0.5),
                                        str=mean(q_0.9-q_0.1),
                                        sstavg=mean(site_current),
                                        nspp=.N),by=c("higherGeography","where")]


# Save figures
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


(p1 <- ggplot(ctihighergeog) +
  geom_abline(slope = 1, intercept = 0, color = "grey90") +
  geom_point(aes(x = sstavg, y = ctiavg, colour = where, size = nspp), alpha = .5) +
  scale_color_manual(values = c("#00c3c9",  "#f58000")) + #"#5151db",
  scale_y_continuous(limits = c(5, 30)) +
  scale_x_continuous(limits = c(5, 30)) +
  scale_size(range = c(2,7)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  ylab("Community Thermal Index (°C)") + 
  xlab("Sea Surface Temperature (°C)") +
  theme_light() +
  theme(panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3, linetype = 2),
        legend.title = element_blank()))

(p2 <- ggplot(ctihighergeog) +
    geom_point(aes(x = sstavg, y = sdcti, colour = where, size = nspp), alpha = .5) +
    scale_color_manual(values = c("#00c3c9",  "#f58000")) + #"#5151db",
    # scale_y_continuous(limits = c(5, 30)) +
    # scale_x_continuous(limits = c(5, 30)) +
    scale_size(range = c(2,7)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    ylab("SD CTI (°C)") + 
    xlab("Sea Surface Temperature (°C)") +
    theme_light() +
    theme(panel.grid = element_blank(),
          panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3, linetype = 2),
          legend.title = element_blank()))

(p3 <- ggplot(ctihighergeog) +
    geom_point(aes(x = sstavg, y = str, colour = where, size = nspp), alpha = .5) +
    scale_color_manual(values = c("#00c3c9",  "#f58000")) + #"#5151db",
    # scale_y_continuous(limits = c(5, 30)) +
    # scale_x_continuous(limits = c(5, 30)) +
    scale_size(range = c(2,7)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    ylab("Species thermal range (average, °C)") + 
    xlab("Sea Surface Temperature (°C)") +
    theme_light() +
    theme(panel.grid = element_blank(),
          panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3, linetype = 2),
          legend.title = element_blank()))
  
library(patchwork)
p1 + p2 + p3 + patchwork::plot_layout(guides = "collect")

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




