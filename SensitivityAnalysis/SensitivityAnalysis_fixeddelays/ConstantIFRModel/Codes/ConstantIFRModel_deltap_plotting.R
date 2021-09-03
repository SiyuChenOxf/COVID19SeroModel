#This code is to plot comparison of exposure & seroprevalence 
#for constant IFR model using differerent values of delta_epsilon.

set.seed(100)

library(rgeos)
library(rgdal)
library(raster)
library(maptools)
library(sf)
library(dplyr)
library(tidyverse)
library(zoo)
library(ggplot2)
library(ggpubr)
library(egg)
library(RColorBrewer)
library(gridExtra)

## common values
ymax_kft <- 0.012
ymax_exposure <- 25
colors_Dark<-brewer.pal(7,"Dark2")
colors_Spectral<-brewer.pal(7,"Spectral")
font_size = 14
font_size_title = 16
lwd = 0.5
pt_size = 0.4
right_margin=1

## Load regional death and seroprevalence data ##
folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

SeroModelConstantIFR <- readRDS(paste(folder,"/SeroModelConstantIFR.rds", sep = ""))#Load posterior estimations for constant IFR model

beta <- rstan::extract(SeroModelConstantIFR)$beta
gamma_NorthEastYorkshireHumber<- rstan::extract(SeroModelConstantIFR)$gamma_NorthEastYorkshireHumber
gamma_London <- rstan::extract(SeroModelConstantIFR)$gamma_London
gamma_NorthWest <- rstan::extract(SeroModelConstantIFR)$gamma_NorthWest
gamma_SouthEast <- rstan::extract(SeroModelConstantIFR)$gamma_SouthEast
gamma_SouthWest <- rstan::extract(SeroModelConstantIFR)$gamma_SouthWest
gamma_Midlands <- rstan::extract(SeroModelConstantIFR)$gamma_Midlands
gamma_EastofEngland <- rstan::extract(SeroModelConstantIFR)$gamma_EastofEngland

London_data<-read.csv(paste(folder,"/London_data.csv", sep = ""),header = TRUE)
NorthWest_data<-read.csv(paste(folder,"/NorthWest_data.csv", sep = ""),header = TRUE)
YorkshireHumber_data<-read.csv(paste(folder,"/YorkshireHumber_data.csv", sep = ""),header = TRUE)
SouthWest_data<-read.csv(paste(folder,"/SouthWest_data.csv", sep = ""),header = TRUE)
EastofEngland_data<-read.csv(paste(folder,"/EastofEngland_data.csv", sep = ""),header = TRUE)
EastMidlands_data<-read.csv(paste(folder,"/EastMidlands_data.csv", sep = ""),header = TRUE)
WestMidlands_data<-read.csv(paste(folder,"/WestMidlands_data.csv", sep = ""),header = TRUE)
NorthEast_data<-read.csv(paste(folder,"/NorthEast_data.csv", sep = ""),header = TRUE)
SouthEast_data<-read.csv(paste(folder,"/SouthEast_data.csv", sep = ""),header = TRUE)

## Load population data ##
pop<-read.csv(paste(folder,"/Pop_data.csv", sep = ""),header = TRUE)

delta_p<-14        #fixed time lag between test for virus and death(seroversion)

#Setting delta_epsilon to different values for sensitivity analysis
delta_epsilon_7<-7  #fixed time lag between exposure and death(seroconversion)
delta_epsilon<-21  #fixed time lag between exposure and death(seroconversion)
delta_epsilon_28<-28  #fixed time lag between exposure and death(seroconversion)
delta_epsilon_14<-14  #fixed time lag between exposure and death(seroconversion)

P0_London<-pop$Pop[which(pop$Region=="London")]      
P0_NorthEast<-pop$Pop[which(pop$Region=="NorthEast")] 
P0_YorkshireHumber<-pop$Pop[which(pop$Region=="YorkshireHumber")] 
P0_NorthEastYorkshireHumber<-P0_NorthEast+P0_YorkshireHumber       #total population in NorthEast and Yorkshire and the Humber
P0_NorthWest<-pop$Pop[which(pop$Region=="NorthWest")] 
P0_EastMidlands<-pop$Pop[which(pop$Region=="EastMidlands")] 
P0_WestMidlands<-pop$Pop[which(pop$Region=="WestMidlands")]  
P0_Midlands<-P0_WestMidlands+P0_EastMidlands                       #total population in East Midlands and West Midlands
P0_EastofEngland<-pop$Pop[which(pop$Region=="EastofEngland")]  
P0_SouthEast<-pop$Pop[which(pop$Region=="SouthEast")] 
P0_SouthWest<-pop$Pop[which(pop$Region=="SouthWest")]  

#London
London_death_end<-which(!is.na(London_data$daily_death))[length(which(!is.na(London_data$daily_death)))] 
daily_death_London<-London_data$daily_death[1:London_death_end]          #Number of daily death in London
sero_London<-round(P0_London*London_data$sero[!is.na(London_data$sero)]) #Number of seropositive individuals at all serological survey time points in London 
n_days_London<-length(daily_death_London)                                #Number of days for daily death in London      
cumul_death_London<-cumsum(daily_death_London)                           #Number of cumulative death in London
n_days2_London<-length(sero_London)                                      #Number of days for serological survey in London
t_London<-seq(1,n_days_London,by=1)                                       
t2_London<-t_London[!is.na(London_data$sero)]                  

#NorthEast & Yorkshire and the Humber 
NorthEast_death_end<-which(!is.na(NorthEast_data$daily_death_NorthEast))[length(which(!is.na(NorthEast_data$daily_death)))]  
YorkshireHumber_death_end<-which(!is.na(YorkshireHumber_data$daily_death_Yorkshire_Humber))[length(which(!is.na(YorkshireHumber_data$daily_death_Yorkshire_Humber)))] 
daily_death_NorthEastYorkshireHumber<-NorthEast_data$daily_death_NorthEast[1:NorthEast_death_end]+YorkshireHumber_data$daily_death_Yorkshire_Humber[1:YorkshireHumber_death_end]
sero_NorthEastYorkshireHumber<-round(P0_NorthEastYorkshireHumber*NorthEast_data$sero[!is.na(NorthEast_data$sero)])
n_days_NorthEastYorkshireHumber<-length(daily_death_NorthEastYorkshireHumber)
cumul_death_NorthEastYorkshireHumber<-cumsum(daily_death_NorthEastYorkshireHumber)
n_days2_NorthEastYorkshireHumber<-length(sero_NorthEastYorkshireHumber)
t_NorthEastYorkshireHumber<-seq(1,n_days_NorthEastYorkshireHumber,by=1)
t2_NorthEastYorkshireHumber<-t_NorthEastYorkshireHumber[!is.na(NorthEast_data$sero)]

#SouthWest
SouthWest_death_end<-which(!is.na(SouthWest_data$daily_death))[length(which(!is.na(SouthWest_data$daily_death)))]#312
daily_death_SouthWest<-SouthWest_data$daily_death[1:SouthWest_death_end]
sero_SouthWest<-round(P0_SouthWest*SouthWest_data$sero[!is.na(SouthWest_data$sero)])
n_days_SouthWest<-length(daily_death_SouthWest)
cumul_death_SouthWest<-cumsum(daily_death_SouthWest)
n_days2_SouthWest<-length(sero_SouthWest)
t_SouthWest<-seq(1,n_days_SouthWest,by=1)
t2_SouthWest<-t_SouthWest[!is.na(SouthWest_data$sero)]

#NorthWest
NorthWest_death_end<-which(!is.na(NorthWest_data$daily_death))[length(which(!is.na(NorthWest_data$daily_death)))]#312
daily_death_NorthWest<-NorthWest_data$daily_death[1:NorthWest_death_end]
sero_NorthWest<-round(P0_NorthWest*NorthWest_data$sero[!is.na(NorthWest_data$sero)])
n_days_NorthWest<-length(daily_death_NorthWest)
cumul_death_NorthWest<-cumsum(daily_death_NorthWest)
n_days2_NorthWest<-length(sero_NorthWest)
t_NorthWest<-seq(1,n_days_NorthWest,by=1)
t2_NorthWest<-t_NorthWest[!is.na(NorthWest_data$sero)]

#SouthEast
SouthEast_death_end<-which(!is.na(SouthEast_data$daily_death))[length(which(!is.na(SouthEast_data$daily_death)))]#312
daily_death_SouthEast<-SouthEast_data$daily_death[1:SouthEast_death_end]
sero_SouthEast<-round(P0_SouthEast*SouthEast_data$sero[!is.na(SouthEast_data$sero)])
n_days_SouthEast<-length(daily_death_SouthEast)
cumul_death_SouthEast<-cumsum(daily_death_SouthEast)
n_days2_SouthEast<-length(sero_SouthEast)
t_SouthEast<-seq(1,n_days_SouthEast,by=1)
t2_SouthEast<-t_SouthEast[!is.na(SouthEast_data$sero)]

#Midlands
EastMidlands_death_end<-which(!is.na(EastMidlands_data$daily_death_EastMidlands))[length(which(!is.na(EastMidlands_data$daily_death_EastMidlands)))]#312
WestMidlands_death_end<-which(!is.na(WestMidlands_data$daily_death_WestMidlands))[length(which(!is.na(WestMidlands_data$daily_death_WestMidlands)))]
daily_death_Midlands<-WestMidlands_data$daily_death_WestMidlands[1:WestMidlands_death_end]+EastMidlands_data$daily_death_EastMidlands[1:EastMidlands_death_end]
sero_Midlands<-round(P0_Midlands*EastMidlands_data$sero[!is.na(EastMidlands_data$sero)])
n_days_Midlands<-length(daily_death_Midlands)
cumul_death_Midlands<-cumsum(daily_death_Midlands)
n_days2_Midlands<-length(sero_Midlands)
t_Midlands<-seq(1,n_days_Midlands,by=1)
t2_Midlands<-t_Midlands[!is.na(EastMidlands_data$sero)]

#EastofEngland
EastofEngland_death_end<-which(!is.na(EastofEngland_data$daily_death))[length(which(!is.na(EastofEngland_data$daily_death)))]#312
daily_death_EastofEngland<-EastofEngland_data$daily_death[1:EastofEngland_death_end]
sero_EastofEngland<-round(P0_EastofEngland*EastofEngland_data$sero[!is.na(EastofEngland_data$sero)])
n_days_EastofEngland<-length(daily_death_EastofEngland)
cumul_death_EastofEngland<-cumsum(daily_death_EastofEngland)
n_days2_EastofEngland<-length(sero_EastofEngland)
t_EastofEngland<-seq(1,n_days_EastofEngland,by=1)
t2_EastofEngland<-t_EastofEngland[!is.na(EastofEngland_data$sero)]

##################
### LONDON
##################
sim<-length(beta)
x_London<-epsilon_London<-kft_London<-matrix(0,sim,n_days_London)

for (i in 1:sim) {
  x_London[i,]<-rnbinom(rep(1,n_days_London), size= 100, mu=cumsum(exp(beta[i]*t_London)*(1-gamma_London[i])/gamma_London[i]*daily_death_London)/(exp(beta[i]*t_London)))/(P0_London-cumul_death_London)
  epsilon_London[i,]<-cumsum((1-gamma_London[i])/gamma_London[i]*daily_death_London)/(P0_London-cumul_death_London)
  kft_London[i,]<-gamma_London[i]
}
epsilon_London_21<-epsilon_London[,(delta_epsilon+1):n_days_London]
epsilon_London_28<-epsilon_London[,(delta_epsilon_28+1):n_days_London]
epsilon_London_14<-epsilon_London[,(delta_epsilon_14+1):n_days_London]
epsilon_London_7<-epsilon_London[,(delta_epsilon_7+1):n_days_London]


data1London = data.frame(output = c(rep("Exposure (Model 2)", (n_days_London-delta_epsilon)), rep("Seroprevalence (Model 1, 2, 3)", n_days_London), rep("Exposure (Model 3)",(n_days_London-delta_epsilon_28)),rep("Exposure (Model 1)",(n_days_London-delta_epsilon_14))), 
                         t=c(as.Date(London_data$Date)[1:(n_days_London-delta_epsilon)],as.Date(London_data$Date)[1:n_days_London],as.Date(London_data$Date)[1:(n_days_London-delta_epsilon_28)],as.Date(London_data$Date)[1:(n_days_London-delta_epsilon_14)]), 
                         median = c(100*apply(epsilon_London_21, 2, function(x) quantile(x, probs = 0.5)), 
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.5)),
                                    100*apply(epsilon_London_28, 2, function(x) quantile(x, probs = 0.5)),
                                    100*apply(epsilon_London_14, 2, function(x) quantile(x, probs = 0.5))), 
                         lower1 = c(100*apply(epsilon_London_21, 2, function(x) quantile(x, probs = 0.5)), 
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.025)),
                                    100*apply(epsilon_London_28, 2, function(x) quantile(x, probs = 0.025)),
                                    100*apply(epsilon_London_14, 2, function(x) quantile(x, probs = 0.025))), 
                         upper1 = c(100*apply(epsilon_London_21, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(epsilon_London_28, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(epsilon_London_14, 2, function(x) quantile(x, probs = 0.975))),
                         lower2 = c(100*apply(epsilon_London_21, 2, function(x) quantile(x, probs = 0.25)), 
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.25)),
                                    100*apply(epsilon_London_28, 2, function(x) quantile(x, probs = 0.25)),
                                    100*apply(epsilon_London_14, 2, function(x) quantile(x, probs = 0.25))), 
                         upper2 = c(100*apply(epsilon_London_21, 2, function(x) quantile(x, probs = 0.75)),
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.75)),
                                    100*apply(epsilon_London_28, 2, function(x) quantile(x, probs = 0.75)),
                                    100*apply(epsilon_London_14, 2, function(x) quantile(x, probs = 0.75))))
data2London = data.frame( t=as.Date(London_data$Date)[1:n_days_London][t2_London], value= 100*London_data$sero[t2_London], upper= 100*London_data$sero_upper[t2_London], lower = 100*London_data$sero_lower[t2_London])

p1London<-ggplot(data1London, aes(x=t, y = median, group = output,colour=output)) +
  scale_fill_manual(values=c(colors_Dark[1], colors_Dark[4], colors_Dark[3],colors_Dark[2]))+
  geom_line(size = lwd) +  ggtitle("London")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2London, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1London <- p1London +
  # scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" Percentage (%) ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "None",
    legend.title = element_blank(),
    axis.title.y = element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.text.y = element_text(size=font_size),
    axis.text.x = element_text(size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
 # styled1London
##################

##############################################
### North East & Yorkshire and the Humber
##############################################
sim<-length(beta)
x_NorthEastYorkshireHumber<-epsilon_NorthEastYorkshireHumber<-kft_NorthEastYorkshireHumber<-matrix(0,sim,n_days_NorthEastYorkshireHumber)

for (i in 1:sim) {
  x_NorthEastYorkshireHumber[i,]<-rnbinom(rep(1,n_days_NorthEastYorkshireHumber), size= 100, mu=cumsum(exp(beta[i]*t_NorthEastYorkshireHumber)*(1-gamma_NorthEastYorkshireHumber[i])/gamma_NorthEastYorkshireHumber[i]*daily_death_NorthEastYorkshireHumber)/(exp(beta[i]*t_NorthEastYorkshireHumber)))/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber)
  epsilon_NorthEastYorkshireHumber[i,]<-cumsum((1-gamma_NorthEastYorkshireHumber[i])/gamma_NorthEastYorkshireHumber[i]*daily_death_NorthEastYorkshireHumber)/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber)
  kft_NorthEastYorkshireHumber[i,]<-gamma_NorthEastYorkshireHumber[i]
}
epsilon_NorthEastYorkshireHumber_21<-epsilon_NorthEastYorkshireHumber[,(delta_epsilon+1):n_days_NorthEastYorkshireHumber]
epsilon_NorthEastYorkshireHumber_28<-epsilon_NorthEastYorkshireHumber[,(delta_epsilon_28+1):n_days_NorthEastYorkshireHumber]
epsilon_NorthEastYorkshireHumber_14<-epsilon_NorthEastYorkshireHumber[,(delta_epsilon_14+1):n_days_NorthEastYorkshireHumber]


data1NorthEastYorkshireHumber = data.frame(output = c(rep("Exposure (Model 2)", (n_days_NorthEastYorkshireHumber-delta_epsilon)), rep("Seroprevalence (Model 1, 2, 3)", n_days_NorthEastYorkshireHumber), rep("Exposure (Model 3)",(n_days_NorthEastYorkshireHumber-delta_epsilon_28)),rep("Exposure (Model 1)",(n_days_NorthEastYorkshireHumber-delta_epsilon_14))), 
                                           t=c(as.Date(NorthEast_data$Date)[1:(n_days_NorthEastYorkshireHumber-delta_epsilon)],as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber],as.Date(NorthEast_data$Date)[1:(n_days_NorthEastYorkshireHumber-delta_epsilon_28)],as.Date(NorthEast_data$Date)[1:(n_days_NorthEastYorkshireHumber-delta_epsilon_14)]), 
                                           median = c(100*apply(epsilon_NorthEastYorkshireHumber_21, 2, function(x) quantile(x, probs = 0.5)), 
                                                      100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.5)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_28, 2, function(x) quantile(x, probs = 0.5)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_14, 2, function(x) quantile(x, probs = 0.5))), 
                                           lower1 = c(100*apply(epsilon_NorthEastYorkshireHumber_21, 2, function(x) quantile(x, probs = 0.5)), 
                                                      100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.025)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_28, 2, function(x) quantile(x, probs = 0.025)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_14, 2, function(x) quantile(x, probs = 0.025))), 
                                           upper1 = c(100*apply(epsilon_NorthEastYorkshireHumber_21, 2, function(x) quantile(x, probs = 0.975)),
                                                      100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.975)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_28, 2, function(x) quantile(x, probs = 0.975)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_14, 2, function(x) quantile(x, probs = 0.975))),
                                           lower2 = c(100*apply(epsilon_NorthEastYorkshireHumber_21, 2, function(x) quantile(x, probs = 0.25)), 
                                                      100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.25)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_28, 2, function(x) quantile(x, probs = 0.25)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_14, 2, function(x) quantile(x, probs = 0.25))), 
                                           upper2 = c(100*apply(epsilon_NorthEastYorkshireHumber_21, 2, function(x) quantile(x, probs = 0.75)),
                                                      100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.75)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_28, 2, function(x) quantile(x, probs = 0.75)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_14, 2, function(x) quantile(x, probs = 0.75))))
data2NorthEastYorkshireHumber = data.frame( t=as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber][t2_NorthEastYorkshireHumber], value= 100*NorthEast_data$sero[t2_NorthEastYorkshireHumber], upper= 100*NorthEast_data$sero_upper[t2_NorthEastYorkshireHumber], lower = NorthEast_data$sero_lower[t2_NorthEastYorkshireHumber])


p1NorthEastYorkshireHumber<-ggplot(data1NorthEastYorkshireHumber, aes(x=t, y = median, group = output, colour = output)) +
  scale_fill_manual(values=c(colors_Dark[1], colors_Dark[4], colors_Dark[3],colors_Dark[2]))+
  geom_line(size = lwd) +  ggtitle("North East")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2NorthEastYorkshireHumber, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1NorthEastYorkshireHumber <- p1NorthEastYorkshireHumber +
  # scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "None",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1NorthEastYorkshireHumber
##################

##################
### NorthWest
##################
sim<-length(beta)
x_NorthWest<-epsilon_NorthWest<-kft_NorthWest<-matrix(0,sim,n_days_NorthWest)

for (i in 1:sim) {
  x_NorthWest[i,]<-rnbinom(rep(1,n_days_NorthWest), size= 100, mu=cumsum(exp(beta[i]*t_NorthWest)*(1-gamma_NorthWest[i])/gamma_NorthWest[i]*daily_death_NorthWest)/(exp(beta[i]*t_NorthWest)))/(P0_NorthWest-cumul_death_NorthWest)
  epsilon_NorthWest[i,]<-cumsum((1-gamma_NorthWest[i])/gamma_NorthWest[i]*daily_death_NorthWest)/(P0_NorthWest-cumul_death_NorthWest)
  kft_NorthWest[i,]<-gamma_NorthWest[i]
}
epsilon_NorthWest_21<-epsilon_NorthWest[,(delta_epsilon+1):n_days_NorthWest]
epsilon_NorthWest_28<-epsilon_NorthWest[,(delta_epsilon_28+1):n_days_NorthWest]
epsilon_NorthWest_14<-epsilon_NorthWest[,(delta_epsilon_14+1):n_days_NorthWest]

data1NorthWest = data.frame(output = c(rep("Exposure (Model 2)", (n_days_NorthWest-delta_epsilon)), rep("Seroprevalence (Model 1, 2, 3)", n_days_NorthWest), rep("Exposure (Model 3)",(n_days_NorthWest-delta_epsilon_28)),rep("Exposure (Model 1)",(n_days_NorthWest-delta_epsilon_14))), 
                            t=c(as.Date(NorthWest_data$Date)[1:(n_days_NorthWest-delta_epsilon)],as.Date(NorthWest_data$Date)[1:n_days_NorthWest],as.Date(NorthWest_data$Date)[1:(n_days_NorthWest-delta_epsilon_28)],as.Date(NorthWest_data$Date)[1:(n_days_NorthWest-delta_epsilon_14)]), 
                            median = c(100*apply(epsilon_NorthWest_21, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_NorthWest_28, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_NorthWest_14, 2, function(x) quantile(x, probs = 0.5))), 
                            lower1 = c(100*apply(epsilon_NorthWest_21, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_NorthWest_28, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_NorthWest_14, 2, function(x) quantile(x, probs = 0.025))), 
                            upper1 = c(100*apply(epsilon_NorthWest_21, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_NorthWest_28, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_NorthWest_14, 2, function(x) quantile(x, probs = 0.975))),
                            lower2 = c(100*apply(epsilon_NorthWest_21, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_NorthWest_28, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_NorthWest_14, 2, function(x) quantile(x, probs = 0.25))), 
                            upper2 = c(100*apply(epsilon_NorthWest_21, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_NorthWest_28, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_NorthWest_14, 2, function(x) quantile(x, probs = 0.75))))
data2NorthWest = data.frame( t=as.Date(NorthWest_data$Date)[1:n_days_NorthWest][t2_NorthWest], value= 100*NorthWest_data$sero[t2_NorthWest], upper= 100*NorthWest_data$sero_upper[t2_NorthWest], lower = 100*NorthWest_data$sero_lower[t2_NorthWest])

p1NorthWest<-ggplot(data1NorthWest, aes(x=t, y = median, group = output, colour = output)) +
  scale_fill_manual(values=c(colors_Dark[1], colors_Dark[4], colors_Dark[3],colors_Dark[2]))+
  geom_line(size = lwd) +  ggtitle("North West")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2NorthWest, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1NorthWest <- p1NorthWest +
  # scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "None",
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1NorthWest
##################

##################
### SouthEast
##################
sim<-length(beta)
x_SouthEast<-epsilon_SouthEast<-kft_SouthEast<-matrix(0,sim,n_days_SouthEast)

for (i in 1:sim) {
  x_SouthEast[i,]<-rnbinom(rep(1,n_days_SouthEast), size= 100, mu=cumsum(exp(beta[i]*t_SouthEast)*(1-gamma_SouthEast[i])/gamma_SouthEast[i]*daily_death_SouthEast)/(exp(beta[i]*t_SouthEast)))/(P0_SouthEast-cumul_death_SouthEast)
  epsilon_SouthEast[i,]<-cumsum((1-gamma_SouthEast[i])/gamma_SouthEast[i]*daily_death_SouthEast)/(P0_SouthEast-cumul_death_SouthEast)
  kft_SouthEast[i,]<-gamma_SouthEast[i]
}
epsilon_SouthEast_21<-epsilon_SouthEast[,(delta_epsilon+1):n_days_SouthEast]
epsilon_SouthEast_28<-epsilon_SouthEast[,(delta_epsilon_28+1):n_days_SouthEast]
epsilon_SouthEast_14<-epsilon_SouthEast[,(delta_epsilon_14+1):n_days_SouthEast]

data1SouthEast = data.frame(output = c(rep("Exposure (Model 2)", (n_days_SouthEast-delta_epsilon)), rep("Seroprevalence (Model 1, 2, 3)", n_days_SouthEast), rep("Exposure (Model 3)",(n_days_SouthEast-delta_epsilon_28)),rep("Exposure (Model 1)",(n_days_SouthEast-delta_epsilon_14))), 
                            t=c(as.Date(SouthEast_data$Date)[1:(n_days_SouthEast-delta_epsilon)],as.Date(SouthEast_data$Date)[1:n_days_SouthEast],as.Date(SouthEast_data$Date)[1:(n_days_SouthEast-delta_epsilon_28)],as.Date(SouthEast_data$Date)[1:(n_days_SouthEast-delta_epsilon_14)]), 
                            median = c(100*apply(epsilon_SouthEast_21, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_SouthEast_28, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_SouthEast_14, 2, function(x) quantile(x, probs = 0.5))), 
                            lower1 = c(100*apply(epsilon_SouthEast_21, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_SouthEast_28, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_SouthEast_14, 2, function(x) quantile(x, probs = 0.025))), 
                            upper1 = c(100*apply(epsilon_SouthEast_21, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_SouthEast_28, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_SouthEast_14, 2, function(x) quantile(x, probs = 0.975))),
                            lower2 = c(100*apply(epsilon_SouthEast_21, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_SouthEast_28, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_SouthEast_14, 2, function(x) quantile(x, probs = 0.25))), 
                            upper2 = c(100*apply(epsilon_SouthEast_21, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_SouthEast_28, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_SouthEast_14, 2, function(x) quantile(x, probs = 0.75))))
data2SouthEast = data.frame( t=as.Date(SouthEast_data$Date)[1:n_days_SouthEast][t2_SouthEast], value= 100*SouthEast_data$sero[t2_SouthEast], upper= 100*SouthEast_data$sero_upper[t2_SouthEast], lower = 100*SouthEast_data$sero_lower[t2_SouthEast])

p1SouthEast<-ggplot(data1SouthEast, aes(x=t, y = median, group = output, colour = output)) +
  scale_fill_manual(values=c(colors_Dark[1], colors_Dark[4], colors_Dark[3],colors_Dark[2]))+
  geom_line(size = lwd) +  ggtitle("South East")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2SouthEast, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1SouthEast <- p1SouthEast +
  # scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(-0.05,-0.3),
    legend.text=element_text(size=font_size),
    legend.title = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1SouthEast
##################

##################
### SouthWest
##################
sim<-length(beta)
x_SouthWest<-epsilon_SouthWest<-kft_SouthWest<-matrix(0,sim,n_days_SouthWest)

for (i in 1:sim) {
  x_SouthWest[i,]<-rnbinom(rep(1,n_days_SouthWest), size= 100, mu=cumsum(exp(beta[i]*t_SouthWest)*(1-gamma_SouthWest[i])/gamma_SouthWest[i]*daily_death_SouthWest)/(exp(beta[i]*t_SouthWest)))/(P0_SouthWest-cumul_death_SouthWest)
  epsilon_SouthWest[i,]<-cumsum((1-gamma_SouthWest[i])/gamma_SouthWest[i]*daily_death_SouthWest)/(P0_SouthWest-cumul_death_SouthWest)
  kft_SouthWest[i,]<-gamma_SouthWest[i]
}
epsilon_SouthWest_21<-epsilon_SouthWest[,(delta_epsilon+1):n_days_SouthWest]
epsilon_SouthWest_28<-epsilon_SouthWest[,(delta_epsilon_28+1):n_days_SouthWest]
epsilon_SouthWest_14<-epsilon_SouthWest[,(delta_epsilon_14+1):n_days_SouthWest]


data1SouthWest = data.frame(output = c(rep("Exposure (Model 2)", (n_days_SouthWest-delta_epsilon)), rep("Seroprevalence (Model 1, 2, 3)", n_days_SouthWest), rep("Exposure (Model 3)",(n_days_SouthWest-delta_epsilon_28)),rep("Exposure (Model 1)",(n_days_SouthWest-delta_epsilon_14))), 
                            t=c(as.Date(SouthWest_data$Date)[1:(n_days_SouthWest-delta_epsilon)],as.Date(SouthWest_data$Date)[1:n_days_SouthWest],as.Date(SouthWest_data$Date)[1:(n_days_SouthWest-delta_epsilon_28)],as.Date(SouthWest_data$Date)[1:(n_days_SouthWest-delta_epsilon_14)]), 
                            median = c(100*apply(epsilon_SouthWest_21, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_SouthWest_28, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_SouthWest_14, 2, function(x) quantile(x, probs = 0.5))), 
                            lower1 = c(100*apply(epsilon_SouthWest_21, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_SouthWest_28, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_SouthWest_14, 2, function(x) quantile(x, probs = 0.025))), 
                            upper1 = c(100*apply(epsilon_SouthWest_21, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_SouthWest_28, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_SouthWest_14, 2, function(x) quantile(x, probs = 0.975))),
                            lower2 = c(100*apply(epsilon_SouthWest_21, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_SouthWest_28, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_SouthWest_14, 2, function(x) quantile(x, probs = 0.25))), 
                            upper2 = c(100*apply(epsilon_SouthWest_21, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_SouthWest_28, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_SouthWest_14, 2, function(x) quantile(x, probs = 0.75))))
data2SouthWest = data.frame( t=as.Date(SouthWest_data$Date)[1:n_days_SouthWest][t2_SouthWest], value= 100*SouthWest_data$sero[t2_SouthWest], upper= 100*SouthWest_data$sero_upper[t2_SouthWest], lower = 100*SouthWest_data$sero_lower[t2_SouthWest])

p1SouthWest<-ggplot(data1SouthWest, aes(x=t, y = median, group = output, colour = output)) +
  scale_fill_manual(values=c(colors_Dark[1], colors_Dark[4], colors_Dark[3],colors_Dark[2]))+
  geom_line(size = lwd) +  ggtitle("South West")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2SouthWest, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1SouthWest <- p1SouthWest +
  # scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" Percentage (%)") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "None",
    legend.title = element_blank(),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    axis.title.y = element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.text.y = element_text(size=font_size),
    axis.text.x = element_text(size=font_size),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1SouthWest
##################

##################
### Midlands
##################
sim<-length(beta)
x_Midlands<-epsilon_Midlands<-kft_Midlands<-matrix(0,sim,n_days_Midlands)

for (i in 1:sim) {
  x_Midlands[i,]<-rnbinom(rep(1,n_days_Midlands), size= 100, mu=cumsum(exp(beta[i]*t_Midlands)*(1-gamma_Midlands[i])/gamma_Midlands[i]*daily_death_Midlands)/(exp(beta[i]*t_Midlands)))/(P0_Midlands-cumul_death_Midlands)
  epsilon_Midlands[i,]<-cumsum((1-gamma_Midlands[i])/gamma_Midlands[i]*daily_death_Midlands)/(P0_Midlands-cumul_death_Midlands)
  kft_Midlands[i,]<-gamma_Midlands[i]
}
epsilon_Midlands_21<-epsilon_Midlands[,(delta_epsilon+1):n_days_Midlands]
epsilon_Midlands_28<-epsilon_Midlands[,(delta_epsilon_28+1):n_days_Midlands]
epsilon_Midlands_14<-epsilon_Midlands[,(delta_epsilon_14+1):n_days_Midlands]


data1Midlands = data.frame(output = c(rep("Exposure (Model 2)", (n_days_Midlands-delta_epsilon)), rep("Seroprevalence (Model 1, 2, 3)", n_days_Midlands), rep("Exposure (Model 3)",(n_days_Midlands-delta_epsilon_28)),rep("Exposure (Model 1)",(n_days_Midlands-delta_epsilon_14))), 
                           t=c(as.Date(EastMidlands_data$Date)[1:(n_days_Midlands-delta_epsilon)],as.Date(EastMidlands_data$Date)[1:n_days_Midlands],as.Date(EastMidlands_data$Date)[1:(n_days_Midlands-delta_epsilon_28)],as.Date(EastMidlands_data$Date)[1:(n_days_Midlands-delta_epsilon_14)]), 
                           median = c(100*apply(epsilon_Midlands_21, 2, function(x) quantile(x, probs = 0.5)), 
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.5)),
                                      100*apply(epsilon_Midlands_28, 2, function(x) quantile(x, probs = 0.5)),
                                      100*apply(epsilon_Midlands_14, 2, function(x) quantile(x, probs = 0.5))), 
                           lower1 = c(100*apply(epsilon_Midlands_21, 2, function(x) quantile(x, probs = 0.5)), 
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.025)),
                                      100*apply(epsilon_Midlands_28, 2, function(x) quantile(x, probs = 0.025)),
                                      100*apply(epsilon_Midlands_14, 2, function(x) quantile(x, probs = 0.025))), 
                           upper1 = c(100*apply(epsilon_Midlands_21, 2, function(x) quantile(x, probs = 0.975)),
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.975)),
                                      100*apply(epsilon_Midlands_28, 2, function(x) quantile(x, probs = 0.975)),
                                      100*apply(epsilon_Midlands_14, 2, function(x) quantile(x, probs = 0.975))),
                           lower2 = c(100*apply(epsilon_Midlands_21, 2, function(x) quantile(x, probs = 0.25)), 
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.25)),
                                      100*apply(epsilon_Midlands_28, 2, function(x) quantile(x, probs = 0.25)),
                                      100*apply(epsilon_Midlands_14, 2, function(x) quantile(x, probs = 0.25))), 
                           upper2 = c(100*apply(epsilon_Midlands_21, 2, function(x) quantile(x, probs = 0.75)),
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.75)),
                                      100*apply(epsilon_Midlands_28, 2, function(x) quantile(x, probs = 0.75)),
                                      100*apply(epsilon_Midlands_14, 2, function(x) quantile(x, probs = 0.75))))
data2Midlands = data.frame( t=as.Date(EastMidlands_data$Date)[1:n_days_Midlands][t2_Midlands], value= 100*EastMidlands_data$sero[t2_Midlands], upper= 100*EastMidlands_data$sero_upper[t2_Midlands], lower = 100*EastMidlands_data$sero_lower[t2_Midlands])

p1Midlands<-ggplot(data1Midlands, aes(x=t, y = median, group = output, colour = output)) +
  scale_fill_manual(values=c(colors_Dark[1], colors_Dark[4], colors_Dark[3],colors_Dark[2]))+
  geom_line(size = lwd) +  ggtitle("Midlands")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2Midlands, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1Midlands <- p1Midlands +
  # scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "None",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1Midlands
##################

##################
### EastofEngland
##################
sim<-length(beta)
x_EastofEngland<-epsilon_EastofEngland<-kft_EastofEngland<-matrix(0,sim,n_days_EastofEngland)

for (i in 1:sim) {
  x_EastofEngland[i,]<-rnbinom(rep(1,n_days_EastofEngland), size= 100, mu=cumsum(exp(beta[i]*t_EastofEngland)*(1-gamma_EastofEngland[i])/gamma_EastofEngland[i]*daily_death_EastofEngland)/(exp(beta[i]*t_EastofEngland)))/(P0_EastofEngland-cumul_death_EastofEngland)
  epsilon_EastofEngland[i,]<-cumsum((1-gamma_EastofEngland[i])/gamma_EastofEngland[i]*daily_death_EastofEngland)/(P0_EastofEngland-cumul_death_EastofEngland)
  kft_EastofEngland[i,]<-gamma_EastofEngland[i]
}
epsilon_EastofEngland_21<-epsilon_EastofEngland[,(delta_epsilon+1):n_days_EastofEngland]
epsilon_EastofEngland_28<-epsilon_EastofEngland[,(delta_epsilon_28+1):n_days_EastofEngland]
epsilon_EastofEngland_14<-epsilon_EastofEngland[,(delta_epsilon_14+1):n_days_EastofEngland]


data1EastofEngland = data.frame(output = c(rep("Exposure (Model 2)", (n_days_EastofEngland-delta_epsilon)), rep("Seroprevalence (Model 1, 2, 3)", n_days_EastofEngland), rep("Exposure (Model 3)",(n_days_EastofEngland-delta_epsilon_28)),rep("Exposure (Model 1)",(n_days_EastofEngland-delta_epsilon_14))), 
                                t=c(as.Date(EastofEngland_data$Date)[1:(n_days_EastofEngland-delta_epsilon)],as.Date(EastofEngland_data$Date)[1:n_days_EastofEngland],as.Date(EastofEngland_data$Date)[1:(n_days_EastofEngland-delta_epsilon_28)],as.Date(EastofEngland_data$Date)[1:(n_days_EastofEngland-delta_epsilon_14)]), 
                                median = c(100*apply(epsilon_EastofEngland_21, 2, function(x) quantile(x, probs = 0.5)), 
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.5)),
                                           100*apply(epsilon_EastofEngland_28, 2, function(x) quantile(x, probs = 0.5)),
                                           100*apply(epsilon_EastofEngland_14, 2, function(x) quantile(x, probs = 0.5))), 
                                lower1 = c(100*apply(epsilon_EastofEngland_21, 2, function(x) quantile(x, probs = 0.5)), 
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.025)),
                                           100*apply(epsilon_EastofEngland_28, 2, function(x) quantile(x, probs = 0.025)),
                                           100*apply(epsilon_EastofEngland_14, 2, function(x) quantile(x, probs = 0.025))), 
                                upper1 = c(100*apply(epsilon_EastofEngland_21, 2, function(x) quantile(x, probs = 0.975)),
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.975)),
                                           100*apply(epsilon_EastofEngland_28, 2, function(x) quantile(x, probs = 0.975)),
                                           100*apply(epsilon_EastofEngland_14, 2, function(x) quantile(x, probs = 0.975))),
                                lower2 = c(100*apply(epsilon_EastofEngland_21, 2, function(x) quantile(x, probs = 0.25)), 
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.25)),
                                           100*apply(epsilon_EastofEngland_28, 2, function(x) quantile(x, probs = 0.25)),
                                           100*apply(epsilon_EastofEngland_14, 2, function(x) quantile(x, probs = 0.25))), 
                                upper2 = c(100*apply(epsilon_EastofEngland_21, 2, function(x) quantile(x, probs = 0.75)),
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.75)),
                                           100*apply(epsilon_EastofEngland_28, 2, function(x) quantile(x, probs = 0.75)),
                                           100*apply(epsilon_EastofEngland_14, 2, function(x) quantile(x, probs = 0.75))))
data2EastofEngland = data.frame( t=as.Date(EastofEngland_data$Date)[1:n_days_EastofEngland][t2_EastofEngland], value= 100*EastofEngland_data$sero[t2_EastofEngland], upper= 100*EastofEngland_data$sero_upper[t2_EastofEngland], lower = 100*EastofEngland_data$sero_lower[t2_EastofEngland])

p1EastofEngland<-ggplot(data1EastofEngland, aes(x=t, y = median, group = output, colour = output)) +
  scale_fill_manual(values=c(colors_Dark[1], colors_Dark[4], colors_Dark[3],colors_Dark[2]))+
  geom_line(size = lwd) +  ggtitle("East")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2EastofEngland, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1EastofEngland <- p1EastofEngland +
  # scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab("Percentage (%) ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "None",
    legend.title = element_blank(),
    axis.title.y = element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.text.y = element_text(size=font_size),
    axis.text.x = element_text(size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1EastofEngland
##################
folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Results"
folder = paste(folder_strings, sep = "", collapse = "/")

tiff(file=paste(folder,"/ConstantIFR_fitting_deltap.tiff", sep = ""),
     width=30, height=16, units="cm", res=300)
grid.arrange(styled1London,styled1NorthEastYorkshireHumber,styled1NorthWest,styled1SouthWest,styled1SouthEast,styled1Midlands,styled1EastofEngland)
dev.off()

