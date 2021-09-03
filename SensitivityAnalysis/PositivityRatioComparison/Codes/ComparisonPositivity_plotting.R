#This code is to plot CFR, HFR and virus positivity ratios
set.seed(100)

library(rgeos)
library(rgdal)
library(raster)
library(maptools)
library(sf)
library(dplyr)
library(tidyverse)
library(zoo)a
library(ggplot2)
library(ggpubr)
library(egg)
library(RColorBrewer)
library(gridExtra)

## common values for plotting
ymax_kft <- 0.012
ymax_exposure <- 0.30
colors_Dark<-brewer.pal(7,"Dark2")
colors_Spectral<-brewer.pal(7,"Spectral")
font_size = 16
font_size_title = 16
lwd = 0.5
pt_size = 0.4
right_margin=1

## Load data ##
folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

## Load virus positivity data from PHE weekly COVID-19 surveillance report##
positivity_PHEC_Pillar1<-read.csv(paste(folder,"/positivity_PHEC_Pillar1_week39.csv", sep = ""),header = TRUE)
positivity_PHEC_Pillar2<-read.csv(paste(folder,"/positivity_PHEC_Pillar2_week39.csv", sep = ""),header = TRUE)
positivity_PHEC_Pillar1_2<-read.csv(paste(folder,"/positivity_PHEC_Pillar1_week44.csv", sep = ""),header = TRUE)
positivity_PHEC_Pillar2_2<-read.csv(paste(folder,"/positivity_PHEC_Pillar2_week44.csv", sep = ""),header = TRUE)

## Load virus positivity data from GOV.UK dashboard##
positivity<-read.csv(paste(folder,"/pcr_positivity.csv", sep = ""),header = TRUE)

## Load population data ##
pop<-read.csv(paste(folder,"/Pop_data.csv", sep = ""),header = TRUE)

P0_London<-pop$Pop[which(pop$Region=="London")]      
P0_NorthEast<-pop$Pop[which(pop$Region=="NorthEast")] 
P0_YorkshireHumber<-pop$Pop[which(pop$Region=="YorkshireHumber")] 
P0_NorthEastYorkshireHumber<-P0_NorthEast+P0_YorkshireHumber       #total population in NorthEast and Yorkshire and the Humber
P0_NorthWest<-pop$Pop[which(pop$Region=="NorthWest")] 
P0_EastMidlands<-pop$Pop[which(pop$Region=="EastMidlands")] 
P0_WestMidlands<-pop$Pop[which(pop$Region=="WestMidlands")]  
P0_Midlands<-P0_WestMidlands+P0_EastMidlands
P0_EastofEngland<-pop$Pop[which(pop$Region=="EastofEngland")]  
P0_SouthEast<-pop$Pop[which(pop$Region=="SouthEast")] 
P0_SouthWest<-pop$Pop[which(pop$Region=="SouthWest")]  

#virus positivity rates in London
positivity_London<-matrix(c(positivity_PHEC_Pillar1$London[1:14],0.5*positivity_PHEC_Pillar1$London[15:22]+0.5*positivity_PHEC_Pillar2$London[15:22],0.25*positivity_PHEC_Pillar1$London[23:35]+0.25*positivity_PHEC_Pillar2$London[23:35]+0.25*positivity_PHEC_Pillar1_2$London[1:13]+0.25*positivity_PHEC_Pillar2_2$London[1:13],0.5*positivity_PHEC_Pillar1_2$London[14:18]+0.5*positivity_PHEC_Pillar2_2$London[14:18]),nrow = 1)
positivity_London<-c(rep(0,26),as.vector(apply(positivity_London, 2,function(x) rep(x,7))))/100
positivity_LondonNEW<-positivity$London/100
data7<-data.frame(date=c(as.Date(positivity$date),as.Date(positivity$date)),positivity_London=c(positivity_London,positivity_LondonNEW),label=c(rep(c("Pillar 1 & Pillar 2 average","GOV.UK dashboard"),each=length(positivity_London))))
p7<-ggplot(data7,aes(x = date, y = positivity_London, colour = label))+ geom_line(lwd=1.5)+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6), limit = c(0, 0.6))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  theme_minimal() + ggtitle("London")+
  ylab(" Test positivity ratio") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

#virus positivity rates in NorthEast & Yorkshire and the Humber
positivity_NorthEast<-matrix(c(positivity_PHEC_Pillar1$North.East[1:14],0.5*positivity_PHEC_Pillar1$North.East[15:22]+0.5*positivity_PHEC_Pillar2$North.East[15:22],0.25*positivity_PHEC_Pillar1$North.East[23:35]+0.25*positivity_PHEC_Pillar2$North.East[23:35]+0.25*positivity_PHEC_Pillar1_2$North.East[1:13]+0.25*positivity_PHEC_Pillar2_2$North.East[1:13],0.5*positivity_PHEC_Pillar1_2$North.East[14:18]+0.5*positivity_PHEC_Pillar2_2$North.East[14:18]),nrow = 1)
positivity_NorthEast<-c(rep(0,26),as.vector(apply(positivity_NorthEast, 2,function(x) rep(x,7))))/100
positivity_YorkshireHumber<-matrix(c(positivity_PHEC_Pillar1$Yorkshire.and.Humber[1:14],0.5*positivity_PHEC_Pillar1$Yorkshire.and.Humber[15:22]+0.5*positivity_PHEC_Pillar2$Yorkshire.and.Humber[15:22],0.25*positivity_PHEC_Pillar1$Yorkshire.and.Humber[23:35]+0.25*positivity_PHEC_Pillar2$Yorkshire.and.Humber[23:35]+0.25*positivity_PHEC_Pillar1_2$Yorkshire.and.Humber[1:13]+0.25*positivity_PHEC_Pillar2_2$Yorkshire.and.Humber[1:13],0.5*positivity_PHEC_Pillar1_2$Yorkshire.and.Humber[14:18]+0.5*positivity_PHEC_Pillar2_2$Yorkshire.and.Humber[14:18]),nrow = 1)
positivity_YorkshireHumber<-c(rep(0,26),as.vector(apply(positivity_YorkshireHumber, 2,function(x) rep(x,7))))/100
positivity_NorthEastYorshireHumber<-positivity_NorthEast*P0_NorthEast/P0_NorthEastYorkshireHumber+positivity_YorkshireHumber*P0_YorkshireHumber/P0_NorthEastYorkshireHumber
positivity_NorthEastYorshireHumberNEW<-positivity$NorthEast*P0_NorthEast/P0_NorthEastYorkshireHumber+positivity$YorkshireandHumber*P0_YorkshireHumber/P0_NorthEastYorkshireHumber
positivity_NorthEastYorshireHumberNEW<-positivity_NorthEastYorshireHumberNEW/100
data6<-data.frame(date=c(as.Date(positivity$date),as.Date(positivity$date)),positivity_NorthEastYorshireHumber=c(positivity_NorthEastYorshireHumber,positivity_NorthEastYorshireHumberNEW),label=c(rep(c("Pillar 1 & Pillar 2 average","GOV.UK dashboard"),each=length(positivity_NorthEastYorshireHumber))))
p6<-ggplot(data6,aes(x = date, y = positivity_NorthEastYorshireHumber, colour = label))+ geom_line(lwd=1.5)+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6), limit = c(0, 0.6))+
  theme_minimal() + ggtitle("North East")+
  ylab(" Test positivity ratio") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )


#virus positivity rates in NorthWest
positivity_NorthWest<-matrix(c(positivity_PHEC_Pillar1$North.West[1:14],0.5*positivity_PHEC_Pillar1$North.West[15:22]+0.5*positivity_PHEC_Pillar2$North.West[15:22],0.25*positivity_PHEC_Pillar1$North.West[23:35]+0.25*positivity_PHEC_Pillar2$North.West[23:35]+0.25*positivity_PHEC_Pillar1_2$North.West[1:13]+0.25*positivity_PHEC_Pillar2_2$North.West[1:13],0.5*positivity_PHEC_Pillar1_2$North.West[14:18]+0.5*positivity_PHEC_Pillar2_2$North.West[14:18]),nrow = 1)
positivity_NorthWest<-c(rep(0,26),as.vector(apply(positivity_NorthWest, 2,function(x) rep(x,7))))/100
positivity_NorthWestNEW<-positivity$NorthWest/100
data5<-data.frame(date=c(as.Date(positivity$date),as.Date(positivity$date)),positivity_NorthWest=c(positivity_NorthWest,positivity_NorthWestNEW),label=c(rep(c("Pillar 1 & Pillar 2 average","GOV.UK dashboard"),each=length(positivity_NorthWest))))
p5<-ggplot(data5,aes(x = date, y = positivity_NorthWest, colour = label))+ geom_line(lwd=1.5)+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6), limit = c(0, 0.6))+
  theme_minimal() + ggtitle("North West")+
  ylab("") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(-0.01,-0.2),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

#virus positivity rates in SouthWest
positivity_SouthWest<-matrix(c(positivity_PHEC_Pillar1$South.West[1:14],0.5*positivity_PHEC_Pillar1$South.West[15:22]+0.5*positivity_PHEC_Pillar2$South.West[15:22],0.25*positivity_PHEC_Pillar1$South.West[23:35]+0.25*positivity_PHEC_Pillar2$South.West[23:35]+0.25*positivity_PHEC_Pillar1_2$South.West[1:13]+0.25*positivity_PHEC_Pillar2_2$South.West[1:13],0.5*positivity_PHEC_Pillar1_2$South.West[14:18]+0.5*positivity_PHEC_Pillar2_2$South.West[14:18]),nrow = 1)
positivity_SouthWest<-c(rep(0,26),as.vector(apply(positivity_SouthWest, 2,function(x) rep(x,7))))/100
positivity_SouthWestNEW<-positivity$SouthWest/100
data4<-data.frame(date=c(as.Date(positivity$date),as.Date(positivity$date)),positivity_SouthWest=c(positivity_SouthWest,positivity_SouthWestNEW),label=c(rep(c("Pillar 1 & Pillar 2 average","GOV.UK dashboard"),each=length(positivity_SouthWest))))
p4<-ggplot(data4,aes(x = date, y = positivity_SouthWest, colour = label))+ geom_line(lwd=1.5)+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6), limit = c(0, 0.6))+
  theme_minimal() + ggtitle("South West")+
  ylab(" Test positivity ratio") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    # axis.title.y = element_blank(),
    # axis.text.y = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    legend.title = element_blank(),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
#virus positivity rates in SouthEast
positivity_SouthEast<-matrix(c(positivity_PHEC_Pillar1$South.East[1:14],0.5*positivity_PHEC_Pillar1$South.East[15:22]+0.5*positivity_PHEC_Pillar2$South.East[15:22],0.25*positivity_PHEC_Pillar1$South.East[23:35]+0.25*positivity_PHEC_Pillar2$South.East[23:35]+0.25*positivity_PHEC_Pillar1_2$South.East[1:13]+0.25*positivity_PHEC_Pillar2_2$South.East[1:13],0.5*positivity_PHEC_Pillar1_2$South.East[14:18]+0.5*positivity_PHEC_Pillar2_2$South.East[14:18]),nrow = 1)
positivity_SouthEast<-c(rep(0,26),as.vector(apply(positivity_SouthEast, 2,function(x) rep(x,7))))/100
positivity_SouthEastNEW<-positivity$SouthEast/100
data3<-data.frame(date=c(as.Date(positivity$date),as.Date(positivity$date)),positivity_SouthEast=c(positivity_SouthEast,positivity_SouthEastNEW),label=c(rep(c("Pillar 1 & Pillar 2 average","GOV.UK dashboard"),each=length(positivity_SouthEast))))
p3<-ggplot(data3,aes(x = date, y = positivity_SouthEast, colour = label))+ geom_line(lwd=1.5)+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6), limit = c(0, 0.6))+
  theme_minimal() + ggtitle("South East")+
  ylab(" Test positivity ratio") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
#virus positivity rates in Midlands
positivity_WestMidlands<-matrix(c(positivity_PHEC_Pillar1$West.Midlands[1:14],0.5*positivity_PHEC_Pillar1$West.Midlands[15:22]+0.5*positivity_PHEC_Pillar2$West.Midlands[15:22],0.25*positivity_PHEC_Pillar1$West.Midlands[23:35]+0.25*positivity_PHEC_Pillar2$West.Midlands[23:35]+0.25*positivity_PHEC_Pillar1_2$West.Midlands[1:13]+0.25*positivity_PHEC_Pillar2_2$West.Midlands[1:13],0.5*positivity_PHEC_Pillar1_2$West.Midlands[14:18]+0.5*positivity_PHEC_Pillar2_2$West.Midlands[14:18]),nrow = 1)
positivity_EastMidlands<-matrix(c(positivity_PHEC_Pillar1$East.Midlands[1:14],0.5*positivity_PHEC_Pillar1$East.Midlands[15:22]+0.5*positivity_PHEC_Pillar2$East.Midlands[15:22],0.25*positivity_PHEC_Pillar1$East.Midlands[23:35]+0.25*positivity_PHEC_Pillar2$East.Midlands[23:35]+0.25*positivity_PHEC_Pillar1_2$East.Midlands[1:13]+0.25*positivity_PHEC_Pillar2_2$East.Midlands[1:13],0.5*positivity_PHEC_Pillar1_2$East.Midlands[14:18]+0.5*positivity_PHEC_Pillar2_2$East.Midlands[14:18]),nrow = 1)
positivity_WestMidlands<-c(rep(0,26),as.vector(apply(positivity_WestMidlands, 2,function(x) rep(x,7))))/100
positivity_EastMidlands<-c(rep(0,26),as.vector(apply(positivity_EastMidlands, 2,function(x) rep(x,7))))/100
positivity_Midlands<-positivity_EastMidlands*P0_EastMidlands/P0_Midlands+positivity_WestMidlands*P0_WestMidlands/P0_Midlands #using population of Westmidlands and Eastmidlands to adjust the positivity in Midlands
positivity_MidlandsNEW<-positivity$EastMidlands*P0_EastMidlands/P0_Midlands+positivity$WestMidlands*P0_WestMidlands/P0_Midlands
positivity_MidlandsNEW<-positivity_MidlandsNEW/100

data2<-data.frame(date=c(as.Date(positivity$date),as.Date(positivity$date)),positivity_Midlands=c(positivity_Midlands,positivity_MidlandsNEW),label=c(rep(c("Pillar 1 & Pillar 2 average","GOV.UK dashboard"),each=length(positivity_MidlandsNEW))))
p2<-ggplot(data2,aes(x = date, y = positivity_Midlands, colour = label))+ geom_line(lwd=1.5)+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6), limit = c(0, 0.6))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  theme_minimal() + ggtitle("Midlands")+
  ylab(" PCR test positivity ratio") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
#virus positivity rates in EastofEngland
positivity_EastofEngland<-matrix(c(positivity_PHEC_Pillar1$East.of.England[1:14],0.5*positivity_PHEC_Pillar1$East.of.England[15:22]+0.5*positivity_PHEC_Pillar2$East.of.England[15:22],0.25*positivity_PHEC_Pillar1$East.of.England[23:35]+0.25*positivity_PHEC_Pillar2$East.of.England[23:35]+0.25*positivity_PHEC_Pillar1_2$East.of.England[1:13]+0.25*positivity_PHEC_Pillar2_2$East.of.England[1:13],0.5*positivity_PHEC_Pillar1_2$East.of.England[14:18]+0.5*positivity_PHEC_Pillar2_2$East.of.England[14:18]),nrow = 1)
positivity_EastofEngland<-c(rep(0,26),as.vector(apply(positivity_EastofEngland, 2,function(x) rep(x,7))))/100
positivity_EastofEnglandNEW<-positivity$EastofEngland/100

data1<-data.frame(date=c(as.Date(positivity$date),as.Date(positivity$date)),positivity_EastofEngland=c(positivity_EastofEngland,positivity_EastofEnglandNEW),label=c(rep(c("Pillar 1 & Pillar 2 average","GOV.UK dashboard"),each=length(positivity_EastofEnglandNEW))))
p1<-ggplot(data1,aes(x = date, y = positivity_EastofEngland, colour = label))+ geom_line(lwd=1.5)+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6), limit = c(0, 0.6))+
  theme_minimal() + ggtitle("East")+
  ylab(" Test positivity ratio") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

#Saving plot
folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Results"
folder = paste(folder_strings, sep = "", collapse = "/")
tiff(file=paste(folder,"/PCR_PositivityComparision.tiff", sep = ""),
     width=36, height=18, units="cm", res=300)
grid.arrange(p1,p2,p3,p4,p5,p6,p7)
dev.off()

