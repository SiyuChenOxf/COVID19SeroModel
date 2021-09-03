#This code is to plot comparison of exposure & seroprevalence 
#for time-varying IFR model using differerent combinations of delta_epsilon & delta_p.

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

#common values for plotting
ymax_kft <- 0.012
ymax_exposure <- 30
colors_Dark<-brewer.pal(7,"Dark2")
colors_Spectral<-brewer.pal(7,"Spectral")
font_size = 13
font_size_title = 15
lwd = 0.5
pt_size = 0.4
right_margin=1

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

SeroModelTimeVaryingIFR_delta_p_7 <- readRDS(paste(folder,"/SeroModelTimeVaryingIFR_delta_p_7.rds", sep = ""))   #Load posterior estimations for parameters using delta_p as 7 days
SeroModelTimeVaryingIFR_delta_p_14 <- readRDS(paste(folder,"/SeroModelTimeVaryingIFR_delta_p_14.rds", sep = "")) #Load posterior estimations for parameters using delta_p as 14 days
SeroModelTimeVaryingIFR_delta_p_21 <- readRDS(paste(folder,"/SeroModelTimeVaryingIFR_delta_p_21.rds", sep = "")) #Load posterior estimations for parameters using delta_p as 21 days

size = 15000

beta_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$beta,size)
eta_NorthEastYorkshireHumber_delta_p_7<- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$eta_NorthEastYorkshireHumber,size)
gamma_NorthEastYorkshireHumber_delta_p_7<- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$gamma_NorthEastYorkshireHumber,size)
eta_London_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$eta_London,size)
gamma_London_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$gamma_London,size)
eta_NorthWest_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$eta_NorthWest,size)
gamma_NorthWest_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$gamma_NorthWest,size)
eta_SouthEast_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$eta_SouthEast,size)
gamma_SouthEast_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$gamma_SouthEast,size)
eta_SouthWest_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$eta_SouthWest,size)
gamma_SouthWest_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$gamma_SouthWest,size)
eta_Midlands_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$eta_Midlands,size)
gamma_Midlands_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$gamma_Midlands,size)
eta_EastofEngland_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$eta_EastofEngland,size)
gamma_EastofEngland_delta_p_7 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_7)$gamma_EastofEngland,size)


beta_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$beta,size)
eta_NorthEastYorkshireHumber_delta_p_14<- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$eta_NorthEastYorkshireHumber,size)
gamma_NorthEastYorkshireHumber_delta_p_14<- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$gamma_NorthEastYorkshireHumber,size)
eta_London_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$eta_London,size)
gamma_London_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$gamma_London,size)
eta_NorthWest_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$eta_NorthWest,size)
gamma_NorthWest_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$gamma_NorthWest,size)
eta_SouthEast_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$eta_SouthEast,size)
gamma_SouthEast_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$gamma_SouthEast,size)
eta_SouthWest_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$eta_SouthWest,size)
gamma_SouthWest_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$gamma_SouthWest,size)
eta_Midlands_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$eta_Midlands,size)
gamma_Midlands_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$gamma_Midlands,size)
eta_EastofEngland_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$eta_EastofEngland,size)
gamma_EastofEngland_delta_p_14 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$gamma_EastofEngland,size)


beta_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$beta,size)
eta_NorthEastYorkshireHumber_delta_p_21<- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$eta_NorthEastYorkshireHumber,size)
gamma_NorthEastYorkshireHumber_delta_p_21<- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$gamma_NorthEastYorkshireHumber,size)
eta_London_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$eta_London,size)
gamma_London_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$gamma_London,size)
eta_NorthWest_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$eta_NorthWest,size)
gamma_NorthWest_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$gamma_NorthWest,size)
eta_SouthEast_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$eta_SouthEast,size)
gamma_SouthEast_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$gamma_SouthEast,size)
eta_SouthWest_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$eta_SouthWest,size)
gamma_SouthWest_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$gamma_SouthWest,size)
eta_Midlands_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_14)$eta_Midlands,size)
gamma_Midlands_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$gamma_Midlands,size)
eta_EastofEngland_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$eta_EastofEngland,size)
gamma_EastofEngland_delta_p_21 <- sample(rstan::extract(SeroModelTimeVaryingIFR_delta_p_21)$gamma_EastofEngland,size)


data<-data.frame(class=factor(rep(c("Model 4, 5, 6","Model 7, 8","Model 9"),each=length(beta_delta_p_21))),
                 para_beta=c(beta_delta_p_7,beta_delta_p_14,beta_delta_p_21),
                 para_London=c(gamma_London_delta_p_7,gamma_London_delta_p_14,gamma_London_delta_p_21),para_London2=c(eta_London_delta_p_7,eta_London_delta_p_14,eta_London_delta_p_21),
                 para_NorthEast=c(gamma_NorthEastYorkshireHumber_delta_p_7,gamma_NorthEastYorkshireHumber_delta_p_14,gamma_NorthEastYorkshireHumber_delta_p_21),para_NorthEast2=c(eta_NorthEastYorkshireHumber_delta_p_7,eta_NorthEastYorkshireHumber_delta_p_14,eta_NorthEastYorkshireHumber_delta_p_21),
                 para_SouthEast=c(gamma_SouthEast_delta_p_7,gamma_SouthEast_delta_p_14,gamma_SouthEast_delta_p_21),para_SouthEast2=c(eta_SouthEast_delta_p_7,eta_SouthEast_delta_p_14,eta_SouthEast_delta_p_21),
                 para_NorthWest=c(gamma_NorthWest_delta_p_7,gamma_NorthWest_delta_p_14,gamma_NorthWest_delta_p_21),para_NorthWest2=c(eta_NorthWest_delta_p_7,eta_NorthWest_delta_p_14,eta_NorthWest_delta_p_21),
                 para_SouthWest=c(gamma_SouthWest_delta_p_7,gamma_SouthWest_delta_p_14,gamma_SouthWest_delta_p_21),para_SouthWest2=c(eta_SouthWest_delta_p_7,eta_SouthWest_delta_p_14,eta_SouthWest_delta_p_21),
                 para_Midlands=c(gamma_Midlands_delta_p_7,gamma_Midlands_delta_p_14,gamma_Midlands_delta_p_21),para_Midlands2=c(eta_Midlands_delta_p_7,eta_Midlands_delta_p_14,eta_Midlands_delta_p_21),
                 para_EastofEngland=c(gamma_EastofEngland_delta_p_7,gamma_EastofEngland_delta_p_14,gamma_EastofEngland_delta_p_21),para_EastofEngland2=c(eta_EastofEngland_delta_p_7,eta_EastofEngland_delta_p_14,eta_EastofEngland_delta_p_21))
them<-theme(
  text = element_text(size=font_size),
  plot.title = element_text(face = "bold", size = font_size_title),
  legend.background = element_rect(fill = "white", size = 0.5, colour = "white"),
  legend.justification = c(0, 1),
  legend.position = "none",
  # legend.title = element_blank(),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank(),
  axis.ticks = element_line(colour = "grey50", size = 0.2),
  plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
  # panel.grid.major = element_line(colour = "grey50", size = 0.2),
  # panel.grid.minor = element_blank(),
)
p1<-ggplot(data, aes(x=para_beta, fill=class))+geom_density(alpha=.3)+
  xlab(" ")+
  ylab(expression(beta))+
  # scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(expression(beta)) +xlab("  ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "None",
    legend.title = element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    # panel.grid.major = element_line(colour = "grey50", size = 0.2),
    # panel.grid.minor = element_blank(),
    # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

p2<-ggplot(data, aes(x=para_London, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[London]))+theme_minimal()+them
p3<-ggplot(data, aes(x=para_NorthEast, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[NorthEast]))+theme_minimal()+them
p4<-ggplot(data, aes(x=para_NorthWest, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[NorthWest]))+theme_minimal()+them
p5<-ggplot(data, aes(x=para_SouthEast, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[SourthEast]))+theme_minimal()+them
p6<-ggplot(data, aes(x=para_SouthWest, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[SouthWest]))+theme_minimal()+them
p7<-ggplot(data, aes(x=para_Midlands, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[Midlands]))+theme_minimal()+them
p8<-ggplot(data, aes(x=para_EastofEngland, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[East]))+theme_minimal()+them

p9<-ggplot(data, aes(x=para_London2, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(eta[London]))+theme_minimal()+them
p10<-ggplot(data, aes(x=para_NorthEast2, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(eta[NorthEast]))+theme_minimal()+them
p11<-ggplot(data, aes(x=para_NorthWest2, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(eta[NorthWest]))+theme_minimal()+them
p12<-ggplot(data, aes(x=para_SouthEast2, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(eta[SouthEast]))+theme_minimal()+them
p13<-ggplot(data, aes(x=para_SouthWest2, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(eta[SouthWest]))+theme_minimal()+them
p14<-ggplot(data, aes(x=para_Midlands2, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(eta[Midlands]))+theme_minimal()+them
p15<-ggplot(data, aes(x=para_EastofEngland2, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(eta[EastofEngland]))+theme_minimal()+them+theme(legend.position = c(1,1),legend.title = element_blank())

London_data<-read.csv(paste(folder,"/London_data.csv", sep = ""),header = TRUE)
NorthWest_data<-read.csv(paste(folder,"/NorthWest_data.csv", sep = ""),header = TRUE)
YorkshireHumber_data<-read.csv(paste(folder,"/YorkshireHumber_data.csv", sep = ""),header = TRUE)
SouthWest_data<-read.csv(paste(folder,"/SouthWest_data.csv", sep = ""),header = TRUE)
EastofEngland_data<-read.csv(paste(folder,"/EastofEngland_data.csv", sep = ""),header = TRUE)
EastMidlands_data<-read.csv(paste(folder,"/EastMidlands_data.csv", sep = ""),header = TRUE)
WestMidlands_data<-read.csv(paste(folder,"/WestMidlands_data.csv", sep = ""),header = TRUE)
NorthEast_data<-read.csv(paste(folder,"/NorthEast_data.csv", sep = ""),header = TRUE)
SouthEast_data<-read.csv(paste(folder,"/SouthEast_data.csv", sep = ""),header = TRUE)

## Load virus positivity data from GOV.UK##
positivity<-read.csv(paste(folder,"/pcr_positivity.csv", sep = ""),header = TRUE)

## Load population data ##
pop<-read.csv(paste(folder,"/Pop_data.csv", sep = ""),header = TRUE)

delta_p_14<-14        #fixed time lag between test for virus and death(seroversion)
delta_p_21<-21        #fixed time lag between test for virus and death(seroversion)
delta_p_7<-7          #fixed time lag between test for virus and death(seroversion)

delta_epsilon_21<-21  #fixed time lag between exposure and death(seroconversion)
delta_epsilon_28<-28  #fixed time lag between exposure and death(seroconversion)
delta_epsilon_14<-14  #fixed time lag between exposure and death(seroconversion)

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
positivity_London<-positivity$London/100

#virus positivity rates in NorthEast & Yorkshire and the Humber
positivity_NorthEastYorshireHumber<-positivity$NorthEast*P0_NorthEast/P0_NorthEastYorkshireHumber+positivity$YorkshireandHumber*P0_YorkshireHumber/P0_NorthEastYorkshireHumber
positivity_NorthEastYorshireHumber<-positivity_NorthEastYorshireHumber/100

#virus positivity rates in NorthWest
positivity_NorthWest<-positivity$NorthWest/100

#virus positivity rates in SouthWest
positivity_SouthWest<-positivity$SouthWest/100

#virus positivity rates in SouthEast
positivity_SouthEast<-positivity$SouthEast/100

#virus positivity rates in Midlands
positivity_Midlands<-positivity$EastMidlands*P0_EastMidlands/P0_Midlands+positivity$WestMidlands*P0_WestMidlands/P0_Midlands
positivity_Midlands<-positivity_Midlands/100

#virus positivity rates in EastofEngland
positivity_EastofEngland<-positivity$EastofEngland/100

#Epidemic stage in London
positivity_start_London<-which(positivity_London!=0)[1] 
London_positivity_end<-which(positivity_London!=0)[length(which(positivity_London!=0))] 
London_death_end<-which(!is.na(London_data$daily_death))[length(which(!is.na(London_data$daily_death)))] 

Epidemic_Stage_London_delta_p_7<-c(rep(0,positivity_start_London+delta_p_7-1),positivity_London[positivity_start_London:London_positivity_end])[1:London_death_end]
Epidemic_Stage_London_delta_p_7<-cumsum(Epidemic_Stage_London_delta_p_7)/sum(Epidemic_Stage_London_delta_p_7)

Epidemic_Stage_London_delta_p_14<-c(rep(0,positivity_start_London+delta_p_14-1),positivity_London[positivity_start_London:London_positivity_end])[1:London_death_end]
Epidemic_Stage_London_delta_p_14<-cumsum(Epidemic_Stage_London_delta_p_14)/sum(Epidemic_Stage_London_delta_p_14)

Epidemic_Stage_London_delta_p_21<-c(rep(0,positivity_start_London+delta_p_21-1),positivity_London[positivity_start_London:London_positivity_end])[1:London_death_end]
Epidemic_Stage_London_delta_p_21<-cumsum(Epidemic_Stage_London_delta_p_21)/sum(Epidemic_Stage_London_delta_p_21)

#Daily death in London
daily_death_London<-London_data$daily_death[1:London_death_end]

#Epidemic stage in NorthEast & Yorkshire and the Humber 
positivity_start_NorthEastYorkshireHumber<-which(positivity_NorthEastYorshireHumber!=0)[1] 
NorthEastYorkshireHumber_positivity_end<-which(positivity_NorthEastYorshireHumber!=0)[length(which(positivity_NorthEastYorshireHumber!=0))] 
NorthEast_death_end<-which(!is.na(NorthEast_data$daily_death_NorthEast))[length(which(!is.na(NorthEast_data$daily_death)))]  
YorkshireHumber_death_end<-which(!is.na(YorkshireHumber_data$daily_death_Yorkshire_Humber))[length(which(!is.na(YorkshireHumber_data$daily_death_Yorkshire_Humber)))] 

Epidemic_Stage_NorthEastYorkshireHumber_delta_p_7<-c(rep(0,positivity_start_NorthEastYorkshireHumber+delta_p_7-1),positivity_NorthEastYorshireHumber[positivity_start_NorthEastYorkshireHumber:NorthEastYorkshireHumber_positivity_end])[1:NorthEast_death_end]
Epidemic_Stage_NorthEastYorkshireHumber_delta_p_7<-cumsum(Epidemic_Stage_NorthEastYorkshireHumber_delta_p_7)/sum(Epidemic_Stage_NorthEastYorkshireHumber_delta_p_7)

Epidemic_Stage_NorthEastYorkshireHumber_delta_p_14<-c(rep(0,positivity_start_NorthEastYorkshireHumber+delta_p_14-1),positivity_NorthEastYorshireHumber[positivity_start_NorthEastYorkshireHumber:NorthEastYorkshireHumber_positivity_end])[1:NorthEast_death_end]
Epidemic_Stage_NorthEastYorkshireHumber_delta_p_14<-cumsum(Epidemic_Stage_NorthEastYorkshireHumber_delta_p_14)/sum(Epidemic_Stage_NorthEastYorkshireHumber_delta_p_14)

Epidemic_Stage_NorthEastYorkshireHumber_delta_p_21<-c(rep(0,positivity_start_NorthEastYorkshireHumber+delta_p_21-1),positivity_NorthEastYorshireHumber[positivity_start_NorthEastYorkshireHumber:NorthEastYorkshireHumber_positivity_end])[1:NorthEast_death_end]
Epidemic_Stage_NorthEastYorkshireHumber_delta_p_21<-cumsum(Epidemic_Stage_NorthEastYorkshireHumber_delta_p_21)/sum(Epidemic_Stage_NorthEastYorkshireHumber_delta_p_21)

#Daily death in NorthEast & Yorkshire and the Humber
daily_death_NorthEastYorkshireHumber<-NorthEast_data$daily_death_NorthEast[1:NorthEast_death_end]+YorkshireHumber_data$daily_death_Yorkshire_Humber[1:YorkshireHumber_death_end]

#Epidemic stage in SouthEast
positivity_start_SouthEast<-which(positivity_SouthEast!=0)[1]#62
SouthEast_positivity_end<-which(positivity_SouthEast!=0)[length(which(positivity_SouthEast!=0))]#306
SouthEast_death_end<-which(!is.na(SouthEast_data$daily_death))[length(which(!is.na(SouthEast_data$daily_death)))]#312

Epidemic_Stage_SouthEast_delta_p_7<-c(rep(0,positivity_start_SouthEast+delta_p_7-1),positivity_SouthEast[positivity_start_SouthEast:SouthEast_positivity_end])[1:SouthEast_death_end]
Epidemic_Stage_SouthEast_delta_p_7<-cumsum(Epidemic_Stage_SouthEast_delta_p_7)/sum(Epidemic_Stage_SouthEast_delta_p_7)

Epidemic_Stage_SouthEast_delta_p_14<-c(rep(0,positivity_start_SouthEast+delta_p_14-1),positivity_SouthEast[positivity_start_SouthEast:SouthEast_positivity_end])[1:SouthEast_death_end]
Epidemic_Stage_SouthEast_delta_p_14<-cumsum(Epidemic_Stage_SouthEast_delta_p_14)/sum(Epidemic_Stage_SouthEast_delta_p_14)

Epidemic_Stage_SouthEast_delta_p_21<-c(rep(0,positivity_start_SouthEast+delta_p_21-1),positivity_SouthEast[positivity_start_SouthEast:SouthEast_positivity_end])[1:SouthEast_death_end]
Epidemic_Stage_SouthEast_delta_p_21<-cumsum(Epidemic_Stage_SouthEast_delta_p_21)/sum(Epidemic_Stage_SouthEast_delta_p_21)

#Daily death in SouthEast
daily_death_SouthEast<-SouthEast_data$daily_death[1:SouthEast_death_end]

#Epidemic stage in SoutWest
positivity_start_SouthWest<-which(positivity_SouthWest!=0)[1]#62
SouthWest_positivity_end<-which(positivity_SouthWest!=0)[length(which(positivity_SouthWest!=0))]#306
SouthWest_death_end<-which(!is.na(SouthWest_data$daily_death))[length(which(!is.na(SouthWest_data$daily_death)))]#312

Epidemic_Stage_SouthWest_delta_p_7<-c(rep(0,positivity_start_SouthWest+delta_p_7-1),positivity_SouthWest[positivity_start_SouthWest:SouthWest_positivity_end])[1:SouthWest_death_end]
Epidemic_Stage_SouthWest_delta_p_7<-cumsum(Epidemic_Stage_SouthWest_delta_p_7)/sum(Epidemic_Stage_SouthWest_delta_p_7)

Epidemic_Stage_SouthWest_delta_p_14<-c(rep(0,positivity_start_SouthWest+delta_p_14-1),positivity_SouthWest[positivity_start_SouthWest:SouthWest_positivity_end])[1:SouthWest_death_end]
Epidemic_Stage_SouthWest_delta_p_14<-cumsum(Epidemic_Stage_SouthWest_delta_p_14)/sum(Epidemic_Stage_SouthWest_delta_p_14)

Epidemic_Stage_SouthWest_delta_p_21<-c(rep(0,positivity_start_SouthWest+delta_p_21-1),positivity_SouthWest[positivity_start_SouthWest:SouthWest_positivity_end])[1:SouthWest_death_end]
Epidemic_Stage_SouthWest_delta_p_21<-cumsum(Epidemic_Stage_SouthWest_delta_p_21)/sum(Epidemic_Stage_SouthWest_delta_p_21)

#Daily death in SouthWest
daily_death_SouthWest<-SouthWest_data$daily_death[1:SouthWest_death_end]

#Epidemic stage in NorthWest
positivity_start_NorthWest<-which(positivity_NorthWest!=0)[1]#62
NorthWest_positivity_end<-which(positivity_NorthWest!=0)[length(which(positivity_NorthWest!=0))]#306
NorthWest_death_end<-which(!is.na(NorthWest_data$daily_death))[length(which(!is.na(NorthWest_data$daily_death)))]#312

Epidemic_Stage_NorthWest_delta_p_7<-c(rep(0,positivity_start_NorthWest+delta_p_7-1),positivity_NorthWest[positivity_start_NorthWest:NorthWest_positivity_end])[1:NorthWest_death_end]
Epidemic_Stage_NorthWest_delta_p_7<-cumsum(Epidemic_Stage_NorthWest_delta_p_7)/sum(Epidemic_Stage_NorthWest_delta_p_7)

Epidemic_Stage_NorthWest_delta_p_14<-c(rep(0,positivity_start_NorthWest+delta_p_14-1),positivity_NorthWest[positivity_start_NorthWest:NorthWest_positivity_end])[1:NorthWest_death_end]
Epidemic_Stage_NorthWest_delta_p_14<-cumsum(Epidemic_Stage_NorthWest_delta_p_14)/sum(Epidemic_Stage_NorthWest_delta_p_14)

Epidemic_Stage_NorthWest_delta_p_21<-c(rep(0,positivity_start_NorthWest+delta_p_21-1),positivity_NorthWest[positivity_start_NorthWest:NorthWest_positivity_end])[1:NorthWest_death_end]
Epidemic_Stage_NorthWest_delta_p_21<-cumsum(Epidemic_Stage_NorthWest_delta_p_21)/sum(Epidemic_Stage_NorthWest_delta_p_21)

#Daily death in NorthWest
daily_death_NorthWest<-NorthWest_data$daily_death[1:NorthWest_death_end]

#Epidemic stage in Midlands
positivity_start_Midlands<-which(positivity_Midlands!=0)[1]#62
Midlands_positivity_end<-which(positivity_Midlands!=0)[length(which(positivity_Midlands!=0))]#306
EastMidlands_death_end<-which(!is.na(EastMidlands_data$daily_death_EastMidlands))[length(which(!is.na(EastMidlands_data$daily_death_EastMidlands)))]#312
WestMidlands_death_end<-which(!is.na(WestMidlands_data$daily_death_WestMidlands))[length(which(!is.na(WestMidlands_data$daily_death_WestMidlands)))]

Epidemic_Stage_Midlands_delta_p_7<-c(rep(0,positivity_start_Midlands+delta_p_7-1),positivity_Midlands[positivity_start_Midlands:Midlands_positivity_end])[1:EastMidlands_death_end]
Epidemic_Stage_Midlands_delta_p_7<-cumsum(Epidemic_Stage_Midlands_delta_p_7)/sum(Epidemic_Stage_Midlands_delta_p_7)

Epidemic_Stage_Midlands_delta_p_14<-c(rep(0,positivity_start_Midlands+delta_p_14-1),positivity_Midlands[positivity_start_Midlands:Midlands_positivity_end])[1:EastMidlands_death_end]
Epidemic_Stage_Midlands_delta_p_14<-cumsum(Epidemic_Stage_Midlands_delta_p_14)/sum(Epidemic_Stage_Midlands_delta_p_14)

Epidemic_Stage_Midlands_delta_p_21<-c(rep(0,positivity_start_Midlands+delta_p_21-1),positivity_Midlands[positivity_start_Midlands:Midlands_positivity_end])[1:EastMidlands_death_end]
Epidemic_Stage_Midlands_delta_p_21<-cumsum(Epidemic_Stage_Midlands_delta_p_21)/sum(Epidemic_Stage_Midlands_delta_p_21)

#Daily death in Midlands
daily_death_Midlands<-WestMidlands_data$daily_death_WestMidlands[1:WestMidlands_death_end]+EastMidlands_data$daily_death_EastMidlands[1:EastMidlands_death_end]

#Epidemic stage in EastofEngland
positivity_start_EastofEngland<-which(positivity_EastofEngland!=0)[1]#62
EastofEngland_positivity_end<-which(positivity_EastofEngland!=0)[length(which(positivity_EastofEngland!=0))]#306
EastofEngland_death_end<-which(!is.na(EastofEngland_data$daily_death))[length(which(!is.na(EastofEngland_data$daily_death)))]#312

Epidemic_Stage_EastofEngland_delta_p_7<-c(rep(0,positivity_start_EastofEngland+delta_p_7-1),positivity_EastofEngland[positivity_start_EastofEngland:EastofEngland_positivity_end])[1:EastofEngland_death_end]
Epidemic_Stage_EastofEngland_delta_p_7<-cumsum(Epidemic_Stage_EastofEngland_delta_p_7)/sum(Epidemic_Stage_EastofEngland_delta_p_7)

Epidemic_Stage_EastofEngland_delta_p_14<-c(rep(0,positivity_start_EastofEngland+delta_p_14-1),positivity_EastofEngland[positivity_start_EastofEngland:EastofEngland_positivity_end])[1:EastofEngland_death_end]
Epidemic_Stage_EastofEngland_delta_p_14<-cumsum(Epidemic_Stage_EastofEngland_delta_p_14)/sum(Epidemic_Stage_EastofEngland_delta_p_14)

Epidemic_Stage_EastofEngland_delta_p_21<-c(rep(0,positivity_start_EastofEngland+delta_p_21-1),positivity_EastofEngland[positivity_start_EastofEngland:EastofEngland_positivity_end])[1:EastofEngland_death_end]
Epidemic_Stage_EastofEngland_delta_p_21<-cumsum(Epidemic_Stage_EastofEngland_delta_p_21)/sum(Epidemic_Stage_EastofEngland_delta_p_21)

#Daily death in EastofEngland
daily_death_EastofEngland<-EastofEngland_data$daily_death[1:EastofEngland_death_end]

#London
sero_London<-round(P0_London*London_data$sero[!is.na(London_data$sero)])
n_days_London<-length(daily_death_London)
cumul_death_London<-cumsum(daily_death_London)
n_days2_London<-length(sero_London)
t_London<-seq(1,n_days_London,by=1)
t2_London<-t_London[!is.na(London_data$sero)]

#NorthEast &Yorkshire and the Humber 
sero_NorthEastYorkshireHumber<-round(P0_NorthEastYorkshireHumber*NorthEast_data$sero[!is.na(NorthEast_data$sero)])
n_days_NorthEastYorkshireHumber<-length(daily_death_NorthEastYorkshireHumber)
cumul_death_NorthEastYorkshireHumber<-cumsum(daily_death_NorthEastYorkshireHumber)
n_days2_NorthEastYorkshireHumber<-length(sero_NorthEastYorkshireHumber)
t_NorthEastYorkshireHumber<-seq(1,n_days_NorthEastYorkshireHumber,by=1)
t2_NorthEastYorkshireHumber<-t_NorthEastYorkshireHumber[!is.na(NorthEast_data$sero)]

#SouthWest
sero_SouthWest<-round(P0_SouthWest*SouthWest_data$sero[!is.na(SouthWest_data$sero)])
n_days_SouthWest<-length(daily_death_SouthWest)
cumul_death_SouthWest<-cumsum(daily_death_SouthWest)
n_days2_SouthWest<-length(sero_SouthWest)
t_SouthWest<-seq(1,n_days_SouthWest,by=1)
t2_SouthWest<-t_SouthWest[!is.na(SouthWest_data$sero)]

#NorthWest
sero_NorthWest<-round(P0_NorthWest*NorthWest_data$sero[!is.na(NorthWest_data$sero)])
n_days_NorthWest<-length(daily_death_NorthWest)
cumul_death_NorthWest<-cumsum(daily_death_NorthWest)
n_days2_NorthWest<-length(sero_NorthWest)
t_NorthWest<-seq(1,n_days_NorthWest,by=1)
t2_NorthWest<-t_NorthWest[!is.na(NorthWest_data$sero)]

#SouthEast
sero_SouthEast<-round(P0_SouthEast*SouthEast_data$sero[!is.na(SouthEast_data$sero)])
n_days_SouthEast<-length(daily_death_SouthEast)
cumul_death_SouthEast<-cumsum(daily_death_SouthEast)
n_days2_SouthEast<-length(sero_SouthEast)
t_SouthEast<-seq(1,n_days_SouthEast,by=1)
t2_SouthEast<-t_SouthEast[!is.na(SouthEast_data$sero)]

#Midlands
sero_Midlands<-round(P0_Midlands*EastMidlands_data$sero[!is.na(EastMidlands_data$sero)])
n_days_Midlands<-length(daily_death_Midlands)
cumul_death_Midlands<-cumsum(daily_death_Midlands)
n_days2_Midlands<-length(sero_Midlands)
t_Midlands<-seq(1,n_days_Midlands,by=1)
t2_Midlands<-t_Midlands[!is.na(EastMidlands_data$sero)]

#EastofEngland
sero_EastofEngland<-round(P0_EastofEngland*EastofEngland_data$sero[!is.na(EastofEngland_data$sero)])
n_days_EastofEngland<-length(daily_death_EastofEngland)
cumul_death_EastofEngland<-cumsum(daily_death_EastofEngland)
n_days2_EastofEngland<-length(sero_EastofEngland)
t_EastofEngland<-seq(1,n_days_EastofEngland,by=1)
t2_EastofEngland<-t_EastofEngland[!is.na(EastofEngland_data$sero)]

##################
### LONDON
##################
sim<-length(beta_delta_p_7)
x_London_delta_p_7<-x_London_delta_p_14<-x_London_delta_p_21<-kft_London_delta_p_7<-kft_London_delta_p_14<-kft_London_delta_p_21<-matrix(0,sim,n_days_London)
epsilon_London_delta_p_7<-epsilon_London_delta_p_14<-epsilon_London_delta_p_21<-kft_London_cert<-matrix(0,sim,n_days_London)

for (i in 1:sim) {
  x_London_delta_p_7[i,]<-rnbinom(rep(1,n_days_London), size= 100, mu=cumsum(exp(beta_delta_p_7[i]*t_London)*(1-(gamma_London_delta_p_7[i]-(eta_London_delta_p_7[i]*gamma_London_delta_p_7[i])*Epidemic_Stage_London_delta_p_7))/(gamma_London_delta_p_7[i]-(eta_London_delta_p_7[i]*gamma_London_delta_p_7[i])*Epidemic_Stage_London_delta_p_7)*daily_death_London)/(exp(beta_delta_p_7[i]*t_London)))/(P0_London-cumul_death_London)
  x_London_delta_p_14[i,]<-rnbinom(rep(1,n_days_London), size= 100, mu=cumsum(exp(beta_delta_p_14[i]*t_London)*(1-(gamma_London_delta_p_14[i]-(eta_London_delta_p_14[i]*gamma_London_delta_p_14[i])*Epidemic_Stage_London_delta_p_14))/(gamma_London_delta_p_14[i]-(eta_London_delta_p_14[i]*gamma_London_delta_p_14[i])*Epidemic_Stage_London_delta_p_14)*daily_death_London)/(exp(beta_delta_p_14[i]*t_London)))/(P0_London-cumul_death_London)
  x_London_delta_p_21[i,]<-rnbinom(rep(1,n_days_London), size= 100, mu=cumsum(exp(beta_delta_p_21[i]*t_London)*(1-(gamma_London_delta_p_21[i]-(eta_London_delta_p_21[i]*gamma_London_delta_p_21[i])*Epidemic_Stage_London_delta_p_21))/(gamma_London_delta_p_21[i]-(eta_London_delta_p_21[i]*gamma_London_delta_p_21[i])*Epidemic_Stage_London_delta_p_21)*daily_death_London)/(exp(beta_delta_p_21[i]*t_London)))/(P0_London-cumul_death_London)
  
  epsilon_London_delta_p_7[i,]<-cumsum((1-(gamma_London_delta_p_7[i]-(eta_London_delta_p_7[i]*gamma_London_delta_p_7[i])*Epidemic_Stage_London_delta_p_7))/(gamma_London_delta_p_7[i]-(eta_London_delta_p_7[i]*gamma_London_delta_p_7[i])*Epidemic_Stage_London_delta_p_7)*daily_death_London)/(P0_London-cumul_death_London)
  epsilon_London_delta_p_14[i,]<-cumsum((1-(gamma_London_delta_p_14[i]-(eta_London_delta_p_14[i]*gamma_London_delta_p_14[i])*Epidemic_Stage_London_delta_p_14))/(gamma_London_delta_p_14[i]-(eta_London_delta_p_7[i]*gamma_London_delta_p_14[i])*Epidemic_Stage_London_delta_p_14)*daily_death_London)/(P0_London-cumul_death_London)
  epsilon_London_delta_p_21[i,]<-cumsum((1-(gamma_London_delta_p_21[i]-(eta_London_delta_p_21[i]*gamma_London_delta_p_21[i])*Epidemic_Stage_London_delta_p_21))/(gamma_London_delta_p_21[i]-(eta_London_delta_p_7[i]*gamma_London_delta_p_21[i])*Epidemic_Stage_London_delta_p_21)*daily_death_London)/(P0_London-cumul_death_London)
 
  kft_London_delta_p_7[i,]<-gamma_London_delta_p_7[i]-(eta_London_delta_p_7[i]*gamma_London_delta_p_7[i])*Epidemic_Stage_London_delta_p_7
  kft_London_delta_p_14[i,]<-gamma_London_delta_p_14[i]-(eta_London_delta_p_14[i]*gamma_London_delta_p_14[i])*Epidemic_Stage_London_delta_p_14
  kft_London_delta_p_21[i,]<-gamma_London_delta_p_21[i]-(eta_London_delta_p_21[i]*gamma_London_delta_p_21[i])*Epidemic_Stage_London_delta_p_21
}

epsilon_London_deltap7_deltaeps14<-epsilon_London_delta_p_7[,(delta_epsilon_14+1):n_days_London]
epsilon_London_deltap7_deltaeps21<-epsilon_London_delta_p_7[,(delta_epsilon_21+1):n_days_London]
epsilon_London_deltap7_deltaeps28<-epsilon_London_delta_p_7[,(delta_epsilon_28+1):n_days_London]
epsilon_London_deltap14_deltaeps21<-epsilon_London_delta_p_14[,(delta_epsilon_21+1):n_days_London]
epsilon_London_deltap14_deltaeps28<-epsilon_London_delta_p_14[,(delta_epsilon_28+1):n_days_London]
epsilon_London_deltap21_deltaeps28<-epsilon_London_delta_p_21[,(delta_epsilon_28+1):n_days_London]

data1London = data.frame(output = c(rep("Exposure (Model 6)", n_days_London-delta_epsilon_14), 
                                    rep("Exposure (Model 5)", n_days_London-delta_epsilon_21),
                                    rep("Exposure (Model 4)", n_days_London-delta_epsilon_28), 
                                    rep("Exposure (Model 7)", n_days_London-delta_epsilon_21), 
                                    rep("Exposure (Model 8)", n_days_London-delta_epsilon_28),
                                    rep("Exposure (Model 9)", n_days_London-delta_epsilon_28),
                                    rep("Seroprevalence (Model 4, 5, 6)", n_days_London),
                                    rep("Seroprevalence (Model 7, 8)", n_days_London), 
                                    rep("Seroprevalence (Model 9)", n_days_London)), 
                         t=c(as.Date(London_data$Date)[1:(n_days_London-delta_epsilon_14)],
                             as.Date(London_data$Date)[1:(n_days_London-delta_epsilon_21)],
                             as.Date(London_data$Date)[1:(n_days_London-delta_epsilon_28)],
                             as.Date(London_data$Date)[1:(n_days_London-delta_epsilon_21)],
                             as.Date(London_data$Date)[1:(n_days_London-delta_epsilon_28)],
                             as.Date(London_data$Date)[1:(n_days_London-delta_epsilon_28)],
                             as.Date(London_data$Date)[1:n_days_London],
                             as.Date(London_data$Date)[1:n_days_London],
                             as.Date(London_data$Date)[1:n_days_London]), 
                         median = c(100*apply(epsilon_London_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.5)), 
                                    100*apply(epsilon_London_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                    100*apply(epsilon_London_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                    100*apply(epsilon_London_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                    100*apply(epsilon_London_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                    100*apply(epsilon_London_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.5)),
                                    100*apply(x_London_delta_p_7, 2, function(x) quantile(x, probs = 0.5)),
                                    100*apply(x_London_delta_p_14, 2, function(x) quantile(x, probs = 0.5)), 
                                    100*apply(x_London_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                         lower1 = c(100*apply(epsilon_London_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.025)), 
                                    100*apply(epsilon_London_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                    100*apply(epsilon_London_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                    100*apply(epsilon_London_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                    100*apply(epsilon_London_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                    100*apply(epsilon_London_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.025)),
                                    100*apply(x_London_delta_p_7, 2, function(x) quantile(x, probs = 0.025)),
                                    100*apply(x_London_delta_p_14, 2, function(x) quantile(x, probs = 0.025)), 
                                    100*apply(x_London_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                         upper1 = c(100*apply(epsilon_London_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.975)), 
                                    100*apply(epsilon_London_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(epsilon_London_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                    100*apply(epsilon_London_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(epsilon_London_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                    100*apply(epsilon_London_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(x_London_delta_p_7, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(x_London_delta_p_14, 2, function(x) quantile(x, probs = 0.975)), 
                                    100*apply(x_London_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                         lower2 = c(100*apply(epsilon_London_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.25)), 
                                    100*apply(epsilon_London_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                    100* apply(epsilon_London_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                    100*apply(epsilon_London_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                    100*apply(epsilon_London_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                    100*apply(epsilon_London_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.25)),
                                    100*apply(x_London_delta_p_7, 2, function(x) quantile(x, probs = 0.25)),
                                    100*apply(x_London_delta_p_14, 2, function(x) quantile(x, probs = 0.25)), 
                                    100*apply(x_London_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                         upper2 = c(100*apply(epsilon_London_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.75)), 
                                    100*apply(epsilon_London_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                    100*apply(epsilon_London_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                    100*apply(epsilon_London_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                    100*apply(epsilon_London_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                    100*apply(epsilon_London_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.75)),
                                    100*apply(x_London_delta_p_7, 2, function(x) quantile(x, probs = 0.75)),
                                    100*apply(x_London_delta_p_14, 2, function(x) quantile(x, probs = 0.75)), 
                                    100*apply(x_London_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))

data2London = data.frame( t=as.Date(London_data$Date)[1:n_days_London][t2_London], value= 100*London_data$sero[t2_London], upper= 100*London_data$sero_upper[t2_London], lower = 100*London_data$sero_lower[t2_London])

p1London<-ggplot(data1London, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("London")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2London, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[4])
styled1London <- p1London +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" Percentage (%) ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.text.x = element_text(size=font_size),
    axis.text.y = element_text(size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1London

##################
### NorthWest
##################
sim<-length(beta_delta_p_7)
x_NorthWest_delta_p_7<-x_NorthWest_delta_p_14<-x_NorthWest_delta_p_21<-kft_NorthWest<-kft_NorthWest_delta_p_7<-kft_NorthWest_delta_p_14<-kft_NorthWest_delta_p_21<-matrix(0,sim,n_days_NorthWest)
epsilon_NorthWest_delta_p_7<-epsilon_NorthWest_delta_p_14<-epsilon_NorthWest_delta_p_21<-kft_NorthWest_cert<-matrix(0,sim,n_days_NorthWest)

for (i in 1:sim) {
  x_NorthWest_delta_p_7[i,]<-rnbinom(rep(1,n_days_NorthWest), size= 100, mu=cumsum(exp(beta_delta_p_7[i]*t_NorthWest)*(1-(gamma_NorthWest_delta_p_7[i]-(eta_NorthWest_delta_p_7[i]*gamma_NorthWest_delta_p_7[i])*Epidemic_Stage_NorthWest_delta_p_7))/(gamma_NorthWest_delta_p_7[i]-(eta_NorthWest_delta_p_7[i]*gamma_NorthWest_delta_p_7[i])*Epidemic_Stage_NorthWest_delta_p_7)*daily_death_NorthWest)/(exp(beta_delta_p_7[i]*t_NorthWest)))/(P0_NorthWest-cumul_death_NorthWest)
  x_NorthWest_delta_p_14[i,]<-rnbinom(rep(1,n_days_NorthWest), size= 100, mu=cumsum(exp(beta_delta_p_14[i]*t_NorthWest)*(1-(gamma_NorthWest_delta_p_14[i]-(eta_NorthWest_delta_p_14[i]*gamma_NorthWest_delta_p_14[i])*Epidemic_Stage_NorthWest_delta_p_14))/(gamma_NorthWest_delta_p_14[i]-(eta_NorthWest_delta_p_14[i]*gamma_NorthWest_delta_p_14[i])*Epidemic_Stage_NorthWest_delta_p_14)*daily_death_NorthWest)/(exp(beta_delta_p_14[i]*t_NorthWest)))/(P0_NorthWest-cumul_death_NorthWest)
  x_NorthWest_delta_p_21[i,]<-rnbinom(rep(1,n_days_NorthWest), size= 100, mu=cumsum(exp(beta_delta_p_21[i]*t_NorthWest)*(1-(gamma_NorthWest_delta_p_21[i]-(eta_NorthWest_delta_p_21[i]*gamma_NorthWest_delta_p_21[i])*Epidemic_Stage_NorthWest_delta_p_21))/(gamma_NorthWest_delta_p_21[i]-(eta_NorthWest_delta_p_21[i]*gamma_NorthWest_delta_p_21[i])*Epidemic_Stage_NorthWest_delta_p_21)*daily_death_NorthWest)/(exp(beta_delta_p_21[i]*t_NorthWest)))/(P0_NorthWest-cumul_death_NorthWest)
  
  epsilon_NorthWest_delta_p_7[i,]<-cumsum((1-(gamma_NorthWest_delta_p_7[i]-(eta_NorthWest_delta_p_7[i]*gamma_NorthWest_delta_p_7[i])*Epidemic_Stage_NorthWest_delta_p_7))/(gamma_NorthWest_delta_p_7[i]-(eta_NorthWest_delta_p_7[i]*gamma_NorthWest_delta_p_7[i])*Epidemic_Stage_NorthWest_delta_p_7)*daily_death_NorthWest)/(P0_NorthWest-cumul_death_NorthWest)
  epsilon_NorthWest_delta_p_14[i,]<-cumsum((1-(gamma_NorthWest_delta_p_14[i]-(eta_NorthWest_delta_p_14[i]*gamma_NorthWest_delta_p_14[i])*Epidemic_Stage_NorthWest_delta_p_14))/(gamma_NorthWest_delta_p_14[i]-(eta_NorthWest_delta_p_7[i]*gamma_NorthWest_delta_p_14[i])*Epidemic_Stage_NorthWest_delta_p_14)*daily_death_NorthWest)/(P0_NorthWest-cumul_death_NorthWest)
  epsilon_NorthWest_delta_p_21[i,]<-cumsum((1-(gamma_NorthWest_delta_p_21[i]-(eta_NorthWest_delta_p_21[i]*gamma_NorthWest_delta_p_21[i])*Epidemic_Stage_NorthWest_delta_p_21))/(gamma_NorthWest_delta_p_21[i]-(eta_NorthWest_delta_p_7[i]*gamma_NorthWest_delta_p_21[i])*Epidemic_Stage_NorthWest_delta_p_21)*daily_death_NorthWest)/(P0_NorthWest-cumul_death_NorthWest)

  kft_NorthWest_delta_p_7[i,]<-gamma_NorthWest_delta_p_7[i]-(eta_NorthWest_delta_p_7[i]*gamma_NorthWest_delta_p_7[i])*Epidemic_Stage_NorthWest_delta_p_7
  kft_NorthWest_delta_p_14[i,]<-gamma_NorthWest_delta_p_14[i]-(eta_NorthWest_delta_p_14[i]*gamma_NorthWest_delta_p_14[i])*Epidemic_Stage_NorthWest_delta_p_14
  kft_NorthWest_delta_p_21[i,]<-gamma_NorthWest_delta_p_21[i]-(eta_NorthWest_delta_p_21[i]*gamma_NorthWest_delta_p_21[i])*Epidemic_Stage_NorthWest_delta_p_21
  }

epsilon_NorthWest_deltap7_deltaeps14<-epsilon_NorthWest_delta_p_7[,(delta_epsilon_14+1):n_days_NorthWest]
epsilon_NorthWest_deltap7_deltaeps21<-epsilon_NorthWest_delta_p_7[,(delta_epsilon_21+1):n_days_NorthWest]
epsilon_NorthWest_deltap7_deltaeps28<-epsilon_NorthWest_delta_p_7[,(delta_epsilon_28+1):n_days_NorthWest]
epsilon_NorthWest_deltap14_deltaeps21<-epsilon_NorthWest_delta_p_14[,(delta_epsilon_21+1):n_days_NorthWest]
epsilon_NorthWest_deltap14_deltaeps28<-epsilon_NorthWest_delta_p_14[,(delta_epsilon_28+1):n_days_NorthWest]
epsilon_NorthWest_deltap21_deltaeps28<-epsilon_NorthWest_delta_p_21[,(delta_epsilon_28+1):n_days_NorthWest]

data1NorthWest = data.frame(output = c(rep("Exposure (Model 6)", n_days_London-delta_epsilon_14), 
                                       rep("Exposure (Model 5)", n_days_London-delta_epsilon_21),
                                       rep("Exposure (Model 4)", n_days_London-delta_epsilon_28), 
                                       rep("Exposure (Model 7)", n_days_London-delta_epsilon_21), 
                                       rep("Exposure (Model 8)", n_days_London-delta_epsilon_28),
                                       rep("Exposure (Model 9)", n_days_London-delta_epsilon_28),
                                       rep("Seroprevalence (Model 4, 5, 6)", n_days_London),
                                       rep("Seroprevalence (Model 7, 8)", n_days_London), 
                                       rep("Seroprevalence (Model 9)", n_days_London)), 
                            t=c(as.Date(NorthWest_data$Date)[1:(n_days_NorthWest-delta_epsilon_14)],
                                as.Date(NorthWest_data$Date)[1:(n_days_NorthWest-delta_epsilon_21)],
                                as.Date(NorthWest_data$Date)[1:(n_days_NorthWest-delta_epsilon_28)],
                                as.Date(NorthWest_data$Date)[1:(n_days_NorthWest-delta_epsilon_21)],
                                as.Date(NorthWest_data$Date)[1:(n_days_NorthWest-delta_epsilon_28)],
                                as.Date(NorthWest_data$Date)[1:(n_days_NorthWest-delta_epsilon_28)],
                                as.Date(NorthWest_data$Date)[1:n_days_NorthWest],
                                as.Date(NorthWest_data$Date)[1:n_days_NorthWest],
                                as.Date(NorthWest_data$Date)[1:n_days_NorthWest]), 
                            median = c(100*apply(epsilon_NorthWest_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(epsilon_NorthWest_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_NorthWest_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(epsilon_NorthWest_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_NorthWest_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(epsilon_NorthWest_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(x_NorthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(x_NorthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_NorthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                            lower1 = c(100*apply(epsilon_NorthWest_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(epsilon_NorthWest_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_NorthWest_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(epsilon_NorthWest_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_NorthWest_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(epsilon_NorthWest_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(x_NorthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(x_NorthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(x_NorthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                            upper1 = c(100*apply(epsilon_NorthWest_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.975)), 
                                       100*apply(epsilon_NorthWest_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_NorthWest_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                       100*apply(epsilon_NorthWest_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_NorthWest_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                       100*apply(epsilon_NorthWest_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_NorthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_NorthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.975)), 
                                       100*apply(x_NorthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                            lower2 = c(100*apply(epsilon_NorthWest_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(epsilon_NorthWest_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_NorthWest_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(epsilon_NorthWest_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_NorthWest_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(epsilon_NorthWest_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(x_NorthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(x_NorthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_NorthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                            upper2 = c(100*apply(epsilon_NorthWest_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.75)), 
                                       100*apply(epsilon_NorthWest_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_NorthWest_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                       100*apply(epsilon_NorthWest_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_NorthWest_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                       100*apply(epsilon_NorthWest_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_NorthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_NorthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.75)), 
                                       100*apply(x_NorthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))

data2NorthWest = data.frame( t=as.Date(NorthWest_data$Date)[1:n_days_NorthWest][t2_NorthWest], value= 100*NorthWest_data$sero[t2_NorthWest], upper= 100*NorthWest_data$sero_upper[t2_NorthWest], lower = 100*NorthWest_data$sero_lower[t2_NorthWest])

p1NorthWest<-ggplot(data1NorthWest, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("North West")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2NorthWest, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[4])
styled1NorthWest <- p1NorthWest +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
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
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.text.x = element_text(size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1NorthWest

##################
### SouthWest
##################
sim<-length(beta_delta_p_7)
x_SouthWest_delta_p_7<-x_SouthWest_delta_p_14<-x_SouthWest_delta_p_21<-kft_SouthWest<-kft_SouthWest_delta_p_7<-kft_SouthWest_delta_p_14<-kft_SouthWest_delta_p_21<-matrix(0,sim,n_days_SouthWest)
epsilon_SouthWest_delta_p_7<-epsilon_SouthWest_delta_p_14<-epsilon_SouthWest_delta_p_21<-kft_SouthWest_cert<-matrix(0,sim,n_days_SouthWest)

for (i in 1:sim) {
  x_SouthWest_delta_p_7[i,]<-rnbinom(rep(1,n_days_SouthWest), size= 100, mu=cumsum(exp(beta_delta_p_7[i]*t_SouthWest)*(1-(gamma_SouthWest_delta_p_7[i]-(eta_SouthWest_delta_p_7[i]*gamma_SouthWest_delta_p_7[i])*Epidemic_Stage_SouthWest_delta_p_7))/(gamma_SouthWest_delta_p_7[i]-(eta_SouthWest_delta_p_7[i]*gamma_SouthWest_delta_p_7[i])*Epidemic_Stage_SouthWest_delta_p_7)*daily_death_SouthWest)/(exp(beta_delta_p_7[i]*t_SouthWest)))/(P0_SouthWest-cumul_death_SouthWest)
  x_SouthWest_delta_p_14[i,]<-rnbinom(rep(1,n_days_SouthWest), size= 100, mu=cumsum(exp(beta_delta_p_14[i]*t_SouthWest)*(1-(gamma_SouthWest_delta_p_14[i]-(eta_SouthWest_delta_p_14[i]*gamma_SouthWest_delta_p_14[i])*Epidemic_Stage_SouthWest_delta_p_14))/(gamma_SouthWest_delta_p_14[i]-(eta_SouthWest_delta_p_14[i]*gamma_SouthWest_delta_p_14[i])*Epidemic_Stage_SouthWest_delta_p_14)*daily_death_SouthWest)/(exp(beta_delta_p_14[i]*t_SouthWest)))/(P0_SouthWest-cumul_death_SouthWest)
  x_SouthWest_delta_p_21[i,]<-rnbinom(rep(1,n_days_SouthWest), size= 100, mu=cumsum(exp(beta_delta_p_21[i]*t_SouthWest)*(1-(gamma_SouthWest_delta_p_21[i]-(eta_SouthWest_delta_p_21[i]*gamma_SouthWest_delta_p_21[i])*Epidemic_Stage_SouthWest_delta_p_21))/(gamma_SouthWest_delta_p_21[i]-(eta_SouthWest_delta_p_21[i]*gamma_SouthWest_delta_p_21[i])*Epidemic_Stage_SouthWest_delta_p_21)*daily_death_SouthWest)/(exp(beta_delta_p_21[i]*t_SouthWest)))/(P0_SouthWest-cumul_death_SouthWest)
  
  epsilon_SouthWest_delta_p_7[i,]<-cumsum((1-(gamma_SouthWest_delta_p_7[i]-(eta_SouthWest_delta_p_7[i]*gamma_SouthWest_delta_p_7[i])*Epidemic_Stage_SouthWest_delta_p_7))/(gamma_SouthWest_delta_p_7[i]-(eta_SouthWest_delta_p_7[i]*gamma_SouthWest_delta_p_7[i])*Epidemic_Stage_SouthWest_delta_p_7)*daily_death_SouthWest)/(P0_SouthWest-cumul_death_SouthWest)
  epsilon_SouthWest_delta_p_14[i,]<-cumsum((1-(gamma_SouthWest_delta_p_14[i]-(eta_SouthWest_delta_p_14[i]*gamma_SouthWest_delta_p_14[i])*Epidemic_Stage_SouthWest_delta_p_14))/(gamma_SouthWest_delta_p_14[i]-(eta_SouthWest_delta_p_7[i]*gamma_SouthWest_delta_p_14[i])*Epidemic_Stage_SouthWest_delta_p_14)*daily_death_SouthWest)/(P0_SouthWest-cumul_death_SouthWest)
  epsilon_SouthWest_delta_p_21[i,]<-cumsum((1-(gamma_SouthWest_delta_p_21[i]-(eta_SouthWest_delta_p_21[i]*gamma_SouthWest_delta_p_21[i])*Epidemic_Stage_SouthWest_delta_p_21))/(gamma_SouthWest_delta_p_21[i]-(eta_SouthWest_delta_p_7[i]*gamma_SouthWest_delta_p_21[i])*Epidemic_Stage_SouthWest_delta_p_21)*daily_death_SouthWest)/(P0_SouthWest-cumul_death_SouthWest)

  kft_SouthWest_delta_p_7[i,]<-gamma_SouthWest_delta_p_7[i]-(eta_SouthWest_delta_p_7[i]*gamma_SouthWest_delta_p_7[i])*Epidemic_Stage_SouthWest_delta_p_7
  kft_SouthWest_delta_p_14[i,]<-gamma_SouthWest_delta_p_14[i]-(eta_SouthWest_delta_p_14[i]*gamma_SouthWest_delta_p_14[i])*Epidemic_Stage_SouthWest_delta_p_14
  kft_SouthWest_delta_p_21[i,]<-gamma_SouthWest_delta_p_21[i]-(eta_SouthWest_delta_p_21[i]*gamma_SouthWest_delta_p_21[i])*Epidemic_Stage_SouthWest_delta_p_21
  }

epsilon_SouthWest_deltap7_deltaeps14<-epsilon_SouthWest_delta_p_7[,(delta_epsilon_14+1):n_days_SouthWest]
epsilon_SouthWest_deltap7_deltaeps21<-epsilon_SouthWest_delta_p_7[,(delta_epsilon_21+1):n_days_SouthWest]
epsilon_SouthWest_deltap7_deltaeps28<-epsilon_SouthWest_delta_p_7[,(delta_epsilon_28+1):n_days_SouthWest]
epsilon_SouthWest_deltap14_deltaeps21<-epsilon_SouthWest_delta_p_14[,(delta_epsilon_21+1):n_days_SouthWest]
epsilon_SouthWest_deltap14_deltaeps28<-epsilon_SouthWest_delta_p_14[,(delta_epsilon_28+1):n_days_SouthWest]
epsilon_SouthWest_deltap21_deltaeps28<-epsilon_SouthWest_delta_p_21[,(delta_epsilon_28+1):n_days_SouthWest]

data1SouthWest = data.frame(output = c(rep("Exposure (Scenario 6)", n_days_London-delta_epsilon_14), 
                                       rep("Exposure (Scenario 5)", n_days_London-delta_epsilon_21),
                                       rep("Exposure (Scenario 4)", n_days_London-delta_epsilon_28), 
                                       rep("Exposure (Scenario 7)", n_days_London-delta_epsilon_21), 
                                       rep("Exposure (Scenario 8)", n_days_London-delta_epsilon_28),
                                       rep("Exposure (Scenario 9)", n_days_London-delta_epsilon_28),
                                       rep("Seroprevalence (Scenario 4, 5, 6)", n_days_London),
                                       rep("Seroprevalence (Scenario 7, 8)", n_days_London), 
                                       rep("Seroprevalence (Scenario 9)", n_days_London)), 
                            t=c(as.Date(SouthWest_data$Date)[1:(n_days_SouthWest-delta_epsilon_14)],
                                as.Date(SouthWest_data$Date)[1:(n_days_SouthWest-delta_epsilon_21)],
                                as.Date(SouthWest_data$Date)[1:(n_days_SouthWest-delta_epsilon_28)],
                                as.Date(SouthWest_data$Date)[1:(n_days_SouthWest-delta_epsilon_21)],
                                as.Date(SouthWest_data$Date)[1:(n_days_SouthWest-delta_epsilon_28)],
                                as.Date(SouthWest_data$Date)[1:(n_days_SouthWest-delta_epsilon_28)],
                                as.Date(SouthWest_data$Date)[1:n_days_SouthWest],
                                as.Date(SouthWest_data$Date)[1:n_days_SouthWest],
                                as.Date(SouthWest_data$Date)[1:n_days_SouthWest]), 
                            median = c(100*apply(epsilon_SouthWest_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(epsilon_SouthWest_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_SouthWest_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(epsilon_SouthWest_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_SouthWest_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(epsilon_SouthWest_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(x_SouthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(x_SouthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_SouthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                            lower1 = c(100*apply(epsilon_SouthWest_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(epsilon_SouthWest_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_SouthWest_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(epsilon_SouthWest_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_SouthWest_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(epsilon_SouthWest_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(x_SouthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(x_SouthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(x_SouthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                            upper1 = c(100*apply(epsilon_SouthWest_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.975)), 
                                       100*apply(epsilon_SouthWest_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_SouthWest_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                       100*apply(epsilon_SouthWest_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_SouthWest_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                       100*apply(epsilon_SouthWest_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_SouthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_SouthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.975)), 
                                       100*apply(x_SouthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                            lower2 = c(100*apply(epsilon_SouthWest_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(epsilon_SouthWest_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_SouthWest_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(epsilon_SouthWest_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_SouthWest_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(epsilon_SouthWest_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(x_SouthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(x_SouthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_SouthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                            upper2 = c(100*apply(epsilon_SouthWest_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.75)), 
                                       100*apply(epsilon_SouthWest_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_SouthWest_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                       100*apply(epsilon_SouthWest_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_SouthWest_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                       100*apply(epsilon_SouthWest_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_SouthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_SouthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.75)), 
                                       100*apply(x_SouthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))

data2SouthWest = data.frame( t=as.Date(SouthWest_data$Date)[1:n_days_SouthWest][t2_SouthWest], value= 100*SouthWest_data$sero[t2_SouthWest], upper= 100*SouthWest_data$sero_upper[t2_SouthWest], lower = 100*SouthWest_data$sero_lower[t2_SouthWest])

p1SouthWest<-ggplot(data1SouthWest, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("South West")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2SouthWest, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[4])
styled1SouthWest <- p1SouthWest +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" Percentage (%) ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    axis.text.x = element_text(size=font_size),
    axis.text.y = element_text(size=font_size),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1SouthWest

##################
### SouthEast
##################
sim<-length(beta_delta_p_7)
x_SouthEast_delta_p_7<-x_SouthEast_delta_p_14<-x_SouthEast_delta_p_21<-kft_SouthEast<-kft_SouthEast_delta_p_7<-kft_SouthEast_delta_p_14<-kft_SouthEast_delta_p_21<-matrix(0,sim,n_days_SouthEast)
epsilon_SouthEast_delta_p_7<-epsilon_SouthEast_delta_p_14<-epsilon_SouthEast_delta_p_21<-matrix(0,sim,n_days_SouthEast)

for (i in 1:sim) {
  x_SouthEast_delta_p_7[i,]<-rnbinom(rep(1,n_days_SouthEast), size= 100, mu=cumsum(exp(beta_delta_p_7[i]*t_SouthEast)*(1-(gamma_SouthEast_delta_p_7[i]-(eta_SouthEast_delta_p_7[i]*gamma_SouthEast_delta_p_7[i])*Epidemic_Stage_SouthEast_delta_p_7))/(gamma_SouthEast_delta_p_7[i]-(eta_SouthEast_delta_p_7[i]*gamma_SouthEast_delta_p_7[i])*Epidemic_Stage_SouthEast_delta_p_7)*daily_death_SouthEast)/(exp(beta_delta_p_7[i]*t_SouthEast)))/(P0_SouthEast-cumul_death_SouthEast)
  x_SouthEast_delta_p_14[i,]<-rnbinom(rep(1,n_days_SouthEast), size= 100, mu=cumsum(exp(beta_delta_p_14[i]*t_SouthEast)*(1-(gamma_SouthEast_delta_p_14[i]-(eta_SouthEast_delta_p_14[i]*gamma_SouthEast_delta_p_14[i])*Epidemic_Stage_SouthEast_delta_p_14))/(gamma_SouthEast_delta_p_14[i]-(eta_SouthEast_delta_p_14[i]*gamma_SouthEast_delta_p_14[i])*Epidemic_Stage_SouthEast_delta_p_14)*daily_death_SouthEast)/(exp(beta_delta_p_14[i]*t_SouthEast)))/(P0_SouthEast-cumul_death_SouthEast)
  x_SouthEast_delta_p_21[i,]<-rnbinom(rep(1,n_days_SouthEast), size= 100, mu=cumsum(exp(beta_delta_p_21[i]*t_SouthEast)*(1-(gamma_SouthEast_delta_p_21[i]-(eta_SouthEast_delta_p_21[i]*gamma_SouthEast_delta_p_21[i])*Epidemic_Stage_SouthEast_delta_p_21))/(gamma_SouthEast_delta_p_21[i]-(eta_SouthEast_delta_p_21[i]*gamma_SouthEast_delta_p_21[i])*Epidemic_Stage_SouthEast_delta_p_21)*daily_death_SouthEast)/(exp(beta_delta_p_21[i]*t_SouthEast)))/(P0_SouthEast-cumul_death_SouthEast)
  
  epsilon_SouthEast_delta_p_7[i,]<-cumsum((1-(gamma_SouthEast_delta_p_7[i]-(eta_SouthEast_delta_p_7[i]*gamma_SouthEast_delta_p_7[i])*Epidemic_Stage_SouthEast_delta_p_7))/(gamma_SouthEast_delta_p_7[i]-(eta_SouthEast_delta_p_7[i]*gamma_SouthEast_delta_p_7[i])*Epidemic_Stage_SouthEast_delta_p_7)*daily_death_SouthEast)/(P0_SouthEast-cumul_death_SouthEast)
  epsilon_SouthEast_delta_p_14[i,]<-cumsum((1-(gamma_SouthEast_delta_p_14[i]-(eta_SouthEast_delta_p_14[i]*gamma_SouthEast_delta_p_14[i])*Epidemic_Stage_SouthEast_delta_p_14))/(gamma_SouthEast_delta_p_14[i]-(eta_SouthEast_delta_p_7[i]*gamma_SouthEast_delta_p_14[i])*Epidemic_Stage_SouthEast_delta_p_14)*daily_death_SouthEast)/(P0_SouthEast-cumul_death_SouthEast)
  epsilon_SouthEast_delta_p_21[i,]<-cumsum((1-(gamma_SouthEast_delta_p_21[i]-(eta_SouthEast_delta_p_21[i]*gamma_SouthEast_delta_p_21[i])*Epidemic_Stage_SouthEast_delta_p_21))/(gamma_SouthEast_delta_p_21[i]-(eta_SouthEast_delta_p_7[i]*gamma_SouthEast_delta_p_21[i])*Epidemic_Stage_SouthEast_delta_p_21)*daily_death_SouthEast)/(P0_SouthEast-cumul_death_SouthEast)
 
  kft_SouthEast_delta_p_7[i,]<-gamma_SouthEast_delta_p_7[i]-(eta_SouthEast_delta_p_7[i]*gamma_SouthEast_delta_p_7[i])*Epidemic_Stage_SouthEast_delta_p_7
  kft_SouthEast_delta_p_14[i,]<-gamma_SouthEast_delta_p_14[i]-(eta_SouthEast_delta_p_14[i]*gamma_SouthEast_delta_p_14[i])*Epidemic_Stage_SouthEast_delta_p_14
  kft_SouthEast_delta_p_21[i,]<-gamma_SouthEast_delta_p_21[i]-(eta_SouthEast_delta_p_21[i]*gamma_SouthEast_delta_p_21[i])*Epidemic_Stage_SouthEast_delta_p_21
  
  }

epsilon_SouthEast_deltap7_deltaeps14<-epsilon_SouthEast_delta_p_7[,(delta_epsilon_14+1):n_days_SouthEast]
epsilon_SouthEast_deltap7_deltaeps21<-epsilon_SouthEast_delta_p_7[,(delta_epsilon_21+1):n_days_SouthEast]
epsilon_SouthEast_deltap7_deltaeps28<-epsilon_SouthEast_delta_p_7[,(delta_epsilon_28+1):n_days_SouthEast]
epsilon_SouthEast_deltap14_deltaeps21<-epsilon_SouthEast_delta_p_14[,(delta_epsilon_21+1):n_days_SouthEast]
epsilon_SouthEast_deltap14_deltaeps28<-epsilon_SouthEast_delta_p_14[,(delta_epsilon_28+1):n_days_SouthEast]
epsilon_SouthEast_deltap21_deltaeps28<-epsilon_SouthEast_delta_p_21[,(delta_epsilon_28+1):n_days_SouthEast]

data1SouthEast = data.frame(output = c(rep("Exposure (Model 6)", n_days_London-delta_epsilon_14), 
                                       rep("Exposure (Model 5)", n_days_London-delta_epsilon_21),
                                       rep("Exposure (Model 4)", n_days_London-delta_epsilon_28), 
                                       rep("Exposure (Model 7)", n_days_London-delta_epsilon_21), 
                                       rep("Exposure (Model 8)", n_days_London-delta_epsilon_28),
                                       rep("Exposure (Model 9)", n_days_London-delta_epsilon_28),
                                       rep("Seroprevalence (Model 4, 5, 6)", n_days_London),
                                       rep("Seroprevalence (Model 7,8)", n_days_London), 
                                       rep("Seroprevalence (Model 9)", n_days_London)), 
                            t=c(as.Date(SouthEast_data$Date)[1:(n_days_SouthEast-delta_epsilon_14)],
                                as.Date(SouthEast_data$Date)[1:(n_days_SouthEast-delta_epsilon_21)],
                                as.Date(SouthEast_data$Date)[1:(n_days_SouthEast-delta_epsilon_28)],
                                as.Date(SouthEast_data$Date)[1:(n_days_SouthEast-delta_epsilon_21)],
                                as.Date(SouthEast_data$Date)[1:(n_days_SouthEast-delta_epsilon_28)],
                                as.Date(SouthEast_data$Date)[1:(n_days_SouthEast-delta_epsilon_28)],
                                as.Date(SouthEast_data$Date)[1:n_days_SouthEast],
                                as.Date(SouthEast_data$Date)[1:n_days_SouthEast],
                                as.Date(SouthEast_data$Date)[1:n_days_SouthEast]), 
                            median = c(100*apply(epsilon_SouthEast_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(epsilon_SouthEast_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_SouthEast_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(epsilon_SouthEast_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_SouthEast_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(epsilon_SouthEast_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(x_SouthEast_delta_p_7, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(x_SouthEast_delta_p_14, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_SouthEast_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                            lower1 = c(100*apply(epsilon_SouthEast_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(epsilon_SouthEast_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_SouthEast_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(epsilon_SouthEast_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_SouthEast_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(epsilon_SouthEast_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(x_SouthEast_delta_p_7, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(x_SouthEast_delta_p_14, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(x_SouthEast_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                            upper1 = c(100*apply(epsilon_SouthEast_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.975)), 
                                       100*apply(epsilon_SouthEast_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_SouthEast_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                       100*apply(epsilon_SouthEast_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_SouthEast_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                       100*apply(epsilon_SouthEast_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_SouthEast_delta_p_7, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_SouthEast_delta_p_14, 2, function(x) quantile(x, probs = 0.975)), 
                                       100*apply(x_SouthEast_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                            lower2 = c(100*apply(epsilon_SouthEast_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(epsilon_SouthEast_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_SouthEast_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(epsilon_SouthEast_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_SouthEast_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(epsilon_SouthEast_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(x_SouthEast_delta_p_7, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(x_SouthEast_delta_p_14, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_SouthEast_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                            upper2 = c(100*apply(epsilon_SouthEast_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.75)), 
                                       100*apply(epsilon_SouthEast_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_SouthEast_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                       100*apply(epsilon_SouthEast_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_SouthEast_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                       100*apply(epsilon_SouthEast_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_SouthEast_delta_p_7, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_SouthEast_delta_p_14, 2, function(x) quantile(x, probs = 0.75)), 
                                       100*apply(x_SouthEast_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))

data2SouthEast = data.frame( t=as.Date(SouthEast_data$Date)[1:n_days_SouthEast][t2_SouthEast], value= 100*SouthEast_data$sero[t2_SouthEast], upper= 100*SouthEast_data$sero_upper[t2_SouthEast], lower = 100*SouthEast_data$sero_lower[t2_SouthEast])

p1SouthEast<-ggplot(data1SouthEast, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("South East")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2SouthEast, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[4])
styled1SouthEast <- p1SouthEast +
   # scale_fill_brewer(palette = "Dark2")+
   # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    # legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(-0.01,-0.2),
    legend.text = element_text(size=12),
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=font_size),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1SouthEast

##################
### EastofEngland
##################
sim<-length(beta_delta_p_7)
x_EastofEngland_delta_p_7<-x_EastofEngland_delta_p_14<-x_EastofEngland_delta_p_21<-kft_EastofEngland<-kft_EastofEngland_delta_p_7<-kft_EastofEngland_delta_p_14<-kft_EastofEngland_delta_p_21<-matrix(0,sim,n_days_EastofEngland)
epsilon_EastofEngland_delta_p_7<-epsilon_EastofEngland_delta_p_14<-epsilon_EastofEngland_delta_p_21<-kft_EastofEngland_cert<-matrix(0,sim,n_days_EastofEngland)

for (i in 1:sim) {
  x_EastofEngland_delta_p_7[i,]<-rnbinom(rep(1,n_days_EastofEngland), size= 100, mu=cumsum(exp(beta_delta_p_7[i]*t_EastofEngland)*(1-(gamma_EastofEngland_delta_p_7[i]-(eta_EastofEngland_delta_p_7[i]*gamma_EastofEngland_delta_p_7[i])*Epidemic_Stage_EastofEngland_delta_p_7))/(gamma_EastofEngland_delta_p_7[i]-(eta_EastofEngland_delta_p_7[i]*gamma_EastofEngland_delta_p_7[i])*Epidemic_Stage_EastofEngland_delta_p_7)*daily_death_EastofEngland)/(exp(beta_delta_p_7[i]*t_EastofEngland)))/(P0_EastofEngland-cumul_death_EastofEngland)
  x_EastofEngland_delta_p_14[i,]<-rnbinom(rep(1,n_days_EastofEngland), size= 100, mu=cumsum(exp(beta_delta_p_14[i]*t_EastofEngland)*(1-(gamma_EastofEngland_delta_p_14[i]-(eta_EastofEngland_delta_p_14[i]*gamma_EastofEngland_delta_p_14[i])*Epidemic_Stage_EastofEngland_delta_p_14))/(gamma_EastofEngland_delta_p_14[i]-(eta_EastofEngland_delta_p_14[i]*gamma_EastofEngland_delta_p_14[i])*Epidemic_Stage_EastofEngland_delta_p_14)*daily_death_EastofEngland)/(exp(beta_delta_p_14[i]*t_EastofEngland)))/(P0_EastofEngland-cumul_death_EastofEngland)
  x_EastofEngland_delta_p_21[i,]<-rnbinom(rep(1,n_days_EastofEngland), size= 100, mu=cumsum(exp(beta_delta_p_21[i]*t_EastofEngland)*(1-(gamma_EastofEngland_delta_p_21[i]-(eta_EastofEngland_delta_p_21[i]*gamma_EastofEngland_delta_p_21[i])*Epidemic_Stage_EastofEngland_delta_p_21))/(gamma_EastofEngland_delta_p_21[i]-(eta_EastofEngland_delta_p_21[i]*gamma_EastofEngland_delta_p_21[i])*Epidemic_Stage_EastofEngland_delta_p_21)*daily_death_EastofEngland)/(exp(beta_delta_p_21[i]*t_EastofEngland)))/(P0_EastofEngland-cumul_death_EastofEngland)
  
  epsilon_EastofEngland_delta_p_7[i,]<-cumsum((1-(gamma_EastofEngland_delta_p_7[i]-(eta_EastofEngland_delta_p_7[i]*gamma_EastofEngland_delta_p_7[i])*Epidemic_Stage_EastofEngland_delta_p_7))/(gamma_EastofEngland_delta_p_7[i]-(eta_EastofEngland_delta_p_7[i]*gamma_EastofEngland_delta_p_7[i])*Epidemic_Stage_EastofEngland_delta_p_7)*daily_death_EastofEngland)/(P0_EastofEngland-cumul_death_EastofEngland)
  epsilon_EastofEngland_delta_p_14[i,]<-cumsum((1-(gamma_EastofEngland_delta_p_14[i]-(eta_EastofEngland_delta_p_14[i]*gamma_EastofEngland_delta_p_14[i])*Epidemic_Stage_EastofEngland_delta_p_14))/(gamma_EastofEngland_delta_p_14[i]-(eta_EastofEngland_delta_p_7[i]*gamma_EastofEngland_delta_p_14[i])*Epidemic_Stage_EastofEngland_delta_p_14)*daily_death_EastofEngland)/(P0_EastofEngland-cumul_death_EastofEngland)
  epsilon_EastofEngland_delta_p_21[i,]<-cumsum((1-(gamma_EastofEngland_delta_p_21[i]-(eta_EastofEngland_delta_p_21[i]*gamma_EastofEngland_delta_p_21[i])*Epidemic_Stage_EastofEngland_delta_p_21))/(gamma_EastofEngland_delta_p_21[i]-(eta_EastofEngland_delta_p_7[i]*gamma_EastofEngland_delta_p_21[i])*Epidemic_Stage_EastofEngland_delta_p_21)*daily_death_EastofEngland)/(P0_EastofEngland-cumul_death_EastofEngland)

  kft_EastofEngland_delta_p_7[i,]<-gamma_EastofEngland_delta_p_7[i]-(eta_EastofEngland_delta_p_7[i]*gamma_EastofEngland_delta_p_7[i])*Epidemic_Stage_EastofEngland_delta_p_7
  kft_EastofEngland_delta_p_14[i,]<-gamma_EastofEngland_delta_p_14[i]-(eta_EastofEngland_delta_p_14[i]*gamma_EastofEngland_delta_p_14[i])*Epidemic_Stage_EastofEngland_delta_p_14
  kft_EastofEngland_delta_p_21[i,]<-gamma_EastofEngland_delta_p_21[i]-(eta_EastofEngland_delta_p_21[i]*gamma_EastofEngland_delta_p_21[i])*Epidemic_Stage_EastofEngland_delta_p_21
  }

epsilon_EastofEngland_deltap7_deltaeps14<-epsilon_EastofEngland_delta_p_7[,(delta_epsilon_14+1):n_days_EastofEngland]
epsilon_EastofEngland_deltap7_deltaeps21<-epsilon_EastofEngland_delta_p_7[,(delta_epsilon_21+1):n_days_EastofEngland]
epsilon_EastofEngland_deltap7_deltaeps28<-epsilon_EastofEngland_delta_p_7[,(delta_epsilon_28+1):n_days_EastofEngland]
epsilon_EastofEngland_deltap14_deltaeps21<-epsilon_EastofEngland_delta_p_14[,(delta_epsilon_21+1):n_days_EastofEngland]
epsilon_EastofEngland_deltap14_deltaeps28<-epsilon_EastofEngland_delta_p_14[,(delta_epsilon_28+1):n_days_EastofEngland]
epsilon_EastofEngland_deltap21_deltaeps28<-epsilon_EastofEngland_delta_p_21[,(delta_epsilon_28+1):n_days_EastofEngland]

data1EastofEngland = data.frame(output = c(rep("Exposure (Scenario 6)", n_days_London-delta_epsilon_14), 
                                           rep("Exposure (Scenario 5)", n_days_London-delta_epsilon_21),
                                           rep("Exposure (Scenario 4)", n_days_London-delta_epsilon_28), 
                                           rep("Exposure (Scenario 7)", n_days_London-delta_epsilon_21), 
                                           rep("Exposure (Scenario 8)", n_days_London-delta_epsilon_28),
                                           rep("Exposure (Scenario 9)", n_days_London-delta_epsilon_28),
                                           rep("Seroprevalence (Scenario 4, 5, 6)", n_days_London),
                                           rep("Seroprevalence (Scenario 7, 8)", n_days_London), 
                                           rep("Seroprevalence (Scenario 9)", n_days_London)), 
                                t=c(as.Date(EastofEngland_data$Date)[1:(n_days_EastofEngland-delta_epsilon_14)],
                                    as.Date(EastofEngland_data$Date)[1:(n_days_EastofEngland-delta_epsilon_21)],
                                    as.Date(EastofEngland_data$Date)[1:(n_days_EastofEngland-delta_epsilon_28)],
                                    as.Date(EastofEngland_data$Date)[1:(n_days_EastofEngland-delta_epsilon_21)],
                                    as.Date(EastofEngland_data$Date)[1:(n_days_EastofEngland-delta_epsilon_28)],
                                    as.Date(EastofEngland_data$Date)[1:(n_days_EastofEngland-delta_epsilon_28)],
                                    as.Date(EastofEngland_data$Date)[1:n_days_EastofEngland],
                                    as.Date(EastofEngland_data$Date)[1:n_days_EastofEngland],
                                    as.Date(EastofEngland_data$Date)[1:n_days_EastofEngland]), 
                                median = c(100*apply(epsilon_EastofEngland_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.5)), 
                                           100*apply(epsilon_EastofEngland_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                           100*apply(epsilon_EastofEngland_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                           100*apply(epsilon_EastofEngland_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                           100*apply(epsilon_EastofEngland_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                           100*apply(epsilon_EastofEngland_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.5)),
                                           100*apply(x_EastofEngland_delta_p_7, 2, function(x) quantile(x, probs = 0.5)),
                                           100*apply(x_EastofEngland_delta_p_14, 2, function(x) quantile(x, probs = 0.5)), 
                                           100*apply(x_EastofEngland_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                                lower1 = c(100*apply(epsilon_EastofEngland_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.025)), 
                                           100*apply(epsilon_EastofEngland_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                           100*apply(epsilon_EastofEngland_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                           100*apply(epsilon_EastofEngland_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                           100*apply(epsilon_EastofEngland_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                           100*apply(epsilon_EastofEngland_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.025)),
                                           100*apply(x_EastofEngland_delta_p_7, 2, function(x) quantile(x, probs = 0.025)),
                                           100*apply(x_EastofEngland_delta_p_14, 2, function(x) quantile(x, probs = 0.025)), 
                                           100*apply(x_EastofEngland_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                                upper1 = c(100*apply(epsilon_EastofEngland_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.975)), 
                                           100*apply(epsilon_EastofEngland_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                           100*apply(epsilon_EastofEngland_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                           100*apply(epsilon_EastofEngland_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                           100*apply(epsilon_EastofEngland_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                           100*apply(epsilon_EastofEngland_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.975)),
                                           100*apply(x_EastofEngland_delta_p_7, 2, function(x) quantile(x, probs = 0.975)),
                                           100*apply(x_EastofEngland_delta_p_14, 2, function(x) quantile(x, probs = 0.975)), 
                                           100*apply(x_EastofEngland_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                                lower2 = c(100*apply(epsilon_EastofEngland_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.25)), 
                                           100*apply(epsilon_EastofEngland_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                           100*apply(epsilon_EastofEngland_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                           100*apply(epsilon_EastofEngland_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                           100*apply(epsilon_EastofEngland_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                           100*apply(epsilon_EastofEngland_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.25)),
                                           100*apply(x_EastofEngland_delta_p_7, 2, function(x) quantile(x, probs = 0.25)),
                                           100*apply(x_EastofEngland_delta_p_14, 2, function(x) quantile(x, probs = 0.25)), 
                                           100*apply(x_EastofEngland_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                                upper2 = c(100*apply(epsilon_EastofEngland_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.75)), 
                                           100*apply(epsilon_EastofEngland_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                           100*apply(epsilon_EastofEngland_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                           100*apply(epsilon_EastofEngland_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                           100*apply(epsilon_EastofEngland_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                           100*apply(epsilon_EastofEngland_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.75)),
                                           100*apply(x_EastofEngland_delta_p_7, 2, function(x) quantile(x, probs = 0.75)),
                                           100*apply(x_EastofEngland_delta_p_14, 2, function(x) quantile(x, probs = 0.75)), 
                                           100*apply(x_EastofEngland_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))

data2EastofEngland = data.frame( t=as.Date(EastofEngland_data$Date)[1:n_days_EastofEngland][t2_EastofEngland], value= 100*EastofEngland_data$sero[t2_EastofEngland], upper= 100*EastofEngland_data$sero_upper[t2_EastofEngland], lower = 100*EastofEngland_data$sero_lower[t2_EastofEngland])

p1EastofEngland<-ggplot(data1EastofEngland, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("East")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2EastofEngland, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[4])
styled1EastofEngland <- p1EastofEngland +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab("Percentage (%) ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.text.x = element_text(size=font_size),
    axis.text.y = element_text(size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1EastofEngland

##################
### Midlands
##################
sim<-length(beta_delta_p_7)
x_Midlands_delta_p_7<-x_Midlands_delta_p_14<-x_Midlands_delta_p_21<-kft_Midlands<-kft_Midlands_delta_p_7<-kft_Midlands_delta_p_14<-kft_Midlands_delta_p_21<-matrix(0,sim,n_days_Midlands)
epsilon_Midlands_delta_p_7<-epsilon_Midlands_delta_p_14<-epsilon_Midlands_delta_p_21<-kft_Midlands_cert<-matrix(0,sim,n_days_Midlands)

for (i in 1:sim) {
  x_Midlands_delta_p_7[i,]<-rnbinom(rep(1,n_days_Midlands), size= 100, mu=cumsum(exp(beta_delta_p_7[i]*t_Midlands)*(1-(gamma_Midlands_delta_p_7[i]-(eta_Midlands_delta_p_7[i]*gamma_Midlands_delta_p_7[i])*Epidemic_Stage_Midlands_delta_p_7))/(gamma_Midlands_delta_p_7[i]-(eta_Midlands_delta_p_7[i]*gamma_Midlands_delta_p_7[i])*Epidemic_Stage_Midlands_delta_p_7)*daily_death_Midlands)/(exp(beta_delta_p_7[i]*t_Midlands)))/(P0_Midlands-cumul_death_Midlands)
  x_Midlands_delta_p_14[i,]<-rnbinom(rep(1,n_days_Midlands), size= 100, mu=cumsum(exp(beta_delta_p_14[i]*t_Midlands)*(1-(gamma_Midlands_delta_p_14[i]-(eta_Midlands_delta_p_14[i]*gamma_Midlands_delta_p_14[i])*Epidemic_Stage_Midlands_delta_p_14))/(gamma_Midlands_delta_p_14[i]-(eta_Midlands_delta_p_14[i]*gamma_Midlands_delta_p_14[i])*Epidemic_Stage_Midlands_delta_p_14)*daily_death_Midlands)/(exp(beta_delta_p_14[i]*t_Midlands)))/(P0_Midlands-cumul_death_Midlands)
  x_Midlands_delta_p_21[i,]<-rnbinom(rep(1,n_days_Midlands), size= 100, mu=cumsum(exp(beta_delta_p_21[i]*t_Midlands)*(1-(gamma_Midlands_delta_p_21[i]-(eta_Midlands_delta_p_21[i]*gamma_Midlands_delta_p_21[i])*Epidemic_Stage_Midlands_delta_p_21))/(gamma_Midlands_delta_p_21[i]-(eta_Midlands_delta_p_21[i]*gamma_Midlands_delta_p_21[i])*Epidemic_Stage_Midlands_delta_p_21)*daily_death_Midlands)/(exp(beta_delta_p_21[i]*t_Midlands)))/(P0_Midlands-cumul_death_Midlands)
  
  epsilon_Midlands_delta_p_7[i,]<-cumsum((1-(gamma_Midlands_delta_p_7[i]-(eta_Midlands_delta_p_7[i]*gamma_Midlands_delta_p_7[i])*Epidemic_Stage_Midlands_delta_p_7))/(gamma_Midlands_delta_p_7[i]-(eta_Midlands_delta_p_7[i]*gamma_Midlands_delta_p_7[i])*Epidemic_Stage_Midlands_delta_p_7)*daily_death_Midlands)/(P0_Midlands-cumul_death_Midlands)
  epsilon_Midlands_delta_p_14[i,]<-cumsum((1-(gamma_Midlands_delta_p_14[i]-(eta_Midlands_delta_p_14[i]*gamma_Midlands_delta_p_14[i])*Epidemic_Stage_Midlands_delta_p_14))/(gamma_Midlands_delta_p_14[i]-(eta_Midlands_delta_p_7[i]*gamma_Midlands_delta_p_14[i])*Epidemic_Stage_Midlands_delta_p_14)*daily_death_Midlands)/(P0_Midlands-cumul_death_Midlands)
  epsilon_Midlands_delta_p_21[i,]<-cumsum((1-(gamma_Midlands_delta_p_21[i]-(eta_Midlands_delta_p_21[i]*gamma_Midlands_delta_p_21[i])*Epidemic_Stage_Midlands_delta_p_21))/(gamma_Midlands_delta_p_21[i]-(eta_Midlands_delta_p_7[i]*gamma_Midlands_delta_p_21[i])*Epidemic_Stage_Midlands_delta_p_21)*daily_death_Midlands)/(P0_Midlands-cumul_death_Midlands)

  kft_Midlands_delta_p_7[i,]<-gamma_Midlands_delta_p_7[i]-(eta_Midlands_delta_p_7[i]*gamma_Midlands_delta_p_7[i])*Epidemic_Stage_Midlands_delta_p_7
  kft_Midlands_delta_p_14[i,]<-gamma_Midlands_delta_p_14[i]-(eta_Midlands_delta_p_14[i]*gamma_Midlands_delta_p_14[i])*Epidemic_Stage_Midlands_delta_p_14
  kft_Midlands_delta_p_21[i,]<-gamma_Midlands_delta_p_21[i]-(eta_Midlands_delta_p_21[i]*gamma_Midlands_delta_p_21[i])*Epidemic_Stage_Midlands_delta_p_21
}

epsilon_Midlands_deltap7_deltaeps14<-epsilon_Midlands_delta_p_7[,(delta_epsilon_14+1):n_days_Midlands]
epsilon_Midlands_deltap7_deltaeps21<-epsilon_Midlands_delta_p_7[,(delta_epsilon_21+1):n_days_Midlands]
epsilon_Midlands_deltap7_deltaeps28<-epsilon_Midlands_delta_p_7[,(delta_epsilon_28+1):n_days_Midlands]
epsilon_Midlands_deltap14_deltaeps21<-epsilon_Midlands_delta_p_14[,(delta_epsilon_21+1):n_days_Midlands]
epsilon_Midlands_deltap14_deltaeps28<-epsilon_Midlands_delta_p_14[,(delta_epsilon_28+1):n_days_Midlands]
epsilon_Midlands_deltap21_deltaeps28<-epsilon_Midlands_delta_p_21[,(delta_epsilon_28+1):n_days_Midlands]

data1Midlands = data.frame(output = c(rep("Exposure (Scenario 6)", n_days_London-delta_epsilon_14), 
                                      rep("Exposure (Scenario 5)", n_days_London-delta_epsilon_21),
                                      rep("Exposure (Scenario 4)", n_days_London-delta_epsilon_28), 
                                      rep("Exposure (Scenario 7)", n_days_London-delta_epsilon_21), 
                                      rep("Exposure (Scenario 8)", n_days_London-delta_epsilon_28),
                                      rep("Exposure (Scenario 9)", n_days_London-delta_epsilon_28),
                                      rep("Seroprevalence (Scenario 4, 5, 6)", n_days_London),
                                      rep("Seroprevalence (Scenario 7, 8)", n_days_London), 
                                      rep("Seroprevalence (Scenario 9)", n_days_London)), 
                           t=c(as.Date(EastMidlands_data$Date)[1:(n_days_Midlands-delta_epsilon_14)],
                               as.Date(EastMidlands_data$Date)[1:(n_days_Midlands-delta_epsilon_21)],
                               as.Date(EastMidlands_data$Date)[1:(n_days_Midlands-delta_epsilon_28)],
                               as.Date(EastMidlands_data$Date)[1:(n_days_Midlands-delta_epsilon_21)],
                               as.Date(EastMidlands_data$Date)[1:(n_days_Midlands-delta_epsilon_28)],
                               as.Date(EastMidlands_data$Date)[1:(n_days_Midlands-delta_epsilon_28)],
                               as.Date(EastMidlands_data$Date)[1:n_days_Midlands],
                               as.Date(EastMidlands_data$Date)[1:n_days_Midlands],
                               as.Date(EastMidlands_data$Date)[1:n_days_Midlands]), 
                           median = c(100*apply(epsilon_Midlands_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.5)), 
                                      100*apply(epsilon_Midlands_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                      100*apply(epsilon_Midlands_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                      100*apply(epsilon_Midlands_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                      100*apply(epsilon_Midlands_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                      100*apply(epsilon_Midlands_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.5)),
                                      100*apply(x_Midlands_delta_p_7, 2, function(x) quantile(x, probs = 0.5)),
                                      100*apply(x_Midlands_delta_p_14, 2, function(x) quantile(x, probs = 0.5)), 
                                      100*apply(x_Midlands_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                           lower1 = c(100*apply(epsilon_Midlands_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.025)), 
                                      100*apply(epsilon_Midlands_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                      100*apply(epsilon_Midlands_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                      100*apply(epsilon_Midlands_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                      100*apply(epsilon_Midlands_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                      100*apply(epsilon_Midlands_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.025)),
                                      100*apply(x_Midlands_delta_p_7, 2, function(x) quantile(x, probs = 0.025)),
                                      100*apply(x_Midlands_delta_p_14, 2, function(x) quantile(x, probs = 0.025)), 
                                      100*apply(x_Midlands_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                           upper1 = c(100*apply(epsilon_Midlands_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.975)), 
                                      100*apply(epsilon_Midlands_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                      100*apply(epsilon_Midlands_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                      100*apply(epsilon_Midlands_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                      100*apply(epsilon_Midlands_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                      100*apply(epsilon_Midlands_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.975)),
                                      100*apply(x_Midlands_delta_p_7, 2, function(x) quantile(x, probs = 0.975)),
                                      100*apply(x_Midlands_delta_p_14, 2, function(x) quantile(x, probs = 0.975)), 
                                      100*apply(x_Midlands_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                           lower2 = c(100*apply(epsilon_Midlands_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.25)), 
                                      100*apply(epsilon_Midlands_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                      100*apply(epsilon_Midlands_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                      100*apply(epsilon_Midlands_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                      100*apply(epsilon_Midlands_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                      100*apply(epsilon_Midlands_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.25)),
                                      100*apply(x_Midlands_delta_p_7, 2, function(x) quantile(x, probs = 0.25)),
                                      100*apply(x_Midlands_delta_p_14, 2, function(x) quantile(x, probs = 0.25)), 
                                      100*apply(x_Midlands_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                           upper2 = c(100*apply(epsilon_Midlands_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.75)), 
                                      100*apply(epsilon_Midlands_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                      100*apply(epsilon_Midlands_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                      100*apply(epsilon_Midlands_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                      100*apply(epsilon_Midlands_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                      100*apply(epsilon_Midlands_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.75)),
                                      100*apply(x_Midlands_delta_p_7, 2, function(x) quantile(x, probs = 0.75)),
                                      100*apply(x_Midlands_delta_p_14, 2, function(x) quantile(x, probs = 0.75)), 
                                      100*apply(x_Midlands_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))

data2Midlands = data.frame( t=as.Date(EastMidlands_data$Date)[1:n_days_Midlands][t2_Midlands], value= 100*EastMidlands_data$sero[t2_Midlands], upper= 100*EastMidlands_data$sero_upper[t2_Midlands], lower = 100*EastMidlands_data$sero_lower[t2_Midlands])

p1Midlands<-ggplot(data1Midlands, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("Midlands")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2Midlands, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[4])
styled1Midlands <- p1Midlands +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
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
    axis.text.x =element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1Midlands

##################
### NorthEastYorkshireHumber
##################
sim<-length(beta_delta_p_7)
x_NorthEastYorkshireHumber_delta_p_7<-x_NorthEastYorkshireHumber_delta_p_14<-x_NorthEastYorkshireHumber_delta_p_21<-kft_NorthEastYorkshireHumber_delta_p_7<-kft_NorthEastYorkshireHumber_delta_p_14<-kft_NorthEastYorkshireHumber_delta_p_21<-matrix(0,sim,n_days_NorthEastYorkshireHumber)
epsilon_NorthEastYorkshireHumber_delta_p_7<-epsilon_NorthEastYorkshireHumber_delta_p_14<-epsilon_NorthEastYorkshireHumber_delta_p_21<-kft_NorthEastYorkshireHumber_cert<-matrix(0,sim,n_days_NorthEastYorkshireHumber)

for (i in 1:sim) {
  x_NorthEastYorkshireHumber_delta_p_7[i,]<-rnbinom(rep(1,n_days_NorthEastYorkshireHumber), size= 100, mu=cumsum(exp(beta_delta_p_7[i]*t_NorthEastYorkshireHumber)*(1-(gamma_NorthEastYorkshireHumber_delta_p_7[i]-(eta_NorthEastYorkshireHumber_delta_p_7[i]*gamma_NorthEastYorkshireHumber_delta_p_7[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_7))/(gamma_NorthEastYorkshireHumber_delta_p_7[i]-(eta_NorthEastYorkshireHumber_delta_p_7[i]*gamma_NorthEastYorkshireHumber_delta_p_7[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_7)*daily_death_NorthEastYorkshireHumber)/(exp(beta_delta_p_7[i]*t_NorthEastYorkshireHumber)))/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber)
  x_NorthEastYorkshireHumber_delta_p_14[i,]<-rnbinom(rep(1,n_days_NorthEastYorkshireHumber), size= 100, mu=cumsum(exp(beta_delta_p_14[i]*t_NorthEastYorkshireHumber)*(1-(gamma_NorthEastYorkshireHumber_delta_p_14[i]-(eta_NorthEastYorkshireHumber_delta_p_14[i]*gamma_NorthEastYorkshireHumber_delta_p_14[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_14))/(gamma_NorthEastYorkshireHumber_delta_p_14[i]-(eta_NorthEastYorkshireHumber_delta_p_14[i]*gamma_NorthEastYorkshireHumber_delta_p_14[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_14)*daily_death_NorthEastYorkshireHumber)/(exp(beta_delta_p_14[i]*t_NorthEastYorkshireHumber)))/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber)
  x_NorthEastYorkshireHumber_delta_p_21[i,]<-rnbinom(rep(1,n_days_NorthEastYorkshireHumber), size= 100, mu=cumsum(exp(beta_delta_p_21[i]*t_NorthEastYorkshireHumber)*(1-(gamma_NorthEastYorkshireHumber_delta_p_21[i]-(eta_NorthEastYorkshireHumber_delta_p_21[i]*gamma_NorthEastYorkshireHumber_delta_p_21[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_21))/(gamma_NorthEastYorkshireHumber_delta_p_21[i]-(eta_NorthEastYorkshireHumber_delta_p_21[i]*gamma_NorthEastYorkshireHumber_delta_p_21[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_21)*daily_death_NorthEastYorkshireHumber)/(exp(beta_delta_p_21[i]*t_NorthEastYorkshireHumber)))/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber)
  
  epsilon_NorthEastYorkshireHumber_delta_p_7[i,]<-cumsum((1-(gamma_NorthEastYorkshireHumber_delta_p_7[i]-(eta_NorthEastYorkshireHumber_delta_p_7[i]*gamma_NorthEastYorkshireHumber_delta_p_7[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_7))/(gamma_NorthEastYorkshireHumber_delta_p_7[i]-(eta_NorthEastYorkshireHumber_delta_p_7[i]*gamma_NorthEastYorkshireHumber_delta_p_7[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_7)*daily_death_NorthEastYorkshireHumber)/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber)
  epsilon_NorthEastYorkshireHumber_delta_p_14[i,]<-cumsum((1-(gamma_NorthEastYorkshireHumber_delta_p_14[i]-(eta_NorthEastYorkshireHumber_delta_p_14[i]*gamma_NorthEastYorkshireHumber_delta_p_14[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_14))/(gamma_NorthEastYorkshireHumber_delta_p_14[i]-(eta_NorthEastYorkshireHumber_delta_p_7[i]*gamma_NorthEastYorkshireHumber_delta_p_14[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_14)*daily_death_NorthEastYorkshireHumber)/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber)
  epsilon_NorthEastYorkshireHumber_delta_p_21[i,]<-cumsum((1-(gamma_NorthEastYorkshireHumber_delta_p_21[i]-(eta_NorthEastYorkshireHumber_delta_p_21[i]*gamma_NorthEastYorkshireHumber_delta_p_21[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_21))/(gamma_NorthEastYorkshireHumber_delta_p_21[i]-(eta_NorthEastYorkshireHumber_delta_p_7[i]*gamma_NorthEastYorkshireHumber_delta_p_21[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_21)*daily_death_NorthEastYorkshireHumber)/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber)

  kft_NorthEastYorkshireHumber_delta_p_7[i,]<-gamma_NorthEastYorkshireHumber_delta_p_7[i]-(eta_NorthEastYorkshireHumber_delta_p_7[i]*gamma_NorthEastYorkshireHumber_delta_p_7[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_7
  kft_NorthEastYorkshireHumber_delta_p_14[i,]<-gamma_NorthEastYorkshireHumber_delta_p_14[i]-(eta_NorthEastYorkshireHumber_delta_p_14[i]*gamma_NorthEastYorkshireHumber_delta_p_14[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_14
  kft_NorthEastYorkshireHumber_delta_p_21[i,]<-gamma_NorthEastYorkshireHumber_delta_p_21[i]-(eta_NorthEastYorkshireHumber_delta_p_21[i]*gamma_NorthEastYorkshireHumber_delta_p_21[i])*Epidemic_Stage_NorthEastYorkshireHumber_delta_p_21
  }

epsilon_NorthEastYorkshireHumber_deltap7_deltaeps14<-epsilon_NorthEastYorkshireHumber_delta_p_7[,(delta_epsilon_14+1):n_days_NorthEastYorkshireHumber]
epsilon_NorthEastYorkshireHumber_deltap7_deltaeps21<-epsilon_NorthEastYorkshireHumber_delta_p_7[,(delta_epsilon_21+1):n_days_NorthEastYorkshireHumber]
epsilon_NorthEastYorkshireHumber_deltap7_deltaeps28<-epsilon_NorthEastYorkshireHumber_delta_p_7[,(delta_epsilon_28+1):n_days_NorthEastYorkshireHumber]
epsilon_NorthEastYorkshireHumber_deltap14_deltaeps21<-epsilon_NorthEastYorkshireHumber_delta_p_14[,(delta_epsilon_21+1):n_days_NorthEastYorkshireHumber]
epsilon_NorthEastYorkshireHumber_deltap14_deltaeps28<-epsilon_NorthEastYorkshireHumber_delta_p_14[,(delta_epsilon_28+1):n_days_NorthEastYorkshireHumber]
epsilon_NorthEastYorkshireHumber_deltap21_deltaeps28<-epsilon_NorthEastYorkshireHumber_delta_p_21[,(delta_epsilon_28+1):n_days_NorthEastYorkshireHumber]

data1NorthEastYorkshireHumber = data.frame(output = c(rep("Exposure (Scenario 6)", n_days_London-delta_epsilon_14), 
                                                      rep("Exposure (Scenario 5)", n_days_London-delta_epsilon_21),
                                                      rep("Exposure (Scenario 4)", n_days_London-delta_epsilon_28), 
                                                      rep("Exposure (Scenario 7)", n_days_London-delta_epsilon_21), 
                                                      rep("Exposure (Scenario 8)", n_days_London-delta_epsilon_28),
                                                      rep("Exposure (Scenario 9)", n_days_London-delta_epsilon_28),
                                                      rep("Seroprevalence (Scenario 4, 5, 6)", n_days_London),
                                                      rep("Seroprevalence (Scenario 7, 8, 9)", n_days_London), 
                                                      rep("Seroprevalence (Scenario 9)", n_days_London)), 
                                           t=c(as.Date(NorthEast_data$Date)[1:(n_days_NorthEastYorkshireHumber-delta_epsilon_14)],
                                               as.Date(NorthEast_data$Date)[1:(n_days_NorthEastYorkshireHumber-delta_epsilon_21)],
                                               as.Date(NorthEast_data$Date)[1:(n_days_NorthEastYorkshireHumber-delta_epsilon_28)],
                                               as.Date(NorthEast_data$Date)[1:(n_days_NorthEastYorkshireHumber-delta_epsilon_21)],
                                               as.Date(NorthEast_data$Date)[1:(n_days_NorthEastYorkshireHumber-delta_epsilon_28)],
                                               as.Date(NorthEast_data$Date)[1:(n_days_NorthEastYorkshireHumber-delta_epsilon_28)],
                                               as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber],
                                               as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber],
                                               as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber]), 
                                           median = c(100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.5)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.5)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.5)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.5)),
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_7, 2, function(x) quantile(x, probs = 0.5)),
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_14, 2, function(x) quantile(x, probs = 0.5)), 
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                                           lower1 = c(100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.025)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.025)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.025)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.025)),
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_7, 2, function(x) quantile(x, probs = 0.025)),
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_14, 2, function(x) quantile(x, probs = 0.025)), 
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                                           upper1 = c(100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.975)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.975)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.975)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.975)),
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_7, 2, function(x) quantile(x, probs = 0.975)),
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_14, 2, function(x) quantile(x, probs = 0.975)), 
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                                           lower2 = c(100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.25)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.25)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.25)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.25)),
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_7, 2, function(x) quantile(x, probs = 0.25)),
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_14, 2, function(x) quantile(x, probs = 0.25)), 
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                                           upper2 = c(100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps14, 2, function(x) quantile(x, probs = 0.75)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap7_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap14_deltaeps21, 2, function(x) quantile(x, probs = 0.75)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap14_deltaeps28, 2, function(x) quantile(x, probs = 0.75)), 
                                                      100*apply(epsilon_NorthEastYorkshireHumber_deltap21_deltaeps28, 2, function(x) quantile(x, probs = 0.75)),
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_7, 2, function(x) quantile(x, probs = 0.75)),
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_14, 2, function(x) quantile(x, probs = 0.75)), 
                                                      100*apply(x_NorthEastYorkshireHumber_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))

data2NorthEastYorkshireHumber = data.frame( t=as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber][t2_NorthEastYorkshireHumber], value= 100*NorthEast_data$sero[t2_NorthEastYorkshireHumber], upper= 100*NorthEast_data$sero_upper[t2_NorthEastYorkshireHumber], lower = 100*NorthEast_data$sero_lower[t2_NorthEastYorkshireHumber])

p1NorthEastYorkshireHumber<-ggplot(data1NorthEastYorkshireHumber, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("North East")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2NorthEastYorkshireHumber, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[4])
styled1NorthEastYorkshireHumber <- p1NorthEastYorkshireHumber +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
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
    axis.text.x= element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
# styled1NorthEastYorkshireHumber

lwd<-1
IFR_London<-data.frame(output = c(rep("Model 4 & 5 & 6", dim(kft_London_delta_p_7)[2]), 
                                                  rep("Model 7 & 8", dim(kft_London_delta_p_14)[2]),
                                                  rep("Model 9", dim(kft_London_delta_p_21)[2])), 
                                       t=c(as.Date(London_data$Date),
                                           as.Date(London_data$Date),
                                           as.Date(London_data$Date)), 
                                       median = c(apply(kft_London_delta_p_7, 2, function(x) quantile(x, probs = 0.5)), 
                                                  apply(kft_London_delta_p_14, 2, function(x) quantile(x, probs = 0.5)),
                                                  apply(kft_London_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                                       lower1 = c(apply(kft_London_delta_p_7, 2, function(x) quantile(x, probs = 0.025)), 
                                                  apply(kft_London_delta_p_14, 2, function(x) quantile(x, probs = 0.025)),
                                                  apply(kft_London_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                                       upper1 = c(apply(kft_London_delta_p_7, 2, function(x) quantile(x, probs = 0.975)), 
                                                  apply(kft_London_delta_p_14, 2, function(x) quantile(x, probs = 0.975)),
                                                  apply(kft_London_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                                       lower2 = c(apply(kft_London_delta_p_7, 2, function(x) quantile(x, probs = 0.25)), 
                                                  apply(kft_London_delta_p_14, 2, function(x) quantile(x, probs = 0.25)),
                                                  apply(kft_London_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                                       upper2 = c(apply(kft_London_delta_p_7, 2, function(x) quantile(x, probs = 0.75)), 
                                                  apply(kft_London_delta_p_14, 2, function(x) quantile(x, probs = 0.75)),
                                                  apply(kft_London_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))

IFRLondon<-ggplot(IFR_London, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("London")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.3, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0.001,0.002,0.003,0.004,0.005,0.006,0.007), limit = c(0.001, 0.007))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))
styledIFR_London <- IFRLondon +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" IFR ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"))


IFR_NorthWest<-data.frame(output = c(rep("Model 4 & 5 & 6", dim(kft_NorthWest_delta_p_7)[2]), 
                                     rep("Model 7 & 8", dim(kft_NorthWest_delta_p_14)[2]),
                                     rep("Model 9", dim(kft_NorthWest_delta_p_21)[2])), 
                          t=c(as.Date(NorthWest_data$Date),
                              as.Date(NorthWest_data$Date),
                              as.Date(NorthWest_data$Date)), 
                          median = c(apply(kft_NorthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.5)), 
                                     apply(kft_NorthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.5)),
                                     apply(kft_NorthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                          lower1 = c(apply(kft_NorthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.025)), 
                                     apply(kft_NorthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.025)),
                                     apply(kft_NorthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                          upper1 = c(apply(kft_NorthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.975)), 
                                     apply(kft_NorthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.975)),
                                     apply(kft_NorthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                          lower2 = c(apply(kft_NorthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.25)), 
                                     apply(kft_NorthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.25)),
                                     apply(kft_NorthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                          upper2 = c(apply(kft_NorthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.75)), 
                                     apply(kft_NorthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.75)),
                                     apply(kft_NorthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))
IFRNorthWest<-ggplot(IFR_NorthWest, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("North West")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.3, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011), limit = c(0.003, 0.011))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))
styledIFR_NorthWest <- IFRNorthWest +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"))

IFR_SouthWest<-data.frame(output = c(rep("Model 4 & 5 & 6", dim(kft_SouthWest_delta_p_7)[2]), 
                                     rep("Model 7 & 8", dim(kft_SouthWest_delta_p_14)[2]),
                                     rep("Model 9", dim(kft_SouthWest_delta_p_21)[2])), 
                          t=c(as.Date(SouthWest_data$Date),
                              as.Date(SouthWest_data$Date),
                              as.Date(SouthWest_data$Date)), 
                          median = c(apply(kft_SouthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.5)), 
                                     apply(kft_SouthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.5)),
                                     apply(kft_SouthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                          lower1 = c(apply(kft_SouthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.025)), 
                                     apply(kft_SouthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.025)),
                                     apply(kft_SouthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                          upper1 = c(apply(kft_SouthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.975)), 
                                     apply(kft_SouthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.975)),
                                     apply(kft_SouthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                          lower2 = c(apply(kft_SouthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.25)), 
                                     apply(kft_SouthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.25)),
                                     apply(kft_SouthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                          upper2 = c(apply(kft_SouthWest_delta_p_7, 2, function(x) quantile(x, probs = 0.75)), 
                                     apply(kft_SouthWest_delta_p_14, 2, function(x) quantile(x, probs = 0.75)),
                                     apply(kft_SouthWest_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))
IFRSouthWest<-ggplot(IFR_SouthWest, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("South West")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.3, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0.006,0.007,0.008,0.009,0.010,0.011,0.012), limit = c(0.006, 0.012))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))
styledIFR_SouthWest <- IFRSouthWest +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" IFR ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"))

IFR_SouthEast<-data.frame(output = c(rep("Model 4, 5, 6", dim(kft_SouthEast_delta_p_7)[2]), 
                                     rep("Model 7, 8", dim(kft_SouthEast_delta_p_14)[2]),
                                     rep("Model 9", dim(kft_SouthEast_delta_p_21)[2])), 
                          t=c(as.Date(SouthEast_data$Date),
                              as.Date(SouthEast_data$Date),
                              as.Date(SouthEast_data$Date)), 
                          median = c(apply(kft_SouthEast_delta_p_7, 2, function(x) quantile(x, probs = 0.5)), 
                                     apply(kft_SouthEast_delta_p_14, 2, function(x) quantile(x, probs = 0.5)),
                                     apply(kft_SouthEast_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                          lower1 = c(apply(kft_SouthEast_delta_p_7, 2, function(x) quantile(x, probs = 0.025)), 
                                     apply(kft_SouthEast_delta_p_14, 2, function(x) quantile(x, probs = 0.025)),
                                     apply(kft_SouthEast_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                          upper1 = c(apply(kft_SouthEast_delta_p_7, 2, function(x) quantile(x, probs = 0.975)), 
                                     apply(kft_SouthEast_delta_p_14, 2, function(x) quantile(x, probs = 0.975)),
                                     apply(kft_SouthEast_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                          lower2 = c(apply(kft_SouthEast_delta_p_7, 2, function(x) quantile(x, probs = 0.25)), 
                                     apply(kft_SouthEast_delta_p_14, 2, function(x) quantile(x, probs = 0.25)),
                                     apply(kft_SouthEast_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                          upper2 = c(apply(kft_SouthEast_delta_p_7, 2, function(x) quantile(x, probs = 0.75)), 
                                     apply(kft_SouthEast_delta_p_14, 2, function(x) quantile(x, probs = 0.75)),
                                     apply(kft_SouthEast_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))
IFRSouthEast<-ggplot(IFR_SouthEast, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("South East")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.3, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0.005,0.008,0.011,0.014,0.017,0.020,0.023), limit = c(0.005, 0.023))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))
styledIFR_SouthEast <- IFRSouthEast +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(-0.01,-0.2),
    legend.title = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"))

IFR_Midlands<-data.frame(output = c(rep("Model 4 & 5 & 6", dim(kft_Midlands_delta_p_7)[2]), 
                                    rep("Model 7 & 8", dim(kft_Midlands_delta_p_14)[2]),
                                    rep("Model 9", dim(kft_Midlands_delta_p_21)[2])), 
                         t=c(as.Date(EastMidlands_data$Date),
                             as.Date(EastMidlands_data$Date),
                             as.Date(EastMidlands_data$Date)), 
                         median = c(apply(kft_Midlands_delta_p_7, 2, function(x) quantile(x, probs = 0.5)), 
                                    apply(kft_Midlands_delta_p_14, 2, function(x) quantile(x, probs = 0.5)),
                                    apply(kft_Midlands_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                         lower1 = c(apply(kft_Midlands_delta_p_7, 2, function(x) quantile(x, probs = 0.025)), 
                                    apply(kft_Midlands_delta_p_14, 2, function(x) quantile(x, probs = 0.025)),
                                    apply(kft_Midlands_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                         upper1 = c(apply(kft_Midlands_delta_p_7, 2, function(x) quantile(x, probs = 0.975)), 
                                    apply(kft_Midlands_delta_p_14, 2, function(x) quantile(x, probs = 0.975)),
                                    apply(kft_Midlands_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                         lower2 = c(apply(kft_Midlands_delta_p_7, 2, function(x) quantile(x, probs = 0.25)), 
                                    apply(kft_Midlands_delta_p_14, 2, function(x) quantile(x, probs = 0.25)),
                                    apply(kft_Midlands_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                         upper2 = c(apply(kft_Midlands_delta_p_7, 2, function(x) quantile(x, probs = 0.75)), 
                                    apply(kft_Midlands_delta_p_14, 2, function(x) quantile(x, probs = 0.75)),
                                    apply(kft_Midlands_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))
IFRMidlands<-ggplot(IFR_Midlands, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("Midlands")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.3, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0.005,0.006,0.007,0.008,0.009,0.010,0.011), limit = c(0.005, 0.011))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))
styledIFR_Midlands <- IFRMidlands +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"))

IFR_EastofEngland<-data.frame(output = c(rep("Model 4 & 5 & 6", dim(kft_EastofEngland_delta_p_7)[2]), 
                                         rep("Model 7 & 8", dim(kft_EastofEngland_delta_p_14)[2]),
                                         rep("Model 9", dim(kft_EastofEngland_delta_p_21)[2])), 
                              t=c(as.Date(EastofEngland_data$Date),
                                  as.Date(EastofEngland_data$Date),
                                  as.Date(EastofEngland_data$Date)), 
                              median = c(apply(kft_EastofEngland_delta_p_7, 2, function(x) quantile(x, probs = 0.5)), 
                                         apply(kft_EastofEngland_delta_p_14, 2, function(x) quantile(x, probs = 0.5)),
                                         apply(kft_EastofEngland_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                              lower1 = c(apply(kft_EastofEngland_delta_p_7, 2, function(x) quantile(x, probs = 0.025)), 
                                         apply(kft_EastofEngland_delta_p_14, 2, function(x) quantile(x, probs = 0.025)),
                                         apply(kft_EastofEngland_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                              upper1 = c(apply(kft_EastofEngland_delta_p_7, 2, function(x) quantile(x, probs = 0.975)), 
                                         apply(kft_EastofEngland_delta_p_14, 2, function(x) quantile(x, probs = 0.975)),
                                         apply(kft_EastofEngland_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                              lower2 = c(apply(kft_EastofEngland_delta_p_7, 2, function(x) quantile(x, probs = 0.25)), 
                                         apply(kft_EastofEngland_delta_p_14, 2, function(x) quantile(x, probs = 0.25)),
                                         apply(kft_EastofEngland_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                              upper2 = c(apply(kft_EastofEngland_delta_p_7, 2, function(x) quantile(x, probs = 0.75)), 
                                         apply(kft_EastofEngland_delta_p_14, 2, function(x) quantile(x, probs = 0.75)),
                                         apply(kft_EastofEngland_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))
IFREastofEngland<-ggplot(IFR_EastofEngland, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("East")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.3, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0.003,0.005,0.007,0.009,0.011,0.013,0.015), limit = c(0.003, 0.015))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))
styledIFR_EastofEngland <- IFREastofEngland +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab("IFR ") +
  xlab("2020")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"))

IFR_NorthEastYorkshireHumber<-data.frame(output = c(rep("Model 4 & 5 & 6", dim(kft_NorthEastYorkshireHumber_delta_p_7)[2]), 
                                                    rep("Model 7 & 8", dim(kft_NorthEastYorkshireHumber_delta_p_14)[2]),
                                                    rep("Model 9", dim(kft_NorthEastYorkshireHumber_delta_p_21)[2])), 
                                         t=c(as.Date(NorthEast_data$Date),
                                             as.Date(NorthEast_data$Date),
                                             as.Date(NorthEast_data$Date)), 
                                         median = c(apply(kft_NorthEastYorkshireHumber_delta_p_7, 2, function(x) quantile(x, probs = 0.5)), 
                                                    apply(kft_NorthEastYorkshireHumber_delta_p_14, 2, function(x) quantile(x, probs = 0.5)),
                                                    apply(kft_NorthEastYorkshireHumber_delta_p_21, 2, function(x) quantile(x, probs = 0.5))), 
                                         lower1 = c(apply(kft_NorthEastYorkshireHumber_delta_p_7, 2, function(x) quantile(x, probs = 0.025)), 
                                                    apply(kft_NorthEastYorkshireHumber_delta_p_14, 2, function(x) quantile(x, probs = 0.025)),
                                                    apply(kft_NorthEastYorkshireHumber_delta_p_21, 2, function(x) quantile(x, probs = 0.025))), 
                                         upper1 = c(apply(kft_NorthEastYorkshireHumber_delta_p_7, 2, function(x) quantile(x, probs = 0.975)), 
                                                    apply(kft_NorthEastYorkshireHumber_delta_p_14, 2, function(x) quantile(x, probs = 0.975)),
                                                    apply(kft_NorthEastYorkshireHumber_delta_p_21, 2, function(x) quantile(x, probs = 0.975))),
                                         lower2 = c(apply(kft_NorthEastYorkshireHumber_delta_p_7, 2, function(x) quantile(x, probs = 0.25)), 
                                                    apply(kft_NorthEastYorkshireHumber_delta_p_14, 2, function(x) quantile(x, probs = 0.25)),
                                                    apply(kft_NorthEastYorkshireHumber_delta_p_21, 2, function(x) quantile(x, probs = 0.25))), 
                                         upper2 = c(apply(kft_NorthEastYorkshireHumber_delta_p_7, 2, function(x) quantile(x, probs = 0.75)), 
                                                    apply(kft_NorthEastYorkshireHumber_delta_p_14, 2, function(x) quantile(x, probs = 0.75)),
                                                    apply(kft_NorthEastYorkshireHumber_delta_p_21, 2, function(x) quantile(x, probs = 0.75))))
IFRNorthEastYorkshireHumber<-ggplot(IFR_NorthEastYorkshireHumber, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  ggtitle("North East")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.3, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0.006,0.007,0.008,0.009,0.010,0.011,0.012,0.013), limit = c(0.006, 0.013))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))
styledIFR_NorthEastYorkshireHumber <- IFRNorthEastYorkshireHumber +
  # scale_fill_brewer(palette = "Dark2")+
  # scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ylab(" ") +
  xlab("2020")+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=14),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"))

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Results"
folder = paste(folder_strings, sep = "", collapse = "/")

tiff(file=paste(folder,"/TimevaryingIFRModel_paraEsti_deltap.tiff", sep = ""),
     width=27, height=17, units="cm", res=300)
ggarrange(p2,p3, p4,p5, p6,p7, p8, p1,p9,p10, p11,p12, p13,p14, p15)
dev.off()


tiff(file=paste(folder,"/TimevaryingIFRModel_fitting_deltap.tiff", sep = ""),
     width=27, height=17, units="cm", res=300)
grid.arrange(styled1London,styled1NorthEastYorkshireHumber,styled1NorthWest,styled1SouthWest,styled1SouthEast,styled1Midlands,styled1EastofEngland)
dev.off()

tiff(file=paste(folder,"/IFR_delta_p.tiff", sep = ""),
     width=30, height=16, units="cm", res=300)
grid.arrange(styledIFR_London,styledIFR_NorthEastYorkshireHumber,styledIFR_NorthWest,styledIFR_SouthWest,styledIFR_SouthEast,styledIFR_Midlands,styledIFR_EastofEngland)
dev.off()
