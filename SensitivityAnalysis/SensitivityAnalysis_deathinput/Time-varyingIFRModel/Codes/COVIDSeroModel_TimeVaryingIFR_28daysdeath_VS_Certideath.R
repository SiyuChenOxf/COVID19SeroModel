#This code is to genertate comparison results using 28 days positive death and certificate death as model inputs
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

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

SeroModelTimeVaryingIFR <- readRDS(paste(folder,"/SeroModelTimeVaryingIFR_28daysdeath.rds", sep = ""))
CertificateSeroModelTimeVaryingIFR <- readRDS(paste(folder,"/SeroModelTimeVaryingIFR_delta_p_14.rds", sep = ""))

#Posterior estimations for parameters using 28 days positive death as model input
beta <- rstan::extract(SeroModelTimeVaryingIFR)$beta
eta_NorthEastYorkshireHumber<- rstan::extract(SeroModelTimeVaryingIFR)$eta_NorthEastYorkshireHumber
gamma_NorthEastYorkshireHumber<- rstan::extract(SeroModelTimeVaryingIFR)$gamma_NorthEastYorkshireHumber
eta_London <- rstan::extract(SeroModelTimeVaryingIFR)$eta_London
gamma_London <- rstan::extract(SeroModelTimeVaryingIFR)$gamma_London
eta_NorthWest <- rstan::extract(SeroModelTimeVaryingIFR)$eta_NorthWest
gamma_NorthWest <- rstan::extract(SeroModelTimeVaryingIFR)$gamma_NorthWest
eta_SouthEast <- rstan::extract(SeroModelTimeVaryingIFR)$eta_SouthEast
gamma_SouthEast <- rstan::extract(SeroModelTimeVaryingIFR)$gamma_SouthEast
eta_SouthWest <- rstan::extract(SeroModelTimeVaryingIFR)$eta_SouthWest
gamma_SouthWest <- rstan::extract(SeroModelTimeVaryingIFR)$gamma_SouthWest
eta_Midlands <- rstan::extract(SeroModelTimeVaryingIFR)$eta_Midlands
gamma_Midlands <- rstan::extract(SeroModelTimeVaryingIFR)$gamma_Midlands
eta_EastofEngland <- rstan::extract(SeroModelTimeVaryingIFR)$eta_EastofEngland
gamma_EastofEngland <- rstan::extract(SeroModelTimeVaryingIFR)$gamma_EastofEngland

#Posterior estimations for parameters using certificate death as model input
cert_beta <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$beta
cert_eta_NorthEastYorkshireHumber<- rstan::extract(CertificateSeroModelTimeVaryingIFR)$eta_NorthEastYorkshireHumber
cert_gamma_NorthEastYorkshireHumber<- rstan::extract(CertificateSeroModelTimeVaryingIFR)$gamma_NorthEastYorkshireHumber
cert_eta_London <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$eta_London
cert_gamma_London <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$gamma_London
cert_eta_NorthWest <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$eta_NorthWest
cert_gamma_NorthWest <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$gamma_NorthWest
cert_eta_SouthEast <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$eta_SouthEast
cert_gamma_SouthEast <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$gamma_SouthEast
cert_eta_SouthWest <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$eta_SouthWest
cert_gamma_SouthWest <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$gamma_SouthWest
cert_eta_Midlands <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$eta_Midlands
cert_gamma_Midlands <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$gamma_Midlands
cert_eta_EastofEngland <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$eta_EastofEngland
cert_gamma_EastofEngland <- rstan::extract(CertificateSeroModelTimeVaryingIFR)$gamma_EastofEngland

data<-data.frame(class=factor(rep(c("Model 11","Model 7"),each=length(beta))),
                 para_beta=c(beta,cert_beta),
                 para_London=c(gamma_London,cert_gamma_London),para_London2=c(eta_London,cert_eta_London),
                 para_NorthEast=c(gamma_NorthEastYorkshireHumber,cert_gamma_NorthEastYorkshireHumber),para_NorthEast2=c(eta_NorthEastYorkshireHumber,cert_eta_NorthEastYorkshireHumber),
                 para_SouthEast=c(gamma_SouthEast,cert_gamma_SouthEast),para_SouthEast2=c(eta_SouthEast,cert_eta_SouthEast),
                 para_NorthWest=c(gamma_NorthWest,cert_gamma_NorthWest),para_NorthWest2=c(eta_NorthWest,cert_eta_NorthWest),
                 para_SouthWest=c(gamma_SouthWest,cert_gamma_SouthWest),para_SouthWest2=c(eta_SouthWest,cert_eta_SouthWest),
                 para_Midlands=c(gamma_Midlands,cert_gamma_Midlands),para_Midlands2=c(eta_Midlands,cert_eta_Midlands),
                 para_EastofEngland=c(gamma_EastofEngland,cert_gamma_EastofEngland),para_EastofEngland2=c(eta_EastofEngland,cert_eta_EastofEngland))
them<-theme(
  text = element_text(size=font_size),
  plot.title = element_text(face = "bold", size = font_size_title),
  legend.background = element_rect(fill = "white", size = 0.5, colour = "white"),
  legend.justification = c(0, 1),
  legend.position = "none",
  legend.title = element_blank(),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank(),
  axis.ticks = element_line(colour = "grey50", size = 0.2),
  # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
  # panel.grid.major = element_line(colour = "grey50"),
  # panel.grid.minor = element_blank()
)
p1<-ggplot(data, aes(x=para_beta, fill=class))    +geom_density(alpha=.3)+xlab(" ")+ylab(expression(beta))+theme_minimal()+them+theme(legend.position = c(1,1))
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
p15<-ggplot(data, aes(x=para_EastofEngland2, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(eta[EastofEngland]))+theme_minimal()+them+theme(legend.position = c(1,1))
###

London_data28<-read.csv(paste(folder,"/28London_data.csv", sep = ""),header = TRUE)
NorthWest_data28<-read.csv(paste(folder,"/28NorthWest_data.csv", sep = ""),header = TRUE)
YorkshireHumber_data28<-read.csv(paste(folder,"/28YorkshireHumber_data.csv", sep = ""),header = TRUE)
SouthWest_data28<-read.csv(paste(folder,"/28SouthWest_data.csv", sep = ""),header = TRUE)
EastofEngland_data28<-read.csv(paste(folder,"/28EastofEngland_data.csv", sep = ""),header = TRUE)
EastMidlands_data28<-read.csv(paste(folder,"/28EastMidlands_data.csv", sep = ""),header = TRUE)
WestMidlands_data28<-read.csv(paste(folder,"/28WestMidlands_data.csv", sep = ""),header = TRUE)
NorthEast_data28<-read.csv(paste(folder,"/28NorthEast_data.csv", sep = ""),header = TRUE)
SouthEast_data28<-read.csv(paste(folder,"/28SouthEast_data.csv", sep = ""),header = TRUE)

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
delta_epsilon<-21  #fixed time lag between exposure and death(seroconversion)

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

## Load virus positivity data from GOV.UK##
positivity<-read.csv(paste(folder,"/pcr_positivity.csv", sep = ""),header = TRUE)

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

Epidemic_Stage_London<-c(rep(0,positivity_start_London+delta_p-1),positivity_London[positivity_start_London:London_positivity_end])[1:London_death_end]
Epidemic_Stage_London<-cumsum(Epidemic_Stage_London)/sum(Epidemic_Stage_London)

#Epidemic stage in NorthEast & Yorkshire and the Humber 
positivity_start_NorthEastYorkshireHumber<-which(positivity_NorthEastYorshireHumber!=0)[1] 
NorthEastYorkshireHumber_positivity_end<-which(positivity_NorthEastYorshireHumber!=0)[length(which(positivity_NorthEastYorshireHumber!=0))] 
NorthEast_death_end<-which(!is.na(NorthEast_data$daily_death_NorthEast))[length(which(!is.na(NorthEast_data$daily_death)))]  
YorkshireHumber_death_end<-which(!is.na(YorkshireHumber_data$daily_death_Yorkshire_Humber))[length(which(!is.na(YorkshireHumber_data$daily_death_Yorkshire_Humber)))] 

Epidemic_Stage_NorthEastYorkshireHumber<-c(rep(0,positivity_start_NorthEastYorkshireHumber+delta_p-1),positivity_NorthEastYorshireHumber[positivity_start_NorthEastYorkshireHumber:NorthEastYorkshireHumber_positivity_end])[1:NorthEast_death_end]
Epidemic_Stage_NorthEastYorkshireHumber<-cumsum(Epidemic_Stage_NorthEastYorkshireHumber)/sum(Epidemic_Stage_NorthEastYorkshireHumber)

#Epidemic stage in SouthEast
positivity_start_SouthEast<-which(positivity_SouthEast!=0)[1]#62
SouthEast_positivity_end<-which(positivity_SouthEast!=0)[length(which(positivity_SouthEast!=0))]#306
SouthEast_death_end<-which(!is.na(SouthEast_data$daily_death))[length(which(!is.na(SouthEast_data$daily_death)))]#312

Epidemic_Stage_SouthEast<-c(rep(0,positivity_start_SouthEast+delta_p-1),positivity_SouthEast[positivity_start_SouthEast:SouthEast_positivity_end])[1:SouthEast_death_end]
Epidemic_Stage_SouthEast<-cumsum(Epidemic_Stage_SouthEast)/sum(Epidemic_Stage_SouthEast)

#Epidemic stage in SoutWest
positivity_start_SouthWest<-which(positivity_SouthWest!=0)[1]#62
SouthWest_positivity_end<-which(positivity_SouthWest!=0)[length(which(positivity_SouthWest!=0))]#306
SouthWest_death_end<-which(!is.na(SouthWest_data$daily_death))[length(which(!is.na(SouthWest_data$daily_death)))]#312

Epidemic_Stage_SouthWest<-c(rep(0,positivity_start_SouthWest+delta_p-1),positivity_SouthWest[positivity_start_SouthWest:SouthWest_positivity_end])[1:SouthWest_death_end]
Epidemic_Stage_SouthWest<-cumsum(Epidemic_Stage_SouthWest)/sum(Epidemic_Stage_SouthWest)

#Epidemic stage in NorthWest
positivity_start_NorthWest<-which(positivity_NorthWest!=0)[1]#62
NorthWest_positivity_end<-which(positivity_NorthWest!=0)[length(which(positivity_NorthWest!=0))]#306
NorthWest_death_end<-which(!is.na(NorthWest_data$daily_death))[length(which(!is.na(NorthWest_data$daily_death)))]#312

Epidemic_Stage_NorthWest<-c(rep(0,positivity_start_NorthWest+delta_p-1),positivity_NorthWest[positivity_start_NorthWest:NorthWest_positivity_end])[1:NorthWest_death_end]
Epidemic_Stage_NorthWest<-cumsum(Epidemic_Stage_NorthWest)/sum(Epidemic_Stage_NorthWest)

#Epidemic stage in Midlands
positivity_start_Midlands<-which(positivity_Midlands!=0)[1]#62
Midlands_positivity_end<-which(positivity_Midlands!=0)[length(which(positivity_Midlands!=0))]#306
EastMidlands_death_end<-which(!is.na(EastMidlands_data$daily_death_EastMidlands))[length(which(!is.na(EastMidlands_data$daily_death_EastMidlands)))]#312
WestMidlands_death_end<-which(!is.na(WestMidlands_data$daily_death_WestMidlands))[length(which(!is.na(WestMidlands_data$daily_death_WestMidlands)))]

Epidemic_Stage_Midlands<-c(rep(0,positivity_start_Midlands+delta_p-1),positivity_Midlands[positivity_start_Midlands:Midlands_positivity_end])[1:EastMidlands_death_end]
Epidemic_Stage_Midlands<-cumsum(Epidemic_Stage_Midlands)/sum(Epidemic_Stage_Midlands)

#Epidemic stage in EastofEngland
positivity_start_EastofEngland<-which(positivity_EastofEngland!=0)[1]#62
EastofEngland_positivity_end<-which(positivity_EastofEngland!=0)[length(which(positivity_EastofEngland!=0))]#306
EastofEngland_death_end<-which(!is.na(EastofEngland_data$daily_death))[length(which(!is.na(EastofEngland_data$daily_death)))]#312

Epidemic_Stage_EastofEngland<-c(rep(0,positivity_start_EastofEngland+delta_p-1),positivity_EastofEngland[positivity_start_EastofEngland:EastofEngland_positivity_end])[1:EastofEngland_death_end]
Epidemic_Stage_EastofEngland<-cumsum(Epidemic_Stage_EastofEngland)/sum(Epidemic_Stage_EastofEngland)

#London
London_death_end<-which(!is.na(London_data$daily_death))[length(which(!is.na(London_data$daily_death)))] 
daily_death_London<-London_data$daily_death[1:London_death_end]          #Number of daily death in London
daily_death_London28<-London_data28$daily_death[1:London_death_end]          #Number of daily death in London

sero_London<-round(P0_London*London_data$sero[!is.na(London_data$sero)]) #Number of seropositive individuals at all serological survey time points in London 
n_days_London<-length(daily_death_London)                                #Number of days for daily death in London      
cumul_death_London<-cumsum(daily_death_London)                           #Number of cumulative death in London
cumul_death_London28<-cumsum(daily_death_London28)                           #Number of cumulative death in London

n_days2_London<-length(sero_London)                                      #Number of days for serological survey in London
t_London<-seq(1,n_days_London,by=1)                                       
t2_London<-t_London[!is.na(London_data$sero)]                  

#NorthEast & Yorkshire and the Humber 
NorthEast_death_end<-which(!is.na(NorthEast_data$daily_death_NorthEast))[length(which(!is.na(NorthEast_data$daily_death)))]  
YorkshireHumber_death_end<-which(!is.na(YorkshireHumber_data$daily_death_Yorkshire_Humber))[length(which(!is.na(YorkshireHumber_data$daily_death_Yorkshire_Humber)))] 
daily_death_NorthEastYorkshireHumber<-NorthEast_data$daily_death_NorthEast[1:NorthEast_death_end]+YorkshireHumber_data$daily_death_Yorkshire_Humber[1:YorkshireHumber_death_end]
daily_death_NorthEastYorkshireHumber28<-NorthEast_data28$daily_death_NorthEast[1:NorthEast_death_end]+YorkshireHumber_data28$daily_death_Yorkshire_Humber[1:YorkshireHumber_death_end]

sero_NorthEastYorkshireHumber<-round(P0_NorthEastYorkshireHumber*NorthEast_data$sero[!is.na(NorthEast_data$sero)])
n_days_NorthEastYorkshireHumber<-length(daily_death_NorthEastYorkshireHumber)
cumul_death_NorthEastYorkshireHumber<-cumsum(daily_death_NorthEastYorkshireHumber)
cumul_death_NorthEastYorkshireHumber28<-cumsum(daily_death_NorthEastYorkshireHumber28)

n_days2_NorthEastYorkshireHumber<-length(sero_NorthEastYorkshireHumber)
t_NorthEastYorkshireHumber<-seq(1,n_days_NorthEastYorkshireHumber,by=1)
t2_NorthEastYorkshireHumber<-t_NorthEastYorkshireHumber[!is.na(NorthEast_data$sero)]

#SouthWest
SouthWest_death_end<-which(!is.na(SouthWest_data$daily_death))[length(which(!is.na(SouthWest_data$daily_death)))]#312
daily_death_SouthWest<-SouthWest_data$daily_death[1:SouthWest_death_end]
daily_death_SouthWest28<-SouthWest_data28$daily_death[1:SouthWest_death_end]

sero_SouthWest<-round(P0_SouthWest*SouthWest_data$sero[!is.na(SouthWest_data$sero)])
n_days_SouthWest<-length(daily_death_SouthWest)
cumul_death_SouthWest<-cumsum(daily_death_SouthWest)
cumul_death_SouthWest28<-cumsum(daily_death_SouthWest28)

n_days2_SouthWest<-length(sero_SouthWest)
t_SouthWest<-seq(1,n_days_SouthWest,by=1)
t2_SouthWest<-t_SouthWest[!is.na(SouthWest_data$sero)]

#NorthWest
NorthWest_death_end<-which(!is.na(NorthWest_data$daily_death))[length(which(!is.na(NorthWest_data$daily_death)))]#312
daily_death_NorthWest<-NorthWest_data$daily_death[1:NorthWest_death_end]
daily_death_NorthWest28<-NorthWest_data28$daily_death[1:NorthWest_death_end]

sero_NorthWest<-round(P0_NorthWest*NorthWest_data$sero[!is.na(NorthWest_data$sero)])
n_days_NorthWest<-length(daily_death_NorthWest)
cumul_death_NorthWest<-cumsum(daily_death_NorthWest)
cumul_death_NorthWest28<-cumsum(daily_death_NorthWest28)

n_days2_NorthWest<-length(sero_NorthWest)
t_NorthWest<-seq(1,n_days_NorthWest,by=1)
t2_NorthWest<-t_NorthWest[!is.na(NorthWest_data$sero)]

#SouthEast
SouthEast_death_end<-which(!is.na(SouthEast_data$daily_death))[length(which(!is.na(SouthEast_data$daily_death)))]#312
daily_death_SouthEast<-SouthEast_data$daily_death[1:SouthEast_death_end]
daily_death_SouthEast28<-SouthEast_data28$daily_death[1:SouthEast_death_end]

sero_SouthEast<-round(P0_SouthEast*SouthEast_data$sero[!is.na(SouthEast_data$sero)])
n_days_SouthEast<-length(daily_death_SouthEast)
cumul_death_SouthEast<-cumsum(daily_death_SouthEast)
cumul_death_SouthEast28<-cumsum(daily_death_SouthEast28)

n_days2_SouthEast<-length(sero_SouthEast)
t_SouthEast<-seq(1,n_days_SouthEast,by=1)
t2_SouthEast<-t_SouthEast[!is.na(SouthEast_data$sero)]

#Midlands
EastMidlands_death_end<-which(!is.na(EastMidlands_data$daily_death_EastMidlands))[length(which(!is.na(EastMidlands_data$daily_death_EastMidlands)))]#312
WestMidlands_death_end<-which(!is.na(WestMidlands_data$daily_death_WestMidlands))[length(which(!is.na(WestMidlands_data$daily_death_WestMidlands)))]
daily_death_Midlands<-WestMidlands_data$daily_death_WestMidlands[1:WestMidlands_death_end]+EastMidlands_data$daily_death_EastMidlands[1:EastMidlands_death_end]
daily_death_Midlands28<-WestMidlands_data28$daily_death_WestMidlands[1:WestMidlands_death_end]+EastMidlands_data28$daily_death_EastMidlands[1:EastMidlands_death_end]

sero_Midlands<-round(P0_Midlands*EastMidlands_data$sero[!is.na(EastMidlands_data$sero)])
n_days_Midlands<-length(daily_death_Midlands)
cumul_death_Midlands<-cumsum(daily_death_Midlands)
cumul_death_Midlands28<-cumsum(daily_death_Midlands28)

n_days2_Midlands<-length(sero_Midlands)
t_Midlands<-seq(1,n_days_Midlands,by=1)
t2_Midlands<-t_Midlands[!is.na(EastMidlands_data$sero)]

#EastofEngland
EastofEngland_death_end<-which(!is.na(EastofEngland_data$daily_death))[length(which(!is.na(EastofEngland_data$daily_death)))]#312
daily_death_EastofEngland<-EastofEngland_data$daily_death[1:EastofEngland_death_end]
daily_death_EastofEngland28<-EastofEngland_data28$daily_death[1:EastofEngland_death_end]

sero_EastofEngland<-round(P0_EastofEngland*EastofEngland_data$sero[!is.na(EastofEngland_data$sero)])
n_days_EastofEngland<-length(daily_death_EastofEngland)
cumul_death_EastofEngland<-cumsum(daily_death_EastofEngland)
cumul_death_EastofEngland28<-cumsum(daily_death_EastofEngland28)

n_days2_EastofEngland<-length(sero_EastofEngland)
t_EastofEngland<-seq(1,n_days_EastofEngland,by=1)
t2_EastofEngland<-t_EastofEngland[!is.na(EastofEngland_data$sero)]

##################
### LONDON
##################
sim<-length(beta)
x_London<-epsilon_London<-kft_London<-matrix(0,sim,n_days_London)
x_London_cert<-epsilon_London_cert<-kft_London_cert<-matrix(0,sim,n_days_London)

for (i in 1:sim) {
  x_London[i,]<-rnbinom(rep(1,n_days_London), size= 100, mu=cumsum(exp(beta[i]*t_London)*(1-(gamma_London[i]-(eta_London[i]*gamma_London[i])*Epidemic_Stage_London))/(gamma_London[i]-(eta_London[i]*gamma_London[i])*Epidemic_Stage_London)*daily_death_London28)/(exp(beta[i]*t_London)))/(P0_London-cumul_death_London28)
  x_London_cert[i,]<-rnbinom(rep(1,n_days_London), size= 100, mu=cumsum(exp(cert_beta[i]*t_London)*(1-(cert_gamma_London[i]-(cert_eta_London[i]*cert_gamma_London[i])*Epidemic_Stage_London))/(cert_gamma_London[i]-(cert_eta_London[i]*cert_gamma_London[i])*Epidemic_Stage_London)*daily_death_London)/(exp(cert_beta[i]*t_London)))/(P0_London-cumul_death_London)
  
  epsilon_London[i,]<-cumsum((1-(gamma_London[i]-(eta_London[i]*gamma_London[i])*Epidemic_Stage_London))/(gamma_London[i]-(eta_London[i]*gamma_London[i])*Epidemic_Stage_London)*daily_death_London28)/(P0_London-cumul_death_London28)
  epsilon_London_cert[i,]<-cumsum((1-(cert_gamma_London[i]-(cert_eta_London[i]*cert_gamma_London[i])*Epidemic_Stage_London))/(cert_gamma_London[i]-(cert_eta_London[i]*cert_gamma_London[i])*Epidemic_Stage_London)*daily_death_London)/(P0_London-cumul_death_London)
}

epsilon_London<-epsilon_London[,(delta_epsilon+1):n_days_London]
epsilon_London_cert<-epsilon_London_cert[,(delta_epsilon+1):n_days_London]

data1London = data.frame(output = c(rep("Exposure_Scenario 11", n_days_London-delta_epsilon), rep("Seroprevalence_Scenario 11", n_days_London),rep("Exposure_Scenario 7", n_days_London-delta_epsilon), rep("Seroprevalence_Scenario 7", n_days_London)), 
                         t=c(as.Date(London_data$Date)[1:(n_days_London-delta_epsilon)],as.Date(London_data$Date)[1:n_days_London],as.Date(London_data$Date)[1:(n_days_London-delta_epsilon)],as.Date(London_data$Date)[1:n_days_London]), 
                         median = c(100*apply(epsilon_London, 2, function(x) quantile(x, probs = 0.5)), 
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.5)),
                                    100*apply(epsilon_London_cert, 2, function(x) quantile(x, probs = 0.5)), 
                                    100*apply(x_London_cert, 2, function(x) quantile(x, probs = 0.5))), 
                         lower1 = c(100*apply(epsilon_London, 2, function(x) quantile(x, probs = 0.025)), 
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.025)),
                                    100*apply(epsilon_London_cert, 2, function(x) quantile(x, probs = 0.025)), 
                                    100*apply(x_London_cert, 2, function(x) quantile(x, probs = 0.025))), 
                         upper1 = c(100*apply(epsilon_London, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(epsilon_London_cert, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(x_London_cert, 2, function(x) quantile(x, probs = 0.975))),
                         lower2 = c(100*apply(epsilon_London, 2, function(x) quantile(x, probs = 0.25)), 
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.25)),
                                    100*apply(epsilon_London_cert, 2, function(x) quantile(x, probs = 0.25)), 
                                    100*apply(x_London_cert, 2, function(x) quantile(x, probs = 0.25))), 
                         upper2 = c(100*apply(epsilon_London, 2, function(x) quantile(x, probs = 0.75)),
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.75)),
                                    100*apply(epsilon_London_cert, 2, function(x) quantile(x, probs = 0.75)),
                                    100*apply(x_London_cert, 2, function(x) quantile(x, probs = 0.75))))

data2London = data.frame( t=as.Date(London_data$Date)[1:n_days_London][t2_London], value= 100*London_data$sero[t2_London], upper= 100*London_data$sero_upper[t2_London], lower = 100*London_data$sero_lower[t2_London])

p1London<-ggplot(data1London, aes(x=t, y = median, group = output, colour = output)) +
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
### North West
##################
sim<-length(beta)
x_NorthWest<-epsilon_NorthWest<-kft_NorthWest<-matrix(0,sim,n_days_NorthWest)
x_NorthWest_cert<-epsilon_NorthWest_cert<-kft_NorthWest_cert<-matrix(0,sim,n_days_NorthWest)

for (i in 1:sim) {
  x_NorthWest[i,]<-rnbinom(rep(1,n_days_NorthWest), size= 100, mu=cumsum(exp(beta[i]*t_NorthWest)*(1-(gamma_NorthWest[i]-(eta_NorthWest[i]*gamma_NorthWest[i])*Epidemic_Stage_NorthWest))/(gamma_NorthWest[i]-(eta_NorthWest[i]*gamma_NorthWest[i])*Epidemic_Stage_NorthWest)*daily_death_NorthWest28)/(exp(beta[i]*t_NorthWest)))/(P0_NorthWest-cumul_death_NorthWest28)
  x_NorthWest_cert[i,]<-rnbinom(rep(1,n_days_NorthWest), size= 100, mu=cumsum(exp(cert_beta[i]*t_NorthWest)*(1-(cert_gamma_NorthWest[i]-(cert_eta_NorthWest[i]*cert_gamma_NorthWest[i])*Epidemic_Stage_NorthWest))/(cert_gamma_NorthWest[i]-(cert_eta_NorthWest[i]*cert_gamma_NorthWest[i])*Epidemic_Stage_NorthWest)*daily_death_NorthWest)/(exp(cert_beta[i]*t_NorthWest)))/(P0_NorthWest-cumul_death_NorthWest)
  
  epsilon_NorthWest[i,]<-cumsum((1-(gamma_NorthWest[i]-(eta_NorthWest[i]*gamma_NorthWest[i])*Epidemic_Stage_NorthWest))/(gamma_NorthWest[i]-(eta_NorthWest[i]*gamma_NorthWest[i])*Epidemic_Stage_NorthWest)*daily_death_NorthWest28)/(P0_NorthWest-cumul_death_NorthWest28)
  epsilon_NorthWest_cert[i,]<-cumsum((1-(cert_gamma_NorthWest[i]-(cert_eta_NorthWest[i]*cert_gamma_NorthWest[i])*Epidemic_Stage_NorthWest))/(cert_gamma_NorthWest[i]-(cert_eta_NorthWest[i]*cert_gamma_NorthWest[i])*Epidemic_Stage_NorthWest)*daily_death_NorthWest)/(P0_NorthWest-cumul_death_NorthWest)
}
epsilon_NorthWest<-epsilon_NorthWest[,(delta_epsilon+1):n_days_NorthWest]
epsilon_NorthWest_cert<-epsilon_NorthWest_cert[,(delta_epsilon+1):n_days_NorthWest]


data1NorthWest = data.frame(output = c(rep("Exposure_Scenario 11", n_days_NorthWest-delta_epsilon), rep("Seroprevalence_Scenario 11", n_days_NorthWest),rep("Exposure_Scenario 7", n_days_NorthWest-delta_epsilon), rep("Seroprevalence_Scenario 7", n_days_NorthWest)), 
                            t=c(as.Date(NorthWest_data$Date)[1:(n_days_NorthWest-delta_epsilon)],as.Date(NorthWest_data$Date)[1:n_days_NorthWest],as.Date(NorthWest_data$Date)[1:(n_days_NorthWest-delta_epsilon)],as.Date(NorthWest_data$Date)[1:n_days_NorthWest]), 
                            median = c(100*apply(epsilon_NorthWest, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_NorthWest_cert, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_NorthWest_cert, 2, function(x) quantile(x, probs = 0.5))), 
                            lower1 = c(100*apply(epsilon_NorthWest, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_NorthWest_cert, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(x_NorthWest_cert, 2, function(x) quantile(x, probs = 0.025))), 
                            upper1 = c(100*apply(epsilon_NorthWest, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_NorthWest_cert, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_NorthWest_cert, 2, function(x) quantile(x, probs = 0.975))),
                            lower2 = c(100*apply(epsilon_NorthWest, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_NorthWest_cert, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_NorthWest_cert, 2, function(x) quantile(x, probs = 0.25))), 
                            upper2 = c(100*apply(epsilon_NorthWest, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_NorthWest_cert, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_NorthWest_cert, 2, function(x) quantile(x, probs = 0.75))))

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
    legend.position = "none",
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

##################
### South West
##################
sim<-length(beta)
x_SouthWest<-epsilon_SouthWest<-kft_SouthWest<-matrix(0,sim,n_days_SouthWest)
x_SouthWest_cert<-epsilon_SouthWest_cert<-kft_SouthWest_cert<-matrix(0,sim,n_days_SouthWest)

for (i in 1:sim) {
  x_SouthWest[i,]<-rnbinom(rep(1,n_days_SouthWest), size= 100, mu=cumsum(exp(beta[i]*t_SouthWest)*(1-(gamma_SouthWest[i]-(eta_SouthWest[i]*gamma_SouthWest[i])*Epidemic_Stage_SouthWest))/(gamma_SouthWest[i]-(eta_SouthWest[i]*gamma_SouthWest[i])*Epidemic_Stage_SouthWest)*daily_death_SouthWest28)/(exp(beta[i]*t_SouthWest)))/(P0_SouthWest-cumul_death_SouthWest28)
  x_SouthWest_cert[i,]<-rnbinom(rep(1,n_days_SouthWest), size= 100, mu=cumsum(exp(cert_beta[i]*t_SouthWest)*(1-(cert_gamma_SouthWest[i]-(cert_eta_SouthWest[i]*cert_gamma_SouthWest[i])*Epidemic_Stage_SouthWest))/(cert_gamma_SouthWest[i]-(cert_eta_SouthWest[i]*cert_gamma_SouthWest[i])*Epidemic_Stage_SouthWest)*daily_death_SouthWest)/(exp(cert_beta[i]*t_SouthWest)))/(P0_SouthWest-cumul_death_SouthWest)
  
  epsilon_SouthWest[i,]<-cumsum((1-(gamma_SouthWest[i]-(eta_SouthWest[i]*gamma_SouthWest[i])*Epidemic_Stage_SouthWest))/(gamma_SouthWest[i]-(eta_SouthWest[i]*gamma_SouthWest[i])*Epidemic_Stage_SouthWest)*daily_death_SouthWest28)/(P0_SouthWest-cumul_death_SouthWest28)
  epsilon_SouthWest_cert[i,]<-cumsum((1-(cert_gamma_SouthWest[i]-(cert_eta_SouthWest[i]*cert_gamma_SouthWest[i])*Epidemic_Stage_SouthWest))/(cert_gamma_SouthWest[i]-(cert_eta_SouthWest[i]*cert_gamma_SouthWest[i])*Epidemic_Stage_SouthWest)*daily_death_SouthWest)/(P0_SouthWest-cumul_death_SouthWest)
  
}
epsilon_SouthWest<-epsilon_SouthWest[,(delta_epsilon+1):n_days_SouthWest]
epsilon_SouthWest_cert<-epsilon_SouthWest_cert[,(delta_epsilon+1):n_days_SouthWest]


data1SouthWest = data.frame(output = c(rep("Exposure_Scenario 11", n_days_SouthWest-delta_epsilon), rep("Seroprevalence_Scenario 11", n_days_SouthWest),rep("Exposure_Scenario 7", n_days_SouthWest-delta_epsilon), rep("Seroprevalence_Scenario 7", n_days_SouthWest)), 
                            t=c(as.Date(SouthWest_data$Date)[1:(n_days_SouthWest-delta_epsilon)],as.Date(SouthWest_data$Date)[1:n_days_SouthWest],as.Date(SouthWest_data$Date)[1:(n_days_SouthWest-delta_epsilon)],as.Date(SouthWest_data$Date)[1:n_days_SouthWest]), 
                            median = c(100*apply(epsilon_SouthWest, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_SouthWest_cert, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_SouthWest_cert, 2, function(x) quantile(x, probs = 0.5))), 
                            lower1 = c(100*apply(epsilon_SouthWest, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_SouthWest_cert, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(x_SouthWest_cert, 2, function(x) quantile(x, probs = 0.025))), 
                            upper1 = c(100*apply(epsilon_SouthWest, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_SouthWest_cert, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_SouthWest_cert, 2, function(x) quantile(x, probs = 0.975))),
                            lower2 = c(100*apply(epsilon_SouthWest, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_SouthWest_cert, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_SouthWest_cert, 2, function(x) quantile(x, probs = 0.25))), 
                            upper2 = c(100*apply(epsilon_SouthWest, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_SouthWest_cert, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_SouthWest_cert, 2, function(x) quantile(x, probs = 0.75))))

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
    axis.text.y = element_text(size=font_size),
    axis.text.x = element_text(size=font_size),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

##################
### South East
##################
sim<-length(beta)
x_SouthEast<-epsilon_SouthEast<-kft_SouthEast<-matrix(0,sim,n_days_SouthEast)
x_SouthEast_cert<-epsilon_SouthEast_cert<-kft_SouthEast_cert<-matrix(0,sim,n_days_SouthEast)

for (i in 1:sim) {
  x_SouthEast[i,]<-rnbinom(rep(1,n_days_SouthEast), size= 100, mu=cumsum(exp(beta[i]*t_SouthEast)*(1-(gamma_SouthEast[i]-(eta_SouthEast[i]*gamma_SouthEast[i])*Epidemic_Stage_SouthEast))/(gamma_SouthEast[i]-(eta_SouthEast[i]*gamma_SouthEast[i])*Epidemic_Stage_SouthEast)*daily_death_SouthEast28)/(exp(beta[i]*t_SouthEast)))/(P0_SouthEast-cumul_death_SouthEast28)
  x_SouthEast_cert[i,]<-rnbinom(rep(1,n_days_SouthEast), size= 100, mu=cumsum(exp(cert_beta[i]*t_SouthEast)*(1-(cert_gamma_SouthEast[i]-(cert_eta_SouthEast[i]*cert_gamma_SouthEast[i])*Epidemic_Stage_SouthEast))/(cert_gamma_SouthEast[i]-(cert_eta_SouthEast[i]*cert_gamma_SouthEast[i])*Epidemic_Stage_SouthEast)*daily_death_SouthEast)/(exp(cert_beta[i]*t_SouthEast)))/(P0_SouthEast-cumul_death_SouthEast)
  
  epsilon_SouthEast[i,]<-cumsum((1-(gamma_SouthEast[i]-(eta_SouthEast[i]*gamma_SouthEast[i])*Epidemic_Stage_SouthEast))/(gamma_SouthEast[i]-(eta_SouthEast[i]*gamma_SouthEast[i])*Epidemic_Stage_SouthEast)*daily_death_SouthEast28)/(P0_SouthEast-cumul_death_SouthEast28)
  epsilon_SouthEast_cert[i,]<-cumsum((1-(cert_gamma_SouthEast[i]-(cert_eta_SouthEast[i]*cert_gamma_SouthEast[i])*Epidemic_Stage_SouthEast))/(cert_gamma_SouthEast[i]-(cert_eta_SouthEast[i]*cert_gamma_SouthEast[i])*Epidemic_Stage_SouthEast)*daily_death_SouthEast)/(P0_SouthEast-cumul_death_SouthEast)
  
}
epsilon_SouthEast<-epsilon_SouthEast[,(delta_epsilon+1):n_days_SouthEast]
epsilon_SouthEast_cert<-epsilon_SouthEast_cert[,(delta_epsilon+1):n_days_SouthEast]


data1SouthEast = data.frame(output = c(rep("Exposure (Model 11)", n_days_SouthEast-delta_epsilon), rep("Seroprevalence (Model 11)", n_days_SouthEast),rep("Exposure (Model 7)", n_days_SouthEast-delta_epsilon), rep("Seroprevalence (Model 7)", n_days_SouthEast)), 
                            t=c(as.Date(SouthEast_data$Date)[1:(n_days_SouthEast-delta_epsilon)],as.Date(SouthEast_data$Date)[1:n_days_SouthEast],as.Date(SouthEast_data$Date)[1:(n_days_SouthEast-delta_epsilon)],as.Date(SouthEast_data$Date)[1:n_days_SouthEast]), 
                            median = c(100*apply(epsilon_SouthEast, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.5)),
                                       100*apply(epsilon_SouthEast_cert, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_SouthEast_cert, 2, function(x) quantile(x, probs = 0.5))), 
                            lower1 = c(100*apply(epsilon_SouthEast, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.025)),
                                       100*apply(epsilon_SouthEast_cert, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(x_SouthEast_cert, 2, function(x) quantile(x, probs = 0.025))), 
                            upper1 = c(100*apply(epsilon_SouthEast, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(epsilon_SouthEast_cert, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_SouthEast_cert, 2, function(x) quantile(x, probs = 0.975))),
                            lower2 = c(100*apply(epsilon_SouthEast, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.25)),
                                       100*apply(epsilon_SouthEast_cert, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_SouthEast_cert, 2, function(x) quantile(x, probs = 0.25))), 
                            upper2 = c(100*apply(epsilon_SouthEast, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(epsilon_SouthEast_cert, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_SouthEast_cert, 2, function(x) quantile(x, probs = 0.75))))

data2SouthEast = data.frame( t=as.Date(SouthEast_data$Date)[1:n_days_SouthEast][t2_SouthEast], value= 100*SouthEast_data$sero[t2_SouthEast], upper= 100*SouthEast_data$sero_upper[t2_SouthEast], lower = 100*SouthEast_data$sero_lower[t2_SouthEast])

p1SouthEast<-ggplot(data1SouthEast, aes(x=t, y = median, group = output, colour = output)) +
  scale_fill_manual(values=c(colors_Dark[1], colors_Dark[4], colors_Dark[3],colors_Dark[2]))+
  geom_line(size = lwd) +  ggtitle("South East")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure),labels = scales::percent_format(accuracy = 1))+
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
    legend.position = c(-0.05,-0.2),
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

##################
### East of England
##################

sim<-length(beta)
x_EastofEngland<-epsilon_EastofEngland<-kft_EastofEngland<-matrix(0,sim,n_days_EastofEngland)
x_EastofEngland_cert<-epsilon_EastofEngland_cert<-kft_EastofEngland_cert<-matrix(0,sim,n_days_EastofEngland)

for (i in 1:sim) {
  x_EastofEngland[i,]<-rnbinom(rep(1,n_days_EastofEngland), size= 100, mu=cumsum(exp(beta[i]*t_EastofEngland)*(1-(gamma_EastofEngland[i]-(eta_EastofEngland[i]*gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland))/(gamma_EastofEngland[i]-(eta_EastofEngland[i]*gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland)*daily_death_EastofEngland28)/(exp(beta[i]*t_EastofEngland)))/(P0_EastofEngland-cumul_death_EastofEngland28)
  x_EastofEngland_cert[i,]<-rnbinom(rep(1,n_days_EastofEngland), size= 100, mu=cumsum(exp(cert_beta[i]*t_EastofEngland)*(1-(cert_gamma_EastofEngland[i]-(cert_eta_EastofEngland[i]*cert_gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland))/(cert_gamma_EastofEngland[i]-(cert_eta_EastofEngland[i]*cert_gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland)*daily_death_EastofEngland)/(exp(cert_beta[i]*t_EastofEngland)))/(P0_EastofEngland-cumul_death_EastofEngland)
  
  epsilon_EastofEngland[i,]<-cumsum((1-(gamma_EastofEngland[i]-(eta_EastofEngland[i]*gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland))/(gamma_EastofEngland[i]-(eta_EastofEngland[i]*gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland)*daily_death_EastofEngland28)/(P0_EastofEngland-cumul_death_EastofEngland28)
  epsilon_EastofEngland_cert[i,]<-cumsum((1-(cert_gamma_EastofEngland[i]-(cert_eta_EastofEngland[i]*cert_gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland))/(cert_gamma_EastofEngland[i]-(cert_eta_EastofEngland[i]*cert_gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland)*daily_death_EastofEngland)/(P0_EastofEngland-cumul_death_EastofEngland)
  
}
epsilon_EastofEngland<-epsilon_EastofEngland[,(delta_epsilon+1):n_days_EastofEngland]
epsilon_EastofEngland_cert<-epsilon_EastofEngland_cert[,(delta_epsilon+1):n_days_EastofEngland]


data1EastofEngland = data.frame(output = c(rep("Exposure_Scenario 11", n_days_EastofEngland-delta_epsilon), rep("Seroprevalence_Scenario 11", n_days_EastofEngland),rep("Exposure_Scenario 7", n_days_EastofEngland-delta_epsilon), rep("Seroprevalence_Scenario 7", n_days_EastofEngland)), 
                                t=c(as.Date(EastofEngland_data$Date)[1:(n_days_EastofEngland-delta_epsilon)],as.Date(EastofEngland_data$Date)[1:n_days_EastofEngland],as.Date(EastofEngland_data$Date)[1:(n_days_EastofEngland-delta_epsilon)],as.Date(EastofEngland_data$Date)[1:n_days_EastofEngland]), 
                                median = c(100*apply(epsilon_EastofEngland, 2, function(x) quantile(x, probs = 0.5)), 
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.5)),
                                           100*apply(epsilon_EastofEngland_cert, 2, function(x) quantile(x, probs = 0.5)), 
                                           100*apply(x_EastofEngland_cert, 2, function(x) quantile(x, probs = 0.5))), 
                                lower1 = c(100*apply(epsilon_EastofEngland, 2, function(x) quantile(x, probs = 0.025)), 
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.025)),
                                           100*apply(epsilon_EastofEngland_cert, 2, function(x) quantile(x, probs = 0.025)), 
                                           100*apply(x_EastofEngland_cert, 2, function(x) quantile(x, probs = 0.025))), 
                                upper1 = c(100*apply(epsilon_EastofEngland, 2, function(x) quantile(x, probs = 0.975)),
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.975)),
                                           100*apply(epsilon_EastofEngland_cert, 2, function(x) quantile(x, probs = 0.975)),
                                           100*apply(x_EastofEngland_cert, 2, function(x) quantile(x, probs = 0.975))),
                                lower2 = c(100*apply(epsilon_EastofEngland, 2, function(x) quantile(x, probs = 0.25)), 
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.25)),
                                           100*apply(epsilon_EastofEngland_cert, 2, function(x) quantile(x, probs = 0.25)), 
                                           100*apply(x_EastofEngland_cert, 2, function(x) quantile(x, probs = 0.25))), 
                                upper2 = c(100*apply(epsilon_EastofEngland, 2, function(x) quantile(x, probs = 0.75)),
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.75)),
                                           100*apply(epsilon_EastofEngland_cert, 2, function(x) quantile(x, probs = 0.75)),
                                           100*apply(x_EastofEngland_cert, 2, function(x) quantile(x, probs = 0.75))))

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
    axis.text.x = element_text(size = font_size),
    axis.text.y = element_text(size = font_size),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

##################
### Midlands
##################
sim<-length(beta)
x_Midlands<-epsilon_Midlands<-kft_Midlands<-matrix(0,sim,n_days_Midlands)
x_Midlands_cert<-epsilon_Midlands_cert<-kft_Midlands_cert<-matrix(0,sim,n_days_Midlands)

for (i in 1:sim) {
  x_Midlands[i,]<-rnbinom(rep(1,n_days_Midlands), size= 100, mu=cumsum(exp(beta[i]*t_Midlands)*(1-(gamma_Midlands[i]-(eta_Midlands[i]*gamma_Midlands[i])*Epidemic_Stage_Midlands))/(gamma_Midlands[i]-(eta_Midlands[i]*gamma_Midlands[i])*Epidemic_Stage_Midlands)*daily_death_Midlands28)/(exp(beta[i]*t_Midlands)))/(P0_Midlands-cumul_death_Midlands28)
  x_Midlands_cert[i,]<-rnbinom(rep(1,n_days_Midlands), size= 100, mu=cumsum(exp(cert_beta[i]*t_Midlands)*(1-(cert_gamma_Midlands[i]-(cert_eta_Midlands[i]*cert_gamma_Midlands[i])*Epidemic_Stage_Midlands))/(cert_gamma_Midlands[i]-(cert_eta_Midlands[i]*cert_gamma_Midlands[i])*Epidemic_Stage_Midlands)*daily_death_Midlands)/(exp(cert_beta[i]*t_Midlands)))/(P0_Midlands-cumul_death_Midlands)
  
  epsilon_Midlands[i,]<-cumsum((1-(gamma_Midlands[i]-(eta_Midlands[i]*gamma_Midlands[i])*Epidemic_Stage_Midlands))/(gamma_Midlands[i]-(eta_Midlands[i]*gamma_Midlands[i])*Epidemic_Stage_Midlands)*daily_death_Midlands28)/(P0_Midlands-cumul_death_Midlands28)
  epsilon_Midlands_cert[i,]<-cumsum((1-(cert_gamma_Midlands[i]-(cert_eta_Midlands[i]*cert_gamma_Midlands[i])*Epidemic_Stage_Midlands))/(cert_gamma_Midlands[i]-(cert_eta_Midlands[i]*cert_gamma_Midlands[i])*Epidemic_Stage_Midlands)*daily_death_Midlands)/(P0_Midlands-cumul_death_Midlands)
}
epsilon_Midlands<-epsilon_Midlands[,(delta_epsilon+1):n_days_Midlands]
epsilon_Midlands_cert<-epsilon_Midlands_cert[,(delta_epsilon+1):n_days_Midlands]


data1Midlands = data.frame(output = c(rep("Exposure_Scenario 11", n_days_Midlands-delta_epsilon), rep("Seroprevalence_Scenario 11", n_days_Midlands),rep("Exposure_Scenario 1", n_days_Midlands-delta_epsilon), rep("Seroprevalence_Scenario 1", n_days_Midlands)), 
                           t=c(as.Date(EastMidlands_data$Date)[1:(n_days_Midlands-delta_epsilon)],as.Date(EastMidlands_data$Date)[1:n_days_Midlands],as.Date(EastMidlands_data$Date)[1:(n_days_Midlands-delta_epsilon)],as.Date(EastMidlands_data$Date)[1:n_days_Midlands]), 
                           median = c(100*apply(epsilon_Midlands, 2, function(x) quantile(x, probs = 0.5)), 
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.5)),
                                      100*apply(epsilon_Midlands_cert, 2, function(x) quantile(x, probs = 0.5)), 
                                      100*apply(x_Midlands_cert, 2, function(x) quantile(x, probs = 0.5))), 
                           lower1 = c(100*apply(epsilon_Midlands, 2, function(x) quantile(x, probs = 0.025)), 
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.025)),
                                      100*apply(epsilon_Midlands_cert, 2, function(x) quantile(x, probs = 0.025)), 
                                      100*apply(x_Midlands_cert, 2, function(x) quantile(x, probs = 0.025))), 
                           upper1 = c(100*apply(epsilon_Midlands, 2, function(x) quantile(x, probs = 0.975)),
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.975)),
                                      100*apply(epsilon_Midlands_cert, 2, function(x) quantile(x, probs = 0.975)),
                                      100*apply(x_Midlands_cert, 2, function(x) quantile(x, probs = 0.975))),
                           lower2 = c(100*apply(epsilon_Midlands, 2, function(x) quantile(x, probs = 0.25)), 
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.25)),
                                      100*apply(epsilon_Midlands_cert, 2, function(x) quantile(x, probs = 0.25)), 
                                      100*apply(x_Midlands_cert, 2, function(x) quantile(x, probs = 0.25))), 
                           upper2 = c(apply(epsilon_Midlands, 2, function(x) quantile(x, probs = 0.75)),
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.75)),
                                      100*apply(epsilon_Midlands_cert, 2, function(x) quantile(x, probs = 0.75)),
                                      100*apply(x_Midlands_cert, 2, function(x) quantile(x, probs = 0.75))))

data2Midlands = data.frame( t=as.Date(EastMidlands_data$Date)[1:n_days_Midlands][t2_Midlands], value= 100*EastMidlands_data$sero[t2_Midlands], upper= 100*EastMidlands_data$sero_upper[t2_Midlands], lower = 100*EastMidlands_data$sero_lower[t2_Midlands])

p1Midlands<-ggplot(data1Midlands, aes(x=t, y = median, group = output, colour = output)) +
  scale_fill_manual(values=c(colors_Dark[1], colors_Dark[4], colors_Dark[3],colors_Dark[2]))+
  geom_line(size = lwd) +  ggtitle("Midlands")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure),labels = scales::percent_format(accuracy = 1))+
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


##################
### NorthEastYorkshireHumber
##################
sim<-length(beta)
x_NorthEastYorkshireHumber<-epsilon_NorthEastYorkshireHumber<-kft_NorthEastYorkshireHumber<-matrix(0,sim,n_days_NorthEastYorkshireHumber)
x_NorthEastYorkshireHumber_cert<-epsilon_NorthEastYorkshireHumber_cert<-kft_NorthEastYorkshireHumber_cert<-matrix(0,sim,n_days_NorthEastYorkshireHumber)

for (i in 1:sim) {
  x_NorthEastYorkshireHumber[i,]<-rnbinom(rep(1,n_days_NorthEastYorkshireHumber), size= 100, mu=cumsum(exp(beta[i]*t_NorthEastYorkshireHumber)*(1-(gamma_NorthEastYorkshireHumber[i]-(eta_NorthEastYorkshireHumber[i]*gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber))/(gamma_NorthEastYorkshireHumber[i]-(eta_NorthEastYorkshireHumber[i]*gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber)*daily_death_NorthEastYorkshireHumber28)/(exp(beta[i]*t_NorthEastYorkshireHumber)))/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber28)
  x_NorthEastYorkshireHumber_cert[i,]<-rnbinom(rep(1,n_days_NorthEastYorkshireHumber), size= 100, mu=cumsum(exp(cert_beta[i]*t_NorthEastYorkshireHumber)*(1-(cert_gamma_NorthEastYorkshireHumber[i]-(cert_eta_NorthEastYorkshireHumber[i]*cert_gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber))/(cert_gamma_NorthEastYorkshireHumber[i]-(cert_eta_NorthEastYorkshireHumber[i]*cert_gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber)*daily_death_NorthEastYorkshireHumber)/(exp(cert_beta[i]*t_NorthEastYorkshireHumber)))/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber)
  
  epsilon_NorthEastYorkshireHumber[i,]<-cumsum((1-(gamma_NorthEastYorkshireHumber[i]-(eta_NorthEastYorkshireHumber[i]*gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber))/(gamma_NorthEastYorkshireHumber[i]-(eta_NorthEastYorkshireHumber[i]*gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber)*daily_death_NorthEastYorkshireHumber28)/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber28)
  epsilon_NorthEastYorkshireHumber_cert[i,]<-cumsum((1-(cert_gamma_NorthEastYorkshireHumber[i]-(cert_eta_NorthEastYorkshireHumber[i]*cert_gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber))/(cert_gamma_NorthEastYorkshireHumber[i]-(cert_eta_NorthEastYorkshireHumber[i]*cert_gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber)*daily_death_NorthEastYorkshireHumber)/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber)
  
}
epsilon_NorthEastYorkshireHumber<-epsilon_NorthEastYorkshireHumber[,(delta_epsilon+1):n_days_NorthEastYorkshireHumber]
epsilon_NorthEastYorkshireHumber_cert<-epsilon_NorthEastYorkshireHumber_cert[,(delta_epsilon+1):n_days_NorthEastYorkshireHumber]


data1NorthEastYorkshireHumber = data.frame(output = c(rep("Exposure_Scenario 11", n_days_NorthEastYorkshireHumber-delta_epsilon), rep("Seroprevalence_Scenario 11", n_days_NorthEastYorkshireHumber),rep("Exposure_Scenario 7", n_days_NorthEastYorkshireHumber-delta_epsilon), rep("Seroprevalence_Scenario 7", n_days_NorthEastYorkshireHumber)), 
                                           t=c(as.Date(NorthEast_data$Date)[1:(n_days_NorthEastYorkshireHumber-delta_epsilon)],as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber],as.Date(NorthEast_data$Date)[1:(n_days_NorthEastYorkshireHumber-delta_epsilon)],as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber]), 
                                           median = c(100*apply(epsilon_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.5)), 
                                                      100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.5)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_cert, 2, function(x) quantile(x, probs = 0.5)), 
                                                      100*apply(x_NorthEastYorkshireHumber_cert, 2, function(x) quantile(x, probs = 0.5))), 
                                           lower1 = c(100*apply(epsilon_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.025)), 
                                                      100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.025)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_cert, 2, function(x) quantile(x, probs = 0.025)), 
                                                      100*apply(x_NorthEastYorkshireHumber_cert, 2, function(x) quantile(x, probs = 0.025))), 
                                           upper1 = c(100*apply(epsilon_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.975)),
                                                      100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.975)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_cert, 2, function(x) quantile(x, probs = 0.975)),
                                                      100*apply(x_NorthEastYorkshireHumber_cert, 2, function(x) quantile(x, probs = 0.975))),
                                           lower2 = c(100*apply(epsilon_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.25)), 
                                                      100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.25)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_cert, 2, function(x) quantile(x, probs = 0.25)), 
                                                      100*apply(x_NorthEastYorkshireHumber_cert, 2, function(x) quantile(x, probs = 0.25))), 
                                           upper2 = c(100*apply(epsilon_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.75)),
                                                      100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.75)),
                                                      100*apply(epsilon_NorthEastYorkshireHumber_cert, 2, function(x) quantile(x, probs = 0.75)),
                                                      100*apply(x_NorthEastYorkshireHumber_cert, 2, function(x) quantile(x, probs = 0.75))))

data2NorthEastYorkshireHumber = data.frame( t=as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber][t2_NorthEastYorkshireHumber], value= 100*NorthEast_data$sero[t2_NorthEastYorkshireHumber], upper= 100*NorthEast_data$sero_upper[t2_NorthEastYorkshireHumber], lower = 100*NorthEast_data$sero_lower[t2_NorthEastYorkshireHumber])

p1NorthEastYorkshireHumber<-ggplot(data1NorthEastYorkshireHumber, aes(x=t, y = median, group = output, colour = output)) +
  scale_fill_manual(values=c(colors_Dark[1], colors_Dark[4], colors_Dark[3],colors_Dark[2]))+
  geom_line(size = lwd) +  ggtitle("North East")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure),labels = scales::percent_format(accuracy = 1))+
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
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=font_size),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x=element_text(size=font_size),
    legend.title = element_blank(),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Results"
folder = paste(folder_strings, sep = "", collapse = "/")

tiff(file=paste(folder,"/ParameterEsti_TimeVarying_deathinput1.tiff", sep = ""),
     width=27, height=17, units="cm", res=300)
ggarrange(p2,p3, p4,p5, p6,p7, p8, p1,p9,p10, p11,p12, p13,p14, p15)
dev.off()

tiff(file=paste(folder,"/Fitting_TimeVarying_deathinput.tiff", sep = ""),
     width=30, height=16, units="cm", res=300)
grid.arrange(styled1London,styled1NorthEastYorkshireHumber,styled1NorthWest,styled1SouthWest,styled1SouthEast,styled1Midlands,styled1EastofEngland)
dev.off()
