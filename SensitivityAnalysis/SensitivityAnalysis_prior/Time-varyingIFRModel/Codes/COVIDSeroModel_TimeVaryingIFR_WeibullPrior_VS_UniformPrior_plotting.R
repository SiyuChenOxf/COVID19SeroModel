#This code is to plot comparison of posterior estimations for time-varying IFR model using Weibull prior and Uniform prior for beta

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

## common values for plotting
ymax_kft <- 0.012
ymax_exposure <- 0.30
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

SeroModelTimeVaryingIFR_Uniform <- readRDS(paste(folder,"/SeroModelTimeVaryingIFR_delta_p_14.rds", sep = ""))          #Load posterior estimations for time-varying IFR model using Uniform prior for beta
SeroModelTimeVaryingIFR_Weibull <- readRDS(paste(folder,"/TimevaryingSeroModelTimeVaryingIFR_Weibull.rds", sep = ""))  #Load posterior estimations for time-varying IFR model using Weibull prior for beta

#Posterior estimations for time-varying IFR model using Uniform prior for beta
beta_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$beta
eta_NorthEastYorkshireHumber_Uniform<- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$eta_NorthEastYorkshireHumber
gamma_NorthEastYorkshireHumber_Uniform<- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$gamma_NorthEastYorkshireHumber
eta_London_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$eta_London
gamma_London_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$gamma_London
eta_NorthWest_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$eta_NorthWest
gamma_NorthWest_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$gamma_NorthWest
eta_SouthEast_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$eta_SouthEast
gamma_SouthEast_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$gamma_SouthEast
eta_SouthWest_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$eta_SouthWest
gamma_SouthWest_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$gamma_SouthWest
eta_Midlands_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$eta_Midlands
gamma_Midlands_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$gamma_Midlands
eta_EastofEngland_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$eta_EastofEngland
gamma_EastofEngland_Uniform <- rstan::extract(SeroModelTimeVaryingIFR_Uniform)$gamma_EastofEngland

#Posterior estimations for time-varying IFR model using Weibull prior for beta
beta_Weibull <- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$beta
eta_NorthEastYorkshireHumber_Weibull<- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$eta_NorthEastYorkshireHumber
gamma_NorthEastYorkshireHumber_Weibull<- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$gamma_NorthEastYorkshireHumber
eta_London_Weibull <- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$eta_London
gamma_London_Weibull<- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$gamma_London
eta_NorthWest_Weibull <- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$eta_NorthWest
gamma_NorthWest_Weibull <- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$gamma_NorthWest
eta_SouthEast_Weibull <- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$eta_SouthEast
gamma_SouthEast_Weibull <- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$gamma_SouthEast
eta_SouthWest_Weibull <- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$eta_SouthWest
gamma_SouthWest_Weibull <- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$gamma_SouthWest
eta_Midlands_Weibull <- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$eta_Midlands
gamma_Midlands_Weibull <- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$gamma_Midlands
eta_EastofEngland_Weibull <- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$eta_EastofEngland
gamma_EastofEngland_Weibull <- rstan::extract(SeroModelTimeVaryingIFR_Weibull)$gamma_EastofEngland

data<-data.frame(class=factor(rep(c("Model 7","Model 13"),each=length(beta_Uniform))),
                 para_beta=c(beta_Uniform,beta_Weibull),
                 para_London=c(gamma_London_Uniform,gamma_London_Weibull),para_London2=c(eta_London_Uniform,eta_London_Weibull),
                 para_NorthEast=c(gamma_NorthEastYorkshireHumber_Uniform,gamma_NorthEastYorkshireHumber_Weibull),para_NorthEast2=c(eta_NorthEastYorkshireHumber_Uniform,eta_NorthEastYorkshireHumber_Weibull),
                 para_SouthEast=c(gamma_SouthEast_Uniform,gamma_SouthEast_Weibull),para_SouthEast2=c(eta_SouthEast_Uniform,eta_SouthEast_Weibull),
                 para_NorthWest=c(gamma_NorthWest_Uniform,gamma_NorthWest_Weibull),para_NorthWest2=c(eta_NorthWest_Uniform,eta_NorthWest_Weibull),
                 para_SouthWest=c(gamma_SouthWest_Uniform,gamma_SouthWest_Weibull),para_SouthWest2=c(eta_SouthWest_Uniform,eta_SouthWest_Weibull),
                 para_Midlands=c(gamma_Midlands_Uniform,gamma_Midlands_Weibull),para_Midlands2=c(eta_Midlands_Uniform,eta_Midlands_Weibull),
                 para_EastofEngland=c(gamma_EastofEngland_Uniform,gamma_EastofEngland_Weibull),para_EastofEngland2=c(eta_EastofEngland_Uniform,eta_EastofEngland_Weibull))
them<-theme(
  text = element_text(size=font_size),
  plot.title = element_text(face = "bold", size = font_size_title),
  legend.background = element_rect(fill = "white", size = 0.5, colour = "white"),
  legend.justification = c(0, 1),
  legend.position = "none",
  # legend.title = element_blank(),
  axis.ticks = element_line(colour = "grey50", size = 0.2),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank(),
  # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
  # panel.grid.major = element_line(colour = "grey50"),
  # panel.grid.minor = element_blank()
)
p1<-ggplot(data, aes(x=para_beta, fill=class))+geom_density(alpha=.3)+xlab(" ")+ylab(expression(beta))+theme_minimal()+theme(
  text = element_text(size=font_size),
  legend.background = element_rect(fill = "white", size = 0.5, colour = "white"),
  legend.justification = c(0, 1),
  legend.position = c(1,1),
  legend.title = element_blank(),
  axis.ticks = element_line(colour = "grey50", size = 0.2),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank(),
  # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
  # panel.grid.major = element_line(colour = "grey50"),
  # panel.grid.minor = element_blank()
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

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Results"
folder = paste(folder_strings, sep = "", collapse = "/")

tiff(file=paste(folder,"/comparison1_timevarying1.tiff", sep = ""),
     width=27, height=17, units="cm", res=300)
ggarrange(p2,p3, p4,p5, p6,p7, p8, p1,p9,p10, p11,p12, p13,p14, p15)
dev.off()

