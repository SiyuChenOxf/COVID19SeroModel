#This code is to plot comparison of posterior estimations for constant IFR model using Weibull prior and Uniform prior for beta
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

ConstantSeroModelConstantIFR_Weibull <- readRDS(paste(folder,"/ConstantSeroModelConstantIFR_Weibull.rds", sep = ""))#Load posterior estimations for constant IFR model using Weibull prior for beta
SeroModelConstantIFR_Uniform <- readRDS(paste(folder,"/SeroModelConstantIFR.rds", sep = ""))                        #Load posterior estimations for constant IFR model using Uniform prior for beta

#Posterior estimations for parameters for constant IFR model using Weibull prior for beta
beta_Weibull <- rstan::extract(ConstantSeroModelConstantIFR_Weibull)$beta
gamma_NorthEastYorkshireHumber_Weibull<- rstan::extract(ConstantSeroModelConstantIFR_Weibull)$gamma_NorthEastYorkshireHumber
gamma_London_Weibull <- rstan::extract(ConstantSeroModelConstantIFR_Weibull)$gamma_London
gamma_NorthWest_Weibull <- rstan::extract(ConstantSeroModelConstantIFR_Weibull)$gamma_NorthWest
gamma_SouthEast_Weibull <- rstan::extract(ConstantSeroModelConstantIFR_Weibull)$gamma_SouthEast
gamma_SouthWest_Weibull <- rstan::extract(ConstantSeroModelConstantIFR_Weibull)$gamma_SouthWest
gamma_Midlands_Weibull <- rstan::extract(ConstantSeroModelConstantIFR_Weibull)$gamma_Midlands
gamma_EastofEngland_Weibull <- rstan::extract(ConstantSeroModelConstantIFR_Weibull)$gamma_EastofEngland

#Posterior estimations for parameters for constant IFR model using Uniform prior for beta
beta_Uniform <- rstan::extract(SeroModelConstantIFR_Uniform)$beta
gamma_NorthEastYorkshireHumber_Uniform<- rstan::extract(SeroModelConstantIFR_Uniform)$gamma_NorthEastYorkshireHumber
gamma_London_Uniform <- rstan::extract(SeroModelConstantIFR_Uniform)$gamma_London
gamma_NorthWest_Uniform <- rstan::extract(SeroModelConstantIFR_Uniform)$gamma_NorthWest
gamma_SouthEast_Uniform <- rstan::extract(SeroModelConstantIFR_Uniform)$gamma_SouthEast
gamma_SouthWest_Uniform <- rstan::extract(SeroModelConstantIFR_Uniform)$gamma_SouthWest
gamma_Midlands_Uniform <- rstan::extract(SeroModelConstantIFR_Uniform)$gamma_Midlands
gamma_EastofEngland_Uniform <- rstan::extract(SeroModelConstantIFR_Uniform)$gamma_EastofEngland
data<-data.frame(class=factor(rep(c("Model 2","Model 12"),each=length(beta_Uniform))),
                 para_beta=c(beta_Weibull,beta_Uniform),
                 para_London=c(gamma_London_Weibull,gamma_London_Uniform),
                 para_NorthEast=c(gamma_NorthEastYorkshireHumber_Weibull,gamma_NorthEastYorkshireHumber_Uniform),
                 para_SouthEast=c(gamma_SouthEast_Weibull,gamma_SouthEast_Uniform),
                 para_NorthWest=c(gamma_NorthWest_Weibull,gamma_NorthWest_Uniform),
                 para_SouthWest=c(gamma_SouthWest_Weibull,gamma_SouthWest_Uniform),
                 para_Midlands=c(gamma_Midlands_Weibull,gamma_Midlands_Uniform),
                 para_EastofEngland=c(gamma_EastofEngland_Weibull,gamma_EastofEngland_Uniform))
them<-theme(
  text = element_text(size=font_size),
  # plot.title = element_text(face = "bold", size = font_size_title),
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
  legend.background = element_rect(fill = "white", size = 0.5, colour = "white"),
  legend.justification = c(0, 1),
  legend.position = c(1,1),
  legend.title = element_blank(),
  axis.text.y=element_blank(),
  axis.ticks.y=element_blank(),
  axis.ticks = element_line(colour = "grey50", size = 0.2),
  # plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
  # panel.grid.major = element_line(colour = "grey50"),
  # panel.grid.minor = element_blank()
)

p2<-ggplot(data, aes(x=para_London, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[London]))+theme_minimal()+them
p3<-ggplot(data, aes(x=para_NorthEast, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[NorthEast]))+theme_minimal()+them+ggtitle("")
p4<-ggplot(data, aes(x=para_NorthWest, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[NorthWest]))+theme_minimal()+them
p5<-ggplot(data, aes(x=para_SouthEast, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[SouthEast]))+theme_minimal()+them
p6<-ggplot(data, aes(x=para_SouthWest, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[SouthWest]))+theme_minimal()+them
p7<-ggplot(data, aes(x=para_Midlands, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[Midlands]))+theme_minimal()+them
p8<-ggplot(data, aes(x=para_EastofEngland, fill=class)) + geom_density(alpha=.3)+xlab(" ")+ylab(expression(gamma[EastofEngland]))+theme_minimal()+them

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Results"
folder = paste(folder_strings, sep = "", collapse = "/")

tiff(file=paste(folder,"/comparison1.tiff", sep = ""),
     width=35, height=18, units="cm", res=300)
ggarrange(p2,p3, p4,p5, p6,p7, p8, p1)
dev.off()
