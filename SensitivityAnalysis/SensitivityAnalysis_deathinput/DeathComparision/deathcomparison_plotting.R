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
ymax_exposure <- 0.30
colors_Dark<-brewer.pal(7,"Dark2")
colors_Spectral<-brewer.pal(7,"Spectral")
font_size = 14
font_size_title = 16
lwd = 0.5
pt_size = 0.4
right_margin=1


England_data<-read.csv("England_deaths.csv")

Death_data<-data.frame(date=as.Date(England_data$date),death=c(England_data$Certificate.Death,England_data$X28days.Positive.Death),lable=rep(c("certificate death","28 days positive death"),each=length(England_data$X28days.Positive.Death)))

p1<-ggplot(data=Death_data,aes(x=date,y=death, group = lable, colour = lable))+geom_line(lwd = 1)+
  scale_y_continuous(breaks = c(0,200,400,600,800,1000,1200,1400), limit = c(0, 1400))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-04-01", "2020-07-01", "2020-10-01","2021-01-01",  "2021-04-01")), labels=c("Jan 2020","Apr 2020", "Jul 2020", "Oct 2020","Jan 2021","Apr 2021"), limit = as.Date(c("2020-01-01","2021-04-03")))+
  xlab(" ")+ylab("Daily deaths")+theme_minimal()+
  theme(
    text = element_text(size=20),
    plot.title = element_text(face = "bold", size = 22,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 0.5, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(0.5,0.9),
    legend.text = element_text(size=18),
    legend.title = element_blank(),
    axis.title.y = element_text(size=18),
    axis.text.y = element_text(size=18),
    axis.text.x = element_text(size=18),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=18),
    axis.ticks = element_line(colour = "grey50", size = 0.2))

tiff(file="DeathComparsion.tiff",
     width=30, height=15, units="cm", res=300)
ggarrange(p1)
dev.off()

