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

SurviCurve_ImperialReport<-1-pweibull(seq(400),shape=3.67,scale=130.4)
SurviCurve_median<-1-pexp(seq(400),rate=0.0057)
SurviCurve_lower<-1-pexp(seq(400),rate=0.0051)
SurviCurve_upper<-1-pexp(seq(400),rate=0.0063)

data=data.frame(output = c(rep("Imperial report", 400), rep("Model 1, 2, 3", 400)),t=c(seq(400),seq(400)), 
                 median = c(SurviCurve_ImperialReport,SurviCurve_median), 
                 lower1 = c(SurviCurve_ImperialReport, SurviCurve_lower),
                 upper1 = c(SurviCurve_ImperialReport, SurviCurve_upper))

ggplot()+geom_line(data2, aes(x=t, y = median, group = output, colour = output))+
  geom_line(data1, aes(x=t, y = median, group = output, colour = output))+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limit = c(0, 1))+
  theme_minimal() +
  ylab(" Probability of seropositivity persistence after seroconversion ") +
  xlab(" Time (days) ")+
  theme(
    text = element_text(size=19),
    plot.title = element_text(face = "bold", size = 20,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(0.6,0.9),
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    axis.title.y = element_text(size=19),
    axis.title.x = element_text(size=19),
    axis.text.y = element_text(size=19),
    axis.text.x = element_text(size=19),
    plot.margin = margin(t=0, r=1, b=0, l=0, "cm")
  )


p<-ggplot(data, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1.5) +  
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  guides(color=guide_legend(override.aes=list(fill=NA)))
styledp<- p +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1.00), limit = c(0, 1))+
  theme_minimal() +
  ylab(" Probability of seropositivity persistence after seroconversion ") +
  xlab(" Time (days) ")+
  theme(
    text = element_text(size=19),
    plot.title = element_text(face = "bold", size = 20,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(0.6,0.9),
    legend.title = element_blank(),
    legend.text = element_text(size=20),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    axis.title.y = element_text(size=19),
    axis.title.x = element_text(size=19),
    axis.text.y = element_text(size=19),
    axis.text.x = element_text(size=19),
    plot.margin = margin(t=0, r=1, b=0, l=0, "cm")
  )

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Results"
folder = paste(folder_strings, sep = "", collapse = "/")

tiff(file=paste(folder,"/WeibullExpoential.tiff", sep = ""),
     width=23, height=20, units="cm", res=300)
ggarrange(styledp)
dev.off()

