#This code is to generate Figure 3 for time-varying IFR model in the paper,
#so please run COVIDSeroModel_TimeVaryingIFR.R firstly
#OR comment out line 38 and run it.

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
library(bayesplot)
library(rstanarm)

## common values for plotting
ymax_kft <- 0.012
ymax_exposure <- 25
library(RColorBrewer)
colors_Dark<-brewer.pal(7,"Dark2")
colors_Spectral<-brewer.pal(7,"Spectral")
font_size = 14
font_size_title = 16
xlab__size_title = 12
lwd = 0.5
pt_size = 0.4
right_margin=1

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Results"
folder = paste(folder_strings, sep = "", collapse = "/")

# SeroModelTimeVaryingIFR <- readRDS(paste(folder,"/SeroModelTimeVaryingIFR_delta_p_14.rds", sep = ""))

####Plotting posterior in the time varying IFR case####
posterior <- as.matrix(SeroModelTimeVaryingIFR)
color_scheme_set("gray")
tiff(file=paste(folder,"/TimeVaryingIFR_0.95_intervals.tiff", sep = ""),
     width=26, height=14, units="cm", res=300)
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 95% intervals")
p1<-mcmc_areas(posterior,
               pars = c("beta","gamma_London", "gamma_EastofEngland", "gamma_Midlands",
                        "gamma_SouthWest","gamma_SouthEast","gamma_NorthWest","gamma_NorthEastYorkshireHumber"),
               prob = 0.95) + plot_title+
  scale_y_discrete(labels=c(
    "beta" = expression(beta),
    "gamma_EastofEngland" =  expression(gamma[East]),
    "gamma_Midlands" =  expression(gamma[Midlands]),
    "gamma_SouthEast" =  expression(gamma[SouthEast]),
    "gamma_SouthWest" =  expression(gamma[SouthWest]),
    "gamma_NorthEastYorkshireHumber" =  expression(gamma[NorthEast]),
    "gamma_NorthWest" =  expression(gamma[NorthWest]),
    "gamma_London" =  expression(gamma[London])
  ))+
  theme_minimal() +
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    #legend.position = c(0.085,0.85),
    legend.position = "none",
    legend.title = element_blank(),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank()
  )

p2<-mcmc_areas(posterior,
               pars = c("eta_London", "eta_EastofEngland", "eta_Midlands",
                        "eta_SouthWest","eta_SouthEast","eta_NorthWest","eta_NorthEastYorkshireHumber"),
               prob = 0.95) + 
  scale_y_discrete(labels=c(
    "eta_EastofEngland" =  expression(eta[East]),
    "eta_Midlands" =  expression(eta[Midlands]),
    "eta_SouthEast" =  expression(eta[SouthEast]),
    "eta_SouthWest" =  expression(eta[SouthWest]),
    "eta_NorthEastYorkshireHumber" =  expression(eta[NorthEast]),
    "eta_NorthWest" =  expression(eta[NorthWest]),
    "eta_London" =  expression(eta[London])
  ))+
  theme_minimal() +
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    #legend.position = c(0.085,0.85),
    legend.position = "none",
    legend.title = element_blank(),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank()
  )

grid.arrange(p1, p2, nrow = 1)
dev.off()

##Plotting seroprevalence and exposure in the time varying IFR##

##################
### LONDON
##################
sim<-length(beta)
x_London<-epsilon_London<-kft_London<-matrix(0,sim,n_days_London)

for (i in 1:sim) {
  x_London[i,]<-rnbinom(rep(1,n_days_London), size= 100, mu=cumsum(exp(beta[i]*t_London)*(1-(gamma_London[i]-(eta_London[i]*gamma_London[i])*Epidemic_Stage_London))/(gamma_London[i]-(eta_London[i]*gamma_London[i])*Epidemic_Stage_London)*daily_death_London)/(exp(beta[i]*t_London)))/(P0_London-cumul_death_London)
  epsilon_London[i,]<-cumsum((1-(gamma_London[i]-(eta_London[i]*gamma_London[i])*Epidemic_Stage_London))/(gamma_London[i]-(eta_London[i]*gamma_London[i])*Epidemic_Stage_London)*daily_death_London)/(P0_London-cumul_death_London)
  kft_London[i,]<-gamma_London[i]-(eta_London[i]*gamma_London[i])*Epidemic_Stage_London
}
epsilon_London<-epsilon_London[,(delta_epsilon+1):n_days_London]


data1London = data.frame(output = c(rep("Exposure", n_days_London-delta_epsilon), rep("Seroprevalence", n_days_London)), 
                         t=c(as.Date(London_data$Date)[1:(n_days_London-delta_epsilon)],as.Date(London_data$Date)[1:n_days_London]), 
                         median = c(100*apply(epsilon_London, 2, function(x) quantile(x, probs = 0.5)), 
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.5))), 
                         lower1 = c(100*apply(epsilon_London, 2, function(x) quantile(x, probs = 0.025)), 
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.025))), 
                         upper1 = c(100*apply(epsilon_London, 2, function(x) quantile(x, probs = 0.975)),
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.975))),
                         lower2 = c(100*apply(epsilon_London, 2, function(x) quantile(x, probs = 0.25)), 
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.25))), 
                         upper2 = c(100*apply(epsilon_London, 2, function(x) quantile(x, probs = 0.75)),
                                    100*apply(x_London, 2, function(x) quantile(x, probs = 0.75))))
data2London = data.frame( t=as.Date(London_data$Date)[1:n_days_London][t2_London], value= 100*London_data$sero[t2_London], upper= 100*London_data$sero_upper[t2_London], lower = 100*London_data$sero_lower[t2_London])

data3London  = data.frame(output = c(rep("Infection fatality rate", n_days_London)), t=as.Date(London_data$Date)[1:n_days_London],
                          median = apply(kft_London, 2, function(x) quantile(x, probs = 0.5)), 
                          lower1 = apply(kft_London, 2, function(x) quantile(x, probs = 0.25)),
                          upper1 = apply(kft_London, 2, function(x) quantile(x, probs = 0.75)),
                          lower2 = apply(kft_London, 2, function(x) quantile(x, probs = 0.025)),
                          upper2 = apply(kft_London, 2, function(x) quantile(x, probs = 0.975)))


p1London<-ggplot(data1London, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2London, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1London <- p1London +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  ggtitle("London")+
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
    axis.title.y = element_text(size=xlab__size_title),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=xlab__size_title),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank()
  )
# styled1London

##################
#####NorthEast & Yorkshire and the Humber####
##################
sim<-length(beta)
x_NorthEastYorkshireHumber<-epsilon_NorthEastYorkshireHumber<-kft_NorthEastYorkshireHumber<-matrix(0,sim,n_days_NorthEastYorkshireHumber)
for (i in 1:sim) {
  x_NorthEastYorkshireHumber[i,]<-rnbinom(rep(1,n_days_NorthEastYorkshireHumber), size= 100, mu=cumsum(exp(beta[i]*t_NorthEastYorkshireHumber)*(1-(gamma_NorthEastYorkshireHumber[i]-(eta_NorthEastYorkshireHumber[i]*gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber))/(gamma_NorthEastYorkshireHumber[i]-(eta_NorthEastYorkshireHumber[i]*gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber)*daily_death_NorthEastYorkshireHumber)/(exp(beta[i]*t_NorthEastYorkshireHumber)))/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber)
  epsilon_NorthEastYorkshireHumber[i,]<-cumsum((1-(gamma_NorthEastYorkshireHumber[i]-(eta_NorthEastYorkshireHumber[i]*gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber))/(gamma_NorthEastYorkshireHumber[i]-(eta_NorthEastYorkshireHumber[i]*gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber)*daily_death_NorthEastYorkshireHumber)/(P0_NorthEastYorkshireHumber-cumul_death_NorthEastYorkshireHumber)
  kft_NorthEastYorkshireHumber[i,]<-gamma_NorthEastYorkshireHumber[i]-(eta_NorthEastYorkshireHumber[i]*gamma_NorthEastYorkshireHumber[i])*Epidemic_Stage_NorthEastYorkshireHumber
}
epsilon_NorthEastYorkshireHumber<-epsilon_NorthEastYorkshireHumber[,(delta_epsilon+1):n_days_NorthEastYorkshireHumber]

data1NorthEastYorkshireHumber  = data.frame(output = c(rep("Exposure", n_days_NorthEastYorkshireHumber-delta_epsilon), rep("Seroprevalence", n_days_NorthEastYorkshireHumber)), 
                                            t=c(as.Date(NorthEast_data$Date)[1:(n_days_NorthEastYorkshireHumber-delta_epsilon)],as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber]), 
                                            median = c(100*apply(epsilon_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.5)), 
                                                       100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.5))), 
                                            lower1 = c(100*apply(epsilon_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.025)), 
                                                       100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.025))), 
                                            upper1 = c(100*apply(epsilon_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.975)),
                                                       100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.975))),
                                            lower2 = c(100*apply(epsilon_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.25)), 
                                                       100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.25))), 
                                            upper2 = c(100*apply(epsilon_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.75)),
                                                       100*apply(x_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.75)))
)
data2NorthEastYorkshireHumber = data.frame( t=as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber][t2_NorthEastYorkshireHumber], value= 100*NorthEast_data$sero[t2_NorthEastYorkshireHumber], upper= 100*NorthEast_data$sero_upper[t2_NorthEastYorkshireHumber], lower = 100*NorthEast_data$sero_lower[t2_NorthEastYorkshireHumber])

data3NorthEastYorkshireHumber  = data.frame(output = c(rep("Infection fatality rate", n_days_NorthEastYorkshireHumber)), t=as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber],
                                            median = apply(kft_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.5)), 
                                            lower1 = apply(kft_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.25)),
                                            upper1 = apply(kft_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.75)),
                                            lower2 = apply(kft_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.025)),
                                            upper2 = apply(kft_NorthEastYorkshireHumber, 2, function(x) quantile(x, probs = 0.975)))

p1NorthEastYorkshireHumber<-ggplot(data1NorthEastYorkshireHumber, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2NorthEastYorkshireHumber, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1NorthEastYorkshireHumber <- p1NorthEastYorkshireHumber +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  ggtitle("North East")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=xlab__size_title),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank()
  )
# styled1NorthEastYorkshireHumber

##################
#####NorthWest####
##################
sim<-length(beta)
x_NorthWest<-epsilon_NorthWest<-kft_NorthWest<-matrix(0,sim,n_days_NorthWest)
for (i in 1:sim) {
  x_NorthWest[i,]<-rnbinom(rep(1,n_days_NorthWest), size= 100, mu=cumsum(exp(beta[i]*t_NorthWest)*(1-(gamma_NorthWest[i]-(eta_NorthWest[i]*gamma_NorthWest[i])*Epidemic_Stage_NorthWest))/(gamma_NorthWest[i]-(eta_NorthWest[i]*gamma_NorthWest[i])*Epidemic_Stage_NorthWest)*daily_death_NorthWest)/(exp(beta[i]*t_NorthWest)))/(P0_NorthWest-cumul_death_NorthWest)
  epsilon_NorthWest[i,]<-cumsum((1-(gamma_NorthWest[i]-(eta_NorthWest[i]*gamma_NorthWest[i])*Epidemic_Stage_NorthWest))/(gamma_NorthWest[i]-(eta_NorthWest[i]*gamma_NorthWest[i])*Epidemic_Stage_NorthWest)*daily_death_NorthWest)/(P0_NorthWest-cumul_death_NorthWest)
  kft_NorthWest[i,]<-gamma_NorthWest[i]-(eta_NorthWest[i]*gamma_NorthWest[i])*Epidemic_Stage_NorthWest
}
epsilon_NorthWest<-epsilon_NorthWest[,(delta_epsilon+1):n_days_NorthWest]

data1NorthWest  = data.frame(output = c(rep("Exposure", n_days_NorthWest-delta_epsilon), rep("Seroprevalence", n_days_NorthWest)), 
                             t=c(as.Date(NorthWest_data$Date)[1:(n_days_NorthWest-delta_epsilon)],as.Date(NorthWest_data$Date)[1:n_days_NorthWest]), 
                             median = c(100*apply(epsilon_NorthWest, 2, function(x) quantile(x, probs = 0.5)), 
                                        100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.5))), 
                             lower1 = c(100*apply(epsilon_NorthWest, 2, function(x) quantile(x, probs = 0.025)), 
                                        100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.025))), 
                             upper1 = c(100*apply(epsilon_NorthWest, 2, function(x) quantile(x, probs = 0.975)),
                                        100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.975))),
                             lower2 = c(100*apply(epsilon_NorthWest, 2, function(x) quantile(x, probs = 0.25)), 
                                        100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.25))), 
                             upper2 = c(100*apply(epsilon_NorthWest, 2, function(x) quantile(x, probs = 0.75)),
                                        100*apply(x_NorthWest, 2, function(x) quantile(x, probs = 0.75))))
data2NorthWest = data.frame( t=as.Date(NorthWest_data$Date)[1:n_days_NorthWest][t2_NorthWest], value= 100*NorthWest_data$sero[t2_NorthWest], upper= 100*NorthWest_data$sero_upper[t2_NorthWest], lower = 100*NorthWest_data$sero_lower[t2_NorthWest])

data3NorthWest  = data.frame(output = c(rep("Infection fatality rate", n_days_NorthWest)), t=as.Date(NorthWest_data$Date)[1:n_days_NorthWest],
                             median = apply(kft_NorthWest, 2, function(x) quantile(x, probs = 0.5)), 
                             lower1 = apply(kft_NorthWest, 2, function(x) quantile(x, probs = 0.25)),
                             upper1 = apply(kft_NorthWest, 2, function(x) quantile(x, probs = 0.75)),
                             lower2 = apply(kft_NorthWest, 2, function(x) quantile(x, probs = 0.025)),
                             upper2 = apply(kft_NorthWest, 2, function(x) quantile(x, probs = 0.975)))

p1NorthWest<-ggplot(data1NorthWest, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2NorthWest, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1NorthWest <- p1NorthWest +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  ggtitle("North West")+
  theme_minimal() +
  ylab(" ") +
  xlab(" 2020 ")+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=xlab__size_title),
    legend.title = element_blank(),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank()
  )
# styled1NorthWest

##################
### SouthWest
##################
sim<-length(beta)
x_SouthWest<-epsilon_SouthWest<-kft_SouthWest<-matrix(0,sim,n_days_SouthWest)
for (i in 1:sim) {
  x_SouthWest[i,]<-rnbinom(rep(1,n_days_SouthWest), size= 100, mu=cumsum(exp(beta[i]*t_SouthWest)*(1-(gamma_SouthWest[i]-(eta_SouthWest[i]*gamma_SouthWest[i])*Epidemic_Stage_SouthWest))/(gamma_SouthWest[i]-(eta_SouthWest[i]*gamma_SouthWest[i])*Epidemic_Stage_SouthWest)*daily_death_SouthWest)/(exp(beta[i]*t_SouthWest)))/(P0_SouthWest-cumul_death_SouthWest)
  epsilon_SouthWest[i,]<-cumsum((1-(gamma_SouthWest[i]-(eta_SouthWest[i]*gamma_SouthWest[i])*Epidemic_Stage_SouthWest))/(gamma_SouthWest[i]-(eta_SouthWest[i]*gamma_SouthWest[i])*Epidemic_Stage_SouthWest)*daily_death_SouthWest)/(P0_SouthWest-cumul_death_SouthWest)
  kft_SouthWest[i,]<-gamma_SouthWest[i]-(eta_SouthWest[i]*gamma_SouthWest[i])*Epidemic_Stage_SouthWest
}
epsilon_SouthWest<-epsilon_SouthWest[,(delta_epsilon+1):n_days_SouthWest]


data1SouthWest = data.frame(output = c(rep("Exposure", n_days_SouthWest-delta_epsilon), rep("Seroprevalence", n_days_SouthWest)), 
                            t=c(as.Date(SouthWest_data$Date)[1:(n_days_SouthWest-delta_epsilon)],as.Date(SouthWest_data$Date)[1:n_days_SouthWest]), 
                            median = c(100*apply(epsilon_SouthWest, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.5))), 
                            lower1 = c(100*apply(epsilon_SouthWest, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.025))), 
                            upper1 = c(100*apply(epsilon_SouthWest, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.975))),
                            lower2 = c(100*apply(epsilon_SouthWest, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.25))), 
                            upper2 = c(100*apply(epsilon_SouthWest, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_SouthWest, 2, function(x) quantile(x, probs = 0.75))))
data2SouthWest = data.frame(t=as.Date(SouthWest_data$Date)[1:n_days_SouthWest][t2_SouthWest], value= 100*SouthWest_data$sero[t2_SouthWest], upper= 100*SouthWest_data$sero_upper[t2_SouthWest], lower = 100*SouthWest_data$sero_lower[t2_SouthWest])

data3SouthWest  = data.frame(output = c(rep("Infection fatality rate", n_days_SouthWest)), t=as.Date(SouthWest_data$Date)[1:n_days_SouthWest],
                             median = apply(kft_SouthWest, 2, function(x) quantile(x, probs = 0.5)), 
                             lower1 = apply(kft_SouthWest, 2, function(x) quantile(x, probs = 0.25)),
                             upper1 = apply(kft_SouthWest, 2, function(x) quantile(x, probs = 0.75)),
                             lower2= apply(kft_SouthWest, 2, function(x) quantile(x, probs = 0.025)),
                             upper2 = apply(kft_SouthWest, 2, function(x) quantile(x, probs = 0.975)))


p1SouthWest<-ggplot(data1SouthWest, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2SouthWest, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1SouthWest <- p1SouthWest +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  ggtitle("South West")+
  theme_minimal() +
  ylab("") +
  xlab("2020 ")+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(0.085,0.85),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=xlab__size_title),
    legend.title = element_blank(),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank()
  )
# styled1SouthWest

##################
####SouthEast#####
##################
sim<-length(beta)
x_SouthEast<-epsilon_SouthEast<-kft_SouthEast<-matrix(0,sim,n_days_SouthEast)
for (i in 1:sim) {
  x_SouthEast[i,]<-rnbinom(rep(1,n_days_SouthEast), size= 100, mu=cumsum(exp(beta[i]*t_SouthEast)*(1-(gamma_SouthEast[i]-(eta_SouthEast[i]*gamma_SouthEast[i])*Epidemic_Stage_SouthEast))/(gamma_SouthEast[i]-(eta_SouthEast[i]*gamma_SouthEast[i])*Epidemic_Stage_SouthEast)*daily_death_SouthEast)/(exp(beta[i]*t_SouthEast)))/(P0_SouthEast-cumul_death_SouthEast)
  epsilon_SouthEast[i,]<-cumsum((1-(gamma_SouthEast[i]-(eta_SouthEast[i]*gamma_SouthEast[i])*Epidemic_Stage_SouthEast))/(gamma_SouthEast[i]-(eta_SouthEast[i]*gamma_SouthEast[i])*Epidemic_Stage_SouthEast)*daily_death_SouthEast)/(P0_SouthEast-cumul_death_SouthEast)
  kft_SouthEast[i,]<-gamma_SouthEast[i]-(eta_SouthEast[i]*gamma_SouthEast[i])*Epidemic_Stage_SouthEast
}
epsilon_SouthEast<-epsilon_SouthEast[,(delta_epsilon+1):n_days_SouthEast]

data1SouthEast = data.frame(output = c(rep("Exposure", n_days_SouthEast-delta_epsilon), rep("Seroprevalence", n_days_SouthEast)), 
                            t=c(as.Date(SouthEast_data$Date)[1:(n_days_SouthEast-delta_epsilon)],as.Date(SouthEast_data$Date)[1:n_days_SouthEast]), 
                            median = c(100*apply(epsilon_SouthEast, 2, function(x) quantile(x, probs = 0.5)), 
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.5))), 
                            lower1 = c(100*apply(epsilon_SouthEast, 2, function(x) quantile(x, probs = 0.025)), 
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.025))), 
                            upper1 = c(100*apply(epsilon_SouthEast, 2, function(x) quantile(x, probs = 0.975)),
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.975))),
                            lower2 = c(100*apply(epsilon_SouthEast, 2, function(x) quantile(x, probs = 0.25)), 
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.25))), 
                            upper2 = c(100*apply(epsilon_SouthEast, 2, function(x) quantile(x, probs = 0.75)),
                                       100*apply(x_SouthEast, 2, function(x) quantile(x, probs = 0.75))))
data2SouthEast = data.frame(t=as.Date(SouthEast_data$Date)[1:n_days_SouthEast][t2_SouthEast], value= 100*SouthEast_data$sero[t2_SouthEast], upper= 100*SouthEast_data$sero_upper[t2_SouthEast], lower = 100*SouthEast_data$sero_lower[t2_SouthEast])

data3SouthEast  = data.frame(output = c(rep("Infection fatality rate", n_days_SouthEast)), t=as.Date(SouthEast_data$Date)[1:n_days_SouthEast],
                             median = apply(kft_SouthEast, 2, function(x) quantile(x, probs = 0.5)), 
                             lower1 = apply(kft_SouthEast, 2, function(x) quantile(x, probs = 0.25)),
                             upper1 = apply(kft_SouthEast, 2, function(x) quantile(x, probs = 0.75)),
                             lower2 = apply(kft_SouthEast, 2, function(x) quantile(x, probs = 0.025)),
                             upper2 = apply(kft_SouthEast, 2, function(x) quantile(x, probs = 0.975)))


p1SouthEast<-ggplot(data1SouthEast, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2SouthEast, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1SouthEast <- p1SouthEast +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  ggtitle("South East")+
  theme_minimal() +
  ylab(" Percentage (%)") +
  xlab(" 2020")+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    #legend.position = c(0.085,0.85),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=xlab__size_title),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=xlab__size_title),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank()
  )
# styled1SouthEast

##################
####Midlands#####
sim<-length(beta)
x_Midlands<-epsilon_Midlands<-kft_Midlands<-matrix(0,sim,n_days_Midlands)
for (i in 1:sim) {
  x_Midlands[i,]<-rnbinom(rep(1,n_days_Midlands), size= 100, mu=cumsum(exp(beta[i]*t_Midlands)*(1-(gamma_Midlands[i]-(eta_Midlands[i]*gamma_Midlands[i])*Epidemic_Stage_Midlands))/(gamma_Midlands[i]-(eta_Midlands[i]*gamma_Midlands[i])*Epidemic_Stage_Midlands)*daily_death_Midlands)/(exp(beta[i]*t_Midlands)))/(P0_Midlands-cumul_death_Midlands)
  epsilon_Midlands[i,]<-cumsum((1-(gamma_Midlands[i]-(eta_Midlands[i]*gamma_Midlands[i])*Epidemic_Stage_Midlands))/(gamma_Midlands[i]-(eta_Midlands[i]*gamma_Midlands[i])*Epidemic_Stage_Midlands)*daily_death_Midlands)/(P0_Midlands-cumul_death_Midlands)
  kft_Midlands[i,]<-gamma_Midlands[i]-(eta_Midlands[i]*gamma_Midlands[i])*Epidemic_Stage_Midlands
}
epsilon_Midlands<-epsilon_Midlands[,(delta_epsilon+1):n_days_Midlands]

data1Midlands = data.frame(output = c(rep("Exposure", n_days_Midlands-delta_epsilon), rep("Seroprevalence", n_days_Midlands)), 
                           t=c(as.Date(EastMidlands_data$Date)[1:(n_days_Midlands-delta_epsilon)],as.Date(EastMidlands_data$Date)[1:n_days_Midlands]), 
                           median = c(100*apply(epsilon_Midlands, 2, function(x) quantile(x, probs = 0.5)), 
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.5))), 
                           lower1 = c(100*apply(epsilon_Midlands, 2, function(x) quantile(x, probs = 0.025)), 
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.025))), 
                           upper1 = c(100*apply(epsilon_Midlands, 2, function(x) quantile(x, probs = 0.975)),
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.975))),
                           lower2 = c(100*apply(epsilon_Midlands, 2, function(x) quantile(x, probs = 0.25)), 
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.25))), 
                           upper2 = c(100*apply(epsilon_Midlands, 2, function(x) quantile(x, probs = 0.75)),
                                      100*apply(x_Midlands, 2, function(x) quantile(x, probs = 0.75))))
data2Midlands = data.frame(t=as.Date(EastMidlands_data$Date)[1:n_days_Midlands][t2_Midlands], value= 100*EastMidlands_data$sero[t2_Midlands], upper= 100*EastMidlands_data$sero_upper[t2_Midlands], lower = 100*EastMidlands_data$sero_lower[t2_Midlands])

data3Midlands  = data.frame(output = c(rep("Infection fatality rate", n_days_Midlands)), t=as.Date(EastMidlands_data$Date)[1:n_days_Midlands],
                            median = apply(kft_Midlands, 2, function(x) quantile(x, probs = 0.5)), 
                            lower1 = apply(kft_Midlands, 2, function(x) quantile(x, probs = 0.25)),
                            upper1 = apply(kft_Midlands, 2, function(x) quantile(x, probs = 0.75)),
                            lower2 = apply(kft_Midlands, 2, function(x) quantile(x, probs = 0.025)),
                            upper2 = apply(kft_Midlands, 2, function(x) quantile(x, probs = 0.975)))

p1Midlands<-ggplot(data1Midlands, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2Midlands, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1Midlands <- p1Midlands +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  ggtitle("Midlands")+
  theme_minimal() +
  ylab("") +
  xlab(" 2020")+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    #legend.position = c(0.085,0.85),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=xlab__size_title),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank()
  )
# styled1Midlands

##################
####EastofEngland#
##################
sim<-length(beta)
x_EastofEngland<-epsilon_EastofEngland<-kft_EastofEngland<-matrix(0,sim,n_days_EastofEngland)
for (i in 1:sim) {
  x_EastofEngland[i,]<-rnbinom(rep(1,n_days_EastofEngland), size= 100, mu=cumsum(exp(beta[i]*t_EastofEngland)*(1-(gamma_EastofEngland[i]-(eta_EastofEngland[i]*gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland))/(gamma_EastofEngland[i]-(eta_EastofEngland[i]*gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland)*daily_death_EastofEngland)/(exp(beta[i]*t_EastofEngland)))/(P0_EastofEngland-cumul_death_EastofEngland)
  epsilon_EastofEngland[i,]<-cumsum((1-(gamma_EastofEngland[i]-(eta_EastofEngland[i]*gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland))/(gamma_EastofEngland[i]-(eta_EastofEngland[i]*gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland)*daily_death_EastofEngland)/(P0_EastofEngland-cumul_death_EastofEngland)
  kft_EastofEngland[i,]<-gamma_EastofEngland[i]-(eta_EastofEngland[i]*gamma_EastofEngland[i])*Epidemic_Stage_EastofEngland
}
epsilon_EastofEngland<-epsilon_EastofEngland[,(delta_epsilon+1):n_days_EastofEngland]


data1EastofEngland = data.frame(output = c(rep("Exposure", n_days_EastofEngland-delta_epsilon), rep("Seroprevalence", n_days_EastofEngland)), 
                                t=c(as.Date(EastofEngland_data$Date)[1:(n_days_EastofEngland-delta_epsilon)],as.Date(EastofEngland_data$Date)[1:n_days_EastofEngland]), 
                                median = c(100*apply(epsilon_EastofEngland, 2, function(x) quantile(x, probs = 0.5)), 
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.5))), 
                                lower1 = c(100*apply(epsilon_EastofEngland, 2, function(x) quantile(x, probs = 0.025)), 
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.025))), 
                                upper1 = c(100*apply(epsilon_EastofEngland, 2, function(x) quantile(x, probs = 0.975)),
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.975))),
                                lower2 = c(100*apply(epsilon_EastofEngland, 2, function(x) quantile(x, probs = 0.25)), 
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.25))), 
                                upper2 = c(100*apply(epsilon_EastofEngland, 2, function(x) quantile(x, probs = 0.75)),
                                           100*apply(x_EastofEngland, 2, function(x) quantile(x, probs = 0.75))))
data2EastofEngland = data.frame(t=as.Date(EastofEngland_data$Date)[1:n_days_EastofEngland][t2_EastofEngland], value= 100*EastofEngland_data$sero[t2_EastofEngland], upper=100*EastofEngland_data$sero_upper[t2_EastofEngland] , lower = 100*EastofEngland_data$sero_lower[t2_EastofEngland])

data3EastofEngland  = data.frame(output = c(rep("Infection fatality rate", n_days_NorthEastYorkshireHumber)), t=as.Date(NorthEast_data$Date)[1:n_days_NorthEastYorkshireHumber],
                                 median = apply(kft_EastofEngland, 2, function(x) quantile(x, probs = 0.5)), 
                                 lower1 = apply(kft_EastofEngland, 2, function(x) quantile(x, probs = 0.25)),
                                 upper1 = apply(kft_EastofEngland, 2, function(x) quantile(x, probs = 0.75)), 
                                 lower2 = apply(kft_EastofEngland, 2, function(x) quantile(x, probs = 0.025)),
                                 upper2 = apply(kft_EastofEngland, 2, function(x) quantile(x, probs = 0.975)))

p1EastofEngland<-ggplot(data1EastofEngland, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  geom_pointrange(data=data2EastofEngland, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
styled1EastofEngland <- p1EastofEngland +
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ggtitle("East")+
  ylab("") +
  xlab("2020 ")+
  scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  theme(
    text = element_text(size=font_size),
    plot.title = element_text(face = "bold", size = font_size_title,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    #legend.position = c(0.085,0.85),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=xlab__size_title),
    axis.ticks = element_line(colour = "grey50", size = 0.2),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
    panel.grid.major = element_line(colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank()
  )
# styled1EastofEngland

#Plotting exposure map
eng <- rgdal::readOGR(paste0("https://opendata.arcgis.com/datasets/",
                             "8d3a9e6e7bd445e2bdcc26cdf007eac7_4.geojson"))
countries <- rgdal::readOGR(paste0("https://opendata.arcgis.com/datasets/",
                                   "92ebeaf3caa8458ea467ec164baeefa4_0.geojson"))
eng <- sf::st_as_sf(eng)
countries <- sf::st_as_sf(countries)
UK <- countries[-1,] 
names(eng)[3] <- "Region"
names(UK)[3] <- "Region"
UK$objectid <- 10:12
eng <- eng[-2]
UK <- UK[c(1, 3, 9:11)]
UK <- rbind(eng, UK)
map_orginal <- ggplot2::ggplot(UK, ggplot2::aes(fill = Region)) + ggplot2::geom_sf()


# make map for our purposes, combine midlands (east and west) and yorkshire to NorthEast
eng <- rgdal::readOGR(paste0("https://opendata.arcgis.com/datasets/",
                             "8d3a9e6e7bd445e2bdcc26cdf007eac7_4.geojson"))
countries <- rgdal::readOGR(paste0("https://opendata.arcgis.com/datasets/",
                                   "92ebeaf3caa8458ea467ec164baeefa4_0.geojson"))
eng <- sf::st_as_sf(eng)
countries <- sf::st_as_sf(countries)
UK <- countries[-1,] 
names(eng)[3] <- "Region"
eng <- eng[-2]
eng$objectid[5] <- 4
eng$objectid[3] <- 1

eng2<-eng %>% 
  group_by(objectid) %>% 
  summarize() %>% 
  ungroup() 
eng2$Exposure <- c( median(epsilon_NorthEastYorkshireHumber[,length(epsilon_NorthEastYorkshireHumber[1,])]),median(epsilon_NorthWest[,length(epsilon_NorthWest[1,])]) , median(epsilon_Midlands[,length(epsilon_Midlands[1,])]), median(epsilon_EastofEngland[,length(epsilon_EastofEngland[1,])]), median(epsilon_London[,length(epsilon_London[1,])]), median(epsilon_SouthEast[,length(epsilon_SouthEast[1,])]), median(epsilon_SouthWest[,length(epsilon_SouthWest[1,])])) # exposure for seven regions on 2020-10-17
pmap <- ggplot2::ggplot(eng2, ggplot2::aes(fill = Exposure)) + ggplot2::geom_sf(colour = "black") + 
  scale_fill_distiller(palette = "Spectral",limits=c(0.05,0.225),breaks = c(0.05, 0.10, 0.15, 0.20),labels = scales::percent_format(accuracy = 1))+
  theme_void()+ 
  theme(
    legend.position = c(0.15,0.5),
    legend.text=element_text(size=14), 
    plot.title = element_text(hjust = 0.5, face = "bold", size = font_size_title),
    legend.title = element_blank(),
    legend.key.width = unit(1, "cm"), legend.key.height = unit(1, "cm"))+
  ggtitle("Predicted exposure \n from time-varying IFR model \n on 17/10/2020")

#Saving figures
lay <- rbind(c(1,1,2,3,4,5),
             c(1,1,6,7,8,9))
tiff(file=paste(folder,"/COVIDSeroModel_TimeVaryingIFR.tiff", sep = ""),
     width=34, height=14, units="cm", res=300)
grid.arrange(pmap,styled1London,styled1NorthEastYorkshireHumber,styled1NorthWest,styled1SouthWest,styled1SouthEast,styled1Midlands,styled1EastofEngland, layout_matrix = lay, widths=c(0.75,0.75,1.35,1.1,1.1,1.1))
dev.off()

df.ifr<-rbind(data3London,data3NorthEastYorkshireHumber,
              data3NorthWest,data3SouthWest,data3SouthEast,
              data3Midlands,data3EastofEngland)
df.ifr$region<-c(rep("London",length(data3London$output)),rep("NorthEast",length(data3NorthEastYorkshireHumber$output)),rep("NorthWest",length(data3NorthWest$output)),
                 rep("SouthWest",length(data3SouthWest$output)),rep("SouthEast",length(data3SouthEast$output)),rep("Midlands",length(data3Midlands$output)), rep("East",length(data3EastofEngland$output)))

tiff(file=paste(folder,"/ifrs.tiff", sep = ""),
     width=20, height=14, units="cm", res=300)
ggplot(df.ifr,aes(x=t,y=median,color=as.factor(region)))+geom_line(size = lwd)+theme_minimal()+ylab("IFR")+xlab("2020")+
  geom_ribbon(aes(ymin=lower2, ymax=upper2,fill=as.factor(region)), alpha=0.1, colour = NA)+
  theme(text = element_text(size=font_size), axis.text = element_text(size = font_size), legend.title = element_blank(),axis.title.x = element_text(angle = 0, vjust = -0.005, hjust = 1,size=xlab__size_title))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))
dev.off()
