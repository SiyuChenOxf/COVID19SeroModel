#This code is to generate posterior estimations for constant IFR model using Weibull prior for beta

set.seed(100)
#options (mc.cores = parallel::detectCores ())

library(rstan)

## Load regional death and seroprevalence data ##
folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

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

data_sir <- list(  n_days2_SouthWest = n_days2_SouthWest,n_days2_NorthWest = n_days2_NorthWest,n_days2_London = n_days2_London ,
                   n_days2_SouthEast = n_days2_SouthEast,n_days2_NorthEastYorkshireHumber = n_days2_NorthEastYorkshireHumber,n_days2_Midlands = n_days2_Midlands ,n_days2_EastofEngland = n_days2_EastofEngland ,
                   
                   n_days_SouthWest=n_days_SouthWest,n_days_NorthWest=n_days_NorthWest,n_days_London=n_days_London ,
                   n_days_SouthEast = n_days_SouthEast,n_days_NorthEastYorkshireHumber = n_days_NorthEastYorkshireHumber,n_days_Midlands = n_days_Midlands ,n_days_EastofEngland = n_days_EastofEngland ,
                   
                   sero_SouthWest= sero_SouthWest,sero_NorthWest= sero_NorthWest, sero_London = sero_London ,
                   sero_SouthEast = sero_SouthEast,sero_NorthEastYorkshireHumber = sero_NorthEastYorkshireHumber,sero_Midlands = sero_Midlands ,sero_EastofEngland = sero_EastofEngland,                   
                                      
                   daily_death_SouthWest =daily_death_SouthWest,daily_death_NorthWest =daily_death_NorthWest,daily_death_London =daily_death_London,
                   daily_death_SouthEast = daily_death_SouthEast,daily_death_NorthEastYorkshireHumber = daily_death_NorthEastYorkshireHumber,daily_death_Midlands = daily_death_Midlands ,daily_death_EastofEngland = daily_death_EastofEngland ,
                   
                   t_SouthWest=t_SouthWest,t_NorthWest=t_NorthWest, t_London=t_London ,
                   t_SouthEast = t_SouthEast,t_NorthEastYorkshireHumber = t_NorthEastYorkshireHumber,t_Midlands = t_Midlands ,t_EastofEngland = t_EastofEngland ,
                   
                   t2_SouthWest=t2_SouthWest,t2_NorthWest=t2_NorthWest, t2_London=t2_London,
                   t2_SouthEast = t2_SouthEast,t2_NorthEastYorkshireHumber = t2_NorthEastYorkshireHumber,t2_Midlands = t2_Midlands ,t2_EastofEngland = t2_EastofEngland
)

niter <- 20000
chains<-4

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Codes"
folder = paste(folder_strings, sep = "", collapse = "/")

SeroModelConstantIFR<-stan(file=paste(folder,"/SeroEstimationConstantIFR_WeibullPrior.stan", sep = ""),
                           data = data_sir,
                           iter = niter,
                           chains = chains)
print(SeroModelConstantIFR)

beta <- rstan::extract(SeroModelConstantIFR)$beta
gamma_NorthEastYorkshireHumber<- rstan::extract(SeroModelConstantIFR)$gamma_NorthEastYorkshireHumber
gamma_London <- rstan::extract(SeroModelConstantIFR)$gamma_London
gamma_NorthWest <- rstan::extract(SeroModelConstantIFR)$gamma_NorthWest
gamma_SouthEast <- rstan::extract(SeroModelConstantIFR)$gamma_SouthEast
gamma_SouthWest <- rstan::extract(SeroModelConstantIFR)$gamma_SouthWest
gamma_Midlands <- rstan::extract(SeroModelConstantIFR)$gamma_Midlands
gamma_EastofEngland <- rstan::extract(SeroModelConstantIFR)$gamma_EastofEngland

quantile(beta, probs = c(0.025,0.5,0.975))
quantile(gamma_NorthEastYorkshireHumber, probs = c(0.025,0.5,0.975))
quantile(gamma_London, probs = c(0.025,0.5,0.975))
quantile(gamma_NorthWest, probs = c(0.025,0.5,0.975))
quantile(gamma_SouthEast, probs = c(0.025,0.5,0.975))
quantile(gamma_SouthWest, probs = c(0.025,0.5,0.975))
quantile(gamma_Midlands, probs = c(0.025,0.5,0.975))
quantile(gamma_EastofEngland, probs = c(0.025,0.5,0.975))

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")
saveRDS(SeroModelConstantIFR, file=paste(folder,"/ConstantSeroModelConstantIFR_Weibull.rds", sep = ""))

