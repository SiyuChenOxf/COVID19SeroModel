#This code is to generate posterior estimations for time-varying IFR model using Weibull prior for beta

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

## Load virus positivity data from GOV.UK##
positivity<-read.csv(paste(folder,"/pcr_positivity.csv", sep = ""),header = TRUE)

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
sero_London<-round(P0_London*London_data$sero[!is.na(London_data$sero)])
daily_death_London<-London_data$daily_death[1:London_death_end]
n_days_London<-length(daily_death_London)
Epidemic_Stage_London<-Epidemic_Stage_London[1:n_days_London]
cumul_death_London<-cumsum(daily_death_London)
n_days2_London<-length(sero_London)
t_London<-seq(1,n_days_London,by=1)
t2_London<-t_London[!is.na(London_data$sero)]

#NorthEast &Yorkshire and the Humber 
sero_NorthEastYorkshireHumber<-round(P0_NorthEastYorkshireHumber*NorthEast_data$sero[!is.na(NorthEast_data$sero)])
daily_death_NorthEastYorkshireHumber<-NorthEast_data$daily_death_NorthEast[1:NorthEast_death_end]+YorkshireHumber_data$daily_death_Yorkshire_Humber[1:YorkshireHumber_death_end]
n_days_NorthEastYorkshireHumber<-length(daily_death_NorthEastYorkshireHumber)
Epidemic_Stage_NorthEastYorkshireHumber<-Epidemic_Stage_NorthEastYorkshireHumber[1:n_days_NorthEastYorkshireHumber]
cumul_death_NorthEastYorkshireHumber<-cumsum(daily_death_NorthEastYorkshireHumber)
n_days2_NorthEastYorkshireHumber<-length(sero_NorthEastYorkshireHumber)
t_NorthEastYorkshireHumber<-seq(1,n_days_NorthEastYorkshireHumber,by=1)
t2_NorthEastYorkshireHumber<-t_NorthEastYorkshireHumber[!is.na(NorthEast_data$sero)]

#SouthWest
sero_SouthWest<-round(P0_SouthWest*SouthWest_data$sero[!is.na(SouthWest_data$sero)])
daily_death_SouthWest<-SouthWest_data$daily_death[1:SouthWest_death_end]
n_days_SouthWest<-length(daily_death_SouthWest)
Epidemic_Stage_SouthWest<-Epidemic_Stage_SouthWest[1:n_days_SouthWest]
cumul_death_SouthWest<-cumsum(daily_death_SouthWest)
n_days2_SouthWest<-length(sero_SouthWest)
t_SouthWest<-seq(1,n_days_SouthWest,by=1)
t2_SouthWest<-t_SouthWest[!is.na(SouthWest_data$sero)]

#NorthWest
sero_NorthWest<-round(P0_NorthWest*NorthWest_data$sero[!is.na(NorthWest_data$sero)])
daily_death_NorthWest<-NorthWest_data$daily_death[1:NorthWest_death_end]
n_days_NorthWest<-length(daily_death_NorthWest)
Epidemic_Stage_NorthWest<-Epidemic_Stage_NorthWest[1:n_days_NorthWest]
cumul_death_NorthWest<-cumsum(daily_death_NorthWest)
n_days2_NorthWest<-length(sero_NorthWest)
t_NorthWest<-seq(1,n_days_NorthWest,by=1)
t2_NorthWest<-t_NorthWest[!is.na(NorthWest_data$sero)]

#SouthEast
sero_SouthEast<-round(P0_SouthEast*SouthEast_data$sero[!is.na(SouthEast_data$sero)])
daily_death_SouthEast<-SouthEast_data$daily_death[1:SouthEast_death_end]
n_days_SouthEast<-length(daily_death_SouthEast)
Epidemic_Stage_SouthEast<-Epidemic_Stage_SouthEast[1:n_days_SouthEast]
cumul_death_SouthEast<-cumsum(daily_death_SouthEast)
n_days2_SouthEast<-length(sero_SouthEast)
t_SouthEast<-seq(1,n_days_SouthEast,by=1)
t2_SouthEast<-t_SouthEast[!is.na(SouthEast_data$sero)]

#Midlands
sero_Midlands<-round(P0_Midlands*EastMidlands_data$sero[!is.na(EastMidlands_data$sero)])
daily_death_Midlands<-WestMidlands_data$daily_death_WestMidlands[1:WestMidlands_death_end]+EastMidlands_data$daily_death_EastMidlands[1:EastMidlands_death_end]
n_days_Midlands<-length(daily_death_Midlands)
Epidemic_Stage_Midlands<-Epidemic_Stage_Midlands[1:n_days_Midlands]
cumul_death_Midlands<-cumsum(daily_death_Midlands)
n_days2_Midlands<-length(sero_Midlands)
t_Midlands<-seq(1,n_days_Midlands,by=1)
t2_Midlands<-t_Midlands[!is.na(EastMidlands_data$sero)]

#EastofEngland
sero_EastofEngland<-round(P0_EastofEngland*EastofEngland_data$sero[!is.na(EastofEngland_data$sero)])
daily_death_EastofEngland<-EastofEngland_data$daily_death[1:EastofEngland_death_end]
n_days_EastofEngland<-length(daily_death_EastofEngland)
Epidemic_Stage_EastofEngland<-Epidemic_Stage_EastofEngland[1:n_days_EastofEngland]
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
                   
                   Epidemic_Stage_SouthWest=Epidemic_Stage_SouthWest,Epidemic_Stage_NorthWest=Epidemic_Stage_NorthWest,Epidemic_Stage_London=Epidemic_Stage_London ,
                   Epidemic_Stage_SouthEast = Epidemic_Stage_SouthEast,Epidemic_Stage_NorthEastYorkshireHumber = Epidemic_Stage_NorthEastYorkshireHumber,Epidemic_Stage_Midlands = Epidemic_Stage_Midlands ,Epidemic_Stage_EastofEngland = Epidemic_Stage_EastofEngland ,
                   
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

SeroModelTimeVaryingIFR<-stan(file=paste(folder,"/SeroEstimationTimeVaryingIFR_WeibullPrior.stan", sep = ""),
                           data = data_sir,
                           iter = niter,
                           chains = chains)

print(SeroModelTimeVaryingIFR)

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


quantile(beta, probs = c(0.025,0.5,0.975))
quantile(eta_NorthEastYorkshireHumber, probs = c(0.025,0.5,0.975))
quantile(gamma_NorthEastYorkshireHumber, probs = c(0.025,0.5,0.975))
quantile(eta_London, probs = c(0.025,0.5,0.975))
quantile(gamma_London, probs = c(0.025,0.5,0.975))
quantile(eta_NorthWest, probs = c(0.025,0.5,0.975))
quantile(gamma_NorthWest, probs = c(0.025,0.5,0.975))
quantile(eta_SouthEast, probs = c(0.025,0.5,0.975))
quantile(gamma_SouthEast, probs = c(0.025,0.5,0.975))
quantile(eta_SouthWest, probs = c(0.025,0.5,0.975))
quantile(gamma_SouthWest, probs = c(0.025,0.5,0.975))
quantile(eta_Midlands, probs = c(0.025,0.5,0.975))
quantile(gamma_Midlands, probs = c(0.025,0.5,0.975))
quantile(eta_EastofEngland, probs = c(0.025,0.5,0.975))
quantile(gamma_EastofEngland, probs = c(0.025,0.5,0.975))

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

saveRDS(SeroModelTimeVaryingIFR, file=paste(folder,"/TimevaryingSeroModelTimeVaryingIFR_Weibull.rds", sep = ""))


