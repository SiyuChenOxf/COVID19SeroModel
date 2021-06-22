library(hrbrthemes) 
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

delta_h<-12        #fixed time lag between hospitalisation and death
delta_p<-14        #fixed time lag between test for virus and death(seroversion)

Maver<-function(x,step) {                                      #This function is to calculate moving average for case, hospitalisation and virus positivity data
  if(step%%2 ==0) {print("please give an odd step"); break}    #step is moving length for averaging
  else if(step==1){return(x)}  
  else {
    M<-rep(0,(length(x)-3))
    for (i in ((step-1)/2+1):(length(x)-(step-1)/2)) {
      M[i]<-1/step*(x[i-1]+x[i-2]+x[i-3]+x[i]+x[i+1]+x[i+2]+x[i+3])
    }
    return(M[-c(1:((step-1)/2))])
  }
}

step<-7

England_data<-read.csv("England_data.csv",header = TRUE)

positivity_England_Pillar1<-read.csv("positivity_England_Pillar1_week39.csv",header = TRUE)
positivity_England_Pillar2<-read.csv("positivity_England_Pillar2_week39.csv",header = TRUE)
positivity_England_Pillar1_2<-read.csv("positivity_England_Pillar1_week44.csv",header = TRUE)
positivity_England_Pillar2_2<-read.csv("positivity_England_Pillar2_week44.csv",header = TRUE)


#virus positivity rates in England
positivity_England<-matrix(c(positivity_England_Pillar1$positivity[1:14],0.5*positivity_England_Pillar1$positivity[15:22]+0.5*positivity_England_Pillar2$positivity[15:22],0.25*positivity_England_Pillar1$positivity[23:35]+0.25*positivity_England_Pillar2$positivity[23:35]+0.25*positivity_England_Pillar1_2$positivity[1:13]+0.25*positivity_England_Pillar2_2$positivity[1:13],0.5*positivity_England_Pillar1_2$positivity[14:18]+0.5*positivity_England_Pillar2_2$positivity[14:18]),nrow = 1)
positivity_England<-c(rep(0,26),as.vector(apply(positivity_England, 2,function(x) rep(x,7))))/100

######PLOTTING START HERE#####
CFR<-Maver(England_data$daily_death[(delta_p+1):length(England_data$daily_death)],step)/Maver(England_data$daily_case[1:(length(England_data$daily_case)-delta_p)],step)
HFR<-Maver(England_data$daily_death[(delta_h+1):length(England_data$daily_death)],step)/Maver(England_data$hospitalization[1:(length(England_data$hospitalization)-delta_h)],step)

normalised_CFR<-c(CFR/max(CFR[!is.na(CFR)]))
normalised_HFR<-c(HFR/max(HFR[!is.na(HFR)]))

data1<-data.frame(t=as.Date(England_data$date), smoothed_case_per_100k= England_data$daily_case*100000/56e6,smoothed_death_per_100k=England_data$daily_death*100000/56e6)

p1<-ggplot(data1, aes(x =t)) +
  geom_line(aes(y = smoothed_case_per_100k, color = "case_per_100k"), size = 2) +
  ggtitle("Case rates") 

p2<-ggplot(data1, aes(x =t)) +
  geom_line(aes(y = smoothed_death_per_100k, color = "death_per_10k"), size = 2) +
  ggtitle("Death rates") 
# p1+p2
cbbPalette <- c( "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

coeff <- 30
p3<-ggplot(data1, aes(x=t)) +
  geom_line( aes(y=smoothed_death_per_100k,color = "Deaths"),size=1.5) +  
  geom_line( aes(y=smoothed_case_per_100k/coeff,color = "Cases"),size=1.5) + 
  scale_y_continuous(
    name = "Daily Deaths per 100,000",
    sec.axis = sec_axis(~.*coeff, name="Daily Cases per 100,000")
  )+ 
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  theme_light() +
  theme(
    axis.title.y = element_text( size=14, vjust = 1.95),
    axis.title.y.right = element_text(size=14, vjust = 1.95)
  )+
  scale_fill_manual( values=cbbPalette)+
  scale_colour_manual( values=cbbPalette)+
  xlab("  ")+
  theme(
    text = element_text(size=15),
    plot.title = element_text(face = "bold", size = 20),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(0.65,0.95),
    legend.title = element_blank(),
    axis.ticks = element_line(colour = "grey90", size = 0.2),
    panel.grid.major = element_line(colour = "grey90", size = 0.2),
    panel.grid.minor = element_blank()
  )  

l<-length(as.Date(England_data$date)[((step-1)/2+1):(length(normalised_CFR)+(step-1)/2)])
data2<-data.frame(t=as.Date(England_data$date)[((step-1)/2+1):(length(normalised_CFR)+(step-1)/2)],
                  normalised_smoothed_virus_positivity=(positivity_England/max(positivity_England))[1:l],
                  normalised_smoothed_CFR=normalised_CFR[1:l],
                  normalised_smoothed_HFR=normalised_HFR[1:l])
p4<-ggplot(data2, aes(x =t)) +
  geom_line(aes(y = normalised_smoothed_virus_positivity, color = "Normalised Positivity Ratio"), size = 1.5) +
  geom_line(aes(y = normalised_smoothed_CFR, color = "Normalised CFR"), size = 1.5) +
  geom_line(aes(y = normalised_smoothed_HFR, color = "Normalised HFR"), size = 1.5) +
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01",  "2020-11-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov"), limit = as.Date(c("2020-01-01","2020-11-07")))+
  scale_fill_manual( values=cbbPalette)+
  scale_colour_manual( values=cbbPalette)+
  theme_light() +
  ylab("") +
  xlab("  ")+
  ylim(0, 1) + 
  theme(
    text = element_text(size=15),
    plot.title = element_text(face = "bold", size = 20),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(0.5,0.95),
    legend.title = element_blank(),
    axis.ticks = element_line(colour = "grey90", size = 0.2),
    panel.grid.major = element_line(colour = "grey90", size = 0.2),
    panel.grid.minor = element_blank()
  )

tiff("CFRHFR.tiff", width = 20, height = 22,units = "cm", res=300)
cowplot::plot_grid(p3, NULL, p4, align = "v", ncol = 1, rel_heights = c(0.4,0.04,0.4))
dev.off()
