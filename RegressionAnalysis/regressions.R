library(ggplot2)
library(gridExtra)
library(viridis)

data<-read.table("regress_data.txt",header = T)

equation = function(x) {
  lm_coef <- list(r2 = round(summary(x)$r.squared, digits = 2));
  lm_eq <- substitute(italic(R)^2~"="~r2,lm_coef)
  as.character(as.expression(lm_eq));                 
}
fit <- lm(IFR ~ Prop_over_60, data = data)
confint(fit)
fit <- lm(IFR ~ Prop_over_45, data = data)
confint(fit)
p11 <- ggplot(data, aes(x=Prop_over_45	, y=IFR)) + geom_point(shape=16) + geom_smooth(method=lm, color="red", fill="#69b3a2",alpha=0.25) +
  scale_x_continuous(name = "Proportion of people over the age of 45",
                     limits = c(0.34,0.5)) +
  scale_y_continuous(name = "Estimated IFR") +
  geom_text(
    label=data$Regions, 
    nudge_x = 0, nudge_y = 0.0002, 
    check_overlap = T)+
  annotate("text", x = 0.375, y = 0.014, label = equation(fit), parse = TRUE) +
  theme(axis.line.x = element_line(size=.5, colour = "black"),
        axis.line.y = element_line(size=.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12),
        text = element_text(size = 14),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())
p11

fit2 <- lm(IFR ~ rate_CH_deaths, data = data)
confint(fit2)
p12 <- ggplot(data, aes(x=rate_CH_deaths	, y=IFR)) + geom_point(shape=16) + geom_smooth(method=lm, color="red", fill="#69b3a2",alpha=0.25) +
  scale_x_continuous(name = "Deaths in the community relative to \n deaths in care homes",
                     limits = c(2.4,10)) +
  scale_y_continuous(name = "Estimated IFR") +
  geom_text(
    label=data$Regions, 
    nudge_x = 0, nudge_y = 0.0002, 
    check_overlap = T)+
  annotate("text", x = 2.86, y = 0.014, label = equation(fit2), parse = TRUE) +
  theme(axis.line.x = element_line(size=.5, colour = "black"),
        axis.line.y = element_line(size=.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12),
        text = element_text(size = 14),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())
p12


fit3 <- lm(IFR ~ diabetes, data = data)
confint(fit3)
p13 <- ggplot(data, aes(x=diabetes	, y=IFR)) + geom_point(shape=16) + geom_smooth(method=lm, color="red", fill="#69b3a2",alpha=0.25) +
  scale_x_continuous(name = "Diabetes prevalence",
                     limits = c(6.35,7.8)) +
  scale_y_continuous(name = "Estimated IFR") +
  geom_text(
    label=data$Regions, 
    nudge_x = 0, nudge_y = 0.0002, 
    check_overlap = T)+
  annotate("text", x = 6.5, y = 0.014, label = equation(fit3), parse = TRUE) +
  theme(axis.line.x = element_line(size=.5, colour = "black"),
        axis.line.y = element_line(size=.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12),
        text = element_text(size = 14),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())
p13

fit4 <- lm(IFR ~ chronic_liver_disease, data = data)
confint(fit4)
p14 <- ggplot(data, aes(x=chronic_liver_disease	, y=IFR)) + geom_point(shape=16) + geom_smooth(method=lm, color="red", fill="#69b3a2",alpha=0.25) +
  scale_x_continuous(name = "Chronic liver disease \n mortality (per 100,000)",
                     limits = c(9.8,17.5)) +
  scale_y_continuous(name = "Estimated IFR") +
  geom_text(
    label=data$Regions, 
    nudge_x = 0, nudge_y = 0.0002, 
    check_overlap = T)+
  annotate("text", x = 10.5, y = 0.014, label = equation(fit4), parse = TRUE) +
  theme(axis.line.x = element_line(size=.5, colour = "black"),
        axis.line.y = element_line(size=.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12),
        text = element_text(size = 14),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())
p14

fit5 <- lm(IFR ~ chronic_obstructive_pulmonary_disease, data = data)
confint(fit5)
p15 <- ggplot(data, aes(x=chronic_obstructive_pulmonary_disease	, y=IFR)) + geom_point(shape=16) + geom_smooth(method=lm, color="red", fill="#69b3a2",alpha=0.25) +
  scale_x_continuous(name = "Chronic obstructive pulmonary disease \n mortality (per 100,000)",
                     limits = c(38,65.5)) +
  scale_y_continuous(name = "Estimated IFR") +
  geom_text(
    label=data$Regions, 
    nudge_x = 0, nudge_y = 0.0002, 
    check_overlap = T)+
  annotate("text", x = 41.5, y = 0.014, label = equation(fit5), parse = TRUE) +
  theme(axis.line.x = element_line(size=.5, colour = "black"),
        axis.line.y = element_line(size=.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12),
        text = element_text(size = 14),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())
p15

fit6 <- lm(IFR ~ CH_beds_per_100_ppl_over_75, data = data)
confint(fit6)
p16 <- ggplot(data, aes(x=CH_beds_per_100_ppl_over_75	, y=IFR)) + geom_point(shape=16) + geom_smooth(method=lm, color="red", fill="#69b3a2",alpha=0.25) +
  scale_x_continuous(name = "Care Home beds per 100 people over the age of 75", 
                     limits = c(6.85,10.75))+
  scale_y_continuous(name = "Estimated IFR") +
  geom_text(
    label=data$Regions, 
    nudge_x = 0, nudge_y = 0.0002, 
    check_overlap = T)+
  annotate("text", x = 7.25, y = 0.014, label = equation(fit6), parse = TRUE) +
  theme(axis.line.x = element_line(size=.5, colour = "black"),
        axis.line.y = element_line(size=.5, colour = "black"),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12),
        text = element_text(size = 14),
        panel.grid.major = element_line(colour = "#d3d3d3"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())
p16

tiff("regress.tiff", width = 40, height = 23, units = 'cm', res = 210)
gt <- arrangeGrob(p11, p12, p16, p13, p14, p15, ncol = 3)
plot(gt)
dev.off()



waves<-read.table("waves.txt",header = T)
attributes(waves)
waves$Age<-factor(waves$Age, ordered=T,levels = c("45_49", "50_54", "55_59", "60_64", "65_69","70_74",
                              "75_79","80_84","85_89","90+","0_59","60+"))

p<-ggplot(waves, aes(fill=Regions, y=Ratio, x=Age)) + 
  geom_bar(stat="identity",width = 0.65,position = position_dodge(0.8))

myPal <- colorRampPalette(c("midnightblue","deepskyblue3",
                            "lightskyblue1","cyan4","darkolivegreen2",
                            "khaki2","sandybrown","red3","rosybrown2"))
tiff("waves2.tiff", width = 24, height = 16, units = 'cm', res = 300)
p+scale_fill_manual(values=myPal(9))+theme_light()+
  theme(
    text = element_text(size=16),
    # axis.title.y = element_text( size=18, vjust = 1.95),
    plot.title = element_text(face = "bold", size = 16),
    legend.background = element_rect(fill = "white", size = 1.15, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(0.15,0.99),
    legend.title = element_blank(),
    axis.ticks = element_line(colour = "grey90", size = 0.2),
    panel.grid.major = element_line(colour = "grey90", size = 0.5),
    panel.grid.minor = element_blank()
  )  
dev.off()


p+scale_fill_brewer(palette="Spectral")+theme_bw()
p+scale_fill_viridis(discrete = T)+theme_bw()

