#Generates figures for manuscript

#Inputs:"data/cycle_period_summary.csv";
#       "data/MTE_EDM_metadata.csv";
#       "data/MTE_EDM_output.csv";
#       "data/TPC_EDM_output.csv";
#       "data/constant_temp_expts/Halbach_1984.csv";
#       "data/constant_temp_expts/Delong_2020.csv"
#Outputs:"figures/*"

library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)
source("code/ggplot_themes.R")

# Cycle periods -----------------------------------------------------------

periods<-read.csv("data/cycle_period_summary.csv",stringsAsFactor=F)

#fit lms
rotifer_lm<-lm(log(Period)~arhen_temp,data=filter(periods, Species=="Rotifer"))
summary(rotifer_lm)
didinium_lm<-lm(log(Period)~arhen_temp,data=filter(periods, Species=="Didinium"))
summary(didinium_lm)
paramecium_lm<-lm(log(Period)~arhen_temp,data=filter(periods, Species=="Paramecium"))
summary(paramecium_lm)

ggplot(data=periods, aes(x = arhen_temp, y = log(Period))) + 
  geom_point(aes(shape=Species), size = 3,stroke=1.1) + 
  geom_smooth(aes(group=Species,linetype=Ref),method="lm", se = F,color="Black",show.legend = F)  + 
  scale_x_continuous(sec.axis=sec_axis(~(.*8.62e-5)^-1 - 273.15, name = "Temperature (°C)"), limits=c(38,40.5)) + 
  scale_y_continuous(limits=c(0.75,3.25)) + 
  scale_shape_manual(values = c(16,17,1)) + 
  xlab("1/(kT)") + ylab("ln cycle period (days)") +
  classic +
  theme(legend.justification = c(1,0), legend.position = c(1,0), #legend.title=element_blank(), 
        legend.background = element_rect(fill=NA,color="black"), legend.key.height = unit(0.2,"in"))
ggsave("figures/cycleperiods.pdf",width = 4, height = 3.5, useDingbats=FALSE)

# Time series with rescaled time
# Time is scaled by the Arrhenius factor, with the MTE activation energy of 0.65.

halbach <- read.csv("data/constant_temp_expts/Halbach_1984.csv",stringsAsFactors=F)
delong<-read.csv("data/constant_temp_expts/Delong_2020.csv",stringsAsFactors=F)

colors1=c("blue2","gold","red2")
colors2=c("blue4", "blue2","gold3","gold","red2","red4")

#Note - in the Halbach paper, there are two replicates for each temperature.
#Halbach calendar time:
calendar<-halbach %>% 
  ggplot(aes(x = Time, y = Abundance, color = factor(Temperature_C))) + 
  geom_line(aes(linetype=substr(Population..temp.replicate.,4,5))) + 
  geom_point() + xlab("Calendar time (days)") + 
  scale_color_manual(values = colors1, name=NULL) + 
  theme_bw() + 
  theme(legend.position=c(1,1),legend.justification = c(1,1),
        plot.title = element_text(size=12),
        legend.background = element_blank(),
        legend.margin = margin(-5, 5, 0, 0),
        legend.spacing.x = unit(0, "pt")) + 
  ggtitle("Rotifers (Halback et al. 1984)") + guides(linetype="none")
#Halbach biological time
bio<-halbach  %>% mutate(bio_time = Time/exp(0.65/(8.62e-5*(Temperature_C+273.15)))) %>% 
  ggplot(aes(x = bio_time, y = Abundance, color = factor(Temperature_C))) + 
  geom_line(aes(linetype=substr(Population..temp.replicate.,4,5))) + 
  geom_point() + xlab("Biological time (days exp(-0.65/kT))") + 
  scale_color_manual(values = colors1) + 
  theme_bw() +
  theme(legend.position="none") + guides(linetype="none")

#Note - the Delong paper had time series for both the predator and prey. 
#These are given separate plots.
#Delong Didinium calendar time:
delong_calendar<-delong %>% filter(Species == "Didinium",!is.na(Density)) %>% 
  ggplot(aes(x = Day, y = Density, color = factor(Temperature.C))) + 
  geom_line()+ geom_point() + xlab("Calendar time (days)") + 
  #scale_color_brewer(palette = "RdYlBu") +
  theme_bw() + 
  scale_color_manual(name=NULL,values = colors2) + 
  theme(legend.position=c(0,1),legend.justification = c(0,1), 
        plot.title = element_text(size=12),
        legend.background = element_blank(),
        legend.margin = margin(-5, 0, 0, 5),
        legend.spacing.x = unit(0, "pt")) + 
  ggtitle("Didinium (Delong et al. 2020)")

#Delong Didinium biological time:
delong_bio<-delong %>% filter(Species == "Didinium",!is.na(Density))%>% 
  mutate(bio_time = Day/exp(0.65/(8.62e-5*(Temperature.C+273.15)))) %>% 
  ggplot(aes(x = bio_time, y = Density, color = factor(Temperature.C))) + 
  geom_line()+ geom_point() + xlab("Biological time (days exp(-0.65/kT))") + 
  scale_color_manual(values = colors2) + 
  theme_bw() + 
  theme(legend.position="none")

#Delong Paramecium calendar time:
delong_paramecium_calendar<-delong %>% 
  filter(Species == "Paramecium",!is.na(Density)) %>% 
  ggplot(aes(x = Day, y = Density, color = factor(Temperature.C))) + 
  geom_line()+ geom_point() + xlab("Calendar time (days)") + 
  scale_color_manual(values = colors2, name = NULL) + 
  theme_bw() + 
  theme(legend.position=c(1,1),legend.justification = c(1,1),
        plot.title = element_text(size=12),
        legend.background = element_blank(),
        legend.margin = margin(-5, 5, 0, 0),
        legend.spacing.x = unit(0, "pt")) + 
  ggtitle("Paramecium (Delong et al. 2020)")

#Delong Paramecium biological time:
delong_paramecium_bio<-delong %>% filter(Species == "Paramecium",!is.na(Density))%>% 
  mutate(bio_time = Day/exp(0.65/(8.62e-5*(Temperature.C+273.15)))) %>% 
  ggplot(aes(x = bio_time, y = Density, color = factor(Temperature.C))) + 
  geom_line()+ geom_point() + xlab("Biological time (days exp(-0.65/kT))") + 
  scale_color_manual(values = colors2) +
  theme_bw() + 
  theme(legend.position="none")

plot_grid(calendar, delong_calendar, delong_paramecium_calendar, bio, 
          delong_bio, delong_paramecium_bio, rel_heights = c(1,0.95))
ggsave("./figures/experimental_ts_rescaling.png", width = 10, height = 6.5)

# Main results ------------------------------------------------------------

meta=read.csv("data/MTE_EDM_metadata.csv",stringsAsFactors = F)
res=read.csv("data/MTE_EDM_output_Simplex_UTD_cov.csv",stringsAsFactors = F)
res2=read.csv("data/MTE_EDM_output_MTE.csv",stringsAsFactors = F)
restpc=read.csv("data/TPC_EDM_output.csv",stringsAsFactors = F)

res=meta %>% left_join(res) %>% left_join(res2) %>% left_join(restpc)

colnames(res)
res$SpeciesAb=factor(res$SpeciesAb,levels=unique(res$SpeciesAb))
res$Location=factor(res$Location,levels=unique(res$Location))

#absolute R2
reslong=gather(res,Model,R2,R2_Simplex,R2_UTD,R2_MTE,R2_TPC,R2_cov)
reslong$Model=sub("R2_","", reslong$Model)
reslong$Model=sub("cov","Covariate", reslong$Model)
reslong$SpeciesAb=factor(reslong$SpeciesAb,levels=unique(reslong$SpeciesAb))
reslong$Model=factor(reslong$Model, levels=c("Simplex","UTD","MTE","Covariate","TPC"))

#change in R2
res$delta_R2_MTE=res$R2_MTE-res$R2_Simplex
res$delta_R2_UTD=res$R2_UTD-res$R2_Simplex
res$delta_R2_Covariate=res$R2_cov-res$R2_Simplex
res$delta_R2_TPC=res$R2_TPC-res$R2_Simplex
res$delta_R2_TPC2=res$R2_TPC-res$R2_MTE

reslong2=gather(res,Model,deltaR2,delta_R2_MTE,delta_R2_UTD,delta_R2_Covariate,delta_R2_TPC)
reslong2$Model=sub("delta_R2_","", reslong2$Model)
reslong2$SpeciesAb=factor(reslong2$SpeciesAb,levels=unique(reslong2$SpeciesAb))
reslong2$Model=factor(reslong2$Model, levels=c("UTD","MTE","Covariate","TPC"))

#likelihood ratios
res$LR_MTE_Simplex=res$logL_MTE-res$logL_Simplex
res$LR_UTD_Simplex=res$logL_UTD-res$logL_Simplex
res$LR_Covariate_Simplex=res$logL_cov-res$logL_Simplex
res$LR_TPC_Simplex=res$logL_TPC-res$logL_Simplex

reslong3=gather(res,Model,LR,LR_MTE_Simplex,LR_UTD_Simplex,LR_Covariate_Simplex,LR_TPC_Simplex)
reslong3$Model=sub("LR_","", reslong3$Model)
reslong3$Model=sub("_Simplex","", reslong3$Model)
reslong3$SpeciesAb=factor(reslong3$SpeciesAb,levels=unique(reslong3$SpeciesAb))
reslong3$Model=factor(reslong3$Model, levels=c("UTD","MTE","Covariate","TPC"))

#L ratio test
res$LR_pval=pchisq(2*res$LR_MTE_Simplex,1,lower.tail = F)

#R2 ranges and means
range(res$R2_Simplex)
mean(res$R2_Simplex)
range(res$R2_MTE)
mean(res$R2_MTE)

#performance improvement (percent)
resecto=filter(res, Site!="Portal") #ectotherms only
length(ectoimp <- which(resecto$delta_R2_MTE>0)) #ectotherms with improvement
mean(resecto$delta_R2_MTE[ectoimp]/resecto$R2_Simplex[ectoimp]*100) # %increase, ones with improvement
mean(resecto$delta_R2_MTE/resecto$R2_Simplex*100) # %increase, all
# % improvement for endotherms
mean(res$delta_R2_MTE[20:22]/res$R2_Simplex[20:22]*100)
# % improvement using tpc vs mte
mean(resecto$delta_R2_TPC2/resecto$R2_MTE*100, na.rm=T)

## Model comparison barplots

#absolute R2
ggplot(reslong, aes(x=SpeciesAb,y=R2, fill=Model)) +
  geom_col(position = position_dodge()) +
  #geom_point(aes(color=Location,y=1),size=5) +
  classic + xlabvert + scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  geom_vline(xintercept = seq(0.5, 23, by = 1), color="gray", size=.5, alpha=.5) +
  scale_fill_brewer(palette = "Paired") +
  labs(y=expression(R^2),x=NULL)
#ggsave("./figures/R2bar_species.png", width = 8, height = 5)

#change in R2
colvals=RColorBrewer::brewer.pal(5,"Paired")
R2plot=ggplot(reslong2, aes(x=SpeciesAb,y=deltaR2, fill=Model)) +
  geom_col(position = position_dodge()) +
  #geom_point(aes(color=Location,y=1),size=5) +
  classic + xlabvert +
  geom_vline(xintercept = seq(0.5, 23, by = 1), color="gray", size=.5, alpha=.5) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=colvals[-1]) +
  labs(y=expression(Delta~R^2~relative~to~Simplex),x=NULL) +
  theme(legend.background = element_rect(fill="white", color = "black"),
        legend.justification = c(1,1),legend.position = c(1,1), legend.title = element_text(size=10),
        legend.text = element_text(size=9),
        legend.key.size = unit(0.15,"in"), axis.title.y.left = element_text(hjust = -0.05, vjust = -0.5))
#ggsave("./figures/delta_R2bar_species.png", R2plot, width = 8, height = 5)

#LR
colvals=RColorBrewer::brewer.pal(5,"Paired")
LRplot=ggplot(reslong3, aes(x=SpeciesAb,y=LR, fill=Model)) +
  geom_col(position = position_dodge()) +
  #geom_point(aes(color=Location,y=1),size=5) +
  classic + xlabvert +
  geom_vline(xintercept = seq(0.5, 23, by = 1), color="gray", size=.5, alpha=.5) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values=colvals[-1]) +
  labs(y="log likelihood ratio\nrelative to Simplex",x=NULL) +
  theme(legend.background = element_rect(fill="white", color = "black"),
        legend.justification = c(1,1),legend.position = c(1,1), legend.title = element_text(size=10),
        legend.text = element_text(size=9),
        legend.key.size = unit(0.15,"in"))

## Activation energies

#points
ggplot(res, aes(x=SpeciesAb,y=OptE0)) +
  geom_hline(yintercept = 0.65, lty=2) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin=OptE0-stdE0, ymax=OptE0+stdE0), width=0) +
  geom_point(size=1.5, aes(color=Location)) +
  classic + xlabvert + 
  geom_vline(xintercept = seq(0.5, 23, by = 1), color="gray", size=.5, alpha=.5) +
  scale_color_brewer(palette = "Dark2") +
  labs(y="Optimized Activation Energy",x=NULL)
ggsave("./figures/E0_species.png", width = 8, height = 5)

#histogram
E0hist=ggplot(res, aes(x=OptE0, fill=Location)) +
  geom_histogram(binwidth = 0.2, boundary=0) +
  classic + 
  geom_vline(xintercept = 0.65, lty=2) +
  scale_fill_brewer(palette = "Dark2") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(x="Optimized Activation Energy",y="Count")

## 4 panel plot

#improvement vs temperature variation
p1=ggplot(filter(res,Location!="Portal, Arizona"), aes(x=std_T,y=delta_R2_MTE,fill=Location)) +
  geom_hline(yintercept = 0) + 
  geom_point(size=3, pch=21, alpha=0.9, color="black") +
  classic + scale_fill_brewer(palette = "Dark2",drop=F) +
  labs(y=expression(R[MTE]^2-R[Simplex]^2), x=expression(SD[Temperature]))

plots=plot_grid(p1 + theme(legend.position="none"), E0hist + theme(legend.position="none"), nrow = 1, align = "h", labels=c("(b)","(c)"), rel_widths = c(1,0.9))
legend <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 3), legend.background = element_blank()))
P3=plot_grid(plots,legend, ncol = 2, rel_widths = c(1,0.3))
#plots2=plot_grid(p0,p0.5, nrow = 1, rel_widths = c(1,0.27),labels=c("a","b"),align = "h")
P4=plot_grid(R2plot,P3,nrow=2,rel_heights = c(1,0.7), labels = "(a)") +
  theme(panel.background=element_rect(fill="white",color = "white"))
ggsave("./figures/multipanel1.png", P4, width = 7, height = 6.5)
ggsave("./figures/Fig3.pdf", P4, width = 7, height = 6.5)

#alternate using LR
p1=ggplot(filter(res,Location!="Portal, Arizona"), aes(x=std_T,y=LR_MTE_Simplex,fill=Location)) +
  geom_hline(yintercept = 0) + 
  geom_point(size=3, pch=21, alpha=0.9, color="black") +
  classic + scale_fill_brewer(palette = "Dark2",drop=F) +
  labs(y="log likelihood ratio\n(MTE/Simplex)", x=expression(SD[Temperature]))

plots=plot_grid(p1 + theme(legend.position="none"), E0hist + theme(legend.position="none"), nrow = 1, align = "h", labels=c("b","c"), rel_widths = c(1,0.9))
legend <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 3), legend.background = element_blank()))
P3=plot_grid(plots,legend, ncol = 2, rel_widths = c(1,0.3))
#plots2=plot_grid(p0,p0.5, nrow = 1, rel_widths = c(1,0.27),labels=c("a","b"),align = "h")
P4=plot_grid(LRplot,P3,nrow=2,rel_heights = c(1,0.7), labels = "a") +
  theme(panel.background=element_rect(fill="white",color = "white"))
ggsave("./figures/multipanel2.png", P4, width = 7, height = 6.5)

# others plots, not used ####

#original 4 panel plot
p0=ggplot(res, aes(x=SpeciesAb,y=delta_R2_MTE,fill=Location, size=R2_MTE)) +
  geom_hline(yintercept = 0) +  
  geom_point(pch=23, alpha=0.9, color="black") +
  geom_point(data=filter(res,Location=="Portal, Arizona"), pch=4,color="black",size=2.5) + 
  classic + xlabvert + scale_fill_brewer(palette = "Dark2") +
  labs(y=expression(R[MTE]^2-R[Simplex]^2),size=expression(R[MTE]^2),x=NULL) +
  guides(fill = "none", size=guide_legend(nrow = 1)) +
  ylim(c(-0.13,0.3)) + scale_size(breaks = seq(0.3,0.8,0.1), limits=c(0.2,0.9)) +
  theme(legend.background = element_rect(fill="gray95", color = "black"),
        legend.justification = c(1,1),legend.position = c(1,1), legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.key.height = unit(0.04,"in"),legend.direction = "horizontal")
p0.5=ggplot(filter(res,!is.na(delta_R2_TPC2)), aes(x=SpeciesAb,y=delta_R2_TPC2, fill=Location)) +
  geom_hline(yintercept = 0) +  
  geom_point(pch=21, alpha=0.9, color="black", size=3, show.legend = F) +
  classic + xlabvert + scale_fill_brewer(palette = "Dark2",drop=F) +
  labs(y=expression(R[TPC]^2-R[MTE]^2),size=expression(R[MTE]),x=NULL) +
  ylim(c(-0.13,0.3)) +
  theme(legend.background = element_rect(fill=NA, color = "black"),
        legend.justification = c(0,0),legend.position = c(0,0))
p1=ggplot(filter(res,Location!="Portal, Arizona"), aes(x=std_T,y=LR_MTE_Simplex,fill=Location)) +
  geom_hline(yintercept = 0) + 
  geom_point(size=3, pch=21, alpha=0.9, color="black") +
  classic + scale_fill_brewer(palette = "Dark2",drop=F) +
  labs(y="log likelihood ratio\n(MTE/Simplex)", x=expression(SD[Temperature]))
p2=ggplot(filter(res,Location!="Portal, Arizona"), aes(x=log10(Mass_mg),y=LR_MTE_Simplex,fill=Location)) +
  geom_hline(yintercept = 0) + 
  geom_point(size=3, pch=21, alpha=0.9, color="black") +
  classic + scale_fill_brewer(palette = "Dark2",drop=F) +
  labs(y="log likelihood ratio\n(MTE/Simplex)", x=expression(log[10]~Mass~(mg)))

plots=plot_grid(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), nrow = 1, labels=c("c","d"))
legend <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 3)))
P3=plot_grid(plots,legend, ncol = 2, rel_widths = c(1,0.3))
plots2=plot_grid(p0,p0.5, nrow = 1, rel_widths = c(1,0.27),labels=c("a","b"),align = "h")
P4=plot_grid(plots2,P3,nrow=2,rel_heights = c(1,0.7))
ggsave("./figures/deltaR2_all.png", P4, width = 7, height = 6.5)

#CV vs mass
ggplot(res, aes(x=log10(Mass_mg),y=cv_T^2,fill=Location)) +
  geom_hline(yintercept = 0) +
  geom_point(size=5, pch=21, alpha=0.9, color="black") +
  classic + scale_fill_brewer(palette = "Dark2")

#relationship with embedding dimension
p1=ggplot(res, aes(x=E_MTE,y=delta_R2_MTE,fill=Location)) +
  geom_hline(yintercept = 0) + 
  geom_point(size=4, pch=21, alpha=0.9, color="black") +
  classic + scale_fill_brewer(palette = "Dark2") +
  labs(y=expression(R[MTE]^2-R[Simplex]^2), x=expression(E[MTE]))
p2=ggplot(res, aes(x=delta_E,y=delta_R2_MTE,fill=Location)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point(size=4, pch=21, alpha=0.9, color="black") +
  classic + scale_fill_brewer(palette = "Dark2") +
  labs(y=expression(R[MTE]^2-R[Simplex]^2), x=expression(E[MTE]-E[Simplex]))
plots=plot_grid(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), nrow = 1, labels=c("A","B"))
legend <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 3)))
plot_grid(plots,legend, ncol = 2, rel_widths = c(1,0.2))
ggsave("./figures/deltaR2_E_deltaE.png", width = 8, height = 3)

#sampling interval
ggplot(res, aes(x=Step_length_days,y=delta_R2_MTE,fill=Location)) +
  geom_hline(yintercept = 0) + 
  geom_point(size=4, pch=21, alpha=0.9, color="black") +
  classic + scale_fill_brewer(palette = "Dark2") +
  labs(y=expression(R[MTE]^2-R[Simplex]^2), x="Step length (days)")
ggsave("./figures/deltaR2_steplength.png", width = 6, height = 4)

