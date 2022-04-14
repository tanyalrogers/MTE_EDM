#Computes dominant cycle period for constant temperature experiments

#Inputs: "data/constant_temp_expts/*.csv" (4 files)
#Outputs: "data/cycle_period_summary.csv"

library(dplyr)
library(ggplot2)
library(cowplot)

# Datasets ----------------------------------------------------------------

halbach <- read.csv("data/constant_temp_expts/Halbach_1984.csv",stringsAsFactors=F)
head(halbach)
unique(halbach$Population..temp.replicate.)

fussman<-read.csv("data/constant_temp_expts/Fussman_2014.csv",stringsAsFactors=F)
head(fussman)
unique(fussman$treatment)

delong<-read.csv("data/constant_temp_expts/Delong_2020.csv",stringsAsFactors=F)
head(delong)
str(delong)

laughton<-read.csv("data/constant_temp_expts/Laughton_Knell_2019.csv")
head(laughton)
unique(laughton$box)

#Plot the data
halbach_abund_plot <- halbach %>% 
  ggplot(aes(x = Time, y = Abundance, linetype=Population..temp.replicate.)) + 
  geom_point(size=1) + geom_line() + scale_linetype_manual(values=rep(c(1,2),3)) +
  facet_wrap(~Temperature_C) + theme(legend.position="none")

fussman %>% 
  ggplot(aes(x=time, y=log(dens),linetype=group)) + 
  geom_point() + geom_line() + facet_wrap(~temp+treatment)
#fussman %>% filter(treatment=="B") %>% ggplot(aes(x=time, y=log(dens),linetype=group)) + geom_point() + geom_line() + facet_wrap(~temp)
#fussman %>% filter(treatment=="C") %>% ggplot(aes(x=time, y=log(dens),linetype=group)) + geom_point() + geom_line() + facet_wrap(~temp)
#fussman %>% filter(treatment=="D") %>% ggplot(aes(x=time, y=log(dens),linetype=group)) + geom_point() + geom_line() + facet_wrap(~temp)

delong_abund_plot <- delong %>% filter(!is.na(Density)) %>% 
  ggplot(aes(x=Day, y = Density, linetype=Species)) + 
  geom_point() + geom_line() + facet_wrap(~Temperature.C) + 
  theme(legend.position=c(0.85,0.9))

laughton %>% 
  ggplot(aes(x=week,y=total,group=box)) + 
  geom_point() + geom_line() + facet_wrap(food~temp)

#Function to compute Fourier spectrum

getspectrum=function(dataframe, ts, lam=0.01, scale=T, trim=T) {
  
  xts=dataframe[,ts]
  
  #approximate spectrum calculated using ridge regression
  if(all(is.na(xts))) {
    return(data.frame(Frequency=NA, Period=NA, Power=NA))
  }
  
  #standardize data
  if(scale==T) {
    xtss=(xts-mean(xts, na.rm=T))/sd(xts, na.rm=T)
  } else {
    xtss=xts-mean(xts, na.rm=T)
  }

  #get spectrum
  if(trim) {  # ts length, exclude NAs on ends
    tsl=length(which(!is.na(xtss))[1]:which(!is.na(xtss))[length(which(!is.na(xtss)))])
  } else { # do not trim NAs on ends
    tsl=length(xtss)
  }
  tsle=floor(tsl/2)*2 #ts length rounded down to even number (if odd)
  fi=(1:(tsle/2))/tsle #frequencies (cycles per timestep)
  per=1/fi #periods (timesteps per cycle)
  wi=2*pi*fi #frequencies in radians per timestep
  
  times=1:length(xtss)
  cosbf=cos(outer(times,wi))
  sinbf=sin(outer(times,wi))
  allbf=cbind(cosbf,sinbf) #all basis functions
  
  y=xtss[complete.cases(xtss)] #remove missing timepoints
  X=allbf[complete.cases(xtss),] #remove missing timepoints
  coefsdir=solve(t(X)%*%X + lam*diag(ncol(X)))%*%t(X)%*%y
  lmr=sqrt(coefsdir[1:(length(coefsdir)/2),]^2+coefsdir[(length(coefsdir)/2+1):length(coefsdir),]^2)
  
  return(data.frame(Frequency=fi, Period=per, Power=lmr))
}

# Compute spectra ---------------------------------------------------------

#### Halbach data ---------------------------------------------------------

halbach_spectra <- data.frame(Frequency=numeric(), Period=numeric(), 
                              Power=numeric(), Temperature_C=numeric(),
                              Replicate=character())

for(i in 1:length(unique(halbach$Population..temp.replicate.))){
	
	data <- filter(halbach, Population..temp.replicate. == unique(halbach$Population..temp.replicate.)[i])
	
	spectrum <- getspectrum(data, 3, lam=0.005, scale=T, trim=T)
	spectrum$Temperature_C <- rep(data$Temperature_C[1],nrow(spectrum))
	spectrum$Replicate <- rep(data$Population..temp.replicate.[1],nrow(spectrum))
	
	halbach_spectra <- rbind(halbach_spectra, spectrum)
	
}
halbach_spectra <- halbach_spectra %>% 
  mutate(arhen_temp = 1/(8.62e-5*(Temperature_C+273.15))) %>% 
  group_by(Replicate) %>% 
  filter(Period != max(Period)) #exclude longest (not well resolved) 
head(halbach_spectra)

#Plot all the specta

halbach_spectra_plot <- halbach_spectra %>% 
  ggplot(aes(x = Period, y = Power, linetype=Replicate)) + 
  geom_line() + geom_point(size = 1) + 
  scale_linetype_manual(values=rep(c(1,2),3)) + 
  facet_wrap(~Temperature_C) + theme(legend.position="none")

#Pull out the one with max power 
max_power_halbach <- halbach_spectra %>% 
  group_by(Replicate) %>% filter(Power==max(Power))

halbach_ols <- lm(log(Period)~arhen_temp, data = max_power_halbach)
summary(halbach_ols)
confint(halbach_ols)

halbach_fig <- 
  ggplot(max_power_halbach, aes(x = arhen_temp, y = log(Period))) + 
  geom_smooth(method = lm, se = F, color = "black", linetype=1) + 
  geom_point(aes(color = Replicate), size = 3) +
  xlab("1/kT") + ylab("ln cycle period (days)") + ylim(c(1.6,3.4)) +
  scale_x_continuous(sec.axis=sec_axis(~(.*8.62e-5)^-1-273.15, name = "Temp in C"), limits=c(38.7,40.5)) 

plot_grid(halbach_abund_plot,halbach_spectra_plot, nrow=2)
halbach_fig

#### Delong data ---------------------------------------------------------

delong <- delong %>% 
  mutate(unique_species_temp = paste0(Species,"_", Temperature.C)) %>%
  data.frame() 

delong_spectra <- data.frame(Frequency=numeric(), Period=numeric(), 
                             Power=numeric(), Temperature_C=numeric(),
                             unique_species_temp=character(), Species=character())

for(i in 1:length(unique(delong$unique_species_temp))){
	
	data <- filter(delong, unique_species_temp == unique(delong$unique_species_temp)[i])
	
	spectrum <- getspectrum(data, 4, lam=0.01, scale=T, trim=T)
	spectrum$Temperature.C <- rep(data$Temperature.C[1], nrow(spectrum))
	spectrum$unique_species_temp <- rep(data$unique_species_temp[1], nrow(spectrum))
	spectrum$Species <- rep(data$Species[1], nrow(spectrum))

	delong_spectra <- rbind(delong_spectra, spectrum)
	
}
delong_spectra <- delong_spectra %>% 
  mutate(arhen_temp = 1/(8.62e-5 * (Temperature.C + 273.15)))

delong_spectra_plot <- delong_spectra %>% 
  ggplot(aes(x = Period, y = Power, linetype=Species)) + 
  geom_line() + geom_point(size=1) + theme(legend.position="none") + 
  facet_wrap(~Temperature.C)

plot_grid(delong_abund_plot, delong_spectra_plot,nrow=2)
#Fourier spectrum does not provide a good estimate of cycle period.

#Just find the day when density goes to zero. Use this as the cycle period. 
delong_period_df <-delong %>% 
  group_by(unique_species_temp) %>% 
  filter(Density == 0) %>% slice(1) %>% ungroup() %>% 
  mutate(arhen_temp = 1/(8.62e-5 * (Temperature.C + 273.15)))

delong_fig <- delong_period_df %>% 
  ggplot(aes(x=arhen_temp, y=log(Day), shape = Species)) + 
  geom_point(size=3) + geom_smooth(method="lm",se=F,color="black") + 
  ylim(c(0.5,3.5)) + xlab("1/kT") + ylab("ln cycle period (days)") + 
  scale_x_continuous(sec.axis=sec_axis(~(.*8.62e-5)^-1-273.15, name = "Temp in C"), limits=c(37.5,40.5)) + 
  theme(legend.position = c(0.05,0.85),legend.title=element_blank())

delong_fig

#### Fussman data ---------------------------------------------------------

#compute the mean density at each time point

fussman <- fussman %>% 
  mutate(unique_species_temp_treat = paste0(treatment,"_", temp,"_",group))
fussman_avg <-fussman %>% 
  group_by(unique_species_temp_treat, time) %>% 
  mutate(mean_dens = mean(dens)) %>% 
  slice(1) %>% ungroup() %>% data.frame()
head(fussman_avg)

ggplot(fussman_avg, aes(x=time, y=dens, color=treatment)) +
  facet_wrap(unique_species_temp_treat~., scales="free") +
  geom_line() + geom_point()

fussman_spectra <- data.frame(Frequency=numeric(), Period=numeric(), 
                              Power=numeric(), treatment=character(),temp=numeric(),
                              unique_species_temp=character(), group=character())

for(i in 1:length(unique(fussman_avg$unique_species_temp_treat))){
	
	data <- filter(fussman_avg, unique_species_temp_treat == unique(fussman_avg$unique_species_temp_treat)[i])
	
	spectrum <- getspectrum(data, 4, lam=0.01, scale=T, trim=F)
	spectrum$treatment <- rep(data$treatment[1],nrow(spectrum))
	spectrum$temp <- rep(data$temp[1], nrow(spectrum))
	spectrum$unique_species_temp_treat <- rep(data$unique_species_temp_treat[1], nrow(spectrum))
	spectrum$group <- rep(data$group[1], nrow(spectrum))

	fussman_spectra <- rbind(fussman_spectra, spectrum)
	
}
fussman_spectra <- fussman_spectra %>% 
  filter(!(treatment == "A" & group=="cil")) %>% 
  mutate(arhen_temp = 1/(8.62e-5 * (temp + 273.15)))

fussman_spectra %>% 
  ggplot(aes(x = Period, y = Power)) + 
  geom_line(aes(linetype=group)) + geom_point(size = 1) + 
  facet_wrap(~temp+treatment, ncol=3)

max_fussman_period <- fussman_spectra %>% 
  group_by(unique_species_temp_treat) %>% 
  filter(Period != max(Period)) %>% 
  filter(Power == max(Power))

fussman_cil_periods <- max_fussman_period %>% filter(group == "cil") %>% 
  ggplot(aes(x = arhen_temp, y = log(Period),  shape = treatment)) + 
  geom_point(size = 3) + ggtitle("ciliates (Tetra)") + theme(legend.position="none")
fussman_bact_periods <- max_fussman_period %>% filter(group == "bact") %>% 
  ggplot(aes(x = arhen_temp, y = log(Period),  shape = treatment)) + 
  geom_point(size = 3) + ggtitle("bacteria (Pseudo)") + theme(legend.position="none")

plot_grid(fussman_cil_periods, fussman_bact_periods)

cil_ols=lm(log(Period)~arhen_temp, filter(max_fussman_period, group == "cil"))
summary(cil_ols)
confint(cil_ols)
bact_ols=lm(log(Period)~arhen_temp, filter(max_fussman_period, group == "bact"))
summary(bact_ols)
confint(bact_ols)

#### Laughton data ---------------------------------------------------------

laughton <- laughton %>%  
  mutate(unique_box_food_temp = paste0(box,"_", food,"_",temp)) %>% data.frame()
laughton_spectra <- data.frame(Frequency=numeric(), Period=numeric(), 
                               Power=numeric(), Temperature_C=numeric(),
                               unique_box_food_temp=character(), food=character(), box=numeric())

for(i in 1:length(unique(laughton$unique_box_food_temp))){
	
	data <- filter(laughton, unique_box_food_temp == unique(laughton$unique_box_food_temp)[i])
	
	spectrum <- getspectrum(data, 9, lam=0.005, scale=T, trim=T)
	spectrum$Temperature_C <- rep(data$temp[1],nrow(spectrum))
	spectrum$unique_box_food_temp <- rep(data$unique_box_food_temp[1],nrow(spectrum))
	spectrum$food<-rep(data$food[1],nrow(spectrum))
	spectrum$box<-rep(data$box[1],nrow(spectrum))

	laughton_spectra <- rbind(laughton_spectra, spectrum)
	
}
laughton_spectra <- laughton_spectra %>% 
  mutate(arhen_temp = 1/(8.62e-5 * (Temperature_C + 273.15))) %>%
  group_by(unique_box_food_temp) %>% filter(Period!=max(Period))
head(laughton_spectra)

laughton_spectra %>% 
  ggplot(aes(x = Period, y = Power,group=box)) + 
  geom_line() + geom_point(size = 1) + 
  facet_wrap(food~Temperature_C)

laughton %>% 
  ggplot(aes(x=week,y=total,group=box)) + 
  geom_point() + geom_line() + 
  facet_wrap(food~temp)

max_power_laughton = laughton_spectra %>% 
  group_by(unique_box_food_temp) %>% filter(Power==max(Power))
 
ggplot(max_power_laughton, aes(x=arhen_temp, y = log(Period))) + 
  geom_point(size=3) + facet_wrap(~food,nrow=1)+theme_bw()

bad_ols=lm(log(Period)~arhen_temp, filter(max_power_laughton, food=="bad"))
summary(bad_ols)
confint(bad_ols)
good_ols=lm(log(Period)~arhen_temp, filter(max_power_laughton, food=="good"))
summary(good_ols)
confint(good_ols)

# Regressions -------------------------------------------------------------

didinium_ols <- lm(log(Day)~arhen_temp, data = filter(delong_period_df,Species=="Didinium"))
summary(didinium_ols)

paramecium_ols <- lm(log(Day)~arhen_temp, data = filter(delong_period_df,Species=="Paramecium"))
summary(paramecium_ols)

halbach_ols <- lm(log(Period)~arhen_temp, data = max_power_halbach)
summary(halbach_ols)

# Export ------------------------------------------------------------------

delong_export <- delong_period_df %>% 
  select(Species, Temperature_C=Temperature.C, arhen_temp, Period=Day) %>% 
  mutate(Ref="Delong et al. 2020", Replicate = NA)
halbach_export <- max_power_halbach %>% 
  select(Temperature_C, arhen_temp, Period, Replicate) %>% 
  mutate(Ref="Halbach et al. 1984", Species="Rotifer")
exportdf <- bind_rows(delong_export,halbach_export) %>% as.data.frame()
write.csv(exportdf, "data/cycle_period_summary.csv", row.names = F)
