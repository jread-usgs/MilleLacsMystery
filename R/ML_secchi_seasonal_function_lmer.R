
library(ggplot2)
library(dplyr)
library(lme4)
library(reshape2)
library(readr)
library(ggrepel)
library(RColorBrewer)


#data from Heidi
mydata=read.table("ML_secchi.txt", sep='\t', header=TRUE)
mydata$date=as.Date(mydata$date, format="%m/%d/%Y")
mydata$doy=as.numeric(format(mydata$date, format="%j"))

year.medians=summarise(group_by(mydata, year), median.secchi=median(secchi_m))

#exclude data outside window used for TOHA
toha.secchi=subset(mydata, doy>100 & doy<320)
year.medians.toha=summarise(group_by(toha.secchi, year), median.secchi=median(secchi_m))
year.medians.toha$date=as.Date(paste(year.medians.toha$year,"-07-01", sep=""))


#create function for seasonal pattern
toha.secchi$week=as.numeric(format(toha.secchi$date, "%U"))
toha.secchi$year.f=factor(toha.secchi$year)

customLogitModel=function(x,max, min){min+((max-min)/(1+exp(.046*(x-192))))}
#k=0.046, midseason=192 as estimated from nls

customLogitModelGradient <- deriv(
  body(customLogitModel)[[2]], 
  namevec = c("max", "min"), 
  function.arg=customLogitModel
)

#starting guesses
#secchi.nls=nls(secchi_m~customLogitModel(doy, max, min, k, mid.season), data=toha.secchi,start = c(max=4, min=1.5, mid.season=200, k=1))

# Fit the model
model <- nlmer(
  secchi_m ~ customLogitModelGradient(x=doy, max, min) ~ 
    # Random effects with a second ~
    (max | year.f) + (min|year.f) , 
  data = toha.secchi, 
  start = c(max=4, min=2.3)
)

summary(model)
year.estimates=coef(model)$year.f
general.estimates=fixef(model)
write.csv(year.estimates, "year_secchi_max_min_estimates.csv", row.names=T)
write.csv(data.frame(general.estimates), "generic_secchi_max_min_estimates.csv",row.names=F)

#create time series - if year is missing or bad estimate, use generic
yearly.secchi=data.frame(year.estimates )
yearly.secchi$year=row.names(yearly.secchi)

secchi.params=data.frame(year=seq(1974,2016))
secchi.params=merge(secchi.params, yearly.secchi, by="year", all=T)
#bad estimates if min is bigger than max
secchi.params$max[secchi.params$max<secchi.params$min]=NA
secchi.params$min[secchi.params$max<secchi.params$min]=NA
secchi.params$max[is.na(secchi.params$max)]=general.estimates[1]
secchi.params$min[is.na(secchi.params$min)]=general.estimates[2]



#create data frame of secchi depths
day.list=data.frame(date= seq(as.Date("1974-03-31"), as.Date("2016-12-31"), by="1 day"))
day.list$doy=as.numeric(format(day.list$date, "%j"))
day.list$year=as.numeric(format(day.list$date, "%Y"))
daily.kd=data.frame("doy"=0, "year"=0, "date"=as.Date(0), "Secchi"=0)

for(i in 1:nrow(secchi.params))
{
  this.year=subset(day.list, year==secchi.params$year[i])
  open.water=data.frame(doy=seq(100,320), Secchi=customLogitModel(x=seq(100,320), max=secchi.params$max[i], min=secchi.params$min[i]), year=as.numeric(secchi.params$year[i]))
  this.year=merge(this.year, open.water, by=c("doy", "year"), all=T)
  daily.kd=rbind(daily.kd, this.year)
}
daily.kd=daily.kd[2:nrow(daily.kd),]

daily.kd$Secchi=na.approx(daily.kd$Secchi, na.rm=F)
daily.kd$kd=1.7/daily.kd$Secchi
write_tsv(daily.kd, "daily_kd_interpolated.txt")

ggplot()+geom_path(data=daily.kd, aes(date, kd))+theme_classic()
ggsave("kd_time_series_interpolated.png", height=5, width=8, units="in")

ggplot()+geom_path(data=daily.kd, aes(doy, Secchi, group=year, colour=factor(year)))+theme_classic()
ggsave("secchi_seasonal_interpolated.png", height=5, width=8, units="in")
