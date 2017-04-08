library(lubridate)
library(ggplot2)
library(gridExtra)
library(zyp)

secchi = read.table('inst/extdata/ML_secchi.txt', sep='\t', header=TRUE, as.is=TRUE)
secchi$date = mdy(secchi$date)
secchi$yday = yday(secchi$date)
secchi$week = week(secchi$date)
secchi$biweek = 2*floor(week(secchi$date)/2)
secchi$year = year(secchi$date)

boxplot(secchi_m~biweek, secchi)


g1 = ggplot(subset(secchi, year < 1997), aes(as.factor(biweek), secchi_m)) + geom_boxplot() + ylim(0,7) + ggtitle('Pre-1997')
g2 = ggplot(subset(secchi, year >= 1997), aes(as.factor(biweek), secchi_m)) + geom_boxplot() + ylim(0,7) + ggtitle('Post-1997')

gridExtra::grid.arrange(g1, g2)


ggplot(secchi, aes(as.factor(biweek), secchi_m)) + geom_boxplot() + ylim(0,7) + ggtitle('All Years')

trends = plyr::ddply(secchi, 'biweek', function(df){zyp.sen(secchi_m~year, df)$coeff})

plot(year~biweek, trends, ylim=c(-0.05,0.1), ylab='Secchi Trend (m/yr)')

boxplot(secchi_m~year, secchi)
