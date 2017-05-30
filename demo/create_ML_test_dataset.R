#create ML datset for modeling

# model 


kd.daily <- readr::read_csv("inst/extdata/ml_model_kd.csv") %>% rename(DateTime = time)
nldas.drivers <- readr::read_csv("inst/extdata/ml_nldas_drivers.csv") %>% rename(DateTime = time)
obs.mod <- readr::read_csv('out/obsmod.csv') %>% 
  left_join(kd.daily) %>% 
  left_join(nldas.drivers)
library(dplyr)

depths <- group_by(obs.mod, Depth) %>% tally %>% arrange(desc(n)) %>% .$Depth

z.names <- c("emotional", "dome", "discovery", "debate", "fragment", "orangutang")
out.df <- obs.mod
out.df$Depth <- NA
for (i in 1:length(z.names)){
  out.df$Depth[obs.mod$Depth == depths[i]] <- z.names[i]
}

z_score <- function(data){
  return((data - mean(data))/sd(data))
}

out.df <- out.df %>% filter(!is.na(Depth), DateTime > as.Date("1990-06-14")) %>% rename(dim_2 = Depth, cactus = Modeled_temp, coffee = Kd) %>% mutate(error = cactus - Observed_temp) %>% 
  mutate(villan = lubridate::yday(DateTime), taco = lubridate::year(DateTime) - min(lubridate::year(DateTime)) + 1, dim_1 = as.numeric(DateTime) - min(as.numeric(DateTime))) %>% 
  mutate(groceries = ShortWave, rabbit = ifelse(AirTemp < 0, 1, 0), donkey = AirTemp, trailmix = LongWave, saddle = RelHum, whisper = WindSpeed) %>% 
  select(dim_1, dim_2, error, everything(), -DateTime, -Observed_temp, -ShortWave, -LongWave, -AirTemp, -RelHum, -WindSpeed, -Rain, -Snow)


# write key!!

saveRDS(out.df, 'inst/extdata/full_training_ML.rds')

challenge.df <- out.df %>% filter(dim_1 < 5600 | dim_1 > 8200)
write.table(x = challenge.df, file = 'out/ML_challenge_single_system.tsv', sep = '\t', row.names = FALSE, quote = TRUE)

TrainData <- challenge.df[, -2]
TrainClasses = as.factor(challenge.df$dim_2)
