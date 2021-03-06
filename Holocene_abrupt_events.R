### Detecting abrupt events in Holocene speleothem records

setwd(".../speleothem_8_2_kyr_signals/")

library(dplyr)
library(strucchange)
library(ggplot2)
library(RMySQL)

# connect to SISAL database
mydb <- dbConnect(MySQL(), user = "root", password = "", dbname = "sisal_v2", 
                  host = "localhost")

# load Holocene data
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING (site_id) JOIN sample USING (entity_id) JOIN original_chronology USING (sample_id) JOIN d18O USING (sample_id)
                       WHERE (interp_age BETWEEN 0 AND 12000);")
Raw_Data <- Raw_Data %>% filter(entity_status != "superseded")

# Select SISAL data with sufficient resolution for abrupt event detection

min_res <- 30

## load function for calculating mean temporal resolution (excluding hiatuses and gaps) - 'get_ent_sampling'
source("entity_sampling_mean_res.R")

res_out <- data.frame()
highres_dat <- data.frame()
for (i in unique(Raw_Data$entity_id)){ # for every entity
  subdat <- Raw_Data %>% filter(entity_id == i) # filter to that entity
  
  subdat2 <- data.frame()
  for (j in seq(300,12000,300)){ # for every 300 year bin within entity
    sub_subdat <- subdat %>% filter(interp_age >= (j-300) & interp_age <= j) #filter to bin
    if (nrow(sub_subdat) <= 1){ next } # if no data, skip
    mean_res <- get_ent_sampling(entity_id = i, age_start = (j-500), age_end = j)$sampling_mean #calc mean sampling res
    
    res_out <- rbind(res_out,
                     data.frame(unique(subdat[c("site_id","site_name","entity_id")]),
                                bin_centre = j-150,
                                mean_res = mean_res)) # save
    if (mean_res >= min_res){ next } # don't include bins with insufficient sampling res
    
    subdat2 <- rbind(subdat2, sub_subdat) #save
  }
  highres_dat <- rbind(highres_dat, subdat2)
}


# breakpoint analysis for every entity and every 1000 year window (with 50% overlap)
n_records <- data.frame()
bp_out <- data.frame()
dtrend_dat <- data.frame()
nentities <- data.frame()
ptm <- proc.time()
for (i in seq(1000,12000,500)){ # for each window
  
  window_dat <- highres_dat %>% filter(interp_age >= (i-1000) & interp_age <= i) # filter to window
  
  ## filter to records > 200 years long (don't want records that are too short)
  dat_length <- data.frame()
  for (j in unique(window_dat$entity_id)){ #for each entity
    sub_length <- window_dat %>% filter(entity_id == j) #filter dat to entity
    length <- max(sub_length$interp_age) - min(sub_length$interp_age) #calculate length
    sub_df <- data.frame(entity_id = j, length = length) #save
    dat_length <- rbind(dat_length, sub_df)
  }
  dat_length <- dat_length %>% filter(length <= 200) # filter
  window_dat <- window_dat %>% filter(!entity_id %in% dat_length$entity_id) # filter to entities of sufficient length
  
  #save number of records within this window
  n_records <- rbind(n_records, 
                     data.frame(win_start = (i-1000), win_end = i, n_entities = length(unique(window_dat$entity_id))))
  
  # breakpoint analysis
  bp_wind <- data.frame()
  dtrend_wind <- data.frame()
  all_entbins <- data.frame()
  for (j in unique(window_dat$entity_id)){ # for each entity
    subdat <- window_dat %>% filter(entity_id == j) #filter data to that entity
    
    ## detrend using linear regression (remove long term trend)
    subdat_lm <- lm(d18O_measurement ~ interp_age, data = subdat)
    lm_predicted <- predict(subdat_lm)
    subdat$detrended_d18O <- residuals(subdat_lm)
    
    if (nrow(subdat) <= 13){ next }  
    bp <- breakpoints(subdat$detrended_d18O ~ 1)# #breakpoint analysis
    
    if (length(bp$breakpoints) == 0){ next } else if (is.na(length(bp$breakpoints))) {
      next} else if (length(bp$breakpoints == 0)){ 
        next } else {
    
    ci_x <- confint(bp, breaks = length(bp$breakpoints)) #get timings of bp's with conf intervals
    
    if (any(ci_x$confint[,c(1,3)] == 0) | any(ci_x$confint[,c(1,3)] < 0) | any(is.na(ci_x$confint[,c(1,3)]))){ next }
    
    # output breakpoints
    ci_ages <- data.frame(site_id = unique(subdat$site_id), entity_id = unique(subdat$entity_id), entity_name = unique(subdat$entity_name
    ), longitude = unique(subdat$longitude), latitude = unique(subdat$latitude),
    bp = subdat$interp_age[ci_x$confint[,2]],
    CI2_5 = subdat$interp_age[ci_x$confint[,1]],
    CI97_5 = subdat$interp_age[ci_x$confint[,3]])
    
    bp_wind <- rbind(bp_wind, ci_ages)
    
    # output detrended data
    sub_dtrend <- data.frame(site_id = unique(subdat$site_id), entity_id = unique(subdat$entity_id), entity_name = unique(subdat$entity_name), longitude = unique(subdat$longitude), latitude = unique(subdat$latitude),
                             sample_id = subdat$sample_id,
                             interp_age = subdat$interp_age,
                             d18O_detrended = subdat$detrended_d18O)
    
    dtrend_wind <- rbind(dtrend_wind, sub_dtrend)
    
    # output which 200 year bins within the 1000yr window the entity records
    subdat$bin <- cut(subdat$interp_age, breaks = seq(0,12000,300), labels = seq(150,11850,300))
    ent_bins <- subdat %>% group_by(site_id, entity_id, entity_name, latitude, longitude, bin) %>% summarise(n())
    
    all_entbins <- rbind(all_entbins, ent_bins)
  }
  
  bp_wind$win_start <- (i-1000); bp_wind$win_end <- i
  dtrend_wind$win_start <- (i-1000); dtrend_wind$win_end <- i
  
  bp_out <- rbind(bp_out, bp_wind) 
  dtrend_dat <- rbind(dtrend_dat, dtrend_wind)
  nentities <- rbind(nentities, all_entbins)
  
}
proc.time() - ptm 



## Remove replicated bp's (due to overlapping windows)
bp_out2 <- data.frame()
for (i in unique(bp_out$entity_id)){
  subdat <- bp_out %>% filter(entity_id == i)
  if (nrow(subdat) == 1){ subdat <- subdat[,-c(9:10)] 
  } else {
    distmat <- dist(subdat$bp) # identify distances between all bp's
    distmat <- as.matrix(distmat)
    #dimnames(distmat) <- list(1:ncol(distmat),1:ncol(distmat))
    xy <- t(combn(colnames(distmat), 2))
    dist_df <- data.frame(xy, dist = distmat[xy])
    
    #identify pairs where distance <= 20 yrs
    small_dist <- dist_df %>% filter(dist <= 20)
    
    if (nrow(small_dist) == 0){ subdat <- subdat[,-c(9:10)] } else {
      bp_sub <- data.frame()
      for (j in 1:nrow(small_dist)){
        subsubdat <- subdat[as.numeric(small_dist[j,1:2]),]
        bp_sub <- rbind(bp_sub, 
                        data.frame(unique(subdat[,1:5]),
                                   bp = mean(subsubdat$bp),
                                   CI2_5 = mean(subsubdat$CI2_5),
                                   CI97_5 = mean(subsubdat$CI97_5)))
      }
      # remove duplicates from subset
      subdat <- subdat[-as.numeric(as.matrix(small_dist[1:nrow(small_dist),1:2])),]
      #add combined duplicates to subdat
      subdat <- rbind(subdat[,-c(9:10)],bp_sub)
      
    }
  }
  
  bp_out2 <- rbind(bp_out2, subdat)
}
}

nentities <- unique(nentities[,-7]) # remove overlap counts


## avg ar for each bin (for randomly generated data for sig. testing)
get_ar_coeff <- function(x){
  ar_x <- arima(x, order = c(1,0,0))
  
  return(as.numeric(ar_x$coef[1]))
}

ar_out <- xx %>% group_by(entity_id, win_start) %>% summarise(ar = get_ar_coeff(d18O_detrended)) %>% group_by(win_start) %>% summarise(mean(ar))


## bp as %

## calculate n-entities per 100 year bin
nentities$bin <- as.numeric(as.character(nentities$bin))
nentities2 <- nentities %>% group_by(bin) %>% summarise(n())

# add bin and window data to bp dataframe
bp_out2$win_start <- cut(bp_out2$bp, breaks = seq(0,12000,1000), labels = seq(0,11000,1000))
bp_out2$bin <- cut(bp_out2$bp, breaks = seq(0,12000,300), labels = seq(150,11850,300))

bp_out2$win_start <- as.numeric(as.character(bp_out2$win_start))
bp_out2$bin <- as.numeric(as.character(bp_out2$bin))

# calculate number of entities with at least 1 bp in each 300yr bin
bp_out3 <- bp_out2 %>% group_by(entity_id, bin) %>% summarise(n_ent = n()) %>%
  filter((n_ent > 1)) %>% group_by(bin) %>%
  summarise(n_ents_w_bp = n())

#calculate % of entities with at least 1 bp in each bin
bp_pcent <- left_join(bp_out3, nentities2)

bp_pcent$pcent <- bp_pcent$n_ents_w_bp/bp_pcent$`n()`

# plot % entities with bp's through Holocene
bp_pcent$bin_start <- bp_pcent$bin - 150; bp_pcent$bin_end <- bp_pcent$bin + 150

##
bp_pcenta <- bp_pcent; bp_pcentb <- bp_pcent
bp_pcenta$bin_age <- bp_pcenta$bin-150
bp_pcentb$bin_age <- bp_pcentb$bin+150
bp_pcent <- rbind(bp_pcenta, bp_pcentb)

bp_pcent <- bp_pcent %>% arrange(bin)

png("Fig2_Holocene_bp.png", width = 18, height = 12, units = "cm", res = 200)
ggplot(data = bp_pcent, aes(x = bin_age, y = pcent*100)) + 
  geom_line(stat = "identity") + 
  geom_segment(aes(x = 8200, y = 80, xend = 8200, yend = 73), arrow = arrow(length = unit(0.1, "cm")), col = "red") +
  geom_text(aes(x = 8200, y = 82, label = "8.2 ka"), col = "red") +
  scale_x_continuous(breaks = seq(0,12000,1000), expand = c(0.01,0.01)) +
  ylab("% entities") + xlab("Age (years BP)") +
  theme_bw()
dev.off()


