### Detecting abrupt events in Holocene speleothem records

setwd("C:/Users/sarah/OneDrive/Documents/PhD/abrupt_Holocene/speleothem_8_2_kyr_signals/")

library(dplyr)
library(strucchange)
library(ggplot2)
library(RMySQL)

# connect to SISAL database
mydb <- dbConnect(MySQL(), user = "root", password = "", dbname = "sisalv2", 
                  host = "localhost")

# load Holocene data
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING (site_id) JOIN sample USING (entity_id) JOIN original_chronology USING (sample_id) JOIN d18O USING (sample_id)
                       WHERE (interp_age BETWEEN 0 AND 12000);")
Raw_Data <- Raw_Data %>% filter(entity_status != "superseded")

# Select SISAL data with sufficient resolution for abrupt event detection

min_res <- 30

## load function for calculating mean temporal resolution (excluding hiatuses and gaps) - 'get_ent_sampling'
load("entity_sampling_mean_res.R")

res_out <- data.frame()
highres_dat <- data.frame()
for (i in unique(Raw_Data$entity_id)){ # for every entity
  subdat <- Raw_Data %>% filter(entity_id == i)
  
  subdat2 <- data.frame()
  for (j in seq(300,12000,300)){ # for every 300 year bin
    sub_subdat <- subdat %>% filter(interp_age >= (j-300) & interp_age <= j)
    if (nrow(sub_subdat) <= 1){ next }
    mean_res <- get_ent_sampling(entity_id = i, age_start = (j-500), age_end = j)$sampling_mean
    
    res_out <- rbind(res_out,
                     data.frame(unique(subdat[c("site_id","site_name","entity_id")]),
                                bin_centre = j-150,
                                mean_res = mean_res))
    if (mean_res >= min_res){ next }#
    
    subdat2 <- rbind(subdat2, sub_subdat)
  }
  highres_dat <- rbind(highres_dat, subdat2)
}

#load function for finding optimal breakpoints
opt_bpts <- function(x) {
  #x = bpts_sum$RSS["BIC",]
  n <- length(x)
  lowest <- vector("logical", length = n-1)
  lowest[1] <- FALSE
  for (i in 2:n) {
    lowest[i] <- x[i] < x[i-1] & x[i] < x[i+1]
  }
  out <- as.integer(names(x)[lowest])
  return(out)
}

# breakpoint analysis for every entity and every 1000 year window (with 50% overlap)
n_records <- data.frame()
bp_out <- data.frame()
dtrend_dat <- data.frame()
nentities <- data.frame()
ptm <- proc.time()
for (i in seq(1000,12000,500)){
  
  window_dat <- highres_dat %>% filter(interp_age >= (i-1000) & interp_age <= i)
  
  ## filter to records > 200 years long
  dat_length <- data.frame()
  for (j in unique(window_dat$entity_id)){
    sub_length <- window_dat %>% filter(entity_id == j)
    length <- max(sub_length$interp_age) - min(sub_length$interp_age)
    sub_df <- data.frame(entity_id = j, length = length)
    dat_length <- rbind(dat_length, sub_df)
  }
  dat_length <- dat_length %>% filter(length <= 200)
  window_dat <- window_dat %>% filter(!entity_id %in% dat_length$entity_id)
  
  n_records <- rbind(n_records, 
                     data.frame(win_start = (i-1000), win_end = i, n_entities = length(unique(window_dat$entity_id))))
  
  # breakpoint analysis
  bp_wind <- data.frame()
  dtrend_wind <- data.frame()
  all_entbins <- data.frame()
  for (j in unique(window_dat$entity_id)){ # each entity
    subdat <- window_dat %>% filter(entity_id == j)
    
    ## detrend
    subdat_lm <- lm(d18O_measurement ~ interp_age, data = subdat)
    lm_predicted <- predict(subdat_lm)
    subdat$detrended_d18O <- residuals(subdat_lm)
    
    if (nrow(subdat) <= 13){ next }  
    bp <- breakpoints(subdat$detrended_d18O ~ 1)#subdat$interp_age) #breakpoint analysis
    bpts_sum <- summary(bp) 
    opt_brks <- opt_bpts(bpts_sum$RSS["BIC",]) #optimal no. breakpoints
    
    if (length(opt_brks) > 1){ opt_brks <- as.numeric(names(which.min(bpts_sum$RSS["BIC",]))) }
    
    if (length(opt_brks) == 0){ next }
    if (is.na(opt_brks)) {next}
    if (opt_brks == 0){ next }
    
    
    ci_x <- confint(bp, breaks = opt_brks) #get timings of bp's with conf intervals
    
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

nentities <- unique(nentities[,-7]) # remove overlap counts


#write.csv(bp_out2, "C:/Users/sarah/OneDrive/Documents/PhD/abrupt_Holocene/Hol_bp.csv", row.names = F)
#write.csv(dtrend_dat, "C:/Users/sarah/OneDrive/Documents/PhD/abrupt_Holocene/Hol_dtrend_dat.csv", row.names = F)
#write.csv(nentities, "C:/Users/sarah/OneDrive/Documents/PhD/abrupt_Holocene/Hol_bp_nentities.csv", row.names = F)
bp_out2 <- read.csv("C:/Users/ph805612/OneDrive - University of Reading/Documents/abrupt_Holocene/Hol_bp.csv")
dtrend_dat <- read.csv("C:/Users/ph805612/OneDrive - University of Reading/Documents/abrupt_Holocene/Hol_dtrend_dat.csv")
nentities <- read.csv("C:/Users/ph805612/OneDrive - University of Reading/Documents/abrupt_Holocene/Hol_bp_nentities.csv")

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
#bp_pcent <- bp_pcent %>% group_by(bin, pcent) %>% summarise(bin_age = c(bin-150, bin+150))

##
bp_pcenta <- bp_pcent; bp_pcentb <- bp_pcent
bp_pcenta$bin_age <- bp_pcenta$bin-150
bp_pcentb$bin_age <- bp_pcentb$bin+150
bp_pcent <- rbind(bp_pcenta, bp_pcentb)

bp_pcent <- bp_pcent %>% arrange(bin)

png("C:/Users/ph805612/OneDrive - University of Reading/Documents/abrupt_Holocene/Fig2_Holocene_bp.png", width = 18, height = 12, units = "cm", res = 200)
ggplot(data = bp_pcent, aes(x = bin_age, y = pcent*100)) + 
  geom_line(stat = "identity") + 
  geom_segment(aes(x = 8200, y = 80, xend = 8200, yend = 73), arrow = arrow(length = unit(0.1, "cm")), col = "red") +
  geom_text(aes(x = 8200, y = 82, label = "8.2 ka"), col = "red") +
  scale_x_continuous(breaks = seq(0,12000,1000), expand = c(0.01,0.01)) +
  ylab("% entities") + xlab("Age (years BP)") +
  theme_bw()
dev.off()


