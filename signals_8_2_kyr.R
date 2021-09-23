## Constrain 8.2 kyr signals 

setwd("C:/Users/sarah/OneDrive/Documents/PhD/abrupt_Holocene/speleothem_8_2_kyr_signals/")

library(dplyr)
library(strucchange)
library(ggplot2)
library(RMySQL)

# connect to SISAL database
mydb <- dbConnect(MySQL(), user = "root", password = "", dbname = "sisalv2", 
                  host = "localhost")

# load 7.4 to 9.0 ka  data
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING (site_id) JOIN sample USING (entity_id) JOIN original_chronology USING (sample_id) JOIN d18O USING (sample_id)
                       WHERE (interp_age BETWEEN 7400 AND 9000);")
Raw_Data <- Raw_Data %>% filter(entity_status != "superseded")

## Select entities with a sufficient resolution between 7800 and 8400 years (of 30 years)
min_res <- 30

# load function for calculating mean temporal resolution (excluding hiatuses and gaps) - 'get_ent_sampling'
source("entity_sampling_mean_res.R")

res_out <- data.frame()
for (i in unique(Raw_Data$entity_id)){ # for every entity
  subdat <- Raw_Data %>% filter(entity_id == i & interp_age >= 7800 & interp_age <= 8400)
  if (nrow(subdat) <= 1){ next }
  mean_res <- get_ent_sampling(entity_id = i, age_start = 7800, age_end = 8400)$sampling_mean
    
  res_out <- rbind(res_out,
                     data.frame(unique(subdat[c("site_id","site_name","entity_id")]),
                                mean_res = mean_res))
}

res_out <- res_out %>% filter(mean_res <= 30)

## breakpoint analysis
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

dtrend_dat <- data.frame()
bp_dat <- data.frame()
# breakpoint analysis
for (i in unique(res_out$entity_id)){
  subdat <- Raw_Data %>% filter(entity_id == i & interp_age >= 7800 & interp_age <= 8400)
  
  if (nrow(subdat) <= 15){ next }
  
  ## detrend
  subdat_lm <- lm(d18O_measurement ~ interp_age, data = subdat)
  lm_predicted <- predict(subdat_lm)
  subdat$detrended_d18O <- residuals(subdat_lm)
  
  bp <- breakpoints(subdat$detrended_d18O ~ 1) #breakpoint analysis
  bpts_sum <- summary(bp) 
  opt_brks <- opt_bpts(bpts_sum$RSS["BIC",]) #optimal no. breakpoints
  
  if (length(opt_brks) > 1){ opt_brks <- as.numeric(names(which.min(bpts_sum$RSS["BIC",]))) }
  if (length(opt_brks) == 0){ next }
  if (is.na(opt_brks)) {next}
  if (opt_brks == 0){ next }
  
  ci_x <- confint(bp, breaks = opt_brks) #get timings of bp's with conf intervals
  
  if (any(ci_x$confint[,c(1,3)] == 0) | any(ci_x$confint[,c(1,3)] < 0) | any(is.na(ci_x$confint[,c(1,3)]))){ print(i) } else {
    # output breakpoints
    ci_ages <- data.frame(site_id = unique(subdat$site_id), site_name = unique(subdat$site_name), entity_id = unique(subdat$entity_id), entity_name = unique(subdat$entity_name),
                          longitude = unique(subdat$longitude), latitude = unique(subdat$latitude),
                          bp = subdat$interp_age[ci_x$confint[,2]],
                          CI2_5 = subdat$interp_age[ci_x$confint[,1]],
                          CI97_5 = subdat$interp_age[ci_x$confint[,3]])
    
    bp_dat <- rbind(bp_dat, ci_ages)
    
    # output detrended data
    sub_dtrend <- data.frame(site_id = unique(subdat$site_id), site_name = unique(subdat$site_name), entity_id = unique(subdat$entity_id), entity_name = unique(subdat$entity_name), longitude = unique(subdat$longitude), latitude = unique(subdat$latitude),
                             sample_id = subdat$sample_id,
                             interp_age = subdat$interp_age,
                             d18O_detrended = subdat$detrended_d18O)
    
    dtrend_dat <- rbind(dtrend_dat, sub_dtrend)
    
  }
}

## summarise results
x <- bp_dat %>% group_by(entity_id) %>% summarise(n_bp = n())

## grps
grps <- read.csv("82_signal_grps.csv")
colnames(grps)[1] <- "entity_id"

## entities with 2 breakpoints
x2 <- x %>% filter(n_bp == 2) # 26

for (i in seq(9,28,9)){
  entities <- x2$entity_id[(i-8):i]
  subdat <- dtrend_dat %>% filter(entity_id %in% entities)
  sub_bp <- bp_dat %>% filter(entity_id %in% entities)
  
  file_no <- ifelse(i == 9, 1, ifelse(i == 18, 2, 3))
  filename <- paste("signals_82_", file_no, ".pdf", sep = "")
  
  p <- ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
    geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
    facet_wrap(.~ entity_id)
  pdf(filename, width = 20/2.54, height = 18/2.54)
  print(p)
  dev.off()
}

# exclude entities: 53, 54, 295, 374, 388, 395, 540, 608, 690
x2 <- x2 %>% filter(!entity_id %in% c(53, 54, 295, 374, 388, 395, 540, 608, 690)) # 20 entities

anom_2bp <- data.frame()
for (i in unique(x2$entity_id)){
  sub_bp <- bp_dat %>% filter(entity_id == i)
  subdat <- dtrend_dat %>% filter(entity_id == i)
  
  d18O_event <- subdat %>% filter(interp_age >= min(sub_bp$bp) & interp_age <= max(sub_bp$bp)) %>% summarise(mean(d18O_detrended))
  event_CI25 <- subdat %>% filter(interp_age >= min(sub_bp$CI2_5) & interp_age <= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
  d18O_base <- subdat %>% filter(interp_age <= min(sub_bp$bp) | interp_age >= max(sub_bp$bp)) %>% summarise(mean(d18O_detrended))
  
  anom <- d18O_event - d18O_base

  sub_df <- data.frame(unique(sub_bp[,c(1:6)]), 
                       min_bp = as.numeric(min(sub_bp$bp)), max_bp = as.numeric(max(sub_bp$bp)),
                       sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-as.numeric(min(sub_bp$bp)))),"sample_id"],
                       sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-as.numeric(max(sub_bp$bp)))),"sample_id"],
                       anom = as.numeric(anom),
                       d18O_base = as.numeric(d18O_base),
                       d18O_event = as.numeric(d18O_event))
  
  anom_2bp <- rbind(anom_2bp, sub_df)
}

ggplot(data = anom_2bp, aes(x = longitude, y = latitude, fill = anom)) +
  geom_point(shape = 21) +
  borders("world") +
  scale_fill_gradient2(high = "red", mid = "white", low = "blue")


## entities with 3 breakpoints
x3 <- x %>% filter(n_bp == 3) # 23

# visualise entities
for (i in seq(9,28,9)){
  entities <- x3$entity_id[(i-8):i]
  subdat <- dtrend_dat %>% filter(entity_id %in% entities)
  sub_bp <- bp_dat %>% filter(entity_id %in% entities)
  
  file_no <- ifelse(i == 9, 1, ifelse(i == 18, 2, 3))
  filename <- paste("signals_82_3bp_", file_no, ".pdf", sep = "")
  
  p <- ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
    geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
    facet_wrap(.~ entity_id)
  pdf(filename, width = 20/2.54, height = 18/2.54)
  print(p)
  dev.off()
}

# exclude entities: 63, 142, 436, 279, 546, 613
x3 <- x3 %>% filter(!entity_id %in% c(63, 142, 279, 436, 546, 613)) # 14 entities



## entities with 4 breakpoints
x4 <- x %>% filter(n_bp == 4) # 9

# visualise
subdat <- dtrend_dat %>% filter(entity_id %in% x4$entity_id)
sub_bp <- bp_dat %>% filter(entity_id %in% x4$entity_id)
  
filename <- "signals_82_4bp.pdf"
  
p <- ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
    geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
    facet_wrap(.~ entity_id)
pdf(filename, width = 20/2.54, height = 18/2.54)
print(p)
dev.off()


# exclude entities: 220, 351, 415 = hiatus
x4 <- x4 %>% filter(!entity_id %in% c(220,351,415))



## 5 bp's
x5 <- x %>% filter(n_bp == 5) # 2
x5 <- x5 %>% filter(entity_id == 591)

sub_bp <- bp_dat %>% filter(entity_id %in% x5$entity_id)
sub_bp <- sub_bp[order(sub_bp$bp),]
subdat <- dtrend_dat %>% filter(entity_id %in% x5$entity_id)

ggplot() + geom_line(data =subdat, aes(x = interp_age, y = d18O_detrended)) + 
  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
  facet_wrap(.~ entity_id)


# calculate anomalies
anom_bp <- data.frame()
x_all <- rbind(x3,x4)
for (i in unique(x_all$entity_id)){
  sub_bp <- bp_dat %>% filter(entity_id == i) %>% arrange(bp)
  sub_bp <- sub_bp[order(sub_bp$bp),]
  subdat <- dtrend_dat %>% filter(interp_age >= 7400 & interp_age <= 9000 & 
                                    entity_id == i)
  
  if (nrow(sub_bp) == 3){ lab_vec = 1:4 } else { lab_vec = 1:5 }
  subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp,max(interp_age)), labels = lab_vec))
  
  event_grp <- as.numeric(grps[which(grps$entity_id == i),2:4])
  event_grp <- na.omit(event_grp)
  
  event_d18O <- subdat %>% filter(grp %in% event_grp) %>% summarise(mean_d18O = mean(d18O_detrended))
  base_d18O <- subdat %>% filter(grp %in% c(1,4)) %>% summarise(mean_d18O = mean(d18O_detrended))
  #base_d18O <- subdat %>% filter(interp_age >= 7400 & interp_age <= 7900 | interp_age >= 8500 & interp_age <= 9000) %>% summarise(mean_d18O = mean(d18O_detrended))
  x_sub <- data.frame(unique(subdat[which(subdat$entity_id == i),1:6]),
                      min_bp = sub_bp$bp[min(event_grp,na.rm = T)-1], max_bp = sub_bp$bp[max(event_grp, na.rm = T)], 
                      sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[min(event_grp,na.rm = T)-1])),"sample_id"],
                      sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[max(event_grp, na.rm = T)])),"sample_id"],
                      anom = as.numeric(event_d18O) - as.numeric(base_d18O),
                      d18O_base = as.numeric(base_d18O), d18O_event = as.numeric(event_d18O))
  anom_bp <- rbind(anom_bp, x_sub)
}


all_dat <- rbind(anom_2bp, anom_bp)
all_dat$duration <- all_dat$max_bp - all_dat$min_bp
all_dat$event_centre <- all_dat$max_bp - (all_dat$duration/2)


## add non-SISAL records
nonSISAL <- read.csv("C:/Users/sarah/OneDrive/Documents/PhD/abrupt_Holocene/nonSISAL_82_signals.csv")

all_dat <- rbind(all_dat, nonSISAL)


## entities with no bp
nosignal <- Raw_Data %>% filter(entity_id %in% c(53, 54, 295, 374, 395, 540, 690, 63, 142, 279, 436, 546, 613, 220,351, 305))
no_signal <- unique(nosignal[c("entity_id","site_id","site_name","latitude","longitude")])

ggplot() + 
  geom_point(data = all_dat, aes(x = longitude, y = latitude, fill = anom), shape = 21, size = 3) + borders("world") +
  geom_point(data = no_signal, aes(x = longitude, y = latitude)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  coord_fixed(ylim = c(-50,60), xlim = c(-120,150))

write.csv(nosignal, "abrupt_Holocene/spel_nosignal.csv", row.names = F)

y <- dtrend_dat %>% filter(entity_id %in% c(295,385,387))
y_bp <- bp_dat %>% filter(entity_id %in% c(295,385,387))
ggplot() +
  geom_vline(data = y_bp, aes(xintercept = bp, col = as.factor(entity_id))) +
  geom_line(data = y, aes(x = interp_age, y = d18O_detrended, col = as.factor(entity_id)))
