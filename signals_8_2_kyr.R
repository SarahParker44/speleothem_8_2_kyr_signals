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
    ci_ages <- data.frame(site_id = unique(subdat$site_id), entity_id = unique(subdat$entity_id), entity_name = unique(subdat$entity_name),
                          longitude = unique(subdat$longitude), latitude = unique(subdat$latitude),
                          bp = subdat$interp_age[ci_x$confint[,2]],
                          CI2_5 = subdat$interp_age[ci_x$confint[,1]],
                          CI97_5 = subdat$interp_age[ci_x$confint[,3]])
    
    bp_dat <- rbind(bp_dat, ci_ages)
    
    # output detrended data
    sub_dtrend <- data.frame(site_id = unique(subdat$site_id), entity_id = unique(subdat$entity_id), entity_name = unique(subdat$entity_name), longitude = unique(subdat$longitude), latitude = unique(subdat$latitude),
                             sample_id = subdat$sample_id,
                             interp_age = subdat$interp_age,
                             d18O_detrended = subdat$detrended_d18O)
    
    dtrend_dat <- rbind(dtrend_dat, sub_dtrend)
    
  }
}

## summarise results
x <- bp_dat %>% group_by(entity_id) %>% summarise(n_bp = n())

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

# exclude entities: 54, 295, 374, 395, 540, 690
x2 <- x2 %>% filter(!entity_id %in% c(54, 295, 374, 395, 540, 690)) # 20 entities

anom_2bp <- data.frame()
for (i in unique(x2$entity_id)){
  sub_bp <- bp_dat %>% filter(entity_id == i)
  subdat <- dtrend_dat %>% filter(entity_id == i)
  
  d18O_event <- subdat %>% filter(interp_age >= min(sub_bp$bp) & interp_age <= max(sub_bp$bp)) %>% summarise(mean(d18O_detrended))
  event_CI25 <- subdat %>% filter(interp_age >= min(sub_bp$CI2_5) & interp_age <= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
  d18O_base <- subdat %>% filter(interp_age <= min(sub_bp$bp) | interp_age >= max(sub_bp$bp)) %>% summarise(mean(d18O_detrended))
  
  anom <- d18O_event - d18O_base

  sub_df <- data.frame(unique(sub_bp[,c(1:5)]), 
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

# exclude entities: 63, 142, 241, 433, 434, 436, 279, 546, 


## entities with 4 breakpoints
x4 <- x %>% filter(n_bp == 4) # 9

subdat <- dtrend_dat %>% filter(entity_id %in% x4$entity_id)
sub_bp <- bp_dat %>% filter(entity_id %in% x4$entity_id)
  
filename <- "signals_82_4bp.pdf"
  
p <- ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
    geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
    facet_wrap(.~ entity_id)
pdf(filename, width = 20/2.54, height = 18/2.54)
print(p)
dev.off()
}

# exclude entities:


## entities with 4 breakpoints
x4 <- x %>% filter(n_bp == 4) # 9

i <- x5$entity_id[2]
sub_bp <- bp_dat %>% filter(entity_id == i)
sub_bp <- sub_bp[order(sub_bp$bp),]
subdat <- dtrend_dat %>% filter(entity_id == i)

ggplot() + geom_line(data =subdat, aes(x = interp_age, y = d18O_detrended)) + 
  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
  ggtitle(unique(filter(Raw_Data, entity_id == i)$site_name))

## 5 bp's
x5 <- x %>% filter(n_bp == 5) # 2




## entities with no bp
nosignal <- left_join(nentities,bp_out2)
nosignal <- nosignal[is.na(nosignal$bp),1:6]
write.csv(nosignal, "abrupt_Holocene/spel_nosignal.csv", row.names = F)

# Focusing on 8.1 to 8.3 bp's
x8_2 <- bp_out2 %>% filter(bin %in% c(7950,8250))
x <- data.frame()
for (i in unique(x8_2$entity_id)){
  sub_bp <- x8_2 %>% filter(entity_id == i)
  subdat <- dtrend_dat %>% filter(interp_age >= 7800 & interp_age <= 8400 & entity_id == i)
  min_age <- min(subdat$interp_age)
  max_age <- max(subdat$interp_age)
  sub_bp <- sub_bp %>% filter(!bp %in% c(min_age, max_age))
  
  x <- rbind(x, sub_bp) 
}

## remove lower res/non-composite entities 
mult <- x %>% group_by(site_id, entity_id) %>% 
  summarise(n()) %>% group_by(site_id) %>% 
  mutate(n()) %>% filter(`n()` > 1)

select_ent <- c()
for (i in unique(mult$site_id)){
  ent <- mult %>% filter(site_id == i)
  subdat <- Raw_Data %>% filter(entity_id %in% ent$entity_id & interp_age >= 7800 & interp_age <= 8400) 
  
  # if site has a composite entity, select this:
  if (any(subdat$speleothem_type == "composite")){
    sub_subdat <- subdat %>% filter(speleothem_type == "composite")
    select_ent <- rbind(select_ent, sub_subdat$entity_id)
  } else { # if not, select higher res entity
    sub_subdat <- subdat %>% group_by(entity_id) %>% summarise(n()) %>% filter(`n()` == max(`n()`))
    select_ent <- rbind(select_ent, sub_subdat$entity_id)
  }
}
remove_ent <- mult %>% filter(!entity_id %in% select_ent)

x <- x %>% filter(!entity_id %in% remove_ent$entity_id)

# no bp
nosignal8_2 <- nosignal %>% filter(bin >= 8100 & bin <= 8300)
nosignal8_2 <- nosignal8_2 %>% group_by(entity_id) %>% filter(n() > 1)
nosignal8_2 <- nosignal8_2 %>% filter(!site_id %in% x$site_id)
nosignal8_2 <- unique(nosignal8_2[,-5]) # 0 entities


# 1 bp
x1 <- x %>% group_by(entity_id) %>% mutate(n_bp = n()) %>% filter(n_bp == 1) # 4 entities
nosignal8_2 <- rbind(nosignal8_2, x1[,1:6]) # 7 entities / 10


# 2 bp's
x2 <- x %>% group_by(entity_id) %>% mutate(n_bp = n()) %>% filter(n_bp == 2) # 20 entities
#sub_bp <- x2 %>% filter(entity_id == 351)
#sub_bp <- sub_bp[order(sub_bp$bp),]
#subdat <- dtrend_dat %>% filter(interp_age >= 7800 & interp_age <= 8400 & entity_id == 351)

#ggplot() + geom_line(data =subdat, aes(x = interp_age, y = d18O_detrended)) + geom_vline(data = sub_bp, aes(xintercept = bp), col = "red")

anom_2bp <- data.frame()
for (i in unique(x2$entity_id)){
  sub_bp <- x2 %>% filter(entity_id == i)
  subdat <- dtrend_dat %>% filter(interp_age >= 7800 & interp_age <= 8400 & entity_id == i)
  
  ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) + 
    geom_vline(data = sub_bp, aes(xintercept = bp), col = "red", lty = 2) + 
    geom_errorbarh(data = sub_bp, aes(xmin = CI2_5, xmax = CI97_5, y = 0.32), col = "red", height = 0.02) +
    ggtitle(paste("entity_id ", i, sep = "= "))
  
  d18O_event <- subdat %>% filter(interp_age >= min(sub_bp$bp) & interp_age <= max(sub_bp$bp)) %>% summarise(mean(d18O_detrended))
  event_CI25 <- subdat %>% filter(interp_age >= min(sub_bp$CI2_5) & interp_age <= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
  d18O_base <- subdat %>% filter(interp_age <= min(sub_bp$bp) | interp_age >= max(sub_bp$bp)) %>% summarise(mean(d18O_detrended))
  
  anom <- d18O_event - d18O_base
  
  
  
  sub_df <- data.frame(unique(sub_bp[,c(1:5)]), 
                       min_bp = as.numeric(min(sub_bp$bp)), max_bp = as.numeric(max(sub_bp$bp)),
                       sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-as.numeric(min(sub_bp$bp)))),"sample_id"],
                       sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-as.numeric(max(sub_bp$bp)))),"sample_id"],
                       anom = as.numeric(anom),
                       d18O_base = as.numeric(d18O_base),
                       d18O_event = as.numeric(d18O_event))
  
  anom_2bp <- rbind(anom_2bp, sub_df)
}

anom_2bp$duration <- anom_2bp$max_bp - anom_2bp$min_bp
anom_2bp$event_centre <- anom_2bp$max_bp - (anom_2bp$duration/2)


# 3 bp's
x3 <- x %>% group_by(entity_id) %>% mutate(n_bp = n()) %>% filter(n_bp == 3) # 8 entities / 17

i = 385
sub_bp <- x3 %>% filter(entity_id == i)
sub_bp <- sub_bp[order(sub_bp$bp),]
subdat <- dtrend_dat %>% filter(interp_age >= 7800 & interp_age <= 8400 & entity_id == i)

ggplot() + geom_line(data =subdat, aes(x = interp_age, y = d18O_detrended)) + geom_vline(data = sub_bp, aes(xintercept = bp), col = "red")

subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp,max(interp_age)), labels = 1:4))
d18O.lm <- lm(d18O_detrended ~ grp, data = subdat)
d18O.av <- aov(d18O.lm)

tukey.test <- TukeyHSD(d18O.av)
tukey_grps <- tukey.test$grp

#anom_3bp <- data.frame()
#for (i in unique(x3$entity_id)){
#  sub_bp <- x3 %>% filter(entity_id == i) %>% arrange(bp)
#  subdat <- dtrend_dat %>% filter(interp_age >= 7800 & interp_age <= 8400 & entity_id == i)

#  ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) + 
#    geom_vline(data = sub_bp, aes(xintercept = bp), col = "red", lty = 2) + 
#    geom_errorbarh(data = sub_bp, aes(xmin = CI2_5, xmax = CI97_5, y = 0.32), col = "red", height = 0.02) +
#    ggtitle(paste("entity_id ", i, sep = "= "))

#  if (max(subdat$interp_age) == max(sub_bp$bp) | min(subdat$interp_age) == min(sub_bp$bp)){
#    print(i)
#  } else {
#    subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp,max(interp_age)), labels = 1:4))
#    

#    d18O.lm <- lm(d18O_detrended ~ grp, data = subdat)
#   d18O.av <- aov(d18O.lm)

#    tukey.test <- TukeyHSD(d18O.av)
#    tukey_grps <- tukey.test$grp

#    if (tukey_grps["3-2","p adj"] <= 0.01){
#      if (tukey_grps["2-1","p adj"] <= 0.01 &
#          tukey_grps["4-2","p adj"] <= 0.01 &
#          tukey_grps["3-1","p adj"] <= 0.01 &
#          tukey_grps["4-3","p adj"] <= 0.01){ # Nuanhe (466), Dongge (446)

#        event_1 <- subdat[which(subdat$interp_age >= sub_bp$bp[1] & subdat$interp_age <= sub_bp$bp[2]),]
#        event_2 <- subdat[which(subdat$interp_age >= sub_bp$bp[2] & subdat$interp_age <= sub_bp$bp[3]),]

#        d18O_base = mean(subdat[which(subdat$interp_age < sub_bp$bp[1] | subdat$interp_age > sub_bp$bp[3]), "d18O_detrended"])
#        d18O_base_25 = mean(subdat[which(subdat$interp_age < sub_bp$CI2_5[1] | subdat$interp_age > sub_bp$CI2_5[3]), "d18O_detrended"])
#        d18O_base_975 = mean(subdat[which(subdat$interp_age < sub_bp$CI97_5[1] | subdat$interp_age > sub_bp$CI97_5[3]), "d18O_detrended"])
#        d18O_base_25_975 = mean(subdat[which(subdat$interp_age < sub_bp$CI2_5[1] | subdat$interp_age > sub_bp$CI97_5[3]), "d18O_detrended"])
#        d18O_base_975_25 = mean(subdat[which(subdat$interp_age < sub_bp$CI97_5[1] | subdat$interp_age > sub_bp$CI2_5[3]), "d18O_detrended"])

#        d18O_longer = mean(subdat[which(subdat$interp_age >= sub_bp$bp[1] & subdat$interp_age <= sub_bp$bp[3]), "d18O_detrended"])
#        d18O_longer_25 = mean(subdat[which(subdat$interp_age >= sub_bp$CI2_5[1] & subdat$interp_age <= sub_bp$CI2_5[3]), "d18O_detrended"])
#        d18O_longer_975 = mean(subdat[which(subdat$interp_age >= sub_bp$CI97_5[1] & subdat$interp_age <= sub_bp$CI97_5[3]), "d18O_detrended"])
#        d18O_longer_25_975 = mean(subdat[which(subdat$interp_age >= sub_bp$CI2_5[1] & subdat$interp_age <= sub_bp$CI97_5[3]), "d18O_detrended"])
#        d18O_longer_975_25 = mean(subdat[which(subdat$interp_age >= sub_bp$CI97_5[1] & subdat$interp_age <= sub_bp$CI2_5[3]), "d18O_detrended"])

#        anom  = d18O_longer - d18O_base
#        anom_25 = d18O_longer_25 - d18O_base_25
#        anom_975 = d18O_longer_975 - d18O_base_975
#        anom_25_975 = d18O_longer_25_975 - d18O_base_25_975
#        anom_975_25 = d18O_longer_975_25 - d18O_base_975_25

#        longer_event <- data.frame(unique(subdat[,1:5]), min_bp = sub_bp$bp[1], max_bp = sub_bp$bp[3], sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[1])),"sample_id"],
#                                   sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[3])),"sample_id"], anom = d18O_longer - d18O_base,  d18O_base = d18O_base, d18O_event = d18O_longer, event ="longer")

#        if (max(abs(event_1$d18O_detrended)) > max(abs(event_2$d18O_detrended))){ # bp1-2 more xtreme values than bp2-3
#          d18O_shorter = mean(subdat[which(subdat$interp_age >= sub_bp$bp[1] & subdat$interp_age <= sub_bp$bp[2]), "d18O_detrended"])
#          d18O_shorter_25 = mean(subdat[which(subdat$interp_age >= sub_bp$CI2_5[1] & subdat$interp_age <= sub_bp$CI2_5[2]), "d18O_detrended"])
#          d18O_shorter_975 = mean(subdat[which(subdat$interp_age >= sub_bp$CI97_5[1] & subdat$interp_age <= sub_bp$CI97_5[2]), "d18O_detrended"])
#          d18O_shorter_25_975 = mean(subdat[which(subdat$interp_age >= sub_bp$CI2_5[1] & subdat$interp_age <= sub_bp$CI97_5[2]), "d18O_detrended"])
#          d18O_shorter_975_25 = mean(subdat[which(subdat$interp_age >= sub_bp$CI97_5[1] & subdat$interp_age <= sub_bp$CI2_5[2]), "d18O_detrended"])

#          anom  = d18O_shorter - d18O_base
#          anom_25 = d18O_shorter_25 - d18O_base_25
#          anom_975 = d18O_shorter_975 - d18O_base_975
#          anom_25_975 = d18O_shorter_25_975 - d18O_base_25_975
#          anom_975_25 = d18O_shorter_975_25 - d18O_base_975_25

#          shorter_event <- data.frame(unique(subdat[,1:5]), min_bp = sub_bp$bp[2], max_bp = sub_bp$bp[3], sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[2])),"sample_id"],
#                                      sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[3])),"sample_id"], anom = d18O_shorter - d18O_base,  d18O_base = d18O_base, d18O_event = d18O_shorter, event = "shorter")
#        } else {
#          d18O_shorter = mean(subdat[which(subdat$interp_age >= sub_bp$bp[2] & subdat$interp_age <= sub_bp$bp[3]), "d18O_detrended"])
#          d18O_shorter_25 = mean(subdat[which(subdat$interp_age >= sub_bp$CI2_5[1] & subdat$interp_age <= sub_bp$CI2_5[2]), "d18O_detrended"])
#          d18O_shorter_975 = mean(subdat[which(subdat$interp_age >= sub_bp$CI97_5[1] & subdat$interp_age <= sub_bp$CI97_5[2]), "d18O_detrended"])
#          d18O_shorter_25_975 = mean(subdat[which(subdat$interp_age >= sub_bp$CI2_5[1] & subdat$interp_age <= sub_bp$CI97_5[2]), "d18O_detrended"])
#          d18O_shorter_975_25 = mean(subdat[which(subdat$interp_age >= sub_bp$CI97_5[1] & subdat$interp_age <= sub_bp$CI2_5[2]), "d18O_detrended"])
#          
#          shorter_event <- data.frame(unique(subdat[,1:5]), 
#                                      min_bp = sub_bp$bp[2], 
#                                      max_bp = sub_bp$bp[3], 
#                                      sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[2])),"sample_id"],
#                                      sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[3])),"sample_id"], 
#                                      anom = d18O_shorter - d18O_base,  
#                                      d18O_base = d18O_base, 
#                                      d18O_event = d18O_shorter, 
#                                      event = "shorter")
#        }
#       out_df <- rbind(longer_event, shorter_event)
#        
#        
#      } else { # group 1
#        d18O_base = mean(subdat[which(subdat$interp_age < sub_bp$bp[1] | subdat$interp_age > sub_bp$bp[3]), "d18O_detrended"])
#        # identify the sig event
#        if (tukey_grps["2-1","p adj"] <0.01 & tukey_grps["4-2","p adj"] <0.1){
#         #bp1-2
#          d18O_shorter <- mean(subdat[which(subdat$interp_age >= sub_bp$bp[1] & subdat$interp_age <= sub_bp$bp[2]), "d18O_detrended"])
#        }  else {
#          #bp2-3
#          d18O_shorter <- mean(subdat[which(subdat$interp_age >= sub_bp$bp[2] & subdat$interp_age <= sub_bp$bp[3]), "d18O_detrended"])
#       }
#        out_df <- data.frame(unique(subdat[,1:5]), min_bp = sub_bp$bp[2], max_bp = sub_bp$bp[3], sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[2])),"sample_id"],
#                            sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[3])),"sample_id"], anom = d18O_shorter - d18O_base, d18O_base = d18O_base, d18O_event = d18O_shorter, event = "shorter")
#     }
#    } else { # group 3
#      event_1 <- subdat[which(subdat$interp_age >= sub_bp$bp[1] & subdat$interp_age <= sub_bp$bp[2]),]
#      event_2 <- subdat[which(subdat$interp_age >= sub_bp$bp[2] & subdat$interp_age <= sub_bp$bp[3]),]
#     
#      d18O_base = mean(subdat[which(subdat$interp_age < sub_bp$bp[1] | subdat$interp_age > sub_bp$bp[3]), "d18O_detrended"])
#      d18O_longer = mean(subdat[which(subdat$interp_age >= sub_bp$bp[1] & subdat$interp_age <= sub_bp$bp[3]), "d18O_detrended"])
#      
#      longer_event <- data.frame(unique(subdat[,1:5]), min_bp = sub_bp$bp[1], max_bp = sub_bp$bp[3],  sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[1])),"sample_id"],
#                                 sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[3])),"sample_id"], anom = d18O_longer - d18O_base, d18O_base = d18O_base, d18O_event = d18O_longer, event ="longer")
#      
#      if (max(abs(event_1$d18O_detrended)) > max(abs(event_2$d18O_detrended))){ # bp1-2 more xtreme values than bp2-3
#        d18O_shorter = mean(subdat[which(subdat$interp_age >= sub_bp$bp[1] & subdat$interp_age <= sub_bp$bp[2]), "d18O_detrended"])
#        shorter_event <- data.frame(unique(subdat[,1:5]), min_bp = sub_bp$bp[2], max_bp = sub_bp$bp[3], anom = d18O_shorter - d18O_base,  d18O_base = d18O_base, d18O_event = d18O_shorter, event = "shorter")
#      } else {
#        d18O_shorter = mean(subdat[which(subdat$interp_age >= sub_bp$bp[2] & subdat$interp_age <= sub_bp$bp[3]), "d18O_detrended"])
#        shorter_event <- data.frame(unique(subdat[,1:5]), min_bp = sub_bp$bp[2], max_bp = sub_bp$bp[3],  sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[2])),"sample_id"],
#                                    sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[3])),"sample_id"], anom = d18O_shorter - d18O_base,  d18O_base = d18O_base, d18O_event = d18O_shorter, event = "shorter")
#      }
#     out_df <- rbind(longer_event, shorter_event)
#      
#    }
#    anom_3bp <- rbind(anom_3bp, out_df)
#    
#  }
#}

anom_3bp$duration <- anom_3bp$max_bp - anom_3bp$min_bp
anom_3bp$event_centre <- anom_3bp$max_bp - (anom_3bp$duration/2)
anom_3bp <- anom_3bp[anom_3bp$event == "shorter",-13];


## 4 bp's
x4 <- x %>% group_by(entity_id) %>% mutate(n_bp = n()) %>% filter(n_bp == 4) # 10 entities / 12 entities

i = 475
sub_bp <- x4 %>% filter(entity_id == i)
sub_bp <- sub_bp[order(sub_bp$bp),]
subdat <- dtrend_dat %>% filter(interp_age >= 7800 & interp_age <= 8400 & entity_id == i)

ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
  geom_vline(data = sub_bp, aes(xintercept = bp), colour = "red") 

subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp,max(interp_age)), labels = 1:5))
sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
sub.av <- aov(sub.lm)
summary(sub.av)

tukey_sub <- as.data.frame(TukeyHSD(sub.av)$grp)
tukey_sub[order(abs(tukey_sub$diff)),] # grp 4 most sig diff

grps <- data.frame(entity_id = c(305,#496,
                                 292,591,61,89,200,#640,
                                 254,#466,529),
                                 590),
                   grp1 = c(2,#2,
                            3,4,3,3,3,#3,
                            2,#2,2),
                            3),
                   grp2 = c(3,#3,
                            NA,NA,NA,4,NA,#4,
                            3,#3,3),
                            NA),
                   grp3 = c(NA,#NA,
                            NA,NA,NA,NA,NA,#NA,
                            4,#4,NA))
                            NA))
xx4 <- data.frame()
for (i in unique(x4$entity_id)){
  sub_bp <- x4 %>% filter(entity_id == i) %>% arrange(bp)
  sub_bp <- sub_bp[order(sub_bp$bp),]
  subdat <- dtrend_dat %>% filter(interp_age >= 7400 & interp_age <= 9000 & 
                                    entity_id == i)
  
  ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
    geom_vline(data = sub_bp, aes(xintercept = bp), col = "red", lty = 2) +
    geom_errorbarh(data = sub_bp, aes(xmin = CI2_5, xmax = CI97_5, y = 1.5), col = "red", height = 0.1, position = position_dodge())
  
  if (i %in% c(61,529)){
    subdat$grp <- cut(subdat$interp_age, breaks = c(sub_bp$bp,max(subdat$interp_age)), labels = 1:4)
  } else {
    subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp,max(interp_age)), labels = 1:5))
    #subdat$grp_25 <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$CI2_5,max(interp_age)), labels = 1:5))
    #subdat$grp_975 <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$CI97_5,max(interp_age)), labels = 1:5))
    
  }
  event_grp <- grps[which(grps$entity_id == i),2:4]
  
  event_d18O <- subdat %>% filter(grp %in% event_grp) %>% summarise(mean_d18O = mean(d18O_detrended))
  base_d18O <- subdat %>% filter(!grp %in% event_grp) %>% summarise(mean_d18O = mean(d18O_detrended))
  #base_d18O <- subdat %>% filter(interp_age >= 7400 & interp_age <= 7900 | interp_age >= 8500 & interp_age <= 9000) %>% summarise(mean_d18O = mean(d18O_detrended))
  x_sub <- data.frame(unique(subdat[which(subdat$entity_id == i),1:5]),
                      min_bp = sub_bp$bp[min(event_grp,na.rm = T)-1], max_bp = sub_bp$bp[max(event_grp, na.rm = T)], 
                      sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[min(event_grp,na.rm = T)-1])),"sample_id"],
                      sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[max(event_grp, na.rm = T)])),"sample_id"],
                      anom = as.numeric(event_d18O) - as.numeric(base_d18O),
                      d18O_base = as.numeric(base_d18O), d18O_event = as.numeric(event_d18O))
  xx4 <- rbind(xx4, x_sub)
}
xx4$duration <- xx4$max_bp - xx4$min_bp
xx4$event_centre <- xx4$max_bp - (xx4$duration/2)


## 5bp
x5 <- x %>% group_by(entity_id) %>% mutate(n_bp = n()) %>% filter(n_bp == 5) # 9 entities / 7
# grps 3,4 for ent_id = 379, grp 4 for ent_id = 190
#grp4 for 129
# grp4 for 150
#ignore 591

grps <- data.frame(entity_id = c(379,190,129,150,#),
                                 51,496,520,640,466),
                   grp1 = c(3,4,4,4,#),
                            4,2,3,3,3),
                   grp2 = c(4,NA,NA,NA,#))
                            NA,3,4,NA,NA))

i <- 466
sub_bp <- x5 %>% filter(entity_id == i)
sub_bp <- sub_bp[order(sub_bp$bp),]
subdat <- dtrend_dat %>% filter(interp_age >= 7800 & interp_age <= 8400 & entity_id == i)

ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
  geom_vline(data = sub_bp, aes(xintercept = bp), colour = "red")

subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp,max(interp_age)), labels = 1:6))
sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
sub.av <- aov(sub.lm)
summary(sub.av)

tukey_sub <- as.data.frame(TukeyHSD(sub.av)$grp)
tukey_sub[order(abs(tukey_sub$diff)),] # grp 4 most sig diff
xx5 <- data.frame()
for (i in unique(x5$entity_id)){
  if(i == 591){ next }
  sub_bp <- x5 %>% filter(entity_id == i)
  sub_bp <- sub_bp[order(sub_bp$bp),]
  subdat <- dtrend_dat %>% filter(interp_age >= 7800 & interp_age <= 8400 & entity_id == i)
  
  if (i %in% c(150)){
    subdat$grp <- cut(subdat$interp_age, breaks = c(sub_bp$bp,max(subdat$interp_age)), labels = 1:5)
  } else {
    subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp,max(interp_age)), labels = 1:6))
  }
  event_grp <- grps[which(grps$entity_id == i),2:3]
  
  event_d18O <- subdat %>% filter(grp %in% event_grp) %>% summarise(mean_d18O = mean(d18O_detrended))
  base_d18O <- subdat %>% filter(!grp %in% event_grp) %>% summarise(mean_d18O = mean(d18O_detrended))
  x_sub <- data.frame(unique(subdat[which(subdat$entity_id == i),1:5]),
                      min_bp = sub_bp$bp[min(event_grp,na.rm = T)-1], max_bp = sub_bp$bp[max(event_grp, na.rm = T)],
                      sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[min(event_grp,na.rm = T)-1])),"sample_id"],
                      sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[max(event_grp, na.rm = T)])),"sample_id"],
                      anom = as.numeric(event_d18O) - as.numeric(base_d18O),
                      d18O_base = as.numeric(base_d18O), d18O_event = as.numeric(event_d18O))
  
  xx5 <- rbind(xx5, x_sub)
}

xx5$duration <- xx5$max_bp - xx5$min_bp
xx5$event_centre <- xx5$max_bp - (xx5$duration/2)



## 6bp
x6 <- x %>% group_by(entity_id) %>% mutate(n_bp = n()) %>% filter(n_bp == 6) # 2 entities

xx6 <- data.frame()
for (i in unique(x6$entity_id)){
  sub_bp <- x6 %>% filter(entity_id == i)
  sub_bp <- sub_bp[order(sub_bp$bp),]
  subdat <- dtrend_dat %>% filter(interp_age >= 7800 & interp_age <= 8400 & entity_id == i)
  
  subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp,max(interp_age)), labels = 1:7))
  
  if (i == 608){ 
    event_d18O <- subdat %>% filter(grp == 3) %>% summarise(mean_d18O = mean(d18O_detrended))
    base_d18O <- subdat %>% filter(grp != 3) %>% summarise(mean_d18O = mean(d18O_detrended))
    x_sub <- data.frame(unique(subdat[which(subdat$entity_id == i),1:5]),
                        min_bp = sub_bp$bp[2], max_bp = sub_bp$bp[3],
                        sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[2])),"sample_id"],
                        sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[3])),"sample_id"],
                        anom = as.numeric(event_d18O) - as.numeric(base_d18O),
                        d18O_base = as.numeric(base_d18O),
                        d18O_event = as.numeric(event_d18O))
  } else {
    event_d18O <- subdat %>% filter(grp %in% c(2:6)) %>% summarise(mean_d18O = mean(d18O_detrended))
    base_d18O <- subdat %>% filter(!grp %in% c(2:6)) %>% summarise(mean_d18O = mean(d18O_detrended))
    x_sub <- data.frame(unique(subdat[which(subdat$entity_id == i),1:5]),
                        min_bp = sub_bp$bp[2], max_bp = sub_bp$bp[6],
                        sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[2])),"sample_id"],
                        sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[6])),"sample_id"],
                        anom = as.numeric(event_d18O) - as.numeric(base_d18O),
                        d18O_base = as.numeric(base_d18O),
                        d18O_event = as.numeric(event_d18O))
  }
  xx6 <- rbind(xx6, x_sub)
}

xx6$duration <- xx6$max_bp - xx6$min_bp
xx6$event_centre <- xx6$max_bp - (xx6$duration/2)


## 7bp
x7 <- x %>% group_by(entity_id) %>% mutate(n_bp = n()) %>% filter(n_bp == 7)

i = unique(x7$entity_id)[1] #1 entities
sub_bp <- x7 %>% filter(entity_id == i)
sub_bp <- sub_bp[order(sub_bp$bp),]
subdat <- dtrend_dat %>% filter(interp_age >= 7800 & interp_age <= 8400 & entity_id == i)

ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
  geom_vline(data = sub_bp, aes(xintercept = bp), colour = "red", lty = 2) +
  scale_y_reverse() +
  ggtitle(paste("entity_id =", i)) 

subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp,max(interp_age)), labels = 1:8))
sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
sub.av <- aov(sub.lm)
summary(sub.av)

tukey_sub <- as.data.frame(TukeyHSD(sub.av)$grp)
tukey_sub[order(abs(tukey_sub$diff)),] # grp 4 most sig diff

event_d18O <- subdat %>% filter(grp == 4) %>% summarise(mean_d18O = mean(d18O_detrended))
base_d18O <- subdat %>% filter(grp != 4) %>% summarise(mean_d18O = mean(d18O_detrended))
xx7 <- data.frame(unique(subdat[,1:5]),
                  min_bp = sub_bp$bp[3], max_bp = sub_bp$bp[4],
                  sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[3])),"sample_id"],
                  sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[4])),"sample_id"],
                  anom = as.numeric(event_d18O) - as.numeric(base_d18O),
                  d18O_base = as.numeric(base_d18O),
                  d18O_event = as.numeric(event_d18O))  #site_id, entity_id, lo, lat, min_bp, max_bp, anom, duration, event_centre
xx7$duration <- xx7$max_bp - xx7$min_bp
xx7$event_centre <- xx7$max_bp - (xx7$duration/2)

## 8bp
x8 <- x %>% group_by(entity_id) %>% mutate(n_bp = n()) %>% filter(n_bp == 8) #0

i = unique(x8$entity_id)[1] #1 entities
sub_bp <- x8 %>% filter(entity_id == i)
sub_bp <- sub_bp[order(sub_bp$bp),]
subdat <- dtrend_dat %>% filter(interp_age >= 7800 & interp_age <= 8400 & entity_id == i)

ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
  geom_vline(data = sub_bp, aes(xintercept = bp), colour = "red", lty = 2) +
  scale_y_reverse() +
  ggtitle(paste("entity_id =", i)) 

subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp,max(interp_age)), labels = 1:9))
sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
sub.av <- aov(sub.lm)
summary(sub.av)

tukey_sub <- as.data.frame(TukeyHSD(sub.av)$grp)
tukey_sub[order(abs(tukey_sub$diff)),] # grps 6 and 7 most sig diff

event_d18O <- subdat %>% filter(grp %in% c(6,7)) %>% summarise(mean_d18O = mean(d18O_detrended))
base_d18O <- subdat %>% filter(!grp %in% c(6,7)) %>% summarise(mean_d18O = mean(d18O_detrended))
xx8 <- data.frame(unique(subdat[,1:4]),
                  min_bp = sub_bp$bp[5], max_bp = sub_bp$bp[7],
                  anom = as.numeric(event_d18O) - as.numeric(base_d18O))  #site_id, entity_id, lo, lat, min_bp, max_bp, anom, duration, event_centre
xx8$duration <- xx8$max_bp - xx8$min_bp
xx8$event_centre <- xx8$max_bp - (xx8$duration/2)


## all entities - 8.2 ka event
all_dat <- rbind(anom_2bp, anom_3bp, xx4, xx5, xx6, xx7)
write.csv(all_dat, "abrupt_Holocene/spel_82_signals.csv", row.names = F)
ggplot(data = all_dat, aes(x = as.factor(site_id), ymin = min_bp, ymax = max_bp)) + geom_errorbar() +
  geom_hline(yintercept = mean(all_dat$event_centre), col = "red") + geom_hline(yintercept = c(mean(all_dat$min_bp), mean(all_dat$max_bp)), col = "red", lty = 2)

#all_dat <- all_dat %>% filter(!site_id %in% c(8,182))

min_bp_dat <- data.frame(grp = c("Greenland ice cores","mean","median"), 
                         val = c(8080,mean(all_dat$min_bp),median(all_dat$min_bp)))
ggplot(data  =all_dat, aes(x = min_bp)) + geom_density() + 
  geom_vline(data = min_bp_dat, aes(xintercept = val, col = grp, lty = grp)) +
  theme(legend.title = element_blank()) +
  xlab("event end (years BP)")

ggplot(data = all_dat, aes(x = longitude, y = latitude, fill = min_bp)) +
  geom_point(shape = 21, size = 4) + borders("world") + coord_fixed(xlim = c(-130,170), ylim = c(-45,55)) +
  colorspace::scale_fill_continuous_sequential(palette = "plasma") +
  labs(fill = "event end", x = "", y = "") +
  theme_bw()

#geom_vline(xintercept = median(all_dat$min_bp),lty = 2, col = "red") +
#geom_vline(xintercept = mean(all_dat$min_bp), lty = 3, col = "blue") +
#geom_vline(xintercept = 8080, col = "green") #+
#geom_vline(xintercept = density(all_dat$min_bp)$x[199])
#which.max(density(all_dat$min_bp)$y)
max_bp_dat <- data.frame(grp = c("Greenland ice cores","mean","median"), 
                         val = c(8245,mean(all_dat$max_bp),median(all_dat$max_bp)))
ggplot(data  =all_dat, aes(x = max_bp)) + geom_density() + 
  geom_vline(data = max_bp_dat, aes(xintercept = val, col = grp, lty = grp)) +
  theme(legend.title = element_blank()) +
  xlab("event start (years BP)")

ggplot(data = all_dat, aes(x = longitude, y = latitude, fill = max_bp)) +
  geom_point(shape = 21, size = 4) + borders("world") + coord_fixed(xlim = c(-130,170), ylim = c(-45,55)) +
  colorspace::scale_fill_continuous_sequential(palette = "plasma") +
  labs(fill = "event start", x = "", y = "") +
  theme_bw()
#geom_vline(xintercept = median(all_dat$max_bp),lty = 2, col = "red") +
#geom_vline(xintercept = mean(all_dat$max_bp), lty = 3, col = "blue") +
#geom_vline(xintercept = 8245, col = "green") #+
#geom_vline(xintercept = density(all_dat$max_bp)$x[342])
#which.max(density(all_dat$max_bp)$y)

duration_dat <- data.frame(grp = c("Greenland ice cores","mean","median"), 
                           val = c(160.5,mean(all_dat$duration),median(all_dat$duration)))
ggplot(data  =all_dat, aes(x = duration)) + geom_density() + 
  geom_vline(data = duration_dat, aes(xintercept = val, col = grp, lty = grp)) +
  theme(legend.title = element_blank())  

ggplot(data = all_dat, aes(x = longitude, y = latitude, fill = duration)) +
  geom_point(shape = 21, size = 4) + borders("world") + coord_fixed(xlim = c(-130,170), ylim = c(-45,55)) +
  colorspace::scale_fill_continuous_sequential(palette = "viridis") +
  labs(fill = "duration", x = "", y = "") +
  theme_bw()
#geom_vline(xintercept = median(all_dat$duration),lty = 2, col = "red") +
#geom_vline(xintercept = mean(all_dat$duration), lty = 3, col = "blue") +
#geom_vline(xintercept = 160.5, col = "green") +
#geom_vline(xintercept = density(all_dat$duration)$x[282])
#which.max(density(all_dat$duration)$y)


## get 8.2 ka hiatus data
hiatus_8_2 <- unique(Raw_Data[1:9]) %>% filter(entity_id %in% c(109,291,415,565))
write.csv(hiatus_8_2, "abrupt_Holocene/spel_hiatus_8_2.csv", row.names = F)

ggplot() + geom_point(data = all_dat, aes(x = longitude, y = latitude, fill = anom), shape = 21, size = 4) +
  borders("world") + scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  geom_point(data = nosignal8_2, aes(x = longitude, y = latitude), shape = 4, size = 2, stroke = 1.5) +
  geom_text(data = hiatus_8_2, aes(x = longitude, y = latitude, label = "H")) +
  coord_fixed(ylim = c(-40,55), xlim = c(-130,170)) + ylab("") + xlab("")



## structure

# Asia
Asia <- all_dat %>% filter(longitude >= 70 & longitude <= 140 & latitude >= -5 & latitude <= 45)
Asia_dat <- dtrend_dat %>% filter(entity_id %in% Asia$entity_id & interp_age >= 7800 & interp_age <= 8400)
ggplot(data = filter(Asia_dat, entity_id %in% Asia$entity_id[7:9]), aes(x = interp_age, y = d18O_detrended, group = as.factor(entity_id), col = as.factor(entity_id))) + geom_line()

# Europe
Europe <- all_dat %>% filter(longitude >= -10 & longitude <= 45 & latitude >= 20 & latitude <= 50)
Europe_dat <- dtrend_dat %>% filter(entity_id %in% Europe$entity_id & interp_age >= 7800 & interp_age <= 8400)
ggplot(data = filter(Europe_dat, entity_id %in% Europe$entity_id[9:10]), aes(x = interp_age, y = d18O_detrended, group = entity_id, col = entity_id)) + geom_line()

#bp_pcent <- data.frame()
#for (i in seq(0,12000,500)){
#  subdat <- bp_out2 %>% filter(win_start == i)

# count for 100 year bins
#  wind_min <- i
#  wind_max <- i+1000
#  subdat$bin <- cut((subdat$bp), breaks = seq(wind_min,wind_max,100), labels = seq((wind_min+50),(wind_max-50),100))

#  subdat <- subdat %>% group_by(win_start, win_end, bin) %>% summarise(n_bp = n())

#  bp_pcent <- rbind(bp_pcent, subdat)
#}  

#bp_pcent <- left_join(bp_pcent, n_records)
#bp_pcent$pcent <- bp_pcent$n_bp/bp_pcent$n_entities

bp_pcent2 <- data.frame()
for (i in seq(1000,12000,500)){
  subdat <- bp_pcent %>% filter(win_end == i)
  subdat <- subdat[-c(1:2,(nrow(subdat)-1),(nrow(subdat))),]
  
  bp_pcent2 <- rbind(bp_pcent2, subdat)
}

ggplot(data = bp_pcent2, aes(x = bin, y = pcent, group = as.factor(win_start), col = as.factor(win_start))) + geom_line()


## focus in on 8.05
bp8_2 <- bp_out %>% filter(bp >= (8050-200) & bp <= (8050+200))

x <- 671
ggplot(data = filter(Raw_Data, entity_id == x & interp_age >= 7500 & interp_age <= 8500), aes(x = interp_age, y = d18O_measurement)) + 
  geom_line() +
  scale_y_reverse() +
  geom_vline(xintercept = filter(bp_out, entity_id == x & bp >= 7500 & bp <= 8500)$bp, col = "red", lty = 2) +
  ggtitle(paste("entity_id =", x))
