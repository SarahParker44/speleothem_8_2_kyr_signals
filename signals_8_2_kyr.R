## Constrain 8.2 kyr signals 

setwd("C:/Users/sarah/OneDrive/Documents/PhD/abrupt_Holocene/speleothem_8_2_kyr_signals/")

library(dplyr)
library(strucchange)
library(ggplot2)
library(RMySQL)

# connect to SISAL database
mydb <- dbConnect(MySQL(), user = "root", password = "BevRed921", dbname = "sisal_v2", 
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
no_signal <- data.frame()
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
  
  if (any(ci_x$confint[,c(1,3)] == 0) | any(ci_x$confint[,c(1,3)] < 0) | any(is.na(ci_x$confint[,c(1,3)]))){ 
    no_signal <- rbind(no_signal, unique(subdat[,c("site_id","site_name","entity_id","entity_name","longitude","latitude")])) 
    } else {
    # output breakpoints
    ci_ages <- data.frame(site_id = unique(subdat$site_id), site_name = unique(subdat$site_name), entity_id = unique(subdat$entity_id), entity_name = unique(subdat$entity_name),
                          longitude = unique(subdat$longitude), latitude = unique(subdat$latitude),
                          bp = subdat$interp_age[ci_x$confint[,2]],
                          CI2_5 = subdat$interp_age[ci_x$confint[,1]],
                          CI97_5 = subdat$interp_age[ci_x$confint[,3]]
                          )
    
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

## entities with 2 breakpoints
x2 <- x %>% filter(n_bp == 2) # 26

# plot each record individually
#for (i in seq(9,27,9)){
#  entities <- x2$entity_id[(i-8):i]
#  subdat <- dtrend_dat %>% filter(entity_id %in% entities)
#  sub_bp <- bp_dat %>% filter(entity_id %in% entities)
  
#  file_no <- ifelse(i == 9, 1, ifelse(i == 18, 2, 3))
#  filename <- paste("signals_82_", file_no, ".pdf", sep = "")

#  pdf(filename, width = 20/2.54, height = 18/2.54)
#  ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
#    geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
#    geom_vline(data = sub_bp, aes(xintercept = CI2_5), col = "red", lty = 2) +
#    geom_vline(data = sub_bp, aes(xintercept = CI97_5), col = "red", lty = 2) +
#    facet_wrap(.~ entity_id)
#  dev.off()
#}

# exclude entities: 53, 54, 295, 374, 388, 395, 540, 608, 690
#x2 <- x2 %>% filter(!entity_id %in% c(53, 54, 295, 374, 388, 395, 540, 608, 690)) # 20 entities
#x2 <- x2 %>% filter(!entity_id %in% c(51, 53, 295, 374, 690))

anom_2bp <- data.frame()
for (i in unique(x2$entity_id)){
  sub_bp <- bp_dat %>% filter(entity_id == i)
  subdat <- dtrend_dat %>% filter(entity_id == i)
  
  # sig diff from base:
  subdat <- subdat %>% arrange(interp_age)
  sub_bp <- sub_bp %>% arrange(bp)
  subdat$grp <- cut(subdat$interp_age, breaks = c(min(subdat$interp_age), sub_bp$bp, max(subdat$interp_age)), labels = 1:3)
  t_test <- t.test(x = subdat[which(subdat$grp == 2),"d18O_detrended"], y = subdat[which(subdat$grp %in% c(1,3)),"d18O_detrended"])
  
  # overall
  d18O_event <- subdat %>% filter(interp_age >= min(sub_bp$bp) & interp_age <= max(sub_bp$bp)) %>% summarise(mean(d18O_detrended))
  d18O_base <- subdat %>% filter(interp_age <= min(sub_bp$bp) | interp_age >= max(sub_bp$bp)) %>% summarise(mean(d18O_detrended))
  d18Osd_base <- subdat %>% filter(interp_age <= min(sub_bp$bp) | interp_age >= max(sub_bp$bp)) %>% summarise(sd(d18O_detrended))
  
  #lower limit uncert
  d18O_event_lower <- subdat %>% filter(interp_age >= min(sub_bp$CI2_5) & interp_age <= max(sub_bp$CI97_5)) %>% summarise(mean(d18O_detrended))
  d18O_base_lower <- subdat %>% filter(interp_age <= min(sub_bp$CI2_5) | interp_age >= max(sub_bp$CI97_5)) %>% summarise(mean(d18O_detrended))
  
  # upper limit uncert
  d18O_event_upper <- subdat %>% filter(interp_age >= min(sub_bp$CI97_5) & interp_age <= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
  d18O_base_upper <- subdat %>% filter(interp_age <= min(sub_bp$CI97_5) | interp_age >= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
  
  anom <- d18O_event - d18O_base
  anom_lower <- d18O_event_lower - d18O_base_lower
  anom_upper <- d18O_event_upper - d18O_base_upper

  sub_df <- data.frame(unique(sub_bp[,c(1:6)]), 
                       min_bp = as.numeric(min(sub_bp$bp)), max_bp = as.numeric(max(sub_bp$bp)),
                       sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-as.numeric(min(sub_bp$bp)))),"sample_id"],
                       sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-as.numeric(max(sub_bp$bp)))),"sample_id"],
                       anom = as.numeric(anom),
                       #anom_lower = as.numeric(anom_lower),
                       #anom_upper = as.numeric(anom_upper),
                       d18O_base = as.numeric(d18O_base),
                       d18O_event = as.numeric(d18O_event),
                       d18Osd_base = as.numeric(d18Osd_base),
                       ttest_Pval = t_test$p.value)
  
  anom_2bp <- rbind(anom_2bp, sub_df)
}

# filter to those with anomalies greater than the base s.d.
anom_2bp <- anom_2bp %>% mutate(diff_from_base_sd = round(abs(anom) - abs(d18Osd_base), digits = 1)) 
no_signal <- rbind(no_signal, filter(anom_2bp, diff_from_base_sd <= 0)[,c("site_id","site_name","entity_id","entity_name","longitude","latitude")])
anom_2bp <- anom_2bp %>% filter(diff_from_base_sd > 0)

# filter to those that are sig diff from base
no_signal <- rbind(no_signal, filter(anom_2bp, ttest_Pval >= 0.001)[,c("site_id","site_name","entity_id","entity_name","longitude","latitude")])
anom_2bp <- anom_2bp %>% filter(ttest_Pval <= 0.001)

#subdat <- dtrend_dat %>% filter(entity_id %in% unique(anom_2bp$entity_id)[c(1,6,7,10,12,15)])
#sub_bp <- bp_dat %>% filter(entity_id %in% unique(anom_2bp$entity_id)[c(1,6,7,10,12,15)]) 

#ggplot() +
#  geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
#  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red", lty = 2) +
#  facet_wrap(.~ as.factor(entity_id), scales = "free_y") +
#  ggtitle("2 breakpoints")

#anom_2bp <- anom_2bp %>% filter(!entity_id %in% c(51, 53, 295, 374, 449, 690))

#xx <- anom_2bp %>% filter(entity_id %in% c(53, 54, 295, 374, 388, 395, 540, 608, 690))

#ggplot(data = anom_2bp, aes(x = longitude, y = latitude, fill = anom)) +
#  geom_point(shape = 21) +
#  borders("world") +
#  scale_fill_gradient2(high = "red", mid = "white", low = "blue")


## grps
#grps <- read.csv("82_signal_grps.csv")
#colnames(grps)[1] <- "entity_id"

## entities with 3 breakpoints
x3 <- x %>% filter(n_bp == 3) # 23

#subdat <- dtrend_dat %>% filter(entity_id %in% unique(x3$entity_id)[c(4,11,17,19,22,23)])
#sub_bp <- bp_dat %>% filter(entity_id %in% unique(x3$entity_id)[c(4,11,17,19,22,23)]) 

#ggplot() +
#  geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
#  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red", lty = 2) +
#  facet_wrap(.~ as.factor(entity_id), scales = "free_y") +
#  ggtitle("3 breakpoints")

# visualise entities
#for (i in seq(9,28,9)){
#  entities <- x3$entity_id[(i-8):i]
#  subdat <- dtrend_dat %>% filter(entity_id %in% entities)
#  sub_bp <- bp_dat %>% filter(entity_id %in% entities)
  
#  file_no <- ifelse(i == 9, 1, ifelse(i == 18, 2, 3))
#  filename <- paste("signals_82_3bp_", file_no, ".pdf", sep = "")
  
#  p <- ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
#    geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
#    facet_wrap(.~ entity_id)
#  pdf(filename, width = 20/2.54, height = 18/2.54)
#  print(p)
#  dev.off()
#}

# 
grps <- data.frame()
for (i in unique(x3$entity_id)){
  subdat <- dtrend_dat %>% filter(entity_id == i) %>% arrange(interp_age)
  sub_bp <- bp_dat %>% filter(entity_id == i) %>% arrange(bp)
  
  subdat$grp <- with(subdat, ifelse(interp_age <= sub_bp$bp[1] | interp_age >= sub_bp$bp[3], 1,
                       ifelse(interp_age >= sub_bp$bp[1] & interp_age <= sub_bp$bp[2], 2, 3)))
  subdat$grp <- as.factor(subdat$grp)
  sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
  sub.av <- aov(sub.lm)
  
  sub_hsd <- TukeyHSD(sub.av)$grp
  #sub_hsd2 <- HSD.test(sub.av, trt = 'grp', alpha = 0.001)
  
  ## QC checks
  sub_hsd2 <- as.data.frame(sub_hsd[which(grepl("1", rownames(sub_hsd))),])
  sub_hsd3 <- t(as.data.frame(sub_hsd[which(!grepl("1", rownames(sub_hsd))),]))
  
  n_sig <- nrow(sub_hsd2[which(sub_hsd2[,4] <= 0.001),])
  
  if (n_sig == 0){
    no_signal <- rbind(no_signal, unique(subdat[,c("site_id","site_name", "entity_id", "entity_name", "longitude", "latitude")]))
  } else if (n_sig == 1){
    if (sub_hsd3[,4] < 0.001){
      event_grp <- rownames(sub_hsd2[which(sub_hsd2[,4] <= 0.001),])
      event_grp <- as.numeric(sub("-1", "", event_grp))
      grps <- rbind(grps, data.frame(entity_id = i, group = event_grp))
    } else {
      grps <- rbind(grps, data.frame(entity_id = rep(i,2), group = c(2,3)))
    }
  } else {
    if (sub_hsd3[,4] < 0.001){
      grp_2 <- mean(subdat[which(subdat$grp == 2),"d18O_detrended"])
      grp_3 <- mean(subdat[which(subdat$grp == 3),"d18O_detrended"])
      
      if (all(c(grp_2, grp_3) > 0) | all(c(grp_2, grp_3) < 0)){
        grps <- rbind(grps, data.frame(entity_id = rep(i, 2), group = c(2,3)))
      } else {
        base_mean <- mean(subdat[which(subdat$grp == 1),"d18O_detrended"])
        anom_grp2 <- round(grp_2 - base_mean, digits = 1)
        anom_grp3 <- round(grp_3 - base_mean, digits = 1)
        
        if (abs(anom_grp2) > abs(anom_grp3)){ event_grp <- 2} else if (abs(anom_grp2) < abs(anom_grp3)) { 
          event_grp <- 3 } else {
            event_grp <- rownames(sub_hsd2[which(sub_hsd2[,4] == min(sub_hsd2$`p adj`)),])
            event_grp <- as.numeric(sub("-1", "", event_grp))
          }
        grps <- rbind(grps, data.frame(entity_id = i, group = event_grp))
      }
    }
  }
  
  base_grp <- sub_hsd2$groups[which(rownames(sub_hsd2$groups) == 1),"groups"]
  # 
  if (length(which(sub_hsd2$groups$groups != base_grp)) == 1){
    event_grp <- as.numeric(rownames(sub_hsd2$groups[which(sub_hsd2$groups$groups != base_grp),]))
  }
}

# exclude entities: 63, 142, 436, 279, 546, 613
#x3 <- x3 %>% filter(!entity_id %in% c(63, 142, 279, 436, 546, 613)) # 14 entities



## entities with 4 breakpoints
x4 <- x %>% filter(n_bp == 4) # 9

#subdat <- dtrend_dat %>% filter(entity_id %in% unique(x4$entity_id))
#sub_bp <- bp_dat %>% filter(entity_id %in% unique(x4$entity_id)) 

#ggplot() +
#  geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
#  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red", lty = 2) +
#  facet_wrap(.~ as.factor(entity_id), scales = "free_y") +
#  ggtitle("4 breakpoints")

# visualise
#subdat <- dtrend_dat %>% filter(entity_id %in% x4$entity_id)
#sub_bp <- bp_dat %>% filter(entity_id %in% x4$entity_id)

#filename <- "signals_82_4bp.pdf"

#p <- ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
#  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
#  facet_wrap(.~ entity_id)
#pdf(filename, width = 20/2.54, height = 18/2.54)
#print(p)
#dev.off()

#
for (i in unique(x4$entity_id)){
  subdat <- dtrend_dat %>% filter(entity_id == i) %>% arrange(interp_age)
  sub_bp <- bp_dat %>% filter(entity_id == i) %>% arrange(bp)
  
  subdat$grp <- with(subdat, ifelse(interp_age <= sub_bp$bp[1] | interp_age >= sub_bp$bp[4], 1,
                                    ifelse(interp_age >= sub_bp$bp[1] & interp_age <= sub_bp$bp[2], 2, 
                                           ifelse(interp_age >= sub_bp$bp[2] & interp_age <= sub_bp$bp[3], 3, 4))))
  subdat$grp <- as.factor(subdat$grp)
  sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
  sub.av <- aov(sub.lm)
  
  sub_hsd <- TukeyHSD(sub.av)$grp
  #sub_hsd2 <- HSD.test(sub.av, trt = 'grp', alpha = 0.001)
  
  ## QC checks
  sub_hsd2 <- as.data.frame(sub_hsd[which(grepl("1", rownames(sub_hsd))),])
  sub_hsd3 <- as.data.frame(sub_hsd[which(!grepl("1", rownames(sub_hsd))),])
  
  n_sig <- nrow(sub_hsd2[which(sub_hsd2[,4] <= 0.001),])
  
  if (n_sig == 0){ # 220
    no_signal <- rbind(no_signal, unique(subdat[,c("site_id","site_name", "entity_id", "entity_name", "longitude", "latitude")]))
    
  } else if (n_sig == 1){ #244
    event_grp <- rownames(sub_hsd2[which(sub_hsd2[,4] <= 0.001),])
    event_grp <- as.numeric(sub("-1", "", event_grp))
    grps <- rbind(grps, data.frame(entity_id = i, group = event_grp))
    
  } else if (n_sig == 2){
    if (nrow(sub_hsd3[which(sub_hsd3$`p adj` <= 0.001),]) == 0){ # all events are insig from one another, ent = 52
      grps <- rbind(grps, data.frame(entity_id = rep(i, 3), group = 2:4))
      } else { #327, 351, 442
        event_grp <- rownames(sub_hsd2[which(sub_hsd2[,4] <= 0.001),])
        event_grp <- as.numeric(sub("-1","", event_grp))
        grps <- rbind(grps, data.frame(entity_id = rep(i, length(event_grp)+1), group = event_grp[1]:event_grp[2]))
      }
    
  } else if (n_sig == 3){
    if (nrow(sub_hsd3[which(sub_hsd3$`p adj` <= 0.001),]) == 3){
      grp_2 <- mean(subdat[which(subdat$grp == 2),"d18O_detrended"])
      grp_3 <- mean(subdat[which(subdat$grp == 3),"d18O_detrended"])
      grp_4 <- mean(subdat[which(subdat$grp == 4),"d18O_detrended"])
      
      if (all(c(grp_2, grp_3, grp_4) > 0) | all(c(grp_2, grp_3, grp_4) < 0)){ #254
        grps <- rbind(grps, data.frame(entity_id = rep(i, 3), group = 2:4))
      }
    } else {
        event_grp <- rownames(sub_hsd3[which(sub_hsd3[,4] >= 0.001),])
        event_grp <- as.numeric(as.vector(strsplit(event_grp, "-"))[[1]])
        grps <- rbind(grps, data.frame(entity_id = rep(i, length(event_grp)+1), group = event_grp[1]:event_grp[2]))
    }
  }
}


# exclude entities: 220, 351, 415 = hiatus
#x4 <- x4 %>% filter(!entity_id %in% c(220,351,415))



## 5 bp's
x5 <- x %>% filter(n_bp == 5) # 2

# visualise
#subdat <- dtrend_dat %>% filter(entity_id %in% x5$entity_id)
#sub_bp <- bp_dat %>% filter(entity_id %in% x5$entity_id)

#filename <- "signals_82_5bp.pdf"

#p <- ggplot() + geom_line(data = subdat, aes(x = interp_age, y = d18O_detrended)) +
#  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
#  facet_wrap(.~ entity_id)
#pdf(filename, width = 20/2.54, height = 18/2.54)
#print(p)
#dev.off()

subdat <- dtrend_dat %>% filter(entity_id == 591)
sub_bp <- bp_dat %>% filter(entity_id == 591)

subdat$grp <- with(subdat, ifelse(interp_age <= sub_bp$bp[1] | interp_age >= sub_bp$bp[5], 1,
                                  ifelse(interp_age >= sub_bp$bp[1] & interp_age <= sub_bp$bp[2], 2, 
                                         ifelse(interp_age >= sub_bp$bp[2] & interp_age <= sub_bp$bp[3], 3,
                                                ifelse(interp_age >= sub_bp$bp[3] & interp_age <= sub_bp$bp[4], 4, 5)))))
subdat$grp <- as.factor(subdat$grp)
sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
sub.av <- aov(sub.lm)

sub_hsd <- TukeyHSD(sub.av)$grp

grps <- rbind(grps, data.frame(entity_id = rep(591,2), group = c(3,4)))

#ggplot() + geom_line(data =subdat, aes(x = interp_age, y = d18O_detrended)) + 
#  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red") +
#  facet_wrap(.~ entity_id)


# calculate anomalies
anom_bp <- data.frame()
x_all <- rbind(x3,x4,x5)
for (i in unique(grps$entity_id)){
  sub_bp <- bp_dat %>% filter(entity_id == i) %>% arrange(bp)
  #sub_bp <- sub_bp[order(sub_bp$bp),]
  subdat <- dtrend_dat %>% filter(interp_age >= 7400 & interp_age <= 9000 & 
                                    entity_id == i)
  
  if (nrow(sub_bp) == 3){ lab_vec = 1:4 } else if (nrow(sub_bp) == 4) { lab_vec = 1:5 } else { lab_vec = 1:6}
  subdat$grp <- with(subdat, cut(interp_age, breaks = c(min(interp_age),sub_bp$bp,max(interp_age)), labels = lab_vec))
  
  event_grp <- as.numeric(grps[which(grps$entity_id == i),2])
  event_grp <- na.omit(event_grp)
  
  event_d18O <- subdat %>% filter(grp %in% event_grp) %>% summarise(mean_d18O = mean(d18O_detrended))
  base_d18O <- subdat %>% filter(grp %in% c(min(lab_vec), max(lab_vec))) %>% summarise(mean_d18O = mean(d18O_detrended))
  d18Osd_base <- subdat %>% filter(grp %in% c(min(lab_vec), max(lab_vec))) %>% summarise(d18Osd_base = sd(d18O_detrended))
  #base_d18O <- subdat %>% filter(interp_age >= 7400 & interp_age <= 7900 | interp_age >= 8500 & interp_age <= 9000) %>% summarise(mean_d18O = mean(d18O_detrended))
  x_sub <- data.frame(unique(subdat[which(subdat$entity_id == i),1:6]),
                      min_bp = sub_bp$bp[min(event_grp,na.rm = T)-1], max_bp = sub_bp$bp[max(event_grp, na.rm = T)], 
                      sample_id_min_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[min(event_grp,na.rm = T)-1])),"sample_id"],
                      sample_id_max_bp = subdat[which.min(abs(subdat$interp_age-sub_bp$bp[max(event_grp, na.rm = T)])),"sample_id"],
                      anom = as.numeric(event_d18O) - as.numeric(base_d18O),
                      d18O_base = as.numeric(base_d18O), d18O_event = as.numeric(event_d18O),
                      d18Osd_base = as.numeric(d18Osd_base$d18Osd_base))
  anom_bp <- rbind(anom_bp, x_sub)
}

# filter to those with anomalies greater than the base s.d.
anom_bp <- anom_bp %>% mutate(diff_from_base_sd = round(abs(anom) - abs(d18Osd_base), digits = 1)) 
no_signal <- rbind(no_signal, filter(anom_bp, diff_from_base_sd <= 0)[,c("site_id","site_name","entity_id","entity_name","longitude","latitude")])
anom_bp <- anom_bp %>% filter(diff_from_base_sd > 0)


# combine
all_dat <- rbind(anom_2bp[,-c(14:17)], anom_bp[,c(1:13)])
all_dat$duration <- all_dat$max_bp - all_dat$min_bp
all_dat$event_centre <- all_dat$max_bp - (all_dat$duration/2)


## add non-SISAL records
nonSISAL <- read.csv("C:/Users/sarah/OneDrive/Documents/PhD/abrupt_Holocene/nonSISAL_82_signals.csv")

all_dat <- rbind(all_dat, nonSISAL)


## save
write.csv(all_dat, "spel_82_signals.csv", row.names = F)
all_dat <- read.csv("spel_82_signals.csv")


### Make table of results for supplement
supp_out <- all_dat %>% select(entity_id, site_name, longitude, latitude, max_bp, min_bp, duration, anom) %>% 
  arrange(entity_id, site_name)
colnames(supp_out)[5:6] <- c("start","end")

write.csv(supp_out, "anom_table.csv", row.names = F)

## entities with no bp
#nosignal <- Raw_Data %>% filter(entity_id %in% c(53, 54, 295, 374, 395, 540, 690, 63, 142, 279, 436, 546, 613, 220,351, 305))
#no_signal <- unique(nosignal[c("entity_id","site_id","site_name","latitude","longitude")])
#no_signal2 <- Raw_Data %>% filter(entity_id %in% c(51, 53, 295, 374, 449, 690, 96,495,563,305))
#no_signal2 <- unique(no_signal2[,c("site_id", "site_name", "entity_id", "entity_name", "longitude", "latitude")])
#no_signal <- rbind(no_signal, no_signal2)

ggplot() + 
  #geom_point(data = no_signal, aes(x = longitude, y = latitude)) +
  geom_point(data = all_dat, aes(x = longitude, y = latitude, fill = anom), shape = 21, size = 3) + borders("world") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  coord_fixed(ylim = c(-50,60), xlim = c(-120,150))

write.csv(no_signal, "C:/Users/sarah/OneDrive/Documents/PhD/abrupt_Holocene/spel_nosignal_82.csv", row.names = F)


## Do sites with multiple entities show agreement amongst records
mult_ent <- all_dat %>% group_by(site_id, site_name) %>% mutate(n_ent = n()) %>% filter(n_ent > 1)

# Q1 are signals is the same direction - e.g. both +ve, both -ve
mult_ent2 <- data.frame()
for (i in unique(mult_ent$site_id)){
  subdat <- mult_ent %>% filter(site_id == i)
  
  if (all(subdat$anom >= 0) | all(subdat$anom <= 0)){
    subdat$signals_agree <- "yes"
  } else {
    subdat$signals_agree <- "no"
  }
  
  mult_ent2 <- rbind(mult_ent2, subdat)
}

# Q2 are signals reasonably contemporaneous (at least overlapping in time)?
mult_ent3 <- data.frame()
for (i in unique(mult_ent2$site_id)){
  subdat <- mult_ent2 %>% filter(site_id == i)
  
  if (min(subdat$max_bp) > max(subdat$min_bp)){
    subdat$overlapping <- "yes"
  } else {
    subdat$overlapping <- "no"
  }
  mult_ent3 <- rbind(mult_ent3, subdat)
}


## timings
ggplot(data = all_dat, aes(x = min_bp)) + geom_density() +
  geom_vline(mapping = aes(xintercept = 8088), col = "red") + 
  geom_vline(mapping = aes(xintercept = median(all_dat$min_bp)))
ggplot(data = all_dat, aes(x = max_bp)) + geom_density() + 
  geom_vline(mapping = aes(xintercept = 8248), col = "red") +
  geom_vline(mapping = aes(xintercept = median(all_dat$max_bp)))


Europe <- all_dat %>% filter(latitude >= 20 & latitude <= 50 & longitude >= -10 & longitude <= 45)
Asia <- all_dat %>% filter(latitude >= 0 & latitude <= 45 & longitude >= 50 & longitude <= 150)
S_America <- all_dat %>% filter(latitude >= -30 & latitude <= 0 & longitude >= -100 & longitude <= -30)
global <- all_dat 

Europe$region <- "Europe"
Asia$region <- "Asia"
S_America$region <- "S_America"
global$region <- "global"

all_dat2 <- rbind(Europe, Asia, S_America, global)

all_dat2$region <- factor(all_dat2$region, levels = c("Europe","S_America","Asia","global"))

ggplot(data = all_dat2, aes(x = region, y = abs(anom))) +
  geom_boxplot() +
  ylab("magnitude") +
  ggtitle("regional differences in event magnitude")

ggplot(data = all_dat2, aes(x = region, y = duration)) +
  geom_boxplot() +
  ylab("duration") +
  ggtitle("regional differences in event duration")

ggplot(data = all_dat2, aes(x = region, y = max_bp)) +
  geom_boxplot() +
  ylab("start timing") +
  ggtitle("regional differences in timing")


t.test(x = Asia$max_bp, y = S_America$max_bp)
t.test(x = Asia$min_bp, y = S_America$min_bp)

t.test(x = abs(Europe$anom), y = abs(Asia$anom))
