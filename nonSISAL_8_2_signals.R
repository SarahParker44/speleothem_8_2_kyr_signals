### breakpoint analysis for non SISAL records

setwd("C:/Users/sarah/OneDrive/Documents/PhD/abrupt_Holocene/")

library(tidyr)
library(dplyr)
library(strucchange)

## Load data
Padre <- read.csv("Padre.csv"); colnames(Padre)[1] <- "Age"; Padre$site_name <- "Padre"
Paixao <- read.csv("Paixao.csv"); colnames(Paixao)[1] <- "Age"; Paixao$site_name <- "Paixao"
Tigre_perdido <- read.csv("tigre-perdido2008.csv"); colnames(Tigre_perdido)[1] <- "Age"; Tigre_perdido$site_name <- "Tigre Perdido"
Kaite <- read.csv("Kaite.csv"); colnames(Kaite)[1] <- "Age"; Kaite$site_name <- "Kaite"
Wuya <- read.csv("Wuya.csv"); colnames(Wuya)[1] <- "Age"; Wuya$site_name <- "Wuya"
Klang <- read.csv("Klang.csv"); colnames(Klang)[1] <- "Age"; Klang$site_name <- "Klang"

dat <- rbind(Padre, Paixao, Tigre_perdido, Kaite[,-2],Wuya, Klang)
dat <- na.omit(dat)
dat <- dat %>% filter(Age <10000)

## Breakpoint analysis
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

# for every 1000 year window (with 50% overlap)

bp_out <- data.frame()
dtrend_dat <- data.frame()


for (i in unique(dat$site_name)){
  subdat <- dat %>% filter(site_name == i)
  
  age_lower <- floor(min(subdat$Age)/1000)*1000
  age_upper <- ceiling(max(subdat$Age)/1000)*1000
  
  for (j in seq(age_lower+1000, age_upper, 500)){
    window_dat <- subdat %>% filter(Age >= (j-1000) & Age <= j)
    
    if (nrow(window_dat) == 0){ next }
    
    ## detrend
    subdat_lm <- lm(d18O ~ Age, data = window_dat)
    lm_predicted <- predict(subdat_lm)
    window_dat$detrended_d18O <- residuals(subdat_lm)
    
    # save detrended
    dtrend_dat <- rbind(dtrend_dat, window_dat)
    
    ## breakpoints
    bp <- breakpoints(window_dat$detrended_d18O ~ 1)#subdat$interp_age) #breakpoint analysis
    bpts_sum <- summary(bp) 
    opt_brks <- opt_bpts(bpts_sum$RSS["BIC",]) #optimal no. breakpoints
    opt_brks <- na.omit(opt_brks)
    
    if (length(opt_brks) > 1){ opt_brks <- as.numeric(names(which.min(bpts_sum$RSS["BIC",]))) }
    if (length(opt_brks) == 0){ next }
    if (is.na(opt_brks)){ next }
    
    ci_x <- confint(bp, breaks = opt_brks) #get timings of bp's with conf intervals
    
    # output breakpoints
    ci_ages <- data.frame(
      site_name = unique(subdat$site_name),
      bp = window_dat$Age[ci_x$confint[,2]],
      CI2_5 = window_dat$Age[ci_x$confint[,1]],
      CI97_5 = window_dat$Age[ci_x$confint[,3]])
    
    # save breakpoints
    bp_out <- rbind(bp_out, ci_ages)
    
  }
}

p <- ggplot() +
  geom_line(data = filter(dtrend_dat, site_name == "Klang"), aes(x = Age, y = detrended_d18O)) +
  geom_vline(data = filter(bp_out, site_name == "Klang"), aes(xintercept = bp), col = "red") +
  scale_y_reverse()
  #scale_x_continuous(breaks = seq(1000,10000,1000))

## Remove replicated bp's (due to overlapping windows)
bp_out2 <- data.frame()
for (i in unique(bp_out$site_name)){
  subdat <- bp_out %>% filter(site_name == i)
  if (nrow(subdat) == 1){ subdat <- subdat
  } else {
    distmat <- dist(subdat$bp) # identify distances between all bp's
    distmat <- as.matrix(distmat)
    
    xy <- t(combn(colnames(distmat), 2))
    dist_df <- data.frame(xy, dist = distmat[xy])
    
    #identify pairs where distance <= 20 yrs
    small_dist <- dist_df %>% filter(dist <= 30)
    
    if (nrow(small_dist) == 0){ subdat <- subdat } else {
      bp_sub <- data.frame()
      for (j in 1:nrow(small_dist)){
        subsubdat <- subdat[as.numeric(small_dist[j,1:2]),]
        bp_sub <- rbind(bp_sub, 
                        data.frame(site_name = i,
                                   bp = mean(subsubdat$bp),
                                   CI2_5 = mean(subsubdat$CI2_5),
                                   CI97_5 = mean(subsubdat$CI97_5)))
      }
      # remove duplicates from subset
      subdat <- subdat[-as.numeric(as.matrix(small_dist[1:nrow(small_dist),1:2])),]
      #add combined duplicates to subdat
      subdat <- rbind(subdat,bp_sub)
      
    }
  }
  
  bp_out2 <- rbind(bp_out2, subdat)
  
}

ggplot() +
  geom_line(data = filter(dtrend_dat, site_name == "Klang"), aes(x = Age, y = detrended_d18O)) +
  geom_vline(data = filter(bp_out2, site_name == "Klang"), aes(xintercept = bp), col = "red") +
  scale_x_continuous(breaks = seq(1000,10000,1000))


## focusing on the 8.2 ka event
bp_8_2 <- bp_out2 %>% filter(bp >= 7800 & bp <= 8400)

# 2 bp's (Tigre Perdido, Paixao)
anom_2bp <- data.frame()
sub_bp <- bp_8_2 %>% filter(site_name == "Tigre Perdido")
subdat <- dtrend_dat %>% filter(Age >= 7800 & Age <= 8400 & site_name == "Tigre Perdido")
  
ggplot() + geom_line(data = subdat, aes(x = Age, y = detrended_d18O)) + 
    geom_vline(data = sub_bp, aes(xintercept = bp), col = "red", lty = 2) + 
    geom_errorbarh(data = sub_bp, aes(xmin = CI2_5, xmax = CI97_5, y = c(0.9, 0.95)), col = "red", height = 0.02) +
    ggtitle("Tigre Perdido")
  
d18O_event <- subdat %>% filter(Age >= min(sub_bp$bp) & Age <= max(sub_bp$bp)) %>% summarise(mean(detrended_d18O))
#event_CI25 <- subdat %>% filter(interp_age >= min(sub_bp$CI2_5) & interp_age <= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
d18O_base <- subdat %>% filter(Age <= min(sub_bp$bp) | Age >= max(sub_bp$bp)) %>% summarise(mean(detrended_d18O))
  
anom <- d18O_event - d18O_base
  
anom_2bp <- data.frame(site_name = "Tigre Perdido", 
                       min_bp = as.numeric(min(sub_bp$bp)), max_bp = as.numeric(max(sub_bp$bp)),
                       anom = as.numeric(anom),
                       d18O_base = as.numeric(d18O_base),
                       d18O_event = as.numeric(d18O_event))

sub_bp <- bp_8_2 %>% filter(site_name == "Paixao")
subdat <- dtrend_dat %>% filter(Age >= 7800 & Age <= 8400 & site_name == "Paixao")

ggplot() + geom_line(data = subdat, aes(x = Age, y = detrended_d18O)) + 
  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red", lty = 2) + 
  geom_errorbarh(data = sub_bp, aes(xmin = CI2_5, xmax = CI97_5, y = c(0.9, 0.95)), col = "red", height = 0.02) +
  ggtitle("Paixao")

d18O_event <- subdat %>% filter(Age >= min(sub_bp$bp) & Age <= max(sub_bp$bp)) %>% summarise(mean(detrended_d18O))
#event_CI25 <- subdat %>% filter(interp_age >= min(sub_bp$CI2_5) & interp_age <= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
d18O_base <- subdat %>% filter(Age <= min(sub_bp$bp) | Age >= max(sub_bp$bp)) %>% summarise(mean(detrended_d18O))

anom <- d18O_event - d18O_base

anom_2bp <- rbind(anom_2bp, data.frame(site_name = "Paixao", 
                       min_bp = as.numeric(min(sub_bp$bp)), max_bp = as.numeric(max(sub_bp$bp)),
                       anom = as.numeric(anom),
                       d18O_base = as.numeric(d18O_base),
                       d18O_event = as.numeric(d18O_event)))

anom_2bp$duration <- anom_2bp$max_bp - anom_2bp$min_bp
anom_2bp$event_centre <- anom_2bp$max_bp - (anom_2bp$duration/2)

# examine sig difference of groups (separated by bp)

# 3 bp's (Kaite, Klang)
anom_3bp <- data.frame()
sub_bp <- bp_8_2 %>% filter(site_name == "Kaite")
subdat <- dtrend_dat %>% filter(Age >= 7800 & Age <= 8400 & site_name == "Kaite")

ggplot() + geom_line(data = subdat, aes(x = Age, y = detrended_d18O)) + 
  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red", lty = 2) + 
  geom_errorbarh(data = sub_bp, aes(xmin = CI2_5, xmax = CI97_5, y = c(0.9, 0.95, 0.97)), col = "red", height = 0.02) +
  ggtitle("Kaite")

sub_bp <- sub_bp %>% arrange(bp)
subdat$grp <- with(subdat, cut(Age, breaks = c(min(Age),sub_bp$bp,max(Age)), labels = 1:4))

d18O.lm <- lm(detrended_d18O ~ grp, data = na.omit(subdat))
d18O.av <- aov(d18O.lm)

tukey_sub <- TukeyHSD(x = d18O.av, conf.level = 0.99)
plot(tukey_sub , las=1 , col="brown")

# grp 2 is only sig diff from background values (grps 1 and 4)
d18O_event <- subdat %>% filter(Age >= sub_bp$bp[1] & Age <= sub_bp$bp[2]) %>% summarise(mean(detrended_d18O))
#event_CI25 <- subdat %>% filter(interp_age >= min(sub_bp$CI2_5) & interp_age <= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
d18O_base <- subdat %>% filter(Age >= min(subdat$Age) & Age <= sub_bp$bp[1] | Age >= sub_bp$bp[3] & Age <= max(subdat$Age)) %>% summarise(mean(detrended_d18O))

anom <- d18O_event - d18O_base

anom_3bp <- data.frame(site_name = "Kaite", 
                       min_bp = as.numeric(sub_bp$bp[1]), max_bp = as.numeric(sub_bp$bp[2]),
                       anom = as.numeric(anom),
                       d18O_base = as.numeric(d18O_base),
                       d18O_event = as.numeric(d18O_event))



sub_bp <- bp_8_2 %>% filter(site_name == "Klang")
subdat <- dtrend_dat %>% filter(Age >= 7800 & Age <= 8400 & site_name == "Klang")

ggplot() + geom_line(data = subdat, aes(x = Age, y = detrended_d18O)) + 
  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red", lty = 2) + 
  geom_errorbarh(data = sub_bp, aes(xmin = CI2_5, xmax = CI97_5, y = c(0.9, 0.95, 0.97)), col = "red", height = 0.02) +
  ggtitle("Klang")

sub_bp <- sub_bp %>% arrange(bp)
subdat$grp <- with(subdat, cut(Age, breaks = c(min(Age),sub_bp$bp,max(Age)), labels = 1:4))

d18O.lm <- lm(detrended_d18O ~ grp, data = na.omit(subdat))
d18O.av <- aov(d18O.lm)

tukey_sub <- TukeyHSD(x = d18O.av, conf.level = 0.99)
plot(tukey_sub , las=1 , col="brown")

# grp 3 is only sig diff from background values (grps 1 and 4)
d18O_event <- subdat %>% filter(Age >= sub_bp$bp[2] & Age <= sub_bp$bp[3]) %>% summarise(mean(detrended_d18O))
#event_CI25 <- subdat %>% filter(interp_age >= min(sub_bp$CI2_5) & interp_age <= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
d18O_base <- subdat %>% filter(Age >= min(subdat$Age) & Age <= sub_bp$bp[1] | Age >= sub_bp$bp[3] & Age <= max(subdat$Age)) %>% summarise(mean(detrended_d18O))

anom <- d18O_event - d18O_base

anom_3bp <- rbind(anom_3bp, data.frame(site_name = "Klang", 
                       min_bp = as.numeric(sub_bp$bp[1]), max_bp = as.numeric(sub_bp$bp[2]),
                       anom = as.numeric(anom),
                       d18O_base = as.numeric(d18O_base),
                       d18O_event = as.numeric(d18O_event)))
anom_3bp$duration <- anom_3bp$max_bp - anom_3bp$min_bp
anom_3bp$event_centre <- anom_3bp$max_bp - (anom_3bp$duration/2)

# 4 bp's (Padre)
anom_4bp <- data.frame()
sub_bp <- bp_8_2 %>% filter(site_name == "Padre")
subdat <- dtrend_dat %>% filter(Age >= 7800 & Age <= 8400 & site_name == "Padre")

ggplot() + geom_line(data = subdat, aes(x = Age, y = detrended_d18O)) + 
  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red", lty = 2) + 
  geom_errorbarh(data = sub_bp, aes(xmin = CI2_5, xmax = CI97_5, y = c(0.9, 0.93, 0.95, 0.97)), col = "red", height = 0.02) +
  ggtitle("Padre")

sub_bp <- sub_bp %>% arrange(bp)
subdat$grp <- with(subdat, cut(Age, breaks = c(min(Age),sub_bp$bp,max(Age)), labels = 1:5))

d18O.lm <- lm(detrended_d18O ~ grp, data = na.omit(subdat))
d18O.av <- aov(d18O.lm)

tukey_sub <- TukeyHSD(x = d18O.av, conf.level = 0.99)
plot(tukey_sub , las=1 , col="brown")

# grp 4 is only sig diff from background values (grps 1 and 4)
d18O_event <- subdat %>% filter(Age >= sub_bp$bp[2] & Age <= sub_bp$bp[3]) %>% summarise(mean(detrended_d18O))
#event_CI25 <- subdat %>% filter(interp_age >= min(sub_bp$CI2_5) & interp_age <= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
d18O_base <- subdat %>% filter(Age >= min(subdat$Age) & Age <= sub_bp$bp[1] | Age >= sub_bp$bp[4] & Age <= max(subdat$Age)) %>% summarise(mean(detrended_d18O))

anom <- d18O_event - d18O_base

anom_4bp <- data.frame(site_name = "Padre", 
                       min_bp = as.numeric(sub_bp$bp[2]), max_bp = as.numeric(sub_bp$bp[3]),
                       anom = as.numeric(anom),
                       d18O_base = as.numeric(d18O_base),
                       d18O_event = as.numeric(d18O_event))
anom_4bp$duration <- anom_4bp$max_bp - anom_4bp$min_bp
anom_4bp$event_centre <- anom_4bp$max_bp - (anom_4bp$duration/2)

# 5 bp's (Wuya)
anom_5bp <- data.frame()
sub_bp <- bp_8_2 %>% filter(site_name == "Wuya")
subdat <- dtrend_dat %>% filter(Age >= 7800 & Age <= 8400 & site_name == "Wuya")

ggplot() + geom_line(data = subdat, aes(x = Age, y = detrended_d18O)) + 
  geom_vline(data = sub_bp, aes(xintercept = bp), col = "red", lty = 2) + 
  geom_errorbarh(data = sub_bp, aes(xmin = CI2_5, xmax = CI97_5, y = c(0.9, 0.93, 0.95, 0.97, 1)), col = "red", height = 0.02) +
  ggtitle("Wuya")

sub_bp <- sub_bp %>% arrange(bp)
subdat$grp <- with(subdat, cut(Age, breaks = c(min(Age),sub_bp$bp,max(Age)), labels = 1:6))

d18O.lm <- lm(detrended_d18O ~ grp, data = na.omit(subdat))
d18O.av <- aov(d18O.lm)

tukey_sub <- TukeyHSD(x = d18O.av, conf.level = 0.99)
plot(tukey_sub , las=1 , col="brown")

# grp 3+4 are sig diff from background values (grps 1 and 6)
d18O_event <- subdat %>% filter(Age >= sub_bp$bp[2] & Age <= sub_bp$bp[4]) %>% summarise(mean(detrended_d18O))
#event_CI25 <- subdat %>% filter(interp_age >= min(sub_bp$CI2_5) & interp_age <= max(sub_bp$CI2_5)) %>% summarise(mean(d18O_detrended))
d18O_base <- subdat %>% filter(Age >= min(subdat$Age) & Age <= sub_bp$bp[2] | Age >= sub_bp$bp[4] & Age <= max(subdat$Age)) %>% summarise(mean(detrended_d18O))

anom <- d18O_event - d18O_base

anom_5bp <- data.frame(site_name = "Wuya", 
                       min_bp = as.numeric(sub_bp$bp[2]), max_bp = as.numeric(sub_bp$bp[4]),
                       anom = as.numeric(anom),
                       d18O_base = as.numeric(d18O_base),
                       d18O_event = as.numeric(d18O_event))
anom_5bp$duration <- anom_5bp$max_bp - anom_5bp$min_bp
anom_5bp$event_centre <- anom_5bp$max_bp - (anom_5bp$duration/2)


anom_all <- rbind(anom_2bp, anom_3bp, anom_4bp, anom_5bp)

# add metadata
anom_all$latitude <- NA; anom_all$longitude <- NA
anom_all[anom_all$site_name == "Kaite",9:10] <- c(43.04,-3.6597)
anom_all[anom_all$site_name == "Padre",9:10] <- c(-13.2167,-44.05)
anom_all[anom_all$site_name == "Paixao",9:10] <- c(-12.39,-41.3)
anom_all[anom_all$site_name == "Tigre Perdido",9:10] <- c(-5.940556,-77.308056)
anom_all[anom_all$site_name == "Wuya",9:10] <- c(33.8167,105.4333)
anom_all[anom_all$site_name == "Klang",9:10] <- c(8.33,98.73)

ggplot(data = anom_all, aes(x = longitude, y = latitude, fill = anom)) + 
  geom_point(shape = 21, size = 3) + borders("world") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  coord_fixed(ylim = c(-50,50), xlim = c(-120,150))

# restructure data
anom_all <- data.frame(site_id = NA,
                       site_name = anom_all$site_name,
                       entity_id = NA, 
                       entity_name = NA,
                       longitude = anom_all$longitude,
                       latitude = anom_all$latitude,
                       min_bp = anom_all$min_bp,
                       max_bp = anom_all$max_bp,
                       sample_id_min_bp = NA,
                       sample_id_max_bp = NA,
                       anom = anom_all$anom,
                       d18O_base = anom_all$d18O_base,
                       d18O_event = anom_all$d18O_event,
                       duration = anom_all$duration,
                       event_centre = anom_all$event_centre)

write.csv(anom_all, "nonSISAL_82_signals.csv", row.names = F)
