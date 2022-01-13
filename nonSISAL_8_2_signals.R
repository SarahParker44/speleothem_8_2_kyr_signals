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
Anjohibe <- read.table("Anjohibe.txt", header = T); Anjohibe$site_name <- "Anjohibe"

dat <- rbind(Padre, Paixao, Tigre_perdido, Kaite[,-2],Wuya, Klang, Anjohibe)
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
  subdat <- dat %>% filter(site_name == i & Age >= 7800 & Age <= 8400)
  
  ## detrend
  subdat_lm <- lm(d18O ~ Age, data = subdat)
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
    ci_ages <- data.frame(site_name = unique(subdat$site_name), 
                          bp = subdat$Age[ci_x$confint[,2]],
                          CI2_5 = subdat$Age[ci_x$confint[,1]],
                          CI97_5 = subdat$Age[ci_x$confint[,3]]
    )
    
    bp_out <- rbind(bp_out, ci_ages)
    
    # output detrended data
    sub_dtrend <- data.frame(site_name = unique(subdat$site_name),
                             Age = subdat$Age,
                             d18O_detrended = subdat$detrended_d18O)
    
    dtrend_dat <- rbind(dtrend_dat, sub_dtrend)
    
  }
}

p <- ggplot() +
  geom_line(data = filter(dtrend_dat, site_name == "Anjohibe"), aes(x = Age, y = d18O_detrended)) +
  geom_vline(data = filter(bp_out, site_name == "Anjohibe"), aes(xintercept = bp), col = "red") +
  scale_y_reverse()
  #scale_x_continuous(breaks = seq(1000,10000,1000))


x <- bp_out %>% group_by(site_name) %>% summarise(n_bp = n())

# 2 bp's (Tigre Perdido, Kaite, Wuya)
x2 <- x %>% filter(n_bp == 2)

anom_2bp <- data.frame()
for (i in x2$site_name){
  subdat <- dtrend_dat %>% filter(site_name == i) %>% arrange(Age)
  sub_bp <- bp_out %>% filter(site_name == i) %>% arrange(bp)
  
  subdat$grp <- cut(subdat$Age, breaks = c(min(subdat$Age), sub_bp$bp, max(subdat$Age)), labels = 1:3)
  t_test <- t.test(x = subdat[which(subdat$grp == 2),"d18O_detrended"], y = subdat[which(subdat$grp %in% c(1,3)),"d18O_detrended"])
  
  d18O_event <- subdat %>% filter(Age >= min(sub_bp$bp) & Age <= max(sub_bp$bp)) %>% summarise(mean(d18O_detrended))
  d18O_base <- subdat %>% filter(Age <= min(sub_bp$bp) | Age >= max(sub_bp$bp)) %>% summarise(mean(d18O_detrended))
  d18Osd_base <- subdat %>% filter(Age <= min(sub_bp$bp) | Age >= max(sub_bp$bp)) %>% summarise(sd(d18O_detrended))
  
  anom <- d18O_event - d18O_base
  
  sub_df <- data.frame(site_name = i, 
                       min_bp = as.numeric(min(sub_bp$bp)), max_bp = as.numeric(max(sub_bp$bp)),
                       anom = as.numeric(anom),
                       d18O_base = as.numeric(d18O_base),
                       d18O_event = as.numeric(d18O_event),
                       d18Osd_base = as.numeric(d18Osd_base),
                       ttest_Pval = t_test$p.value)
  
  anom_2bp <- rbind(anom_2bp, sub_df)
}


no_signal <- data.frame()
grps <- data.frame()
# 3 bp's (Anjohibe, Klang)
x3 <- x %>% filter(n_bp == 3)

subdat <- dtrend_dat %>% filter(site_name == "Anjohibe") %>% arrange(Age)
sub_bp <- bp_out %>% filter(site_name == "Anjohibe") %>% arrange(bp)

subdat$grp <- with(subdat, cut(Age, breaks = c(min(Age),sub_bp$bp,max(Age)), labels = 1:4))

sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
sub.av <- aov(sub.lm)

sub_hsd <- TukeyHSD(sub.av)$grp

grps <- rbind(grps, data.frame(site_name = rep("Anjohibe", 2), group = c(2,3)))

##
subdat <- dtrend_dat %>% filter(site_name == "Klang") %>% arrange(Age)
sub_bp <- bp_out %>% filter(site_name == "Klang") %>% arrange(bp)

subdat$grp <- with(subdat, cut(Age, breaks = c(min(Age),sub_bp$bp,max(Age)), labels = 1:4))

sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
sub.av <- aov(sub.lm)

sub_hsd <- TukeyHSD(sub.av)$grp

grps <- rbind(grps, data.frame(site_name = rep("Klang", 2), group = 3))



# 4 bp's (Padre)
x4 <- x %>% filter(n_bp == 4)

subdat <- dtrend_dat %>% filter(site_name == "Padre") %>% arrange(Age)
sub_bp <- bp_out %>% filter(site_name == "Padre") %>% arrange(bp)

subdat$grp <- with(subdat, cut(Age, breaks = c(min(Age),sub_bp$bp,max(Age)), labels = 1:5))

sub.lm <- lm(d18O_detrended ~ grp, data = subdat)
sub.av <- aov(sub.lm)

sub_hsd <- TukeyHSD(sub.av)$grp

grps <- rbind(grps, data.frame(site_name = "Padre", group = 3))


# calculate anoms

anom_bp <- data.frame()
for (i in unique(grps$site_name)){
  sub_bp <- bp_out %>% filter(site_name == i) %>% arrange(bp)
  #sub_bp <- sub_bp[order(sub_bp$bp),]
  subdat <- dtrend_dat %>% filter(site_name == i)
  
  if (nrow(sub_bp) == 3){ lab_vec = 1:4 } else { lab_vec = 1:5 } 
  subdat$grp <- with(subdat, cut(Age, breaks = c(min(Age),sub_bp$bp,max(Age)), labels = lab_vec))
  
  event_grp <- as.numeric(grps[which(grps$site_name == i),2])
  event_grp <- na.omit(event_grp)
  
  event_d18O <- subdat %>% filter(grp %in% event_grp) %>% summarise(mean_d18O = mean(d18O_detrended))
  base_d18O <- subdat %>% filter(grp %in% c(min(lab_vec), max(lab_vec))) %>% summarise(mean_d18O = mean(d18O_detrended))
  d18Osd_base <- subdat %>% filter(grp %in% c(min(lab_vec), max(lab_vec))) %>% summarise(d18Osd_base = sd(d18O_detrended))
  #base_d18O <- subdat %>% filter(interp_age >= 7400 & interp_age <= 7900 | interp_age >= 8500 & interp_age <= 9000) %>% summarise(mean_d18O = mean(d18O_detrended))
  t_test <- t.test(x = subdat[which(subdat$grp %in% c(min(lab_vec), max(lab_vec))),"d18O_detrended"], y = subdat[which(subdat$grp %in% event_grp),"d18O_detrended"])
  
  
  x_sub <- data.frame(site_name = i,
                      min_bp = sub_bp$bp[min(event_grp,na.rm = T)-1], max_bp = sub_bp$bp[max(event_grp, na.rm = T)], 
                      anom = as.numeric(event_d18O) - as.numeric(base_d18O),
                      d18O_base = as.numeric(base_d18O), d18O_event = as.numeric(event_d18O),
                      d18Osd_base = as.numeric(d18Osd_base$d18Osd_base),
                      ttest_Pval = as.numeric(t_test$p.value))
  anom_bp <- rbind(anom_bp, x_sub)
}

anom_all <- rbind(anom_bp, anom_2bp)


# add metadata
anom_all$latitude <- NA; anom_all$longitude <- NA
anom_all[anom_all$site_name == "Anjohibe",9:10] <- c(-15.54, 46.89)
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

anom_all$duration <- anom_all$max_bp - anom_all$min_bp
anom_all$event_centre <- anom_all$max_bp - (anom_all$duration/2)

# restructure data
anom_all <- data.frame(site_id = rep(NA,6),
                       site_name = anom_all$site_name,
                       entity_id = rep(NA,6), 
                       entity_name = rep(NA,6),
                       longitude = anom_all$longitude,
                       latitude = anom_all$latitude,
                       min_bp = anom_all$min_bp,
                       max_bp = anom_all$max_bp,
                       sample_id_min_bp = rep(NA,6),
                       sample_id_max_bp = rep(NA,6),
                       anom = anom_all$anom,
                       d18O_base = anom_all$d18O_base,
                       d18O_event = anom_all$d18O_event,
                       duration = anom_all$duration,
                       event_centre = rep(NA,6))

write.csv(anom_all, "nonSISAL_82_signals.csv", row.names = F)
