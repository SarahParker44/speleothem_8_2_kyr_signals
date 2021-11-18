setwd("C:/Users/sarah/OneDrive - University of Reading/Documents/SISAL/R_programming/Data/")

library(dplyr)
library(tidyr)
library(ggplot2)

## Load 8.2ka entities 
entities_8_2 <- read.csv("speleothem_8_2_kyr_signals/spel_82_signals.csv")

# connect to SISAL database
mydb <- dbConnect(MySQL(), user = "root", password = "BevRed921", dbname = "sisal_v2", 
                  host = "localhost")

# load Holocene data
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING (site_id) JOIN sample USING (entity_id) JOIN original_chronology USING (sample_id) JOIN sisal_chronology USING (sample_id) JOIN d18O USING (sample_id);")
Raw_Data <- Raw_Data %>% filter(entity_id %in% entities_8_2$entity_id)

files_vec <- list.files(path = "C:/Users/sarah/OneDrive - University of Reading/Documents/SISAL/R_programming/Data/age_model_ensembles/")

no_SISAL_chron <- c()
age_uncert_df <- data.frame()
for (i in 1:nrow(entities_8_2)){
  ent_out <- data.frame()
  filename <- paste(entities_8_2[i,"entity_id"], entities_8_2[i,"entity_name"], sep = "-")
  sub_files <- grep(filename, files_vec, value = T)
  if (length(sub_files) == 0){ no_SISAL_chron <- c(no_SISAL_chron, i) } else {
    
    subdat <- Raw_Data %>% filter(entity_id == entities_8_2$entity_id[i])
    
    ent_sub <- data.frame()
    for (j in 1:length(sub_files)){ 
      load(paste("age_model_ensembles/", sub_files[j], sep = ""))
      dat <- get(gsub("\\..*","",sub_files[j]))
      dat$depth_sample <- dat$depth_sample*10 #mm to cm
      dat2 <- left_join(subdat[,c(1:11,37)], dat)
      
      sub_out <- data.frame()
      for (k in 1:nrow(dat2)){
        submed <- median(as.numeric(dat2[k,14:ncol(dat)]))
        subQ1 <- quantile(as.numeric(dat2[k,14:ncol(dat)]), probs = 0.2, na.rm = T)
        subQ3 <- quantile(as.numeric(dat2[k,14:ncol(dat)]), probs = 0.8, na.rm = T)
        
        age_model <- gsub("\\..*","",sub_files[j])
        age_model <- gsub(paste(unique(subdat$entity_id),unique(subdat$entity_name),sep = "-"),"",age_model)
        age_model <- gsub("-","",age_model)
        
        sub_out <- rbind(sub_out, data.frame(dat2[k,1:13], age_model = age_model,  med = submed, Q1 = subQ1, Q3 = subQ3))
      }
      ent_sub <- rbind(ent_sub, sub_out)
    }
    ent_out <- rbind(ent_out, ent_sub)  
  }
  age_uncert_df <- rbind(age_uncert_df, ent_out)
}
write.csv(age_uncert_df, "C:/Users/sarah/OneDrive/Documents/PhD/abrupt_Holocene/speleothem_8_2_kyr_signals/age_uncert_med_Qs.csv", row.names = F)
#age_uncert_df <- read.csv("speleothem_8_2_kyr_signals/age_uncert_med_Qs.csv")


orig_only <- entities_8_2 %>% filter(entity_id %in% entities_8_2[no_SISAL_chron, "entity_id"])
write.csv(orig_only, "C:/Users/sarah/OneDrive/Documents/PhD/abrupt_Holocene/speleothem_8_2_kyr_signals/82_orig_chron_only.csv", row.names = T)
#orig_only <- read.csv("82_orig_chron_only.csv")

rm(list = ls()[!ls() %in% c("entities_8_2", "age_uncert_df", "orig_only")])


# reshape
entities_8_2 <- entities_8_2[-38,]

age_uncert_82 <- data.frame()
for (i in unique(age_uncert_df$entity_id)){
  if (i == 327){ next }
  sub_ent <- entities_8_2 %>% filter(entity_id == i)
  
  x <- age_uncert_df %>% filter(entity_id == i & sample_id %in% sub_ent[,c("sample_id_min_bp","sample_id_max_bp")])
  
  #if (i == 608){ x <- unique(x)}
  
  x <- x  %>% arrange(age_model, med)
  
  x$event <- rep(c("end","start"),length(unique(x$age_model)))
  
  ## combine SISAL and original chronologies into 1 df
  x <- rbind(subset(x, select = -c(interp_age,depth_sample)),
             data.frame(unique(subset(x, select = c(sample_id:entity_name))),
                        age_model = "original",
                        med = as.numeric(as.matrix(unique(subset(x, select = interp_age)))),
                        Q1 = NA,
                        Q3 = NA,
                        event = c("end","start")))
  
  x <- x %>% gather(key = "grp", value = "age", med:Q3)
  x <- subset(x, select = -sample_id) %>% spread(key = event, value = age)
  
  x <- x %>% group_by(age_model, grp) %>% mutate(duration = start - end)
  
  x$anom <- sub_ent$anom
  
  age_uncert_82 <- rbind(age_uncert_82, x)
  
}

#age_uncert_82$duration <- age_uncert_82$start - age_uncert_82$end  
#age_uncert_82 <- age_uncert_82[-which(age_uncert_82$age_model == "original"),]
age_uncert_82 <- na.omit(age_uncert_82)




## Variation within record age-model combos
age_uncert_var <- age_uncert_82 %>% group_by(site_name, entity_id, age_model) %>%
  summarise(start_variation = max(start) - min(start),
            end_variation = max(end)-min(end),
            duration_variation = max(duration) - min(duration))
age_uncert_var <- age_uncert_var %>% gather(key = "grp", value = "variation", 4:6)

#p <- ggplot(data = age_uncert_var, aes(x = grp, y = variation)) + geom_boxplot()

# reshape
#duration <- age_uncert_var %>% filter(grp == "duration_variation") %>%
#  spread(key = age_model, value = variation) 

#duration$mu <- rowMeans(duration[,4:8], na.rm = T)
#duration <- arrange(duration, mu)

#duration <- duration[,-c(3,9)]
#duration[,3:7] <- round(duration[3:7])

#library(formattable)
#Test2 <- formatter("span",
#                  style = x ~ style(display = "block", font.weight = "bold",color = "black","border-radius" = "4px",
#                                    "padding-right" = "4px",
#                                    "background-color" = ifelse(x >= 120, "red", ifelse(x > 80 & x <= 120, "tomato", ifelse(x > 40 & x <= 80, "white", ifelse(x >20 & x <= 40, "palegreen",ifelse(x < 20, "green",NA)))))))

#tab <- formattable(duration,list(
#  Bacon = Test2,
#  Bchron = Test2, 
#  copRa = Test2,
#  linInterp = Test2,
#  linReg = Test2))

#library("htmltools")
#library("webshot")    

#export_formattable <- function(f, file, width = "100%", height = NULL, 
#                               background = "white", delay = 0.2)
#{
#  w <- as.htmlwidget(f, width = width, height = height)
#  path <- html_print(w, background = background, viewer = NULL)
#  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
#  webshot(url,
#          file = file,
#          selector = ".formattable_widget",
#          delay = delay)
#}
#export_formattable(tab, "duration_uncert.png")
#

#start_uncert <- age_uncert_var %>% filter(grp == "start_variation") %>%
#  spread(key = age_model, value = variation) 

#start_uncert$mu <- rowMeans(start_uncert[,4:8], na.rm = T)
#start_uncert <- arrange(start_uncert, mu)

#start_uncert <- start_uncert[,-c(3,9)]
#start_uncert[,3:7] <- round(start_uncert[3:7])

#Test <- formatter("span",
#                  style = x ~ style(display = "block", font.weight = "bold",color = "black","border-radius" = "4px",
#                                    "padding-right" = "4px",
#                                    "background-color" = ifelse(x >= 200, "red", ifelse(x > 150 & x <= 200, "tomato", ifelse(x > 100 & x <= 150, "white", ifelse(x >50 & x <= 100, "palegreen",ifelse(x < 50, "green",NA)))))))

#tab <- formattable(start_uncert,list(
#  Bacon = Test,
#  Bchron = Test, 
#  copRa = Test,
#  linInterp = Test,
#  linReg = Test))

#export_formattable(tab, "start_uncert.png")


#end_uncert <- age_uncert_var %>% filter(grp == "end_variation") %>%
#  spread(key = age_model, value = variation) 

#end_uncert$mu <- rowMeans(end_uncert[,4:8], na.rm = T)
#end_uncert <- arrange(end_uncert, mu)

#end_uncert <- end_uncert[,-c(3,9)]
#end_uncert[,3:7] <- round(end_uncert[3:7])

#tab <- formattable(end_uncert,list(
#  Bacon = Test,
#  Bchron = Test, 
#  copRa = Test,
#  linInterp = Test,
#  linReg = Test))

#export_formattable(tab, "end_uncert.png")


## record uncert
#record_uncert <- age_uncert_82 %>% group_by(entity_id, site_name) %>%
#  summarise(uncert_start = max(start) - min(start),
#            uncert_end = max(end) - min(end),
#            uncert_duration = max(duration) - min(duration))

#record_uncert$mu <- rowMeans(record_uncert[,3:5], na.rm = T)
#record_uncert <- arrange(record_uncert, mu)

#record_uncert <- record_uncert[,-6]
#record_uncert[,3:5] <- round(record_uncert[3:5])

#tab <- formattable(record_uncert,list(
#  uncert_start = Test,
#  uncert_end = Test, 
#  uncert_duration = Test2))

#export_formattable(tab, "overall_uncert.png")


## overall visualisation
#x <- age_uncert_var %>% spread(key = grp, value = variation)
#x$mu <- rowMeans(x[,4:6])
#x <- x %>% arrange(mu, site_name)

# "good" age uncertainty
y <- age_uncert_82 #%>% filter(entity_id %in% as.numeric(as.matrix(x[1:52,"entity_id"])))
y <- subset(y,select = -c(elevation,geology,rock_age,monitoring))

mydb <- dbConnect(MySQL(), user = "root", password = "BevRed921", dbname = "sisal_v2", 
                  host = "localhost")

# load Holocene metadata
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING (site_id);")
Raw_Data <- Raw_Data %>% filter(entity_id %in% orig_only$entity_id) %>% select(entity_id, site_name)

orig_only <- left_join(orig_only, Raw_Data)
orig_only <- na.omit(orig_only)

orig_only <- with(orig_only, data.frame(entity_id = entity_id, 
                                        site_id = site_id, 
                                        site_name = site_name, 
                                        latitude = latitude, 
                                        longitude = longitude, 
                                        entity_name = entity_name, 
                                        age_model = "original", 
                                        grp = "med", end = min_bp, start = max_bp, 
                                        duration = duration, anom = anom))
y <- rbind(y, orig_only)

y$region <- with(y, ifelse(latitude >= 20 & longitude >= -10 & longitude <= 50, "Europe", ifelse(
  latitude >= -10 & longitude >= 50, "Asia", ifelse(
    latitude <= 0 & longitude <= -50, "South America", ifelse(
      latitude >0 & longitude <= -50, "North America", "other"
    )
  )
)))

y <- arrange(y, region)
y <- subset(y, select = c(site_name, entity_name, region, grp, age_model, end, start, duration))

y$site_name2 <- with(y, paste(region, site_name, entity_name, sep = "_"))

pdf("SISAL_chron_82_timing.pdf", width = 15/2.54, height = 20/2.54)
p <- ggplot(data = y, mapping = aes(xmin = end, xmax = start, y = site_name2, col = age_model, lty = grp)) + 
  geom_errorbarh(position = position_dodge()) +
  geom_vline(xintercept = 8150, lty = 2) +
  geom_vline(xintercept = c(8090,8245), lty = 3)
dev.off()

## identify entities that do not overlap in time with Greenland record

y_out <- data.frame()
for (i in unique(y$site_name)){
  subdat <- y %>% filter(site_name == i)
  
  for (j in unique(subdat$age_model)){
    subsubdat <- subdat %>% filter(age_model == j)
    
    overlap <- c()
    if (subsubdat$start < 8090 | subsubdat$end > 8245){
      overlap = "no"
    } else { overlap = "yes"}
    subsubdat$overlap <- overlap
    
    y_out <- rbind(y_out, subsubdat)
  }
  
}


y_out <- na.omit(y_out)
#y_out <- y_out %>% group_by(site_name) %>% mutate(n_grps = n())

#y_keep <- data.frame()
#y_discard <- data.frame()
y_out2 <- data_frame()
for (i in unique(y_out$entity_name)){
  subdat <- y_out %>% filter(entity_name == i)
  
  subdatout <- data.frame()
  for (j in unique(subdat$age_model)){
    subsubdat <- subdat %>% filter(age_model == j)
    
    if (any(subsubdat$overlap == "yes")){
      subsubdat_out <- data.frame(unique(subsubdat[c(1:3,5)]),
                                  overlap = "yes")
    } else {
      subsubdat_out <- data.frame(unique(subsubdat[c(1:3,5)]),
                                  overlap = "no")
    }
    subdatout <- rbind(subdatout, subsubdat_out)
  }
  
  n_overlap = length(which(subdatout$overlap == "yes"))
  #length(which(subdat$overlap == "yes"))/nrow(subdat)
  
  #if (prop > 0.5){
  y_out2 <- rbind(y_out2, data.frame(unique(subdatout[,1:3]), n_overlap = n_overlap, n_models = nrow(subdatout)))
  #} else {
  #  y_discard <- rbind(y_discard, data.frame(subdat, prop = prop))
  #}
}


## entities with big diffs between Q1 and Q3, e.g. Tonnel'naya
x <- age_uncert_var %>% group_by(site_name, entity_id) %>% summarise(mean_start_uncert = mean(start_variation), mean_end_variation = mean(end_variation))
large_uncert <- x %>% filter(mean_start_uncert > 100)

y_keep <- y_keep %>% filter(!site_name %in% large_uncert$site_name)

ggplot(data = y_keep, mapping = aes(xmin = end, xmax = start, y = site_name2, col = age_model, lty = grp)) + 
  geom_errorbarh(position = position_dodge(), height = 0.5) +
  geom_vline(xintercept = 8150, lty = 2) +
  geom_vline(xintercept = c(8090,8245), lty = 3)

## summary of timings and duration

pdf("Asia_Europe_end.pdf", width = 15/2.54, height = 10/2.54)
ggplot(data = xx, aes(x = region, y = end)) + geom_violin() + 
  geom_hline(yintercept = 8090, lty = 2) + geom_hline(yintercept = c(median(filter(xx, region == "Asia")$end),median(filter(xx, region == "Europe")$end)), col = "red")
dev.off()

# summary data
age_uncert_82_summ <- with(y_keep, data.frame(median_start = median(start, na.rm = T),
                                              Q1_start = quantile(start, na.rm = T)[2],
                                              Q3_start = quantile(start, na.rm = T)[4],
                                              
                                              median_end = median(end, na.rm = T),
                                              Q1_end = quantile(end, na.rm = T)[2],
                                              Q3_end = quantile(end, na.rm = T)[4],
                                              
                                              median_duration = median(duration, na.rm = T),
                                              Q1_duration = quantile(duration, na.rm = T)[2],
                                              Q3_duration = quantile(duration, na.rm = T)[4]))


y_keep2 <- y_keep %>% gather(key = "start_end", value = "val", end:start)
ggplot() + geom_boxplot(data = y_keep2, aes(x = val, y = start_end, fill = start_end), alpha = 0.5) +
  geom_vline(xintercept =  8090, lty = 2, col = "#F8766D") +
  geom_vline(xintercept = 8245, lty = 2, col = "#00BFC4") +
  theme(legend.position = "none") +
  labs(y = "", x = "Years BP")

ggplot() + geom_density(data = y_keep, aes(x = start)) +
  theme_bw() +
  
  geom_vline(xintercept = age_uncert_82_summ$median_start, col = "red") +
  geom_text(aes(label = round(age_uncert_82_summ$median_start), x = age_uncert_82_summ$median_start, y = 0.0001), col = "red") +
  
  geom_vline(xintercept = as.matrix(age_uncert_82_summ[c("Q1_start","Q3_start")]), col = "red", lty = 2) +
  geom_text(aes(label = round(as.numeric(as.matrix(age_uncert_82_summ[c("Q1_start","Q3_start")]))), x = as.matrix(age_uncert_82_summ[c("Q1_start","Q3_start")]), y = 0.0028), col = "red")

ggplot() + geom_density(data = y_keep, aes(x = end)) +
  theme_bw() +
  
  geom_vline(xintercept = age_uncert_82_summ$median_end, col = "red") +
  geom_text(aes(label = round(age_uncert_82_summ$median_end), x = age_uncert_82_summ$median_end, y = 0.0001), col = "red") +
  
  geom_vline(xintercept = as.matrix(age_uncert_82_summ[c("Q1_end","Q3_end")]), col = "red", lty = 2) +
  geom_text(aes(label = round(as.numeric(as.matrix(age_uncert_82_summ[c("Q1_end","Q3_end")]))), x = as.matrix(age_uncert_82_summ[c("Q1_end","Q3_end")]), y = 0.0028), col = "red")

ggplot() + geom_density(data = y_keep, aes(x = duration)) +
  theme_bw() +
  
  geom_vline(xintercept = age_uncert_82_summ$median_duration, col = "red") +
  geom_text(aes(label = round(age_uncert_82_summ$median_duration), x = age_uncert_82_summ$median_duration, y = 0.0001), col = "red") +
  
  geom_vline(xintercept = as.matrix(age_uncert_82_summ[c("Q1_duration","Q3_duration")]), col = "red", lty = 2) +
  geom_text(aes(label = round(as.numeric(as.matrix(age_uncert_82_summ[c("Q1_duration","Q3_duration")]))), x = as.matrix(age_uncert_82_summ[c("Q1_duration","Q3_duration")]), y = 0.0028), col = "red")
