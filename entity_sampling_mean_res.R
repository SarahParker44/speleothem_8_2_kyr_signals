## function that gives mean and stdev of sampling intervals for a record

get_ent_sampling <- function(entity_id, age_start, age_end){
  
  # Load data for entity
  query <- paste("SELECT site.site_id, site.site_name, site.latitude, site.longitude, site.elevation, 
  entity.entity_name, entity.entity_status, entity.corresponding_current, entity.depth_ref,
  sample.*, gap, hiatus, d18O_measurement, d18O_precision, 
  interp_age 
  FROM site JOIN entity USING(site_id) 
  JOIN sample USING(entity_id) 
  LEFT JOIN hiatus USING(sample_id) 
  LEFT JOIN gap USING(sample_id) 
  LEFT JOIN original_chronology USING(sample_id) 
  LEFT JOIN d18O USING(sample_id) 
                 WHERE entity_id =", entity_id, "AND (interp_age BETWEEN", age_start, "AND", age_end, ");", sep = " ")
  
  dt <- dbGetQuery(mydb, query)
  
  #samples and hiatuses/gaps in chronological order
  dt <- dt[order(dt$depth_sample, dt$sample_id),]
  row.names(dt) <- NULL # reset the index
  
  #x <- which(dt$interp_age <= age_end) # select samples <= age_end
  #y <- which(is.na(dt[1:x[length(x)],"interp_age"])) # select any hiatuses/gaps that are <= age_end
  #xy <- sort(c(x,y)) # order both
  #dt <- dt[xy,] # filter to all samples and hiatuses/gaps that are <= age_end
  
  if(nrow(dt) == 1){ 
    dt_site <- data.frame(entity_id = entity_id, sampling_mean = NA, sampling_sd = NA)}
  
  # Scenario 1: entity has gaps and/or hiatuses
   else if (any((dt[['gap']] == 'G') | (dt[['hiatus']] == 'H'), na.rm = T)){
    # group the samples 
    grp_ctr <- 1
    for (k in 1:dim(dt)[1]){
      if (is.na(dt[k,'gap']) & is.na(dt[k,'hiatus'])){
        dt$grp[k] = toString(grp_ctr)
      } else {
        dt$grp[k] = NA
        grp_ctr = grp_ctr + 1
      }
    }
    res_vec_out <- c()
    for (j in 1:max(dt$grp, na.rm = T)){
      dt_sub <- dt %>% filter(grp == j)
      res_vec <- numeric()
      for (i in nrow(dt_sub):2){
        x <- dt_sub$interp_age[i] - dt_sub$interp_age[i-1]
        res_vec[i-1] <- x
      }
      res_vec_out <- c(res_vec_out, res_vec)
    }
    mean_sam <- mean(res_vec_out)
    sd_sam <- sd(res_vec_out)
    
    #df for storing output
    dt_site <- data.frame(entity_id = entity_id, sampling_mean = mean_sam, sampling_sd = sd_sam)
    
  } else {
    # Scenario 2: no hiatuses or gaps
    res_vec <- numeric()
    if (unique(dt$depth_ref) == "from base"){
      for (i in 1:(nrow(dt)-1)){
        x <- dt$interp_age[i] - dt$interp_age[i+1]
        res_vec[i] <- x
      }
    } else {
      for (i in nrow(dt):2){
        x <- dt$interp_age[i] - dt$interp_age[i-1]
        res_vec[i-1] <- x
      }
    }
    if (all(res_vec <0)){ res_vec <- res_vec*-1}
    mean_sam <- mean(res_vec)
    sd_sam <- sd(res_vec)
    
    #df for storing output
    dt_site <- data.frame(entity_id = entity_id, sampling_mean = mean_sam, sampling_sd = sd_sam)
    
  }
  
  return(dt_site)
}


