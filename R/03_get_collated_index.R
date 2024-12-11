# calculates collated index for CBMS and uBMS separately
# but using phenologies calculated with both datasets together

# NOTE: this is included for completeness: raw data is not included, 
# can be obtained after a signed agreement.

# INPUTS
# - ubms/cbms counts and visits
# - list of phenologies
# - s-index dataframe

# OUTPUTS
# - collated index list (dataframe?)

# -------------------------------------------------------------------------
# rm(list=ls())
library(tidyverse)
library(rbms)
library(data.table)
library(lubridate)

# -------------------------------------------------------------------------
all_counts<-read.delim("data/all_counts.txt", sep=";")
all_visits<-read.delim("data/all_visits.txt", sep=";")

u_visits<- data.frame(all_visits %>% filter(SCHEME == "uBMS"))
u_counts<- data.frame(all_counts %>% filter(SCHEME == "uBMS"))
u_visits_unique <- unique(u_visits[,c("DATE","SITE_ID")])

c_visits<- data.frame(all_visits %>% filter(SCHEME == "CBMS"))
c_counts<- data.frame(all_counts %>% filter(SCHEME == "CBMS"))
c_visits_unique <- unique(c_visits[,c("DATE","SITE_ID")])

# -------------------------------------------------------------------------
load("data/phenology_list_CBMS_UBMS.Rdata")
# load("s_index_CBMS_UBMS.Rdata") # do er need this?
s_index <- read.delim("data/s_index_2018_2023_CBMS_uBMS.txt",sep=";")

# -------------------------------------------------------------------------

# Create datasets
counts_datasets <- list(u_counts, c_counts)

# Give name to the datasets
names(counts_datasets) <- c("uBMS", "CBMS")
schemes <- c("uBMS", "CBMS")

visits_datasets <- list(u_visits_unique,c_visits_unique)
# names(visits_datasets) <- c("uBMS","CBMS")

c.sp <- unique(c_counts$SPECIES)
u.sp <- unique(u_counts$SPECIES)

# COLLATED-INDEX for TRANSECTS_WALKS:

collated_index <- list()
for(i.data in 1:length(counts_datasets)){
  if(i.data == 1){
    collated_index[[i.data]] <- vector(mode = "list", length = length(u.sp))
    names(collated_index[[i.data]]) <- u.sp
  }else{
    collated_index[[i.data]] <- vector(mode = "list", length = length(c.sp))
    names(collated_index[[i.data]]) <- c.sp
  }
}
names(collated_index) <- c("uBMS","CBMS")

# -------------------------------------------------------------------------
# go through both datasets separately using a loop - UBMS and CBMS

for(i.data in 1:length(schemes)){
  
  if(names(counts_datasets)[i.data] == "CBMS"){
    my.sp <- c.sp
  }else{
    my.sp <- u.sp
  }
  
  my.counts <- counts_datasets[[i.data]]
  my.visits <- visits_datasets[[i.data]]
  
  sindex_scheme <- subset(s_index, SCHEME == schemes[i.data])
  
  # -------------------------------------------------------------------------
  
  init.year <- min(unique(my.counts$YEAR))
  end.year <- max(unique(my.counts$YEAR))
  
  ts_date <- rbms::ts_dwmy_table(InitYear = init.year, 
                                 LastYear = end.year, 
                                 WeekDay1 = 'monday')
  
  ## add monitoring time 
  ts_season <- rbms::ts_monit_season(ts_date, StartMonth = 3, 
                                     EndMonth = 10, 
                                     StartDay = 1, 
                                     EndDay = NULL, 
                                     CompltSeason = TRUE, 
                                     Anchor = TRUE, 
                                     AnchorLength = 2, 
                                     AnchorLag = 2, 
                                     TimeUnit = 'w')
  
  #Add site visits to the time-series
  ts_season_visit <- rbms::ts_monit_site(ts_season, my.visits)
  
  # -------------------------------------------------------------------------
  # get collated-index for every species
  
  for(i.sp in 1:length(my.sp)){
    
    species_sindex_data <- subset(sindex_scheme, SPECIES == my.sp[i.sp])
    
    cat("scheme:",i.data,",sp:",i.sp,", observations:",nrow(species_sindex_data),"\n")
    
    if(nrow(species_sindex_data) > 0){
      
      # only observed counts (not absences) 
      ts_season_count <- rbms::ts_monit_count_site(ts_season_visit,
                                                   my.counts,
                                                   sp = my.sp[i.sp])
      
      my.coindex <- collated_index(data = species_sindex_data,
                                   s_sp = my.sp[i.sp],
                                   sindex_value = "SINDEX",
                                   glm_weights = TRUE,
                                   rm_zero = TRUE)
      
      my.coindex <- my.coindex$col_index
      my.coindex_b <- my.coindex[COL_INDEX > 0.0001 & COL_INDEX < 100000, ]
      my.coindex_logInd <- my.coindex_b[BOOTi == 0, .(M_YEAR, COL_INDEX)][, log(COL_INDEX)/log(10), by = M_YEAR][, mean_logInd := mean(V1)]
      
      # merge the mean log index
      data.table::setnames(my.coindex_logInd, "V1", "logInd")
      setkey(my.coindex_logInd, M_YEAR)
      setkey(my.coindex_b, M_YEAR)
      
      my.coindex_b <- merge(my.coindex_b, my.coindex_logInd, all.x = TRUE)
      my.coindex_b$SPECIES <- my.sp[i.sp]
      
      # store the observed collated index 
      collated_index[[i.data]][[i.sp]] <- my.coindex_b
    }else{
      collated_index[[i.data]][[i.sp]] <- NA
    }
  }# for i.sp
}# for i.data

# -------------------------------------------------------------------------

save(collated_index,file = "data/collated_index_CBMS_UBMS.Rdata")
