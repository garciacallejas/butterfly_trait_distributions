# calculate s-index separately for UBMS and CBMS data, 
# but using phenologies calculated with both datasets together

# NOTE: this is included for completeness: raw data is not included, 
# can be obtained after a signed agreement.

# INPUTS
# - ubms/cbms counts and visits
# - list of phenologies

# OUTPUTS
# - s-index dataframe

# -------------------------------------------------------------------------

library(tidyverse)
library(rbms)
library(data.table)
library(lubridate)

# -------------------------------------------------------------------------
all_counts<-read.delim("data/all_counts.txt", sep=";")
all_visits<-read.delim("data/all_visits.txt", sep=";")
load("data/phenology_list_CBMS_UBMS.Rdata")

# -------------------------------------------------------------------------
 
# vector with all species observed
all.sp <- unique(all_counts$SPECIES)
# small clean up
all.sp <- all.sp[which(all.sp != "")]

#initialise a time-series with day-week-month-year information.

init.year <- min(unique(all_counts$YEAR))

end.year <- max(unique(all_counts$YEAR))

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
ts_season_visit <- rbms::ts_monit_site(m_visit = all_visits,ts_season = ts_season)
my.counts <-all_counts

# get s-index for every species
s_index <- list()
for(i.sp in 1:length(all.sp)){
  
  # only observed counts (not absences) 
  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, #
                                               my.counts, 
                                               sp = all.sp[i.sp])
  
  impt_counts <- rbms::impute_count(ts_season_count=ts_season_count, 
                                    ts_flight_curve=phenology.list[[all.sp[i.sp]]], 
                                    YearLimit= NULL, 
                                    TimeUnit='w')
  
  # S-INDEX
  s_index[[i.sp]] <- rbms::site_index(butterfly_count = impt_counts, 
                                      MinFC = 0.10) # MinFC only keep sites that have been monitored at least a certain proportion of the flight curve. In this  we set the threshold to 10%,
}# for i.sp

save(s_index,file = "data/s_index_CBMS_UBMS.Rdata")

# -------------
# The following part is to clean sp that should not be in CBMS or uBMS if they've nver have been seen. 
# Then it will be saved as DB

ubms.sp <- all_counts %>% filter(SCHEME == "uBMS") %>%
  select(SPECIES) %>% unique()
cbms.sp <- all_counts %>% filter(SCHEME == "CBMS") %>%
  select(SPECIES) %>% unique() 

# convert list of s-index to DB
s.index.db<-do.call(rbind.data.frame, s_index)

# Add Scheme and separate to then filter by each sp scheme and merge again (YM: I do nto know how to do it shoerter/clenaer)
s.index.db<- data.frame(s.index.db %>% mutate(SCHEME =  case_when(
  str_detect(SITE_ID, pattern = "_uBMS") ~ "uBMS", TRUE ~ "CBMS" )))

uBMS.s.index.db <- s.index.db %>% filter(SCHEME == "uBMS") %>% filter(SPECIES %in% ubms.sp$SPECIES)
CBMS.s.index.db <- s.index.db %>% filter(SCHEME == "CBMS") %>% filter(SPECIES %in% cbms.sp$SPECIES)

s.index.db2 <- rbind(CBMS.s.index.db, uBMS.s.index.db)
write.table(s.index.db2, "data/s_index_2018_2023_CBMS_uBMS.txt", sep=";")

