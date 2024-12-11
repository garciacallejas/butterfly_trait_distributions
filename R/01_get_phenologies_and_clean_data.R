# this script calculates the phenology (flight curves)
# for butterfly species accounting for UBMS and CBMS data together

# NOTE: it is included for completeness: the raw data can be shared 
# after a signed agreement.

# INPUTS
# - CBMS transects and visits
# - UBMS transects and visits

# OUTPUTS
# - list of phenologies

# -------------------------------------------------------------------------

# remotes::install_github("RetoSchmucki/rbms")
#instruccions in: https://retoschmucki.github.io/rbms/articles/Get_Started_1.html
library(rbms)
library(data.table)
library(tidyverse)

set.seed(218795)
# -------------------------------------------------------------------------
# CBMS as needed in rbms

# CBMS transects
c_transects <- read.csv2("data/medicy_transect.csv",sep=",") # this is the CBMS data
c_transects <- c_transects %>%
  rename(SITE_ID = IDitin, altitude = altura)
c_transects$altitude <- as.numeric(c_transects$altitude)
max(c_transects$altitude)
# subset < 600m
c_transects <- c_transects %>%
  filter(altitude <= 600)

# CBMS counts
c_counts <- read.csv2("data/medicy_counts.csv",sep=";") # CBMS data
c_counts <- c_counts %>%
  rename(SITE_ID = IDitin, DATE = Datam, COUNT = Nindiv, YEAR = Any)
c_counts$COUNT <-as.numeric(c_counts$COUNT)
c_counts <-c_counts %>% filter(SITE_ID %in% c_transects$SITE_ID)

# needs species names - code:
sp.list<- read.delim("data/Traits_DB_YM_21032024.txt")[,1:2] # trait data
# duplicated(sp.list$sp.id); head(sp.list)
c_counts<- merge(c_counts, sp.list, by.x = "IDesp", by.y = "sp.id", all.x = T, all.y = F)
missing<- c_counts %>% filter(is.na(c_counts$SPECIES))
# list(unique(missing$IDesp)); unique(sp.list$SPECIES)

# N species that at family or genus level
Nmorphp <- data.frame(
  one_word_count = sum(sapply(strsplit(unique(c_counts$SPECIES), " "), length) == 1),
  two_words_with_sp = sum(grepl("^\\S+ sp\\.$", unique(c_counts$SPECIES)))
)

# transform dates to columns day-month-year-year.week, works better
c_counts<- c_counts %>%
  mutate(DATE = dmy(DATE)) %>%
  mutate(DATE = format(DATE, "%Y-%m-%d"))
c_counts$DAY <- day(ymd(c_counts$DATE))
c_counts$MONTH <- month(ymd(c_counts$DATE))
c_counts$week.of.year <- week(ymd(c_counts$DATE))

#order it 
# names(c_counts)
c_counts <- c_counts[, c(3, 4, 7,8 ,2, 6, 5,9)]

# CBMS visits 
c_visits <- read.csv2("data/medicy_visit.csv",sep=";")[, 1:2]
c_visits <- c_visits %>% filter(SITE_ID %in% c_transects$SITE_ID)

# correct format date for rbms
c_visits<- c_visits %>%
  mutate(DATE = dmy(DATE)) %>%
  mutate(DATE = format(DATE, "%Y-%m-%d"))

# -------------------------------------------------------------------------
# uBMS as needed in rbms
u_counts<-read.csv2("data/output_count_table.csv", sep=",") # uBMS data
u_counts <- u_counts %>%
  rename(SITE_ID = transect_id, DATE = visit_date, COUNT = count, YEAR = year, MONTH = month, DAY = day, SPECIES= species_name)
u_counts$COUNT <-as.numeric(u_counts$COUNT)

# IMP: We will use T and A (walk) as an unique site. Hence, we need to merge them
u_counts <- u_counts %>% mutate(SITE_ID = str_sub(SITE_ID, end = -3))
u_counts<- data.frame(u_counts %>% group_by(SITE_ID, DATE, DAY, MONTH, YEAR, SPECIES,) %>% 
                        summarise(COUNT = sum(COUNT)))
u_counts$week.of.year <- week(ymd(u_counts$DATE))

# select only those in BCN
parks<-read.delim("data/CorrectedDB_2018_2023.txt", sep=";")[, 6:7] # uBMS data
parks<- parks %>% filter(!duplicated(parks), ciudad =="BARCELONA")
parks$SITE_ID <- paste("ES_uBMS_", parks$shortcode, sep="")
u_counts <-u_counts %>% filter(SITE_ID %in% parks$SITE_ID)

# uBMS visits
u_visits <- read.csv2("data/raw_ubms_ebms_visit.csv", sep=",")[, 2:3]
u_visits <- u_visits %>%
  rename(SITE_ID = transect_id, DATE= date_of_visit)

# IMP: We will use T and A (walk) as an unique site. Hence, we need to merge them
u_visits <- u_visits %>% mutate(SITE_ID = str_sub(SITE_ID, end = -3))
u_visits<- u_visits %>% filter(!duplicated(u_visits))
u_visits<- u_visits %>%
  mutate(DATE = ymd(DATE)) %>%
  mutate(DATE = format(DATE, "%Y-%m-%d"))
u_visits <-u_visits %>% filter(SITE_ID %in% parks$SITE_ID)         
# -------------------------------------------------------------------------
# join in a unique DB

# rename COUNT
c_counts$SCHEME <- rep("CBMS")
c_visits$SCHEME <- rep("CBMS")
u_counts$SCHEME <- rep("uBMS")
u_visits$SCHEME <- rep("uBMS")

# join ubms and cbms, counts and visits
all_visits <- rbind(u_visits,c_visits)
all_visits <- unique(all_visits)

all_counts <- rbind(u_counts, c_counts)
all_counts <- unique(all_counts)
sort(unique(all_counts$SPECIES))
# duplicated(unique(all_counts$SPECIES))

# -------------------------------------------------------------------------

# START PHENOLOGIES CODE

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

# -------------------------------------------------------------------------
# get phenology for every species

# store index in a list, one element per species

phenology.list <- list()

# Phenology Transects_Walks:


for(i.sp in 1:length(all.sp)){
  
  # only observed counts (not absences) 
  ts_season_count <- rbms::ts_monit_count_site(ts_season_visit, 
                                                               all_counts, 
                                                               sp = all.sp[i.sp])
  
  #yearly flight curve 
  ts_flight_curve <- rbms::flight_curve(ts_season_count, 
                                                        NbrSample = 200, 
                                                        MinVisit = 3, 
                                                        MinOccur = 2, 
                                                        MinNbrSite = 2, 
                                                        MaxTrial = 4, 
                                                        GamFamily = 'nb', 
                                                        SpeedGam = FALSE, 
                                                        CompltSeason = TRUE, 
                                                        SelectYear = NULL, 
                                                        TimeUnit = 'w')
  
  # -------------------------------------------------------------------------
  # store phenology
  
  phenology.list[[i.sp]] <- ts_flight_curve$pheno
  
}# for i.sp

names(phenology.list) <- all.sp

# -------------------------------------------------------------------------
# store in disk
save(phenology.list,file = "data/phenology_list_CBMS_UBMS.Rdata")

write.table(all_visits, "data/all_visits.txt", sep=";", row.names = F)
write.table(all_counts, "data/all_counts.txt", sep=";", row.names = F)

# the following is just the list of sp and their N per year to check which ones did not work for the pheno etc.
sp<-data.frame(all_counts %>% group_by(SPECIES, YEAR) %>% 
                 summarise(RawCOUNT = sum(COUNT)))
write.table(sp, "data/sp_years_working.txt", sep=";", row.names = F)
