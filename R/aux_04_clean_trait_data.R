
# input trait values for morphospecies

# -------------------------------------------------------------------------

library(tidyverse)
# library(pracma)
# -------------------------------------------------------------------------

traits.raw <- read.csv2("data/traits_raw.csv",dec = ".") %>%
  separate_wider_delim(cols = SPECIES2,delim = " ",names = c("genus","species"),
                       too_few = "align_start",cols_remove = F)

# ad SVTI (SD of mean TEMP) to the list of traits
SVTI<-read.delim("data/SVTI.txt", sep=";")
SVTI2<-read.csv2("data/ZooKeys-367-065-s001.csv",dec = ",")
SVTI2<-SVTI2[,c(1,10)]
SVTI<-rbind(SVTI, SVTI2)
SVTI<- SVTI %>% 
  distinct(SPECIES, temp.sd, .keep_all = TRUE)
SVTI$temp.sd<-as.numeric(SVTI$temp.sd)

traits.raw2 <- merge(traits.raw, SVTI, by="SPECIES", all.x = T, all.y = F)
names(traits.raw2)[names(traits.raw2) == "temp.sd"] <- "SVTI"

traits.subset <- traits.raw2 %>% dplyr::select(SPECIES,genus,Family,sp.id,TAO,SSI,
                                              WI,HSI,FMo_Average,STI,SVTI,
                                              OS_ordinal,Max.Voltinism,Morphosp)

traits.subset  <- as_tibble(traits.subset)

traits.with.na <- traits.subset %>% filter(if_any(everything(), is.na))

# input missing trait values
traits.modified <- traits.subset

for(i.sp in 1:nrow(traits.with.na)){
  # which traits are NA
  my.na.traits <- names(traits.with.na)[which(is.na(traits.with.na[i.sp,]))]
  
  # shortcuts
  my.morphosp <- traits.with.na$Morphosp[i.sp]
  my.genus <- traits.with.na$genus[i.sp]
  my.fam <- traits.with.na$Family[i.sp]
  my.id <- traits.with.na$sp.id[i.sp]
  
  # position of this sp in the original dataframe
  my.index <- which(traits.subset$genus == my.genus & 
                      traits.subset$sp.id == my.id)
  
  # the species pool removing this one
  other.sp <- traits.subset[-my.index,]
  
  # is there an associated morphosp?
  # and if so, is it a genus or a family?
  # select the pool of rows from which to average
  if(my.morphosp != ""){
    if(my.morphosp %in% other.sp$genus){
      my.pool <- which(other.sp$genus == my.morphosp)
    }else{
      my.pool <- which(other.sp$Family == my.morphosp)
    }
  }else{
    my.sp.index <- which(other.sp$genus == my.genus & other.sp$sp.id == my.id)
    if(my.genus %in% other.sp$genus){
      my.pool <- which(other.sp$genus == my.genus)
    }else{
      my.pool <- which(other.sp$Family == my.fam)
    }
  }
  # cat(i.sp,"-",other.sp$genus[my.pool],"-",other.sp$Family[my.pool],"\n")
  
  # obtain trait averages from the species pool and update the data
  for(i.trait in 1:length(my.na.traits)){
    
    if(my.na.traits[i.trait] %in% c("TAO","SSI","WI","HSI","STI","SVTI")){
      my.avg <- colMeans(other.sp[my.pool,my.na.traits[i.trait]],na.rm=T)
    }else{
      my.avg <- pracma::Mode(dplyr::pull(other.sp,my.na.traits[i.trait])[my.pool])
    }
    traits.modified[my.index,my.na.traits[i.trait]] <- my.avg
  }# for i.trait
}# for i.sp

# -------------------------------------------------------------------------
# NOTE
# Lysandra sp. takes its traits from the other Lysandra species, 
# but these sp don't have SVTI values, so it gives NaN.
# The SVTI from Lysandra sp. should be calculated from the other Lycaenidae species,
# so it's easier to just input it here
sum(is.nan(traits.modified$SVTI))
my.lysandra <- which(is.nan(traits.modified$SVTI))
my.value <- traits.modified$SVTI[traits.modified$Morphosp == "Lycaenidae"][1]
traits.modified$SVTI[my.lysandra] <- my.value

# -------------------------------------------------------------------------

write.csv2(traits.modified,"results/trait_data_complete.csv",row.names = F)




