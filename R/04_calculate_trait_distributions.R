
# obtain trait distributions per zone, using the collated index

# INPUTS
# - corrected collated index
# - trait data

# OUTPUTS
# - distribution of traits in ubms/cbms every year
# - trait distribution metrics 

# -------------------------------------------------------------------------

library(tidyverse)
library(traitstrap)
library(moments)

# -------------------------------------------------------------------------

sp.traits <- read.csv2("results/trait_data_complete.csv")
collated.index <- read.csv2("data/collated_index_2018-2023_CBMS_uBMS_corrected.txt",dec = ".")

# rename schemes
collated.index$SCHEME[collated.index$SCHEME == "CBMS"] <- "Regional"
collated.index$SCHEME[collated.index$SCHEME == "uBMS"] <- "Filtered"

my.traits <- c("HSI","SSI","STI","SVTI","TAO","WI")
my.years <- sort(unique(collated.index$YEAR))
my.pool <- unique(collated.index$SCHEME)

# -------------------------------------------------------------------------
names(sp.traits)[1] <- "species"

c.index.clean <- collated.index %>%
  dplyr::select(SCHEME,YEAR,SPECIES,COL_INDEX)
names(c.index.clean) <- c("pool","year","species","c_index")

c.index.clean.2 <- c.index.clean %>% drop_na(c_index,pool)  

traits.clean <- sp.traits %>%
  dplyr::select(species,TAO,SSI,WI,HSI,FMo_Average,STI,SVTI,OS_ordinal,Max.Voltinism) %>%
  pivot_longer(cols = TAO:Max.Voltinism, names_to = "trait", values_to = "value") %>%
  filter(trait %in% my.traits)
traits.clean$trait <- factor(traits.clean$trait,levels = my.traits)
names(traits.clean) <- c("species","trait","value")

treatments <- data.frame(pool = sort(unique(c.index.clean$pool)),
                         year = sort(unique(c.index.clean$year)))

# apparently the trait_fill function needs that the "treatments", that in our case
# are the regional/filtered communities and the year, are also columns in the trait data
# since we have no variation in trait values across pool or year, simply repeat
# the trait dataframe for each combination
traits.repeated <- traits.clean %>%
  expand_grid(pool = sort(unique(c.index.clean$pool)),
              year = sort(unique(c.index.clean$year)))

# -------------------------------------------------------------------------
# compute trait distributions and metrics with associated confidence intervals
nrep <- 1000

trait_filling <- trait_fill(comm = c.index.clean.2,
                            traits = traits.repeated,
                            taxon_col = "species",
                            trait_col = "trait",
                            value_col = "value",
                            abundance_col = "c_index",
                            scale_hierarchy = c("pool","year"))

# non-parametric bootstrap
trait_moments <- trait_np_bootstrap(trait_filling,nrep = nrep)

# range is done manually, and it has no confidence interval associated
trait.dist.range <- trait_filling %>%
  group_by(pool,year,trait) %>%
  summarise(min = min(value),max = max(value)) %>%
  mutate(range = abs(max-min)) %>%
  pivot_longer(cols = range,names_to = "metric",values_to = "value") %>%
  mutate(low = NA, high = NA) %>%
  select(pool,year,trait,metric,value,low,high)

# and the original distribution
raw_moments <- trait_np_bootstrap(filled_traits = trait_filling,raw = T)

# get the statistical moments of the bootstrap data
summarised_trait_dists <- trait_summarise_boot_moments(trait_moments,parametric = F)

# get the full kurtosis, not the excess one
# and the standard dev instead of variance
summarised_trait_dists$standard_deviation <- sqrt(summarised_trait_dists$var)
summarised_trait_dists$ci_low_standard_deviation <- sqrt(summarised_trait_dists$ci_low_var)
summarised_trait_dists$ci_high_standard_deviation <- sqrt(summarised_trait_dists$ci_high_var)

summarised_trait_dists$kurt <- summarised_trait_dists$kurt + 3
summarised_trait_dists$ci_low_kurt <- summarised_trait_dists$ci_low_kurt + 3
summarised_trait_dists$ci_high_kurt <- summarised_trait_dists$ci_high_kurt + 3

# and transform to long format for plotting
summarised_trait_dists_long <- summarised_trait_dists %>%
  pivot_longer(
    cols = c("mean","standard_deviation","skew","kurt"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    low = case_when(
      metric == "mean" ~ ci_low_mean,
      metric == "standard_deviation" ~ ci_low_standard_deviation,
      metric == "skew" ~ ci_low_skew,
      metric == "kurt" ~ ci_low_kurt
    ),
    high = case_when(
      metric == "mean" ~ ci_high_mean,
      metric == "standard_deviation" ~ ci_high_standard_deviation,
      metric == "skew" ~ ci_high_skew,
      metric == "kurt" ~ ci_high_kurt
    )
  ) %>%
  select(pool, year, trait, metric, value, low, high)

summarised_trait_dists_long$metric[summarised_trait_dists_long$metric == "mean"] <- "average"
summarised_trait_dists_long$metric[summarised_trait_dists_long$metric == "skew"] <- "skewness"
summarised_trait_dists_long$metric[summarised_trait_dists_long$metric == "kurt"] <- "kurtosis"
summarised_trait_dists_long$global <- NULL

# -------------------------------------------------------------------------
# NOTE:
# above I obtain the main metrics with the traitstrap package, to get 
# bootstrap confidence intervals as well
# however, I am not sure how to gather the raw distribution 
# so I keep the "manual" way because i need that raw distribution 
# to calculate the DIP multimodality metric
# this serves also as a sanity check, because the figure of trait metrics across years
# is equivalent with the traitstrap method and with the manually obtained distributions

cindex.traits <- c.index.clean.2 %>% left_join(sp.traits) %>%
  dplyr::select(pool,species,year,c_index,TAO,SSI,WI,HSI,FMo_Average,STI,SVTI,OS_ordinal,Max.Voltinism) %>%
  pivot_longer(cols = TAO:Max.Voltinism, names_to = "trait", values_to = "value") %>%
  filter(trait %in% my.traits)
cindex.traits$trait <- factor(cindex.traits$trait,levels = my.traits)

# -------------------------------------------------------------------------
# i need a vector of trait values for every scheme, year, site, and trait
schemes <- sort(unique(cindex.traits$pool))
years <- sort(unique(cindex.traits$year))
traits <- my.traits

# i.scheme <- i.year <- i.trait <- i.sp <-   1

trait.list <- list()
for(i.scheme in 1:length(schemes)){
  for(i.year in 1:length(years)){
    # for(i.site in 1:length(sites)){
      for(i.trait in 1:length(traits)){
        
        trait.data <- subset(cindex.traits,pool == schemes[i.scheme] & 
                               year == years[i.year] &
                               trait == traits[i.trait])
        # only observations with positive c_index values
        trait.data <- subset(trait.data, c_index > 0)
        
        my.trait.dist <- NA
        
        if(nrow(trait.data)>0){
          my.sp <- sort(unique(trait.data$species))
          for(i.sp in 1:length(my.sp)){
            my.trait.value <- trait.data$value[trait.data$species == my.sp[i.sp]]
            my.cindex.value <- trait.data$c_index[trait.data$species == my.sp[i.sp]]
            if(round(my.cindex.value) == 0){
              my.cindex.value <- 1
            }else{
              my.cindex.value <- round(my.cindex.value)
            }
            
            my.trait.sp.dist <- rep(my.trait.value,my.cindex.value)
            
            my.trait.dist <- c(my.trait.dist,my.trait.sp.dist)
            
          }# for i.sp
          my.trait.dist <- my.trait.dist[which(!is.na(my.trait.dist))]
          
          trait.list[[length(trait.list)+1]] <- data.frame(pool = schemes[i.scheme],
                                                                     year = years[i.year],
                                                                     trait = traits[i.trait],
                                                                     value = my.trait.dist)
          
        }# if data
      }# for i.trait
    # }# for i.site
  }# for i.year
}# for i.scheme

trait.dist.df <- bind_rows(trait.list)

# -------------------------------------------------------------------------
# UPDATE
# I also want the confidence interval of the multimodality metric,
# i.e. apply it to bootstrapped distributions
# so let's redo manually the whole workflow of the traitstrap package... sigh...
# at least this ensures that the package works as it should

# trait.dists.manual.bootstrap.list <- list()
trait.dists.bootstrap.metric.list <- list()

# use nrep from above
sample_size <- 200 # as the default in the trait_nt_bootstrap function
for(i.pool in 1:length(my.pool)){
  for(i.year in 1:length(my.years)){
    for(i.trait in 1:length(my.traits)){
      my.subset.data <- subset(trait_filling,pool == my.pool[i.pool] &
                                 year == my.years[i.year] & 
                                 trait == my.traits[i.trait])
      for(i.rep in 1:nrep){
        # this is copied verbatim from the trait_np_bootstrap function,
        # so should be equivalent
        my.dist <- slice_sample(my.subset.data, n = sample_size, 
                                replace = TRUE, weight_by = weight)
        my.dist$rep <- i.rep
        # trait.dists.manual.bootstrap.list[[length(trait.dists.manual.bootstrap.list)+1]] <- my.dist
  
        # obtain all metrics
        my.dist.range <- c(min(my.dist$value,na.rm = T),max(my.dist$value,na.rm = T))
        my.dist.dip <- diptest::dip.test(my.dist$value)
        
        trait.dists.bootstrap.metric.list[[length(trait.dists.bootstrap.metric.list)+1]] <- 
          data.frame(pool = my.pool[i.pool],
                     year = my.years[i.year],
                     trait = my.traits[i.trait],
                     rep = i.rep,
                     average = mean(my.dist$value),
                     standard_deviation = sd(my.dist$value),
                     skewness = moments::skewness(my.dist$value),
                     kurtosis = moments::kurtosis(my.dist$value),
                     range = abs(my.dist.range[2]-my.dist.range[1]),
                     multimodality = my.dist.dip$statistic)
        #cat(i.rep,"\n")
      }# for i.rep
    }# for i.trait
  }# for i.yer
}# for i.pool

bootstrap.metrics.manual <- bind_rows(trait.dists.bootstrap.metric.list)
# bootstrap.distributions <- bind_rows(trait.dists.manual.bootstrap.list)

# -------------------------------------------------------------------------
# replicate the workflow of traitstrap to obtain CI
# copy-paste their function
sd_mult <- 1
ci <- 0.95
parametric <- F

get_ci <- function(data, sd_mult = 1, ci = 0.95, which, parametric = TRUE) {
  if (isTRUE(parametric)) {
    if (which == "high") {
      return(mean(data) + sd(data) * sd_mult)
    }
    if (which == "low") {
      return(mean(data) - sd(data) * sd_mult)
    }
  } else {
    if (which == "high") {
      return(quantile(data, probs = (1 + ci) / 2, type = 1))
    }
    if (which == "low") {
      return(quantile(data, probs = (1 - ci) / 2, type = 1))
    }
  }
}

# copy-paste from their function 
# summarise_bootstrap_moments
# with some adjustments (i.e. sd instead of var)
bootstrap.metrics.ci <- bootstrap.metrics.manual %>% 
  group_by(pool,year,trait) %>%
  summarise(
    mean = mean(.data$average),
    ci_low_mean = get_ci(
      data = .data$average, sd_mult = sd_mult, ci = ci,
      which = "low", parametric = parametric
    ),
    ci_high_mean = get_ci(
      data = .data$average, sd_mult = sd_mult, ci = ci,
      which = "high", parametric = parametric
    ),
    sd = mean(.data$standard_deviation),
    ci_low_sd = get_ci(
      data = .data$standard_deviation, sd_mult = sd_mult, ci = ci,
      which = "low", parametric = parametric
    ),
    ci_high_sd = get_ci(
      data = .data$standard_deviation, sd_mult = sd_mult, ci = ci,
      which = "high", parametric = parametric
    ),
    skew = mean(.data$skewness),
    ci_low_skew = get_ci(
      data = .data$skewness, sd_mult = sd_mult, ci = ci,
      which = "low", parametric = parametric
    ),
    ci_high_skew = get_ci(
      data = .data$skewness, sd_mult = sd_mult, ci = ci,
      which = "high", parametric = parametric
    ),
    kurt = mean(.data$kurtosis),
    ci_low_kurt = get_ci(
      data = .data$kurtosis, sd_mult = sd_mult, ci = ci,
      which = "low", parametric = parametric
    ),
    ci_high_kurt = get_ci(
      data = .data$kurtosis, sd_mult = sd_mult, ci = ci,
      which = "high", parametric = parametric
    ),
    rg = mean(.data$range),
    ci_low_rg = get_ci(
      data = .data$range, sd_mult = sd_mult, ci = ci,
      which = "low", parametric = parametric
    ),
    ci_high_rg = get_ci(
      data = .data$range, sd_mult = sd_mult, ci = ci,
      which = "high", parametric = parametric
    ),
    md = mean(.data$multimodality),
    ci_low_md = get_ci(
      data = .data$multimodality, sd_mult = sd_mult, ci = ci,
      which = "low", parametric = parametric
    ),
    ci_high_md = get_ci(
      data = .data$multimodality, sd_mult = sd_mult, ci = ci,
      which = "high", parametric = parametric
    )
    # cvar = mean(.data$cv),
    # ci_low_cvar = get_ci(
    #   data = .data$cv, sd_mult = sd_mult, ci = ci,
    #   which = "low", parametric = parametric
    # ),
    # ci_high_cvar = get_ci(
    #   data = .data$cv, sd_mult = sd_mult, ci = ci,
    #   which = "high", parametric = parametric
    # )
  )# summarise

# put everything together
# and transform to long format for plotting
summarised_trait_dists_long_2 <- bootstrap.metrics.ci %>%
  pivot_longer(
    cols = c("mean","sd","skew","kurt","rg","md"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    low = case_when(
      metric == "mean" ~ ci_low_mean,
      metric == "sd" ~ ci_low_sd,
      metric == "skew" ~ ci_low_skew,
      metric == "kurt" ~ ci_low_kurt,
      metric == "rg" ~ ci_low_rg,
      metric == "md" ~ ci_low_md,
      # metric == "cvar" ~ ci_low_cvar
    ),
    high = case_when(
      metric == "mean" ~ ci_high_mean,
      metric == "sd" ~ ci_high_sd,
      metric == "skew" ~ ci_high_skew,
      metric == "kurt" ~ ci_high_kurt,
      metric == "rg" ~ ci_high_rg,
      metric == "md" ~ ci_high_md,
      # metric == "cvar" ~ ci_high_cvar
    )
  ) %>%
  select(pool, year, trait, metric, value, low, high)

summarised_trait_dists_long_2$metric[summarised_trait_dists_long_2$metric == "mean"] <- "average"
summarised_trait_dists_long_2$metric[summarised_trait_dists_long_2$metric == "sd"] <- "standard_deviation"
summarised_trait_dists_long_2$metric[summarised_trait_dists_long_2$metric == "skew"] <- "skewness"
summarised_trait_dists_long_2$metric[summarised_trait_dists_long_2$metric == "kurt"] <- "kurtosis"
summarised_trait_dists_long_2$metric[summarised_trait_dists_long_2$metric == "rg"] <- "range"
summarised_trait_dists_long_2$metric[summarised_trait_dists_long_2$metric == "md"] <- "multimodality"
# summarised_trait_dists_long_2$metric[summarised_trait_dists_long_2$metric == "cvar"] <- "coef_variation"

# coefficient of variation
cv.df <- summarised_trait_dists_long_2 %>% filter(metric %in% c("average","standard_deviation","skewness","kurtosis","range","multimodality")) %>%
  group_by(pool,trait,metric) %>%
  # summarise(avg = mean(value), sdev = sd(value), cv = sd(value)/abs(mean(value)))
  summarise(avg = mean(value), sdev = sd(value), cv = sd(value)/mean(abs(value)))

# avg.sd.plot <- ggplot(cv.df, aes(x = abs(avg), y = sdev, color = pool)) + 
#   geom_point() + 
#   geom_abline(slope = 1, intercept = 0) +
#   facet_wrap(vars(trait,metric)) + 
#   theme_bw() +
#   NULL
# svti.sk.data <- subset(summarised_trait_dists_long_2, trait == "SVTI" & metric == "skewness")
# svti.sk.data %>% group_by(pool) %>% summarise(avg = mean(value), sd = sd(value))
# -------------------------------------------------------------------------

write.csv2(trait.dist.df,"results/trait_distributions.csv",row.names = F)
write.csv2(summarised_trait_dists_long_2, "results/trait_distributions_metrics.csv",row.names = F)
write.csv2(cv.df,"results/trait_distributions_coefficient_of_variation.csv",row.names = F)

# -------------------------------------------------------------------------
# OLD CODE and tests/checks to be deleted eventually
# -------------------------------------------------------------------------

# 1 - gather the moments and ci from the traitstrap functions, I can check that they 
# are basically equal to the ones I obtained, as they should
# but in this table the range and multimod have no CI

# summarised_trait_dists_long <- bind_rows(list(summarised_trait_dists_long,trait.dist.range,trait.dist.multimod.long))

# 2 - old-school way to get multimodality, it's equivalent - I did not fully trust the very low p-values in every dist
# test.multimod.list <- list()
# for(i.scheme in 1:length(schemes)){
#   for(i.year in 1:length(years)){
#     for(i.trait in 1:length(traits)){
#       my.data <- subset(trait.dist.df, SCHEME == schemes[i.scheme] &
#                           YEAR == years[i.year] & 
#                           trait == traits[i.trait])
#       my.multimod <- diptest::dip.test(my.data$value)
#       test.multimod.list[[length(test.multimod.list)+1]] <- data.frame(pool = schemes[i.scheme],
#                                                                        year = years[i.year],
#                                                                        trait = traits[i.trait],
#                                                                        dip_statistic = my.multimod$statistic,
#                                                                        dip_pvalue = my.multimod$p.value)
#     }# for i.trait
#   }# for i.year
# }# for i.scheme
# test.multimod.df <- bind_rows(test.multimod.list)

# compute metrics for all trait distributions
# if I afterwards plot this, the figure is equivalent to the metrics
# obtained with traitstrap

# trait.dist.metrics <- trait.dist.df %>%
#   group_by(SCHEME,YEAR,trait) %>%
#   summarise(resu = list(distribution_descriptors(value))) %>%
#   unnest_wider(resu) %>%
#   dplyr::select(-distr_multimodality_pvalue) %>%
#   rename(c(obs = num_obs, 
#            min_value = distr_min, 
#            max_value = distr_max,
#            range = distr_range,
#            average = distr_average,
#            median = distr_median,
#            st_dev = distr_sd,
#            skewness = distr_skewness,
#            kurtosis = distr_kurtosis,
#            distance_to_boundary = distr_dist_to_sk_boundary,
#            multimodality = distr_multimodality_statistic)) %>%
#   pivot_longer(cols = obs:multimodality,names_to = "metric",values_to = "value")
# 





