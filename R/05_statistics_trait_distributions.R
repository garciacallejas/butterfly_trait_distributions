
library(tidyverse)
library(highmean)
library(colorblindr)
# -------------------------------------------------------------------------

trait.dist <- read.csv2("results/trait_distributions.csv")
trait.metrics.dist <- read.csv2("results/trait_distributions_metrics.csv")
cv.df <- read.csv2("results/trait_distributions_coefficient_of_variation.csv")
dist.data.wide <- read.csv2("results/skewness_kurtosis_relationships.csv")

# -------------------------------------------------------------------------
my.pool <- sort(unique(trait.dist$pool),decreasing = T)
my.traits <- c("HSI","SSI","STI","SVTI","TAO","WI")
my.years <- sort(unique(trait.dist$year))
my.metrics <- factor(unique(trait.metrics.dist$metric),levels = c("average","standard_deviation",
                                                                  "skewness","kurtosis","range",
                                                                  "multimodality"))

# -------------------------------------------------------------------------
# 1) test differences between distributions regional-filtered
# mann-whitney test for every pair in trait.dist

dif.dist <- trait.dist %>% group_by(year,trait) %>% 
  summarise(mann_whitney_pvalue = wilcox.test(value ~ pool)$p.value) %>%
  mutate(significant_diff = ifelse(mann_whitney_pvalue<0.05,TRUE,FALSE))
dif.dist$sig.level <- ""
dif.dist$sig.level[dif.dist$mann_whitney_pvalue < 0.05 & dif.dist$mann_whitney_pvalue >= 0.01] <- "*"
dif.dist$sig.level[dif.dist$mann_whitney_pvalue < 0.01 & dif.dist$mann_whitney_pvalue >= 0.001] <- "**"
dif.dist$sig.level[dif.dist$mann_whitney_pvalue < 0.001] <- "***"

# write.csv2(dif.dist,"results/distribution_differences.csv",row.names = F)

# just to check these tests are doing what they should,
# uncomment for a more explicit approach
# for(i.year in my.years){
#   for(i.trait in my.traits){
#     my.filtered.data <- subset(trait.dist,pool == "Filtered" & 
#                                  year == i.year &
#                                  trait == i.trait)
#     my.regional.data <- subset(trait.dist,pool == "Regional" & 
#                                  year == i.year &
#                                  trait == i.trait)
#     my.dif <- wilcox.test(my.filtered.data$value,my.regional.data$value)$p.value
#     sig.level <- "none"
#     if(my.dif < 0.05 & my.dif >= 0.01){
#       sig.level <- "*"
#     }else if(my.dif < 0.01 & my.dif >= 0.001){
#       sig.level <- "**"
#     }else if(my.dif < 0.001){
#       sig.level <- "***"
#     }
#     cat(i.year,"-",i.trait,":",round(my.dif,4)," ",sig.level,"\n")
#   }
# }

# -------------------------------------------------------------------------
# 2) test variation in temporal trends
# for each metric, pool, and trait, how many pairs of "different" years are there?
# use a z-test for that, that computes differences based on mean, sd, and n.

two_sample_z_test <- function(mean1, sd1, n1, mean2, sd2, n2) {
  # Compute standard error
  se1 <- sd1 / sqrt(n1)
  se2 <- sd2 / sqrt(n2)
  
  # Compute z-score
  z <- (mean1 - mean2) / sqrt(se1^2 + se2^2)
  
  # Compute p-value (two-tailed test)
  p_value <- 2 * (1 - pnorm(abs(z)))
  
  # Return results
  return(list(z_score = z, p_value = p_value))
}

generate_year_combinations <- function(years) {
  combinations <- as.data.frame(t(combn(years, 2)))
  colnames(combinations) <- c("Year1", "Year2")
  return(combinations)
}

year_trait_differences <- expand_grid(trait = my.traits, pool = my.pool, metric = my.metrics,
                                      generate_year_combinations(my.years),z_score = NA, z_pvalue = NA)
for(i.diff in 1:nrow(year_trait_differences)){
  my.pool <- year_trait_differences$pool[i.diff]
  my.metric <- year_trait_differences$metric[i.diff]
  my.trait <- year_trait_differences$trait[i.diff]
  my.year1 <- year_trait_differences$Year1[i.diff]
  my.year2 <- year_trait_differences$Year2[i.diff]
  
  my.mean1 <- trait.metrics.dist$value[trait.metrics.dist$pool == my.pool &
                                         trait.metrics.dist$metric == my.metric &
                                         trait.metrics.dist$year == my.year1 &
                                         trait.metrics.dist$trait == my.trait]
  my.sd1 <- trait.metrics.dist$sd[trait.metrics.dist$pool == my.pool &
                                    trait.metrics.dist$metric == my.metric &
                                    trait.metrics.dist$year == my.year1 &
                                    trait.metrics.dist$trait == my.trait]
  my.n1 <- trait.metrics.dist$num_bootstrap_replicates[trait.metrics.dist$pool == my.pool &
                                                         trait.metrics.dist$metric == my.metric &
                                                         trait.metrics.dist$year == my.year1 &
                                                         trait.metrics.dist$trait == my.trait]
  
  my.mean2 <- trait.metrics.dist$value[trait.metrics.dist$pool == my.pool &
                                         trait.metrics.dist$metric == my.metric &
                                         trait.metrics.dist$year == my.year2 &
                                         trait.metrics.dist$trait == my.trait]
  my.sd2 <- trait.metrics.dist$sd[trait.metrics.dist$pool == my.pool &
                                    trait.metrics.dist$metric == my.metric &
                                    trait.metrics.dist$year == my.year2 &
                                    trait.metrics.dist$trait == my.trait]
  my.n2 <- trait.metrics.dist$num_bootstrap_replicates[trait.metrics.dist$pool == my.pool &
                                                         trait.metrics.dist$metric == my.metric &
                                                         trait.metrics.dist$year == my.year2 &
                                                         trait.metrics.dist$trait == my.trait]
  
  my.test <- two_sample_z_test(mean1 = my.mean1,sd1 = my.sd1,n1 = my.n1,
                               mean2 = my.mean2,sd2 = my.sd2,n2 = my.n2)
  
  year_trait_differences$z_score[i.diff] <- my.test$z_score
  year_trait_differences$z_pvalue[i.diff] <- my.test$p_value
}

year_trait_differences$signif <- F
year_trait_differences$signif[abs(year_trait_differences$z_score) > 1.96 & year_trait_differences$z_pvalue < 0.05] <- T

year_trait_summary_differences <- year_trait_differences %>%
  group_by(trait,pool,metric) %>%
  summarise(num.pairs = n(),num.signif = sum(signif),
            proportion_diff = num.signif/num.pairs)
year_trait_summary_differences$trait <- factor(year_trait_summary_differences$trait, 
                                               levels = c("HSI","SSI","STI","SVTI","TAO","WI"))

# -------------------------------------------------------------------------
year_trait_diff_plot <- ggplot(year_trait_summary_differences) +
  geom_point(aes(x = proportion_diff, y = trait, fill = pool), 
             # position=position_dodge(.05),
             shape = 21, 
             size = 3
  ) +
  scale_y_discrete(limits = rev(levels(year_trait_summary_differences$trait))) + 
  facet_grid(rows = vars(metric), cols = vars(pool)) + 
  labs(y = "",x = "freq. of significantly different\npairs of years") +
  # xlim(0,1) +
  scale_fill_OkabeIto(guide = "none") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  NULL
# year_trait_diff_plot

# this is supplementary figure S3
# ggsave("results/images/temporal_metrics_variation.pdf",year_trait_diff_plot,
#        device = cairo_pdf,
#        width = 7, height = 7,dpi = 300)

# -------------------------------------------------------------------------
# 3) differences in the S2/K space. 
# for each pair of trait sets of fig 3 (same trait across years, different schemes)
# compute the statistical test on the multivariate distributions

sk.data <- dist.data.wide %>% 
  select(pool,year,trait,skewness_squared,kurtosis) %>%
  group_by(year,trait) %>% 
  pivot_wider(names_from = pool,values_from = c(skewness_squared,kurtosis)) 

adaptive.diff.data <- list()
sk.data.nested <- sk.data %>% group_by(trait) %>% nest()
for(i.trait in 1:nrow(sk.data.nested)){
  my.filtered.data <- as.matrix(sk.data.nested$data[i.trait][[1]][,c("skewness_squared_Filtered","kurtosis_Filtered")])
  my.regional.data <- as.matrix(sk.data.nested$data[i.trait][[1]][,c("skewness_squared_Regional","kurtosis_Regional")])
  adaptive.diff.data[[length(adaptive.diff.data)+1]] <- tibble(trait = sk.data.nested$trait[i.trait],
                   aSPU_pvalue = cpval_aSPU(sam1 = my.filtered.data,sam2 = my.regional.data)$pval["aSPU"])
}
adaptive.diff.df <- bind_rows(adaptive.diff.data)

# write.csv2(adaptive.diff.df,"results/sk_differences.csv",row.names = F)

