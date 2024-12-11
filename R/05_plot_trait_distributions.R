
# plot trait distributions

# INPUTS
# - trait distributions
# - trait distribution metrics

# OUTPUTS
# - several plots

# -------------------------------------------------------------------------

library(tidyverse)
library(colorblindr)
library(patchwork)
library(ggridges)
library(ggh4x)

# -------------------------------------------------------------------------

trait.dist <- read.csv2("results/trait_distributions.csv")
trait.metrics.dist <- read.csv2("results/trait_distributions_metrics.csv")
cv.df <- read.csv2("results/trait_distributions_coefficient_of_variation.csv")

# -------------------------------------------------------------------------
my.pool <- sort(unique(trait.dist$pool),decreasing = T)
my.traits <- c("HSI","SSI","STI","SVTI","TAO","WI")
my.years <- sort(unique(trait.dist$year))

# -------------------------------------------------------------------------
# which metrics do I want to plot?

my.metrics <- c("average","standard_deviation","skewness","kurtosis","range","multimodality")
trait.metrics.dist <- subset(trait.metrics.dist, metric %in% my.metrics)

trait.metrics.dist$year <- as.factor(trait.metrics.dist$year)
trait.metrics.dist$metric <- factor(trait.metrics.dist$metric, levels = my.metrics)
trait.metrics.dist$trait <- factor(trait.metrics.dist$trait, levels = my.traits)
trait.metrics.dist$pool <- factor(trait.metrics.dist$pool, levels = my.pool)

trait.dist$year <- as.factor(trait.dist$year)
trait.dist$trait <- factor(trait.dist$trait, levels = my.traits)
trait.dist$pool <- factor(trait.dist$pool, levels = my.pool)

# -------------------------------------------------------------------------
# example

my.data <- subset(subset(trait.dist, year == 2019 & trait %in% c("STI","SVTI")))
my.means <- my.data %>% group_by(pool, trait) %>%
  summarise(avg = mean(value, na.rm = T))

# I also highlight these distrib in the main data for the next plots
trait.metrics.dist$highlight <- F
trait.metrics.dist$highlight[which(trait.metrics.dist$year == 2019 &
                                     trait.metrics.dist$trait %in% c("STI","SVTI"))] <- T

example.dist.plot <- ggplot(my.data) + 
  # geom_density(aes(x = value, fill = pool), alpha = .5) +
  geom_histogram(aes(x = value, fill = pool), alpha = .9) +
  geom_vline(data = my.means, aes(xintercept = avg, 
                                  color = pool), 
             linetype = "dashed", linewidth = 1.2) +
  ggh4x::facet_grid2(.~trait, scales = "free", independent = "all") +
  labs(x = "Trait value", y = "Count") +
  scale_fill_OkabeIto() +
  scale_color_OkabeIto(guide = "none", darken = 0.1) +
  theme_bw() + 
  # theme(legend.position = "top") +
  theme(strip.background = element_blank()) +
  theme(legend.title=element_blank()) +
  guides(fill = guide_legend(position = "inside", title = NULL)) + # place legend inside plot
  theme(legend.text = element_text(size = 7.5),
        legend.key.spacing.y = unit(.2, "pt"),
        legend.box.margin = unit(0,"pt"),
        legend.justification.inside = c(.39, .994)) + # top right
  
  NULL

# this is Table 1
my.example.metrics <- subset(trait.metrics.dist, year == 2019 & trait %in% c("STI","SVTI"))

# -------------------------------------------------------------------------

pd <- 0.5
metrics.plot <- ggplot(trait.metrics.dist, aes(x = year, y = value, group = pool)) + 
  # geom_errorbar(aes(x = year,
  #                   ymin = low,
  #                   ymax = high, color = pool),
  #               width = .5,
  #               position = position_dodge(pd)) +
  geom_ribbon(aes(x = year, ymin = low, ymax = high, fill = pool),
              position = position_dodge(pd),
              alpha = .3) +
  geom_line(aes(color = pool), position = position_dodge(pd)) +
  # geom_point(aes(fill = pool), size = 2.5, shape = 21, position = position_dodge(pd)) +
  geom_point(aes(fill = pool, shape = highlight), size = 2.5, position = position_dodge(pd)) + 
  # facet_grid(trait~metric, scales = "free_y", nrow = length(my.traits), ncol = length(my.metrics)) + 
  ggh4x::facet_grid2(trait~metric, scales = "free_y", independent = "y") +
  # scale_fill_OkabeIto(guide = "none") + 
  # scale_color_OkabeIto(guide = "none") +
  scale_shape_manual(values = c(21,24), guide = "none") +
  scale_fill_OkabeIto(guide = "none") + 
  scale_color_OkabeIto() +
  labs(x="Year",y = "Value") +
  theme_bw() + 
  theme(strip.background = element_blank()) +
  theme(legend.title=element_blank()) +
  theme(legend.position="top") +
  # labs(y = "", x = "") +
  # ggtitle(label = my.traits[i.trait]) +
  NULL
# metrics.plot

# -------------------------------------------------------------------------
# trait full distributions
trait.dist$year <- as.factor(trait.dist$year)
trait.dist$pool <- factor(trait.dist$pool, levels = my.pool)
trait.dist$trait <- factor(trait.dist$trait, levels = my.traits)

trait.dist.plot.full <- ggplot(trait.dist) + 
  # geom_density(aes(x = value, fill = pool), alpha = .5) +
  geom_histogram(aes(x = value, fill = pool), alpha = .9) +
  ggh4x::facet_grid2(trait~year, scales = "free", independent = "all") +
  labs(x = "Trait value", y = "Count") +
  scale_fill_OkabeIto() +
  theme_bw() + 
  theme(legend.position = "top") +
  theme(strip.background = element_blank()) +
  theme(legend.title=element_blank()) +
  NULL

# -------------------------------------------------------------------------
# Skewness-kurtosis relationships

dist.data.wide <- subset(trait.metrics.dist, trait %in% my.traits & metric %in% c("skewness","kurtosis","range")) %>%
  pivot_wider(names_from = metric,values_from = c(value,low,high)) %>%
  mutate(skewness_squared = value_skewness^2, low_skewness_squared = low_skewness^2, high_skewness_squared = high_skewness^2) %>%
  rename(skewness = value_skewness, kurtosis = value_kurtosis, range = value_range) %>%
  select(pool,year,trait,skewness,skewness_squared,kurtosis,range,low_skewness,low_skewness_squared,low_kurtosis,low_range,
         high_skewness,high_skewness_squared,high_kurtosis,high_range)

# -------------------------------------------------------------------------
# obtain euclidean distance to boundary of each point

# the boundary is given by S^2 + 1 = K (Gross 2021)
dist.formula <- function(x,y){
  abs(-x + y - 1)/sqrt(2)
}

dist.data.wide$dist.to.boundary <- dist.formula(dist.data.wide$skewness_squared,dist.data.wide$kurtosis)
dist.data.wide <- dist.data.wide %>%
  group_by(trait) %>%
  mutate(norm_range = scales::rescale(range,to = c(0.01,0.99)))

# highlight distributions from Fig 1 here as well
dist.data.wide$highlight <- F
dist.data.wide$highlight[which(dist.data.wide$year == 2019 &
                                 dist.data.wide$trait %in% c("STI","SVTI"))] <- T
# ------------------------------------------

dist.data.plot <- ggplot(dist.data.wide, aes(x = skewness_squared, y = kurtosis, label = year)) + 
  # geom_point(aes(fill = SCHEME), shape = 21, size = 3) +
  geom_errorbar(aes(xmin = low_skewness_squared, xmax = high_skewness_squared, color = pool), width = .025) + 
  geom_errorbar(aes(ymin = low_kurtosis, ymax = high_kurtosis, color = pool), width = .025) + 
  # geom_point(aes(fill = pool, size = norm_range), shape = 21) +
  geom_point(aes(fill = pool, size = norm_range, shape = highlight)) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "grey50", linewidth = 1.5) +
  # geom_label() +
  # geom_text(hjust = 0, nudge_x = 0.001) +
  # geom_text(aes(color = SCHEME), vjust = 0, nudge_y = 0.05) +
  # scale_fill_OkabeIto(name="") +
  scale_shape_manual(values = c(21,24), guide = "none") +
  scale_fill_OkabeIto(name="", guide = "none") +
  scale_color_OkabeIto(name="", guide = "none") +
  scale_size_continuous(range = c(2,6), guide = "none") +
  labs(x = expression(S^{2}), y = "K") +
  facet_wrap(vars(trait), scales = "free",nrow = 1) +
  # ggh4x::force_panelsizes(cols = rep(2,6),total_width = unit(12,"cm")) +
  # ggh4x::force_panelsizes(cols = c(1,1,1,1,1,1)) +
  # facet_grid(cols = vars(trait)) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  NULL
# dist.data.plot

# -------------------------------------------------------------------------
dist.to.bound.plot <- ggplot(dist.data.wide, aes(x = pool, 
                                                 y = dist.to.boundary)) + 
  geom_boxplot(aes(fill = pool)) + 
  scale_fill_OkabeIto(guide = "none") +
  facet_wrap(vars(trait), scales = "free_y",nrow = 1) +
  # ggh4x::force_panelsizes(cols = rep(2,6),total_width = unit(12,"cm")) +
  theme_bw() +
  theme(strip.background = element_blank(),strip.text.x = element_blank()) +
  labs(x = "", y = "Distance to kurtosis\n boundary") +
  NULL
# dist.to.bound.plot

full.sk.plot <- dist.data.plot/dist.to.bound.plot

# -------------------------------------------------------------------------
# coefficient of variation
cv.df$metric <- factor(cv.df$metric,levels = c("average","standard_deviation","skewness","kurtosis","range","multimodality"))
cv.df$trait <- factor(cv.df$trait, levels = c("HSI","SSI","STI","SVTI","TAO","WI"))
cv.df$pool <- factor(cv.df$pool, levels = my.pool)

pd <- 0.5
cv.plot <- ggplot(cv.df) + 
  geom_point(aes(x = cv, y = trait, fill = pool), 
             position=position_dodge(.2),
             shape = 21, 
             size = 3
             ) + 
  # geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_y_discrete(limits = rev(levels(cv.df$trait))) + 
  xlim(0,1) +
  xlab("coefficient of variation") +
  scale_fill_OkabeIto() +
  facet_grid(rows = vars(metric)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.title=element_blank()) +
  guides(fill = guide_legend(position = "inside", title = NULL)) + # place legend inside plot
  theme(#legend.text = element_text(size = 7.5),
        legend.key.spacing.y = unit(.2, "pt"),
        legend.box.margin = unit(0,"pt"),
        legend.justification.inside = c(.99, .994)) + # top right
  NULL
# cv.plot

# -------------------------------------------------------------------------
ggsave("results/images/trait_metrics.pdf",metrics.plot,
       device = cairo_pdf,
       width = 15, height = 11,dpi = 300)

ggsave("results/images/trait_distributions.pdf",trait.dist.plot.full,
       device = cairo_pdf,
       width = 13, height = 9,dpi = 300)

ggsave("results/images/trait_example_distributions.pdf",example.dist.plot,
       device = cairo_pdf,
       width = 8, height = 4,dpi = 300)

ggsave("results/images/skewness_kurtosis_relationships_v3.pdf",full.sk.plot,
       device = cairo_pdf,
       width = 12, height = 7,dpi = 300)

ggsave("results/images/coeffcient_of_variation.pdf",cv.plot,
       device = cairo_pdf,
       width = 4, height = 7,dpi = 300)
