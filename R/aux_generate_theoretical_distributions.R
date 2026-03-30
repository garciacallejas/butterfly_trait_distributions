
# generate probability density functions with varying moments

library(tidyverse)
library(moments)
library(patchwork)
# -------------------------------------------------------------------------

mean.value <- 3
sd.value <- 1.2
variance.value <- sd.value^2

num.observations <- 100000

# A convenient skewed distribution to work with is the gamma. 
# It has two parameters, shape and scale. 
# The mean is shape*scale and the variance is shape*scale*scale. 
# So to match the mean and variance, 
# set the scale of the gamma equal to the ratio of the variance to the mean. 
# Then once I have the scale, set the shape to the mean divided by the scale.

gamma.scale <- variance.value/mean.value
gamma.shape <- mean.value/gamma.scale

# another distribution would be a gaussian with same mean and wider sd
sd.value.wide <- 2

# yet another, a gamma with different scale and/or shape parameters
gamma.shape.skewed <- gamma.shape - 3
gamma.mean.skewed <- gamma.shape.skewed*gamma.scale

# yet another, a gamma with high kurtosis. In this case, sd needs to be larger
# to accomodate kurtosis > 5
gamma.shape.kurt <- (mean.value / 1.75)^2
gamma.scale.kurt <- 1.75 / sqrt(gamma.shape.kurt)

# -------------------------------------------------------------------------
gaussian.1 <- rnorm(num.observations,mean = mean.value,sd = sd.value)
gamma.1 <- rgamma(n = num.observations,shape = gamma.shape,scale = gamma.scale)
gaussian.2 <- rnorm(num.observations,mean = mean.value,sd = sd.value.wide)
# gamma.2.orig <- rgamma(n = num.observations,shape = gamma.shape.skewed,scale = gamma.scale)
gamma.2 <- rgamma(n = num.observations,shape = gamma.shape.kurt,scale = gamma.scale.kurt)

gaussian.1.moments <- data.frame(dist = "gaussian.1", mean = mean(gaussian.1), sd = sd(gaussian.1), skewness = moments::skewness(gaussian.1), kurtosis = moments::kurtosis(gaussian.1))
gaussian.2.moments <- data.frame(dist = "gaussian.2", mean = mean(gaussian.2), sd = sd(gaussian.2), skewness = moments::skewness(gaussian.2), kurtosis = moments::kurtosis(gaussian.2))
gamma.1.moments <- data.frame(dist = "gamma.1", mean = mean(gamma.1), sd = sd(gamma.1), skewness = moments::skewness(gamma.1), kurtosis = moments::kurtosis(gamma.1))
gamma.2.moments <- data.frame(dist = "gamma.2", mean = mean(gamma.2), sd = sd(gamma.2), skewness = moments::skewness(gamma.2), kurtosis = moments::kurtosis(gamma.2))

dist.moments <- bind_rows(gaussian.1.moments,gaussian.2.moments,gamma.1.moments, gamma.2.moments)

dist.data <- data.frame(g1 = gaussian.1, g2 = gaussian.2, gamma1 = gamma.1, gamma2 = gamma.2)
x.axis.lims <- c(-2,9)
y.axis.lims <- c(0,0.36)
# -------------------------------------------------------------------------

g1.plot <- ggplot(dist.data, aes(x = g1)) + 
  geom_density() + 
  geom_vline(aes(xintercept = dist.moments$mean[1]), linetype = "dashed", color = "darkgreen") + 
  # scale_y_continuous(breaks = NULL) +
  theme_bw() + 
  labs(x = "trait value", y = "density") + 
  # xlim(x.axis.lims) +
  scale_x_continuous(limits = x.axis.lims, breaks = c(-2,0,2,4,6,8)) +
  scale_y_continuous(limits = y.axis.lims, breaks = c(0,0.1,0.2,0.3)) +
  theme(panel.grid.minor=element_blank()) +
  # theme(axis.ticks.y = element_blank()) +
  NULL
# g1.plot

g2.plot <- ggplot(dist.data, aes(x = g2)) + 
  geom_density() + 
  geom_vline(aes(xintercept = dist.moments$mean[2]), linetype = "dashed", color = "darkgreen") + 
  # scale_y_continuous(breaks = NULL) +
  theme_bw() + 
  labs(x = "", y = "") +   
  # xlim(x.axis.lims) +
  scale_x_continuous(limits = x.axis.lims, breaks = c(-2,0,2,4,6,8)) +
  scale_y_continuous(limits = y.axis.lims, breaks = c(0,0.1,0.2,0.3)) +
  theme(panel.grid.minor=element_blank()) +
  # theme(axis.ticks.y = element_blank()) +
  NULL
# g2.plot

gamma1.plot <- ggplot(dist.data, aes(x = gamma1)) + 
  geom_density() + 
  geom_vline(aes(xintercept = dist.moments$mean[3]), linetype = "dashed", color = "darkgreen") + 
  # scale_y_continuous(breaks = NULL) +
  theme_bw() + 
  labs(x = "", y = "") + 
  scale_x_continuous(limits = x.axis.lims, breaks = c(-2,0,2,4,6,8)) +
  scale_y_continuous(limits = y.axis.lims, breaks = c(0,0.1,0.2,0.3)) +
  theme(panel.grid.minor=element_blank()) +
  # xlim(x.axis.lims) +
  # theme(axis.ticks.y = element_blank()) +
  NULL
# gamma1.plot

gamma2.plot <- ggplot(dist.data, aes(x = gamma2)) + 
  geom_density() + 
  geom_vline(aes(xintercept = dist.moments$mean[4]), linetype = "dashed", color = "darkgreen") + 
  # scale_y_continuous(breaks = NULL) +
  theme_bw() + 
  labs(x = "trait value", y = "") +  
  scale_x_continuous(limits = x.axis.lims, breaks = c(-2,0,2,4,6,8)) +
  scale_y_continuous(limits = y.axis.lims, breaks = c(0,0.1,0.2,0.3)) +
  theme(panel.grid.minor=element_blank()) +
  # xlim(x.axis.lims) +
  # theme(axis.ticks.y = element_blank()) +
  NULL
# gamma2.plot

# -------------------------------------------------------------------------

# dist.plot <- g1.plot + (g2.plot/gamma1.plot/gamma2.plot)

layout <- "
##BB
AACC
##DD
"
dist.plot <- g1.plot + g2.plot + gamma1.plot + gamma2.plot + plot_layout(design = layout) 

# dist.plot
# -------------------------------------------------------------------------

# ggsave("results/images/theoretical_distributions_raw.pdf",dist.plot,
#        device = cairo_pdf,
#        width = 8, height = 5,dpi = 300)
