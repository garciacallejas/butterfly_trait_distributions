
# correct the c-index of species with very few observations

# NOTE: this is included for completeness: raw data is not included, 
# can be obtained after a signed agreement.

# -------------------------------------------------------------------------

library(tidyverse)
library(lmtest)
library(ggeffects)

# -------------------------------------------------------------------------

# read worked DB with sp list of failed and non-failed phenologies
sp.w<-read.delim("data/sp_years_working.txt",sep = ";")
# head(sp.w)

# read counts DBs and separate per scheme
all_counts<-read.delim("data/all_counts.txt", sep=";")
all_counts.year<- data.frame(all_counts %>% group_by(YEAR, SPECIES, SCHEME) %>% 
                             summarise(COUNT = sum(COUNT)))

u_counts<- data.frame(all_counts.year %>% filter(SCHEME == "uBMS"))
c_counts<- data.frame(all_counts.year %>% filter(SCHEME == "CBMS"))

# prepare a DB for uBMS and CBMS (differ in their species * year)
sp.w.u <- merge(sp.w[,-4],u_counts, by = c("SPECIES","YEAR"), all.x = F, all.y = T)
sp.w.c <- merge(sp.w[,-4],c_counts, by = c("SPECIES","YEAR"), all.x = F, all.y = T)

# read collated DB (it will have uBMS and CBMS):
load("data/collated_index_CBMS_UBMS.Rdata")

# convert list of s-index to DB (per scheme)
c.index.ubms<-do.call(rbind.data.frame, collated_index[["uBMS"]]) %>% drop_na()
c.index.cbms<-do.call(rbind.data.frame, collated_index[["CBMS"]]) %>% drop_na()
c.index.cbms<- c.index.cbms[,c(1,3,5,6,9)] %>% rename(YEAR = M_YEAR)
c.index.ubms<- c.index.ubms[,c(1,3,5,6,9)] %>% rename(YEAR = M_YEAR)

# merge sp.w with the sheme DBs (may need to select columns from c.index DBs). FYI: Datasets do not c zero-values
# uBMS
sp.w.u<- merge(sp.w.u,c.index.ubms, by = c("SPECIES","YEAR"), all.x = T, all.y = T)
# collated index created values for some that are not seen and should be removed
# sp.w.u<- sp.w.u %>% filter(!is.na(PHENO)) #remove those that are never seen (x2 checked in u_counts). 
na.u<- sp.w.u %>% filter(is.na(COL_INDEX))
# list(unique(na.u$SPECIES))

# CBMS
sp.w.c<- merge(sp.w.c,c.index.cbms, by = c("SPECIES","YEAR"), all.x = T, all.y = T)
# sp.w.c<- sp.w.c %>% filter(!is.na(PHENO)) 
na.c<- sp.w.c %>% filter(is.na(COL_INDEX))
# list(unique(na.c$SPECIES))

# ---------------- MODEL AND PREDICT MISSING COL_INDEX ----------------

# get linear regression
##### CBMS #####
lm <- lm(COL_INDEX ~ COUNT, data = sp.w.c)
#plot(lm) # not very good. Could we only use low COUNTS?
sp.w.c.low<- sp.w.c %>% filter(COUNT < 500)
# plot(sp.w.c.low$COL_INDEX ~ sp.w.c.low$COUNT)
# abline(lm(sp.w.c.low$COL_INDEX ~ sp.w.c.low$COUNT), col ="orange")
lm <- lm(COL_INDEX ~ COUNT, data = sp.w.c.low)
#plot(lm)
# summary(lm)

# Breusch-Pagan test - test for heteroskedasticity
# bptest(lm) # p-valor < 8.579e-05

# resolve heteroskedasticity
model <- lm(COL_INDEX ~ COUNT, data = sp.w.c.low)

robustsummary <- function(model) {
  library(sandwich)
  library(lmtest)
  coeftest <- coeftest(model, vcov = vcovHC(model, type ="HC1"))
  summ <- summary(model)
  summ$coefficients[,2] <- coeftest[,2]
  summ$coefficients[,3] <- coeftest[,3]
  summ$coefficients[,4] <- coeftest[,4]
  summ
}
# summary(model)
# robustsummary(model)

# predict values for COUNT = [0 : max COUNT]
# max(sp.w.c$COUNT,na.rm = T)
my.new.db <- ggpredict(model, terms = c("COUNT[0:500]"), type = "re")
# head(my.new.db)
# put the intercept to start from zero:
my.new.db$predicted2 <-round(my.new.db$predicted - 0.1123090, 4)

#merge only with the missing data to add to the plot:
c.na<-sp.w.c.low %>% filter(is.na(COL_INDEX))
c.na <- merge(c.na, my.new.db[,c(1,7,4,5)], by.x = "COUNT", by.y = "x", all.x = T, all.y = F)

# plot simple linear model witht the predicted data
plot_model <- ggplot(sp.w.c.low, aes(x = COUNT, y = COL_INDEX)) +
  geom_point(color = "lightgrey") +
    geom_smooth(method = "lm", se = TRUE, color = "orange") +  # Model estimate and confidence interval
  geom_point(data = c.na, aes(x = COUNT, y = predicted2), color = "salmon", size=3) +  # Predicted values for missing data
  theme_minimal() +
  labs(x = "Number of observations", y = "Species Annual Abundance") +
  annotate("text", x = -Inf, y = Inf, label = "(a)", hjust = - 1, vjust = 2.2, size = 4, color= "black")+
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(colour = "darkgrey"))  # Dark grey axis lines
# plot_model


# merge predicted values for sp with NA col index to fill the CBMS dataset.
sp.w.c <- merge(sp.w.c, my.new.db[,c(1,7,4,5)], by.x = "COUNT", by.y = "x", all.x = T, all.y = F)
sp.w.c2 <- sp.w.c %>% mutate(COL_INDEX =  case_when(
  is.na(COL_INDEX) ~ predicted2, TRUE ~ COL_INDEX ))
# sp.w.c2 %>% filter(!is.na(NSITE))
# head(sp.w.c)



##### uBMS #####
lm <- lm(COL_INDEX ~ COUNT, data = sp.w.u)
#plot(lm) # not very good. Could we only use low COUNTS?
sp.w.u.low<- sp.w.u %>% filter(COUNT < 500)
# plot(sp.w.u.low$COL_INDEX ~ sp.w.u.low$COUNT)
# abline(lm(sp.w.u.low$COL_INDEX ~ sp.w.u.low$COUNT), col ="orange")
lm <- lm(COL_INDEX ~ COUNT, data = sp.w.u.low)
#plot(lm)
# summary(lm)

# Breusch-Pagan test - test for heteroskedasticity
# bptest(lm) # p-valor < 8.579e-05

# resolve heteroskedasticity
model <- lm(COL_INDEX ~ COUNT, data = sp.w.u.low)

robustsummary <- function(model) {
  library(sandwich)
  library(lmtest)
  coeftest <- coeftest(model, vcov = vcovHC(model, type ="HC1"))
  summ <- summary(model)
  summ$coefficients[,2] <- coeftest[,2]
  summ$coefficients[,3] <- coeftest[,3]
  summ$coefficients[,4] <- coeftest[,4]
  summ
}
# summary(model)
# robustsummary(model)

# predict values for COUNT = [0 : max COUNT]
my.new.db <- ggpredict(model, terms = c("COUNT[0:500]"), type = "re")
# head(my.new.db)
# put the intercept to start from zero:
my.new.db$predicted2 <-round(my.new.db$predicted - 0.1123090, 4)

#merge only with the missing data to add to the plot:
u.na<-sp.w.u.low %>% filter(is.na(COL_INDEX))
u.na <- merge(u.na, my.new.db[,c(1,7,4,5)], by.x = "COUNT", by.y = "x", all.x = T, all.y = F)

# plot simple linear model witht the predicted data
plot_model.u <- ggplot(sp.w.u.low, aes(x = COUNT, y = COL_INDEX)) +
  geom_point(color = "lightgrey") +
  geom_smooth(method = "lm", se = TRUE, color = "lightblue") +  # Model estimate and confidence interval
  geom_point(data = u.na, aes(x = COUNT, y = predicted2), color = "blue", size=3) +  # Predicted values for missing data
  theme_minimal() + xlim(0,500)+
  labs(x = "Number of observations", y = "") +
  annotate("text", x = -Inf, y = Inf, label = "(b)", hjust = - 1, vjust = 2.2, size = 4)+
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.line = element_line(colour = "darkgrey"))  # Dark grey axis lines
# plot_model.u

# merge predicted values for sp with NA col index
sp.w.u %>% filter(is.na(COL_INDEX))
sp.w.u <- merge(sp.w.u, my.new.db[,c(1,7,4,5)], by.x = "COUNT", by.y = "x", all.x = T, all.y = F)
sp.w.u <- sp.w.u %>% mutate(COL_INDEX =  case_when(
  is.na(COL_INDEX) ~ predicted2, TRUE ~ COL_INDEX ))


##### plot both #####
# library(gridExtra)
# grid.arrange(plot_model, plot_model.u, ncol=2)

##### merge both DB  #####
collated_index_corrected <- rbind(sp.w.c, sp.w.u)

# remove SPECIES == Cydalima perspectalis == POLILLA DEL BOIX
collated_index_corrected <- collated_index_corrected %>% filter(!(SPECIES %in% c("Cydalima perspectalis","Satyridae")))

write.table(collated_index_corrected, "data/collated_index_2018-2023_CBMS_uBMS_corrected.txt", sep=";", row.names = F)
