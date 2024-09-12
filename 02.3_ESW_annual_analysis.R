# libraries --------------------------------------------------------------------

library(tidyverse)
library(broom)
library(lubridate)
library(patchwork)
#devtools::install_github('cttobin/ggthemr')
library(ggthemr)
ggthemr("fresh") 
library(ggExtra)
library(corrplot)
library(DescTools)
library(ggrepel)

source('00_functions.R', echo=TRUE)
#source('01_processing.R', echo=TRUE)

# species summary ----------------------------------------------------------------

#speciesSummary <- data %>%
#                    group_by(Species, Year) %>%
#                    summarise(nuObs = length(Species)) %>%
#                    ungroup()


# get ESW models for each year  --------------------------------------------------

#I used the data for each year to separately estimate the ESWs for each year

nullDF2014 <- readRDS("Data/distance_passerines_constant_2014.rds") %>% add_column(Year=2014)
nullDF2015 <- readRDS("Data/distance_passerines_constant_2015.rds") %>% add_column(Year=2015)
nullDF2016 <- readRDS("Data/distance_passerines_constant_2016.rds") %>% add_column(Year=2016)
nullDF2017 <- readRDS("Data/distance_passerines_constant_2017.rds") %>% add_column(Year=2017)

# combine all -------------------------------------------------------------------

nullDFannual <- bind_rows(nullDF2014,
                          nullDF2015,
                          nullDF2016,
                          nullDF2017) %>%
                dplyr::select(ESW,Species,Year) %>%
                pivot_wider(names_from = Year,
                            values_from = ESW) %>%
                janitor::clean_names() %>%
                filter(complete.cases(.))

# correlation ------------------------------------------------------------------------

#how correlated are species ESWs across years?
pairs(nullDFannual[,-1])
#pretty good!!
#tidy up this plot for the SI
#remove the X from the year labels

cor(nullDFannual[,-1])

#the concordance correlation coefficient is more apt here
CCC(nullDFannual$x2014, nullDFannual$x2015)$rho.c
#est    lwr.ci   upr.ci
#1 0.8016382 0.7033917 0.869829
CCC(nullDFannual$x2015, nullDFannual$x2016)$rho.c
#est    lwr.ci    upr.ci
#1 0.8509094 0.7727503 0.9036479
CCC(nullDFannual$x2016, nullDFannual$x2017)$rho.c
#est    lwr.ci    upr.ci
#1 0.843113 0.7631393 0.8976518

# overlay --------------------------------------------------------------------------

#swatch()
panela <- ggplot(nullDFannual) +
  geom_point(aes(x = x2014, y = x2015, colour = "2014 vs 2015")) +
  geom_point(aes(x = x2015, y = x2016, colour = "2015 vs 2016")) +
  geom_point(aes(x = x2016, y = x2017, colour = "2016 vs 2017")) +
  xlab("Previous Effective Strip Width (m)") +
  ylab("Effective Strip Width (m)") +
  geom_abline(linetype = "dashed") +
  scale_colour_manual(
    values = c(
      "2014 vs 2015" = swatch()[1],
      "2015 vs 2016" = swatch()[5],
      "2016 vs 2017" = swatch()[9]
    )
  ) +
  labs(colour = "Year Comparison") + 
  theme(
    legend.position = c(0.1, 0.9))  
panela

#tidy and save

# stability of trait associations --------------------------------------------------

#get the traits data frame from the other script

traitsDF <- nullDFannual %>%
  rename(Species = species) %>%
  inner_join(., traits, by="Species") 

lm2014 <- lm(x2014 ~ log(Mass) + ForStrat.ground + Habitat + 
            Trophic.Level, data=traitsDF)
summary(lm2014)

lm2015 <- lm(x2015 ~ log(Mass) + ForStrat.ground + Habitat + 
             Trophic.Level, data=traitsDF)
summary(lm2015)

lm2016 <- lm(x2016 ~ log(Mass) + ForStrat.ground + Habitat + 
             Trophic.Level, data=traitsDF)
summary(lm2016)

lm2017 <- lm(x2017 ~ log(Mass) + ForStrat.ground + Habitat + 
             Trophic.Level, data=traitsDF)
summary(lm2017)

# can you arrange all these outputs into a single table 
# just pull the estimate and confidence intervals of the estimates (get these from the confint() instead of summary()) from each year
# the table could have rows for each coefficient (ignore the intercept)
# and columns for each year
library(modelsummary)
modelsummary(lm2014, statistic = c("conf.low", "conf.high"), 
             shape = term ~ model + statistic,
             output="Figures/Linear_Model_Ouput_2014.docx")

modelsummary(lm2015, statistic = c("conf.low", "conf.high"), 
             shape = term ~ model + statistic,
             output="Figures/Linear_Model_Ouput_2015.docx")

modelsummary(lm2016, statistic = c("conf.low", "conf.high"), 
             shape = term ~ model + statistic,
             output="Figures/Linear_Model_Ouput_2016.docx")

modelsummary(lm2017, statistic = c("conf.low", "conf.high"), 
             shape = term ~ model + statistic,
             output="Figures/Linear_Model_Ouput_2017.docx")

# prediction -----------------------------------------------------------------------

# lets start playing with scenarios of not having the ESW data for a given year

# if we want to predict the ESW in a year (assuming we dont have any direct data for that year), 
# do we use:
# a previous estimate of ESW from the previous year - 'previous estimate'
# the ESW predicted by traits - 'trait-based estimate'
# which is better? lets ask this question:

#lets imagine we are in 2015, but we only have data for 2014


## 2014 to 2015 ------------------------------------------------------------

#we'll use the 2014 data to make trait predictions 
lm1 <- lm(x2014 ~ log(Mass) + ForStrat.ground + Habitat + 
            Trophic.Level, data=traitsDF)
traitsDF$preds2014 <- predict(lm1)

#we'll then compare it to what was actually observed in 2015
summary(traitsDF$x2015 - traitsDF$x2014) #compare to past estimate
summary(traitsDF$x2015 - traitsDF$preds2014) #compare to trait-based estimate
CCC(traitsDF$x2015, traitsDF$x2014)$rho.c
CCC(traitsDF$x2015, traitsDF$preds2014)$rho.c
#looks like the past estimate does a bit better than the predicted estimate based on the traits

#plot individual errors
errors1415 <- data.frame(error = c(CCC(traitsDF$x2015, traitsDF$x2014)$blalt$delta,
                               CCC(traitsDF$x2015, traitsDF$preds2014)$blalt$delta),
                     type = c(rep("past estimate",72), rep("trait-based estimate",72)),
                     species = rep(traitsDF$Species, 2))


ggplot(errors1415, aes(x=type, y = error, groups = species))+
  geom_point() +
  geom_path() + # Brittany can you figure out how to draw lines between the estimates of the same species
  xlab("Approach") + ylab ("Error in ESW (m)")

#maybe we can label species with large error
ggplot(errors1415, aes(x = type, y = error, group = species)) +
  geom_point() +
  geom_line() +  # Draw lines between estimates of the same species
  geom_label_repel(data = subset(errors, abs(error) > 25), aes(label = species), size = 3.5) +  # Use ggrepel for better label placement
  xlab("Approach") +
  ylab("Error in ESW (m)")

#repeat this for the remaining years
# i.e., use 2015 data to predict for 2016, and 2016 data to predict 2017

## 2015 to 2016 ------------------------------------------------------------

#we'll use the 2015 data to make trait predictions 
lm1 <- lm(x2015 ~ log(Mass) + ForStrat.ground + Habitat + 
            Trophic.Level, data=traitsDF)
traitsDF$preds2015 <- predict(lm1)

#we'll then compare it to what was actually observed in 2015
summary(traitsDF$x2016 - traitsDF$x2015) #compare to past estimate
summary(traitsDF$x2016 - traitsDF$preds2015) #compare to trait-based estimate
CCC(traitsDF$x2016, traitsDF$x2015)$rho.c
CCC(traitsDF$x2016, traitsDF$preds2015)$rho.c
#looks like the past estimate does a bit better than the predicted estimate based on the traits

#plot individual errors
errors1516 <- data.frame(error = c(CCC(traitsDF$x2016, traitsDF$x2015)$blalt$delta,
                               CCC(traitsDF$x2016, traitsDF$preds2015)$blalt$delta),
                     type = c(rep("past estimate",72), rep("trait-based estimate",72)),
                     species = rep(traitsDF$Species, 2))

ggplot(errors1516, aes(x=type, y = error, groups = species))+
  geom_point() +
  geom_path() + # Brittany can you figure out how to draw lines between the estimates of the same species
  xlab("Approach") + ylab ("Error in ESW (m)")

#maybe we can label species with large error
ggplot(errors1516, aes(x = type, y = error, group = species)) +
  geom_point() +
  geom_line() +  # Draw lines between estimates of the same species
  geom_label_repel(data = subset(errors, abs(error) > 25), aes(label = species), size = 3.5) +  # Use ggrepel for better label placement
  xlab("Approach") +
  ylab("Error in ESW (m)")

## 2016 to 2017 ------------------------------------------------------------

#we'll use the 2016 data to make trait predictions 
lm1 <- lm(x2016 ~ log(Mass) + ForStrat.ground + Habitat + 
            Trophic.Level, data=traitsDF)
traitsDF$preds2016 <- predict(lm1)

#we'll then compare it to what was actually observed in 2015
summary(traitsDF$x2017 - traitsDF$x2016) #compare to past estimate
summary(traitsDF$x2017 - traitsDF$preds2016) #compare to trait-based estimate
CCC(traitsDF$x2017, traitsDF$x2016)$rho.c
CCC(traitsDF$x2017, traitsDF$preds2016)$rho.c
#looks like the past estimate does a bit better than the predicted estimate based on the traits

#plot individual errors
errors1617 <- data.frame(error = c(CCC(traitsDF$x2017, traitsDF$x2016)$blalt$delta,
                               CCC(traitsDF$x2017, traitsDF$preds2016)$blalt$delta),
                     type = c(rep("past estimate",72), rep("trait-based estimate",72)),
                     species = rep(traitsDF$Species, 2))

ggplot(errors1617, aes(x=type, y = error, groups = species))+
  geom_point() +
  geom_path() + # Brittany can you figure out how to draw lines between the estimates of the same species
  xlab("Approach") + ylab ("Error in ESW (m)")

# Draw lines between estimates of the same species
ggplot(errors1617, aes(x = type, y = error, group = species)) +
  geom_point() +
  geom_line() +  
  xlab("Approach") +
  ylab("Error in ESW (m)")

#maybe we can label species with large error
ggplot(errors1617, aes(x = type, y = error, group = species)) +
  geom_point() +
  geom_line() +  # Draw lines between estimates of the same species
  geom_label_repel(data = subset(errors, abs(error) > 25), aes(label = species), size = 3.5) +  # Use ggrepel for better label placement
  xlab("Approach") +
  ylab("Error in ESW (m)")

#we'll need to think about how we summarise this across years

## Combined plot ------------------------------------------------------

# combine data
errors1415 <- errors1415 %>%
  mutate(years="2014 to 2015")

errors1516 <- errors1516 %>%
  mutate(years="2015 to 2016")

errors1617 <- errors1617 %>%
  mutate(years="2016 to 2017")

errors_all <- rbind(errors1415, errors1516, errors1617)

# plot the data
panelb <- ggplot(errors_all, aes(x = type, y = error, group = species)) +
  geom_point(aes(color=years)) +
  geom_line(aes(color=years)) +  # Draw lines between estimates of the same species
  xlab("Approach") +
  ylab("Error in ESW (m)") +
  facet_wrap(~ years) +
  scale_colour_manual(
    values = c(
      "2014 to 2015" = swatch()[1],
      "2015 to 2016" = swatch()[5],
      "2016 to 2017" = swatch()[9]
    )) +
  scale_x_discrete(labels = c("past estimate" = "Past estimate", "trait-based estimate" = "Trait-based estimate")) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12))
panelb

# population size estimate --------------------------------------------

# lets see the implications of using the predicted/previous ESW instead of observed ESW for a year

# lets continue with the scenario that was want to predict for 2015 but only have data for 2014

#get site-level relative abundance predictions (not corrected for detection)
sitePreds <- readRDS("Data/xgsitepredsData.rds") %>% 
  dplyr::select(Species, preds) %>%
  rename(Rel_Abund = preds)

#merge
sitePreds <- sitePreds %>% inner_join(., traitsDF, by = "Species")

## 2014/2015 ---------------------------------------------------------------

#correct all the site level predictions by detection and sum abundances across all sites (= total pop size)
popOutput1415 <- lapply(sort(unique(sitePreds$Species)), function(myspecies){
  
  sitePreds_sp <- sitePreds %>% filter(Species == myspecies) 
  
  #correct the predictions by the two estimates of ESW (previous estimate or trait-based estimate)
  sitePreds_sp <- sitePreds_sp %>%
    mutate(True_Abs_Abund = Rel_Abund * 100/x2015, # using truth
           Prev_Abs_Abund = Rel_Abund * 100/x2014, # using previous estimate
           Trait_Abs_Abund = Rel_Abund * 100/preds2014) # using trait-based estimate
  
  #calculate the total population size using the observed ESW and predicted ESW
  sitePreds_sp %>%
    summarise(True_total_pop = sum(True_Abs_Abund),
              Prev_total_pop = sum(Prev_Abs_Abund),
              Trait_total_pop = sum(Trait_Abs_Abund)) %>%
    mutate(Pop_error_w_previous = (True_total_pop - Prev_total_pop)/True_total_pop,
           Pop_error_w_traits = (True_total_pop - Trait_total_pop)/True_total_pop) %>%
    add_column(Species = myspecies)
  
}) %>% bind_rows()

popOutput %>%
  dplyr::select(Species, Pop_error_w_previous, Pop_error_w_traits) %>%
  pivot_longer(cols=2:3, names_to ="Type", values_to ="error") %>%
  ggplot(aes(x=Type, y = abs(error)))+
  geom_boxplot()+
  geom_jitter(colour="black", alpha = 0.5) +
  xlab("Prediction approach") + ylab("Population size error (%)")


# repeat the above for 2016/2015 and 2017/2016

## 2015/2016 ---------------------------------------------------------------

#correct all the site level predictions by detection and sum abundances across all sites (= total pop size)
popOutput1516 <- lapply(sort(unique(sitePreds$Species)), function(myspecies){
  
  sitePreds_sp <- sitePreds %>% filter(Species == myspecies) 
  
  #correct the predictions by the two estimates of ESW (previous estimate or trait-based estimate)
  sitePreds_sp <- sitePreds_sp %>%
    mutate(True_Abs_Abund = Rel_Abund * 100/x2016, # using truth
           Prev_Abs_Abund = Rel_Abund * 100/x2015, # using previous estimate
           Trait_Abs_Abund = Rel_Abund * 100/preds2015) # using trait-based estimate
  
  #calculate the total population size using the observed ESW and predicted ESW
  sitePreds_sp %>%
    summarise(True_total_pop = sum(True_Abs_Abund),
              Prev_total_pop = sum(Prev_Abs_Abund),
              Trait_total_pop = sum(Trait_Abs_Abund)) %>%
    mutate(Pop_error_w_previous = (True_total_pop - Prev_total_pop)/True_total_pop,
           Pop_error_w_traits = (True_total_pop - Trait_total_pop)/True_total_pop) %>%
    add_column(Species = myspecies)
  
}) %>% bind_rows()

popOutput %>%
  dplyr::select(Species, Pop_error_w_previous, Pop_error_w_traits) %>%
  pivot_longer(cols=2:3, names_to ="Type", values_to ="error") %>%
  ggplot(aes(x=Type, y = abs(error)))+
  geom_boxplot()+
  geom_jitter(colour="black", alpha = 0.5) +
  xlab("Prediction approach") + ylab("Population size error (%)")


## 2016/2017 ---------------------------------------------------------------

#correct all the site level predictions by detection and sum abundances across all sites (= total pop size)
popOutput1617 <- lapply(sort(unique(sitePreds$Species)), function(myspecies){
  
  sitePreds_sp <- sitePreds %>% filter(Species == myspecies) 
  
  #correct the predictions by the two estimates of ESW (previous estimate or trait-based estimate)
  sitePreds_sp <- sitePreds_sp %>%
    mutate(True_Abs_Abund = Rel_Abund * 100/x2017, # using truth
           Prev_Abs_Abund = Rel_Abund * 100/x2016, # using previous estimate
           Trait_Abs_Abund = Rel_Abund * 100/preds2016) # using trait-based estimate
  
  #calculate the total population size using the observed ESW and predicted ESW
  sitePreds_sp %>%
    summarise(True_total_pop = sum(True_Abs_Abund),
              Prev_total_pop = sum(Prev_Abs_Abund),
              Trait_total_pop = sum(Trait_Abs_Abund)) %>%
    mutate(Pop_error_w_previous = (True_total_pop - Prev_total_pop)/True_total_pop,
           Pop_error_w_traits = (True_total_pop - Trait_total_pop)/True_total_pop) %>%
    add_column(Species = myspecies)
  
}) %>% bind_rows()

popOutput %>%
  dplyr::select(Species, Pop_error_w_previous, Pop_error_w_traits) %>%
  pivot_longer(cols=2:3, names_to ="Type", values_to ="error") %>%
  ggplot(aes(x=Type, y = abs(error)))+
  geom_boxplot()+
  geom_jitter(colour="black", alpha = 0.5) +
  xlab("Prediction approach") + ylab("Population size error (%)")

## Combine plots ----------------------------------------------------------------------

# clean and combine the data 
popOutput1415 <- popOutput1415 %>%
  mutate(year="2014 to 2015") %>%
  dplyr::select(Species, Pop_error_w_previous, Pop_error_w_traits, year) %>%
  pivot_longer(cols=2:3, names_to ="Type", values_to ="error")

popOutput1516 <- popOutput1516 %>%
  mutate(year="2015 to 2016") %>%
  dplyr::select(Species, Pop_error_w_previous, Pop_error_w_traits, year) %>%
  pivot_longer(cols=2:3, names_to ="Type", values_to ="error")

popOutput1617 <- popOutput1617 %>%
  mutate(year="2016 to 2017") %>%
  dplyr::select(Species, Pop_error_w_previous, Pop_error_w_traits, year) %>%
  pivot_longer(cols=2:3, names_to ="Type", values_to ="error")

popOutput_all <- rbind(popOutput1415, popOutput1516, popOutput1617)

# combine all the results on the final plot - use different colours for different years, so the ggplot becomes
ggplot(popOutput_all, aes(x=Type, y = abs(error)))+
  geom_boxplot(aes(fill=year))+
  geom_jitter(aes(colour=year), alpha = 0.5) +
  xlab("Prediction approach") + 
  ylab("Population size error (%)") +
  labs(fill="Year Comparison")

# Create the faceted plot
panelc <- ggplot(popOutput_all, aes(x = Type, y = abs(error))) +
  geom_boxplot(aes(fill = year), alpha = 0.5) +  # Boxplots colored by year
  geom_jitter(aes(colour = year), alpha = 0.5, width = 0.2) +  # Points colored by year
  xlab("Prediction approach") + 
  ylab("Population size error (%)") +
  scale_fill_manual(
    values = c(
      "2014 to 2015" = swatch()[1],
      "2015 to 2016" = swatch()[5],
      "2016 to 2017" = swatch()[9]
    )) +
  scale_colour_manual(
    values = c(
      "2014 to 2015" = swatch()[1],
      "2015 to 2016" = swatch()[5],
      "2016 to 2017" = swatch()[9]
    )) +
  labs(fill = "Year Comparison", colour = "Year Comparison") +  # Legend titles
  facet_wrap(~ year, scales = "free_y") +  # Facet by year
  scale_x_discrete(labels = c("Pop_error_w_previous" = "Past estimate", "Pop_error_w_traits" = "Trait-based estimate")) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12))  # Remove the legend
panelc

# Combine all three plots ------------------------------------------------------

#Figure 2c and Figure 2d
#combine with Figure 2a and 2b above
combined_fig <- (patchwork::wrap_elements(panela) / 
                   patchwork::wrap_elements(panelb) /
                   patchwork::wrap_elements(panelc)) +
  plot_annotation(tag_levels = list(c("a", "b", "c")))
combined_fig
#save in plots folder

ggsave("Figures/annual_analysis_fig.jpeg", height=8, width=8, units="in")

# end --------------------------------------------------------------------------
