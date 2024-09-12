# libraries --------------------------------------------------------------------
library(tidyverse)
library(broom)
library(lubridate)
library(patchwork)
#devtools::install_github('cttobin/ggthemr')
library(ggthemr)
ggthemr("fresh") 
library(ggExtra)

set.seed(1039)

# DOF data ---------------------------------------------------------------------

source('00_functions.R', echo=TRUE)

# Get estimated ESW ------------------------------------------------------------

#these are calculated in the script 02.1_ESW_analysis.R
detDF <- readRDS("Data/distance_passerines_constant.rds")
#for each species, we have an estimate of its estimated strip width (ESW)
#(ESW is essentially a measure of species detectability)

# get traits file --------------------------------------------------------------

#already compiled ecological and morphological traits for each species
traits <- readRDS("Data/traits.rds")
head(traits)

#some tweaks
table(traits$Habitat)
traits$Habitat[traits$Habitat=="Woodland"] <- "Forest"

#for black redstart, change to:
traits$Habitat[traits$Habitat=="Rock"] <- "Human Modified"

#remove L Canus
traits <- traits %>% filter(Species!= "Larus canus")

traitsDF <- detDF %>%
  inner_join(., traits, by="Species") 


# traits selection --------------------------------------------------------

#The brt allows us to include more traits too.
#but there are still a few too many, so lets remove the below
traitsDF <- traitsDF %>%
  select_if(~ !any(is.na(.)))  %>%
  mutate_if(., is.character, as.factor) %>%
  dplyr::select(-ID, -Order, -Family, 
                -Association.outside.the.breeding.season,-Data.source) %>%
  dplyr::select(-Association.during.nesting, -Family1, -Order1,-Sequence,-Avibase.ID1,-Fish_Y,
                -Desert, -Fish_B,-Savanna, -Facultative.migrant,-Marine, -Nocturnal, -Reed,
                -Sexual.dimorphism,-Sedentary,-Short.distance.migrant, -Long.distance.migrant, -EggL_MEAN,
                -EggW_MEAN, -BodyMass.Value, -WeightU_MEAN, -WeightM_MEAN,-WeightF_MEAN,
                -WingM_MEAN, -WingF_MEAN, -TailM_MEAN,-TailF_MEAN, -BillM_MEAN, -BillF_MEAN,
                -TarsusM_MEAN,-TarsusF_MEAN, -Centroid.Latitude) 

# BRT ------------------------------------------------------------------------

#Brittany, this is currently a copy of the lm script
# we are going to repeat that analysis except using a brt model instead of lm
# you can read about brt here
# https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2656.2008.01390.x

#the only difference will be this (instead of lm) 
library(dismo)
library(gbm)

gbm1 <- gbm.step(data=traitsDF, gbm.x =9:ncol(traitsDF), gbm.y = "ESW", family = "gaussian",
                 tree.complexity = 3, learning.rate = 0.0001, step.size = 10,
                 bag.fraction = 0.75, n.folds = 5, n.trees = 2000)
# I tried different variables for learning rate, but 0.0001 was necessary for the 
# leave one analysis later on

# we can look at their relative like this:
summary(gbm1)
gbm.plot(gbm1, n.plots=12)
#save these in the plots folder

# now let's look at the raw points and distribution of data
# use the v parameter to specify which plot you want to view.
gbm.plot.fits(gbm1, v=1)

# let's look at pairwise interactions in the data
interactions <- gbm.interactions(gbm1)
interactions

#you can use 'predict' in the same way as for lm

# modify the sections of the script as needed below to use this brt model instead of the lm

## predictions -------------------------------------------------------------------

# single species 

# do a leave-one-out type thing
# where each species is left out at a time and
# then their ESW is predicted

leave_one_out_function <- function(species_name){
  
  dat_filtered <- traitsDF %>%
    dplyr::filter(Species != species_name)
  
  gbm1 <- gbm.step(data=traitsDF, gbm.x =9:ncol(traitsDF), gbm.y = "ESW", family = "gaussian",
                   tree.complexity = 3, learning.rate = 0.0001, step.size = 10,
                   bag.fraction = 0.75, n.folds = 5, n.trees = 2000)
  
  new_dat <- traitsDF %>%
    dplyr::filter(Species==species_name)
  
  predicted_ESW <- predict.gbm(gbm1, new_dat,
                               n.trees=gbm1$gbm.call$best.trees, 
                               type="response")
  
  loo_summary <- data.frame(Species=species_name,
                            observed_ESW=new_dat$ESW,
                            predicted_ESW=predicted_ESW)
  
  return(loo_summary)
  
}

# the following code will take over an hour to run
leave_one_out_results <- bind_rows(lapply(unique(traitsDF$Species), 
                                          leave_one_out_function))

saveRDS(leave_one_out_results, file="Data/leave_one_out_results_brt.rds")
leave_one_out_results <- readRDS("Data/leave_one_out_results_brt.rds")

#Figure 2
Fig2a <- ggplot(leave_one_out_results, aes(x=observed_ESW, y=predicted_ESW))+
  geom_point()+
  theme_minimal()+
  theme(axis.text=element_text(color="black"),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5))+
  xlab("Observed ESW")+
  ylab("Predicted ESW")+
  geom_smooth(method="lm")
Fig2a

ggsave("Figures/Figure2_brt.jpeg", height=5, width=5, units="in")

#save this plot in the plots folder

cor.test(leave_one_out_results$observed_ESW,
         leave_one_out_results$predicted_ESW)


## simplified model regression -------------------------------------------------

# in the linear model code, we tried running a linear model with just mass
# in this code, we can try running the BRT with only variables that the gbm.simplify 
# function identifies as relevant

#Figure 2a can be the plot above using the full model and 
# Figure 2b can be same but for the model using only body size

# now, let's try simplying the model to see which variables can be dropped
gbm1_simp <- gbm.simplify(gbm1, n.drops = 20)
# this code found 9 variables to drop

# let's try runnnig the model again with these rows dropped.
gbm1_reduced <- gbm.step(data=traitsDF, gbm.x =gbm1_simp$pred.list[[1]], gbm.y = "ESW", family = "gaussian",
                         tree.complexity = 3, learning.rate = 0.0001, step.size = 10,
                         bag.fraction = 0.75, n.folds = 5, n.trees = 2000)

summary(gbm1_reduced)

leave_one_out_function <- function(species_name){
  
  dat_filtered <- traitsDF %>%
    dplyr::filter(Species != species_name)
  
  gbm1 <- gbm.step(data=traitsDF, gbm.x =gbm1_simp$pred.list[[1]], gbm.y = "ESW", family = "gaussian",
                   tree.complexity = 3, learning.rate = 0.0001, step.size = 10,
                   bag.fraction = 0.75, n.folds = 5, n.trees = 2000)
  
  new_dat <- traitsDF %>%
    dplyr::filter(Species==species_name)
  
  predicted_ESW <- predict.gbm(gbm1, new_dat,
                               n.trees=gbm1$gbm.call$best.trees, 
                               type="response")
  
  loo_summary <- data.frame(Species=species_name,
                            observed_ESW=new_dat$ESW,
                            predicted_ESW=predicted_ESW)
  
  return(loo_summary)
  
}

# the following code will take over an hour to run
leave_one_out_results_reduced <- bind_rows(lapply(unique(traitsDF$Species), 
                                          leave_one_out_function))

saveRDS(leave_one_out_results_reduced, file="Data/leave_one_out_results_brt_reduced.rds")
leave_one_out_results_reduced <- readRDS("Data/leave_one_out_results_brt_reduced.rds")

#Figure 2
Fig2b <- ggplot(leave_one_out_results_reduced, aes(x=observed_ESW, y=predicted_ESW))+
  geom_point()+
  theme_minimal()+
  theme(axis.text=element_text(color="black"),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5))+
  xlab("Observed ESW")+
  ylab("Predicted ESW")+
  geom_smooth(method="lm")
Fig2b

cor.test(leave_one_out_results_reduced$observed_ESW,
         leave_one_out_results_reduced$predicted_ESW)

# How many species do we need -------------------------------------------

# # Brittany, ignore this section for the moment
# # These plots are for the SI only
# 
# # multiple species
# # how many species can we leave out at a time and still get a good predictive accuracy of the model?
# 
# #couple of functions:
# random_sampling_function <- function(sample_size){
# 
#   boot_fun <- function(draw_number, sample_size){
# 
#     dat_filtered <- traitsDF %>%
#       sample_n(nrow(traitsDF)-sample_size)
# 
#     mod <- lm(ESW ~ log(Mass) + ForStrat.ground + Habitat + 
#                 Trophic.Level + flockSize,
#               data=dat_filtered)
# 
#     new_dat <- traitsDF %>%
#       dplyr::filter(!Species %in% dat_filtered$Species)
# 
#     predicted_ESW <- predict(mod, new_dat)
# 
#     loo_summary <- data.frame(Species=new_dat$Species,
#                               observed_ESW=new_dat$ESW,
#                               predicted_ESW=predicted_ESW) %>%
#       mutate(draw=draw_number)
# 
#     return(loo_summary)
# 
#   }
# 
#   boot_results <- bind_rows(lapply(c(1:100), function(i){ boot_fun(draw=i)}
#                                    )) %>%
#     mutate(sample=sample_size) %>%
#     mutate(percent_missing=(sample/nrow(traitsDF))*100)
# 
# 
# }
# 
# lapply_with_error <- function(X,FUN,...){
#   lapply(X, function(x, ...) tryCatch(FUN(x, ...),
#                                       error=function(e) NULL))
# }
# 
# # adding the lapply with error
# # allows to 'skip'
# # the possibility when something is left out
# # that has a level not in the model
# # e.g., primary lifestyle this happens
# 
# #Lets run it for 20, 30, 40 etc numbers of species
# random_sampling_results <- bind_rows(lapply_with_error(c(20,30,40,50,60,70),
#                                                        random_sampling_function))
# 
# # plotting
# random_sampling_results %>%
#   mutate(percent_missing=as.numeric(percent_missing)) %>%
#   nest(data = -percent_missing) %>%
#   mutate(
#     fit = map(data, ~ lm(predicted_ESW ~ observed_ESW, data = .x)),
#     tidied = map(fit, tidy),
#     glanced=map(fit, glance)
#   ) %>%
#   unnest(glanced) %>%
#   ggplot(., aes(x=percent_missing, y=r.squared))+
#   geom_jitter(aes(colour=sample))
# 
# # so that kind of makes sense
# # but with the lapply with error....
# # by the time you have the majority of species being 'left-out'
# # then the levels are gone so it isn't working.
# # just as a test case I'm gonna copy and paste the above
# # but just use body mass (as we know we have this for all species and don't have to worry about levels)
# # but how many species can we leave out at a time and still get a good predictive accuracy of the model?
# 
# # body size only
# random_sampling_function_body <- function(sample_size){
# 
#   boot_fun <- function(draw_number){
# 
#     dat_filtered <- traitsDF %>%
#       sample_n(nrow(traitsDF)-sample_size)
# 
#     mod <- lm(ESW ~ log(Mass),
#               data=dat_filtered)
# 
#     new_dat <- traitsDF %>%
#       dplyr::filter(!Species %in% dat_filtered$Species)
# 
#     predicted_ESW <- predict(mod, new_dat)
# 
#     loo_summary <- data.frame(Species=new_dat$Species,
#                               observed_ESW=new_dat$ESW,
#                               predicted_ESW=predicted_ESW) %>%
#       mutate(draw=draw_number)
# 
#     return(loo_summary)
# 
#   }
# 
#   boot_results <- bind_rows(lapply(c(1:50), boot_fun)) %>%
#     mutate(sample=sample_size) %>%
#     mutate(percent_missing=(sample/nrow(traitsDF))*100)
# 
# 
# }
# 
# random_sampling_results <- bind_rows(lapply_with_error(c(1:60), random_sampling_function_body))
# 
# # try again?
# random_sampling_results %>%
#   mutate(percent_missing=as.numeric(percent_missing)) %>%
#   nest(data = -percent_missing) %>%
#   mutate(
#     fit = map(data, ~ lm(predicted_ESW ~ observed_ESW, data = .x)),
#     tidied = map(fit, tidy),
#     glanced=map(fit, glance)
#   ) %>%
#   unnest(glanced) %>%
#   ggplot(., aes(x=percent_missing, y=r.squared))+
#   geom_jitter()

# population size estimate --------------------------------------------

#We now want to understand the implications of using the predicted ESW instead 
# of observed ESW 

#get site-level relative abundance predictions (not corrected for detection)
sitePreds <- readRDS("Data/xgsitepredsData.rds") %>% 
              dplyr::select(Species, preds) %>%
              rename(Rel_Abund = preds)

#also get estimates of detection rates (ESW) calculated above
detectionRates <- readRDS("Data/leave_one_out_results_brt.rds")

#merge
sitePreds <- sitePreds %>% inner_join(., detectionRates, by = "Species")

#correct all the site level predictions by detection and sum abundances across all sites (= total pop size)
popOutput <- lapply(sort(unique(sitePreds$Species)), function(myspecies){
  
sitePreds_sp <- sitePreds %>% filter(Species == myspecies) 

#correct the predictions by the estimated detction rates (ESW)
sitePreds_sp <- sitePreds_sp %>%
                  mutate(Obs_Abs_Abund = Rel_Abund * 100/observed_ESW,
                         Pred_Abs_Abund = Rel_Abund * 100/predicted_ESW)

#calculate the total population size using the observed ESW and predicted ESW
sitePreds_sp %>%
            summarise(Obs_total_pop = sum(Obs_Abs_Abund),
                      Pred_total_pop = sum(Pred_Abs_Abund),
                      Obs_ESW = mean(observed_ESW),
                      Pred_ESW = mean(predicted_ESW)) %>%
            mutate(Pop_error = (Obs_total_pop - Pred_total_pop)/Obs_total_pop,
                   ESW_error = Obs_ESW - Pred_ESW) %>%
            add_column(Species = myspecies)

}) %>% bind_rows()

# plot the relationship and include this as Fig 2c and d (combine with Figur 2a and 2b above)
g_error <- ggplot(popOutput) +
  geom_point(aes(x = ESW_error, y = Pop_error)) +
  geom_hline(yintercept = 0, linetype ="dashed") +
  geom_vline(xintercept = 0, linetype ="dashed") +
  xlab("Effective strip width error (m)") +
  ylab("Population size error (%)")

g_Error <- ggMarginal(g_error, type="density", fill = "lightblue")

g_abs <- ggplot(popOutput) +
  geom_point(aes(x = abs(ESW_error), y = abs(Pop_error))) +
  geom_hline(yintercept = 0, linetype ="dashed") +
  geom_vline(xintercept = 0, linetype ="dashed") +
  xlab("Absolute effective strip width error (m)") +
  ylab("Absolute population size error (%)")

g_Abs <- ggMarginal(g_abs, type="boxplot", fill = "lightblue")
  

#Figure 2c and Figure 2d
#combine with Figure 2a and 2b above
combined_fig2 <- (patchwork::wrap_elements(Fig2a) + patchwork::wrap_elements(Fig2b))/
  (patchwork::wrap_elements(g_Error) + patchwork::wrap_elements(g_Abs)) + 
  plot_annotation(tag_levels = list(c("a", "b", "c", "d")))
combined_fig2
#save in plots folder

ggsave("Figures/Fig2_brt.jpeg", height=7, width=7, units="in")

# end --------------------------------------------------------------------------