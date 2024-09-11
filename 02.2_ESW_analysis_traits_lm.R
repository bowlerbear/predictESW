# libraries --------------------------------------------------------------------
library(tidyverse)
library(broom)
library(lubridate)
library(patchwork)
# devtools::install_github('cttobin/ggthemr')
library(ggthemr)
ggthemr("fresh") 
library(ggExtra)

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

# plot relationships -----------------------------------------------------------

#We want to look at relationships between species traits and their ESW

#Figure 1 of the paper

#rewrite these using ggplot, instead of qplot
#and generally tidy up

g1 <- qplot(Mass, ESW, data=traitsDF) + 
  stat_smooth(method="lm") +
  ylab("Effective strip width") +
  scale_x_log10()

g2 <- qplot(Habitat, ESW, data=traitsDF, geom="boxplot") + 
  geom_jitter(width = 0.2,colour="black", alpha = 0.5)+
  coord_flip() 

g3 <- qplot(Trophic.Level, ESW, data=traitsDF, geom="boxplot") +
  geom_jitter(width = 0.2,colour="black", alpha = 0.5)+
  xlab("Trophic niche") +
  coord_flip()

g4 <- qplot(ForStrat.ground, ESW, data=traitsDF) + 
  stat_smooth(method="lm")+
  xlab("Ground foraging  %")

g5 <- qplot(as.factor(flockSize), ESW, data=traitsDF, geom="boxplot") +
  geom_jitter(width = 0.2,colour="black", alpha = 0.5)+
  xlab("Flocking")

#please tidy this up a bit, maybe use patchwork instead of cowplot?
# Update your plots to reduce white space between the y-axis label and the axis
gRest <- cowplot::plot_grid(g4, g5,
          g2, g3,
          nrow=2,
          labels=c("b","c","d","e"))

cowplot::plot_grid(g1,
          gRest, nrow=2,
          labels = c("a", ""))


# Arrange g4, g5, g2, and g3 in a 2x2 grid with labels
gRest <- (g4 | g5) / 
  (g2 | g3) + 
  plot_annotation(tag_levels = list(c("b", "c", "d", "e")))

# Combine g1 with gRest, with g1 on top and gRest below
final_plot <- g1 / gRest + 
  plot_layout(nrow = 2) +  # Ensure the layout has 2 rows
  plot_annotation(tag_levels = list(c("a", "b", "c", "d", "e")))  # Add labels for the combined plot

# Print the final plot
final_plot
ggsave("Figures/Figure1.jpeg", width=7, height=7, units="in")

#save plot in plots folder

# all traits regression ------------------------------------------------------------------

hist(traitsDF$ESW)

lm1 <- lm(ESW ~ log(Mass) + ForStrat.ground + Habitat + 
             Trophic.Level + flockSize, data=traitsDF)
summary(lm1)

#pull out coefficients and plot them
coefDF <- data.frame(Param = names(coefficients(lm1)),
                     estimate = as.numeric(coefficients(lm1)),
                     lower = as.numeric(confint(lm1)[,1]),
                     upper = as.numeric(confint(lm1)[,2]))

ggplot(coefDF)+
  geom_pointrange(aes(x=Param, y=estimate, ymin=lower, ymax=upper))+
  geom_hline(yintercept=0, linetype="dashed")+
  coord_flip()

# drop to significant variables
lm1 <- lm(ESW ~ log(Mass) + ForStrat.ground + Habitat + 
            Trophic.Level, data=traitsDF)
summary(lm1)

#maybe organise this into a nicely formatted table?
#https://modelsummary.com/

library(modelsummary)
modelsummary(lm1, output="Figures/Linear_Model_Ouput.docx")

## predictions -------------------------------------------------------------------

# single species 

# do a leave-one-out type thing
# where each species is left out at a time and
# then their ESW is predicted

leave_one_out_function <- function(species_name){
  
  dat_filtered <- traitsDF %>%
    dplyr::filter(Species != species_name)
  
  mod <- lm(ESW ~ log(Mass) + ForStrat.ground + Habitat + 
              Trophic.Level,
            data=dat_filtered)
  
  new_dat <- traitsDF %>%
    dplyr::filter(Species==species_name)
  
  predicted_ESW <- predict(mod, new_dat)
  
  loo_summary <- data.frame(Species=species_name,
                            observed_ESW=new_dat$ESW,
                            predicted_ESW=predicted_ESW)
  
  return(loo_summary)
  
}

leave_one_out_results <- bind_rows(lapply(unique(traitsDF$Species), 
                                          leave_one_out_function))

saveRDS(leave_one_out_results, file="Data/leave_one_out_results.rds")

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

ggsave("Figures/Figure2.jpeg", height=5, width=5, units="in")
#save this plot in the plots folder

cor.test(leave_one_out_results$observed_ESW,
         leave_one_out_results$predicted_ESW)


## body size regression -------------------------------------------------

#repeat the previous section of 'all traits regression'
# but only using body size
# so mod <- lm(ESW ~ log(Mass), data=dat_filtered)
#everything else the same

#Figure 2a can be the plot above using the full model and 
# Figure 2b can be same but for the model using only body size

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
sitePreds <- readRDS("outputs/xgsitepredsData.rds") %>% 
              dplyr::select(Species, preds) %>%
              rename(Rel_Abund = preds)

#also get estimates of detection rates (ESW) calculated above
detectionRates <- readRDS("outputs/leave_one_out_results.rds")

#merge
sitePreds <- sitePreds %>% inner_join(., detectionRates, by = "Species")

#also get original effective strip estimates with error on them (we've ignored the original error till now)
detDF$ESW_plus_error <- detDF$ESW + (2 * detDF$ESW_se)
sitePreds$ESW_plus_error <- detDF$ESW_plus_error[match(sitePreds$Species, detDF$Species)]
#we'll use this as our null error of the original observed estimate    

#correct all the site level predictions by detection and sum abundances across all sites (= total pop size)
popOutput <- lapply(sort(unique(sitePreds$Species)), function(myspecies){
  
sitePreds_sp <- sitePreds %>% filter(Species == myspecies) 

#correct the predictions by the estimated detction rates (ESW)
sitePreds_sp <- sitePreds_sp %>%
                  mutate(Obs_Abs_Abund = Rel_Abund * 100/observed_ESW,
                         Pred_Abs_Abund = Rel_Abund * 100/predicted_ESW,
                         Obs_Error_Abs_Abund = Rel_Abund * 100/ESW_plus_error)

#calculate the total population size using the observed ESW and predicted ESW
sitePreds_sp %>%
            summarise(Obs_total_pop = sum(Obs_Abs_Abund),
                      Pred_total_pop = sum(Pred_Abs_Abund),
                      Obs_Error_total_pop = sum(Obs_Error_Abs_Abund),
                      Obs_ESW = mean(observed_ESW),
                      Pred_ESW = mean(predicted_ESW)) %>%
            mutate(Pop_error = (Obs_total_pop - Pred_total_pop)/Obs_total_pop,
                   Null_error = (Obs_total_pop - Obs_Error_total_pop)/Obs_total_pop,
                   ESW_error = Obs_ESW - Pred_ESW) %>%
            add_column(Species = myspecies)

}) %>% bind_rows()

# plot the relationship and include this as Fig 2c and d (combine with Figur 2a and 2b above)
g_error <- ggplot(popOutput) +
  geom_point(aes(x = ESW_error, y = Pop_error)) +
  geom_hline(yintercept =  mean(abs(popOutput$Null_error)), linetype ="dashed", color = "red") + # null error
  geom_hline(yintercept =  mean(-1*abs(popOutput$Null_error)), linetype ="dashed", color = "red") + # null error
  geom_hline(yintercept = 0, linetype ="dashed") +
  geom_vline(xintercept = 0, linetype ="dashed") +
  xlab("Effective strip width error (m)") +
  ylab("Population size error (%)")

g_Error <- ggMarginal(g_error, type="density", fill = "lightblue")

g_abs <- ggplot(popOutput) +
  geom_point(aes(x = abs(ESW_error), y = abs(Pop_error))) +
  geom_hline(yintercept =  mean(abs(popOutput$Null_error)), linetype ="dashed", color = "red") + # null error
  geom_hline(yintercept = 0, linetype ="dashed") +
  geom_vline(xintercept = 0, linetype ="dashed") +
  xlab("Absolute effective strip width error (m)") +
  ylab("Absolute population size error (%)")

g_Abs <- ggMarginal(g_abs, type="boxplot", fill = "lightblue")
  

#Figure 2c and Figure 2d
#combine with Figure 2a and 2b above
patchwork::wrap_elements(g_Error) + patchwork::wrap_elements(g_Abs)
#save in plots folder

# covariate effects ------------------------------------------------------------

#finallly we can ask weather we can use traits to predict environmental effects on the ESW
#this output includes the effect of two covariates tested
#forest cover - we expect higher forest cover leads to lower detections
#paths/roads - we expect some species to be attracted to roads and others to be deterred

detDF <- readRDS("outputs/distance_passerines_covs.rds")

traitsDF <- detDF %>%
  inner_join(., traits, by="Species") 

#save regression output in format for paper
lm1 <- lm(Z_forest ~ Habitat, data=traitsDF)
summary(lm1)
lm1 <- lm(Z_pathroad ~ Habitat, data=traitsDF)
summary(lm1)

g1 <- qplot(Habitat, Z_forest, data=traitsDF, geom="boxplot") + 
  geom_jitter(width = 0.2,colour="black", alpha = 0.5)+
  coord_flip() +
  ylab("Effect of forest cover on ESW")

g2 <- qplot(Habitat, Z_pathroad, data=traitsDF, geom="boxplot") + 
  geom_jitter(width = 0.2,colour="black", alpha = 0.5)+
  coord_flip() +
  ylab("Effect of path/road on ESW")

g1 + g2
#save in plots folder

# end --------------------------------------------------------------------------