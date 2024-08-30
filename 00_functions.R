lapply_with_error <- function(X,FUN,...){    
  lapply(X, function(x, ...) tryCatch(FUN(x, ...),
                                      error=function(e) NULL))
}


plotMap <- function(myspecies, data){
  
  #and just one species
  dataS <- data %>% 
    filter(Species==myspecies) 
  
  infoS <- info %>%
    filter(type !="vinter")
  
  #get all transects for the transe
  allData <- left_join(infoS, dataS) %>%
    mutate(total = ifelse(is.na(total), 0, total),
           X.0 = ifelse(is.na(X.0), 0, X.0),
           X.1 = ifelse(is.na(X.1), 0, X.1),
           X.2 = ifelse(is.na(X.2), 0, X.2)) %>%
    mutate(Species = ifelse(is.na(Species), myspecies, Species))
  
  #add on environmental data
  allData <- inner_join(allData, environData)
  
  #add coordinates
  allData <- addCoords(allData)
  
  #cap the maximum number observation at the 99% quantile
  #top99 <- as.numeric(round(quantile(allData$total,0.99)))
  #allData$total <- ifelse(allData$total > top99,
  #                        top99, allData$total)
  
  #sum counts per X and Y
  summaryData <- allData %>%
    group_by(X,Y) %>%
    summarise(total = mean(total,na.rm=T))
  
  ggplot(summaryData) +
    geom_point(aes(x=X, y=Y, colour=(total+0.01)),size=rel(0.9))+
    scale_color_viridis_c(direction=-1, trans="log10")
  
}

addCoords <- function(df){
  
  coordsDF <- readRDS("environ-data/line_coords.rds") %>%
                as_tibble() %>%
                rename(X = x, Y = y)
  
  df <- inner_join(df, coordsDF, by="kvadratnr")
  
  return(df)
  
}

getSpeciesData <- function(myspecies, data){
  
  #and just one species
  dataS <- data %>% 
    filter(Species==myspecies) 
  
  infoS <- info %>%
    filter(type !="vinter")
  
  #get all transects for the transe
  allData <- left_join(infoS, dataS) %>%
    mutate(total = ifelse(is.na(total), 0, total),
           X.0 = ifelse(is.na(X.0), 0, X.0),
           X.1 = ifelse(is.na(X.1), 0, X.1),
           X.2 = ifelse(is.na(X.2), 0, X.2)) %>%
    mutate(Species = ifelse(is.na(Species), myspecies, Species))
  
  #cap each species count at 99%
  allData <- allData %>%
              group_by(Species) %>%
              mutate(total = ifelse(total>quantile(total,0.99),
                                    round(as.numeric(quantile(total, 0.99))),
                                               total)) %>%
              ungroup()
  
  #add on environmental data
  allData <- inner_join(allData, environData)
  
  #add coordinates
  allData <- addCoords(allData)
  
  return(allData)
  
}

### distance sampling ####

getDistanceData <- function(myspecies){
  
  #and just one species
  dataS <- data %>% 
    filter(Species==myspecies) %>%
    filter(type!="vinter")
  
  #join data - all years
  allData <- left_join(dataS,info)
  
  #add on environmental data
  allData <- inner_join(allData, environData)
  
  #get number of birds seen at each distance
  allData %>% 
    dplyr::select(Year, X.0, X.1, X.2) %>%
    pivot_longer(starts_with("X."), names_to = "distance", values_to = "nu") %>%
    dplyr::filter(!is.na(nu)) %>%
    dplyr::mutate(distance = recode(distance, X.0 = "25", X.1 = "50", X.2 = "100")) %>%
    dplyr::group_by(Year, distance) %>%
    dplyr::summarise(total = sum(nu)) %>%
    dplyr::mutate(distance = factor(distance, levels = c("25", "50", "100"))) %>%
    ungroup() %>%
    add_column(Species = myspecies)
  
}


getSiteDistanceData <- function(myspecies){
  
  #and just one species
  dataS <- data %>% 
    filter(Species==myspecies) %>%
    filter(type!="vinter")
  
  #join data - all years
  allData <- left_join(dataS,info)
  
  #get number of birds seen at each distance
  allData %>% 
    dplyr::select(kvadratnr, X.0, X.1, X.2) %>%
    pivot_longer(starts_with("X."), names_to = "distance", values_to = "nu") %>%
    dplyr::filter(!is.na(nu)) %>%
    dplyr::mutate(distance = recode(distance, X.0 = 1, X.1 = 2, X.2 = 3)) %>%
    dplyr::group_by(kvadratnr, distance) %>%
    dplyr::summarise(total = sum(nu)) %>%
    ungroup() %>%
    group_by(kvadratnr) %>%
    summarise(meanDistance = weighted.mean(distance,total)) %>%
    ungroup() %>%
    filter(!is.nan(meanDistance)) %>%
    add_column(Species = myspecies) %>%
    inner_join(., environData)
  
}


distanceModel <- function(myspecies, data, model="null"){

  #and just one species
  dataS <- data %>% 
    filter(Species==myspecies) %>%
    filter(type!="vinter") 
  
  infoS <- info %>%
    filter(type != "vinter")
  
  #get all transects for the obs
  allData <- inner_join(infoS, dataS)
  
  #make long
  detectionData <- allData %>%
    dplyr::select(kvadratnr,Year, X.0, X.1, X.2) %>%
    pivot_longer(starts_with("X."), names_to = "distance", values_to = "nu") %>%
    dplyr::filter(!is.na(nu)) %>%
    dplyr::mutate(distance = recode(distance, X.0 = "12.5", X.1 = "37.5", X.2 = "75")) %>%
    dplyr::mutate(distance = as.numeric(as.character(distance))) %>%
    rename(size = nu) %>%
    filter(size>0) %>%
    filter(!is.na(size))
  
  # #aggregate
  # detectionData <- detectionData %>%
  #                   group_by(kvadratnr, distance) %>%
  #                   summarise(size = sum(size))
  #   
  #add on environmental data
  detectionData <- inner_join(detectionData, environData)
  
  if(model=="null"){
    
  hr.model <- ds(detectionData, truncation=100,
                 transect="line",
                 formula = ~1,
                 cutpoints = c(0,25,50,100),
                 key = "hn")
  
  temp <- summary(hr.model)
  
  return(data.frame(Species = myspecies,
                    P = predict(hr.model)$fitted[1],
                    P_se = predict(hr.model,  se.fit=TRUE)$se.fit[1],
                    ESW = predict(hr.model,  esw=TRUE)$fitted[1],
                    ESW_se = predict(hr.model, esw=TRUE, se.fit=TRUE)$se.fit[1],
                    AIC = AIC(hr.model)$AIC,
                    N = nrow(detectionData)))
  
  }else if(model=="subset"){
    
  #subset to sites without paths
  detectionData <- detectionData %>%
                       filter(lines_pathroad_orig/lines_mapped < 0.5)
  
  hr.model <- ds(detectionData, truncation=100,
                 transect="line",
                 formula = ~1,
                 cutpoints = c(0,25,50,100),
                 key = "hn")
  
  temp <- summary(hr.model)
  
  return(data.frame(Species = myspecies,
                    P = predict(hr.model)$fitted[1],
                    P_se = predict(hr.model,  se.fit=TRUE)$se.fit[1],
                    ESW = predict(hr.model,  esw=TRUE)$fitted[1],
                    ESW_se = predict(hr.model, esw=TRUE, se.fit=TRUE)$se.fit[1],
                    AIC = AIC(hr.model)$AIC,
                    N = nrow(detectionData)))
  
  }else if(model=="covariates"){
  
  hr.model <- ds(detectionData, truncation=100,
                    transect="line",
                    formula = ~lines_pathroad + lines_forest,
                    cutpoints = c(0,25,50,100),
                    key = "hn")
  
  temp <- summary(hr.model)
  
  return(data.frame(Species = myspecies,
                    P = temp$ds$average.p,
                    P_se = temp$ds$average.p.se,
                    ESW = temp$ds$average.p*100,
                    ESW_se = temp$ds$average.p.se*100,
                    AIC = AIC(hr.model)$AIC,
                    N = nrow(detectionData),
                    Coef_pathroad = temp$ds$coeff$key.scale$estimate[2],
                    Coef_forest = temp$ds$coeff$key.scale$estimate[3],
                    Z_pathroad =  temp$ds$coeff$key.scale$estimate[2]/temp$ds$coeff$key.scale$se[2],
                    Z_forest =  temp$ds$coeff$key.scale$estimate[3]/temp$ds$coeff$key.scale$se[3]))

  }
  
}

resampleddistanceModel <- function(myspecies, data){
  
  #and just one species
  dataS <- data %>% 
    filter(Species==myspecies) %>%
    filter(type!="vinter") 
  
  infoS <- info %>%
    filter(type != "vinter")
  
  #get all transects for the obs
  allData <- inner_join(infoS, dataS)
  
  #make long
  detectionData <- allData %>%
    dplyr::select(kvadratnr,Year, X.0, X.1, X.2) %>%
    pivot_longer(starts_with("X."), names_to = "distance", values_to = "nu") %>%
    dplyr::filter(!is.na(nu)) %>%
    dplyr::mutate(distance = recode(distance, X.0 = "12.5", X.1 = "37.5", X.2 = "75")) %>%
    dplyr::mutate(distance = as.numeric(as.character(distance))) %>%
    rename(size = nu) %>%
    filter(size>0) %>%
    filter(!is.na(size))
  
    vals <- as.numeric()
      
    for(i in 1:1000){
  
    d <-   detectionData[sample(nrow(detectionData), replace=TRUE),]
      
    hr.model <- ds(d, truncation=100,
                   transect="line",
                   formula = ~1,
                   cutpoints = c(0,25,50,100),
                   key = "hn")
    
    temp <- summary(hr.model)
    vals[i] <- temp$ds$average.p * 100
    
    }
  
  
    return(data.frame(Species = myspecies,
                      ESW_samples = vals))
    
}


sitedistanceModel <- function(myspecies, data){
  
  #and just one species
  dataS <- data %>% 
    filter(Species==myspecies) %>%
    filter(type!="vinter") 
  
  infoS <- info %>%
    filter(type != "vinter")
  
  #get all transects for the obs
  allData <- inner_join(infoS, dataS)
  
  #make long
  detectionData <- allData %>%
    dplyr::select(kvadratnr,Year, X.0, X.1, X.2) %>%
    pivot_longer(starts_with("X."), names_to = "distance", values_to = "nu") %>%
    dplyr::filter(!is.na(nu)) %>%
    dplyr::mutate(distance = recode(distance, X.0 = "12.5", X.1 = "37.5", X.2 = "75")) %>%
    dplyr::mutate(distance = as.numeric(as.character(distance))) %>%
    rename(size = nu) %>%
    filter(size>0) %>%
    filter(!is.na(size))
  
  # #aggregate
  # detectionData <- detectionData %>%
  #                   group_by(kvadratnr, distance) %>%
  #                   summarise(size = sum(size))
  #  
  
  #add on environmental data
  detectionData <- inner_join(detectionData, environData)
  
  #fit model with covariates
  
  hr.model <- ds(detectionData, truncation=100,
                   transect="line",
                   formula = ~lines_pathroad + lines_forest,
                   cutpoints = c(0,25,50,100),
                   key = "hn")
  
  predGrid <- environData %>%
                select(kvadratnr, lines_pathroad, lines_forest) %>%
                add_column(distance=100) %>%
                add_column(Species = myspecies)
  
  predGrid$ESWpreds <- predict(hr.model, newdata = predGrid, esw=TRUE)$fitted
  
  return(predGrid)
  
}


comparedistanceModel <- function(myspecies, data){
  
  #and just one species
  dataS <- data %>% 
    filter(Species==myspecies) %>%
    select(-id)
  
  #and use the right season for the species
  dataS <- dataS %>%
    mutate(type = as.character(type)) %>%
    mutate(bestSeason=as.character(bestSeason)) %>%
    dplyr::filter(type==bestSeason)
  
  #subset transects for season of selected bird
  infoS <- info %>%
    filter(type == unique(dataS$bestSeason))
  
  #get all transects for the obs
  allData <- inner_join(infoS, dataS)
  
  #make long
  detectionData <- allData %>%
    dplyr::select(kvadratnr,Year, X.0,X.1,X.2) %>%
    pivot_longer(starts_with("X."),names_to = "distance", values_to = "nu") %>%
    dplyr::filter(!is.na(nu)) %>%
    dplyr::mutate(distance = recode(distance, X.0 = "12.5", X.1 = "37.5", X.2 = "75")) %>%
    dplyr::mutate(distance = as.numeric(as.character(distance))) %>%
    rename(size = nu) %>%
    filter(size>0)
  
  #add on environmental data
  detectionData <- inner_join(detectionData, environData)
  
  #different distance models
  hn.model <- ds(detectionData, truncation=100,
                 transect="line",
                 formula = ~ 1,
                 key = "hn",
                 cutpoints = c(0,25,50,100),
                 adjustment = NULL)
  
  #hr.model <- ds(detectionData, truncation=100,
  #               transect="line",
  #               formula = ~ 1,
  #               key = "hr",
  #               cutpoints = c(0,25,50,100),
  #               adjustment = NULL)
  
  unif.model <- ds(detectionData, truncation=100,
                 transect="line",
                 formula = ~ 1,
                 key = "unif",
                 cutpoints = c(0,25,50,100))
  
  return(data.frame(Species = myspecies,
                    Model = c("hn", "unif"),
                    AIC =  c(AIC(hn.model)$AIC,AIC(unif.model)$AIC
                    )))
}



  
# dsmModel <- function(myspecies, data, type="detection"){
#   
#   require(Distance)
#   require(dsm)
#   
#   #and just one species
#   dataS <- data %>% 
#     filter(Species==myspecies) %>%
#     select(-id)
#   
#   #and use the right season for the species
#   dataS <- dataS %>%
#     mutate(type = as.character(type)) %>%
#     mutate(bestSeason=as.character(bestSeason)) %>%
#     dplyr::filter(type==bestSeason)
#   
#   #subset transects for season of selected bird
#   infoS <- info %>%
#     filter(type == unique(dataS$bestSeason))
#   
#   #get all transects for the obs
#   allData <- inner_join(infoS, dataS)
#   
#   #make long
#   detectionData <- allData %>%
#     dplyr::select(kvadratnr,Year, X.0,X.1,X.2) %>%
#     pivot_longer(starts_with("X."),names_to = "distance", values_to = "nu") %>%
#     dplyr::filter(!is.na(nu)) %>%
#     dplyr::mutate(distance = recode(distance, X.0 = "12.5", X.1 = "37.5", X.2 = "75")) %>%
#     dplyr::mutate(distance = as.numeric(as.character(distance))) %>%
#     rename(size = nu) %>%
#     filter(size>0)
#   
#   #add on environmental data
#   detectionData <- inner_join(detectionData, environData)
#   
#   hr.model <- ds(detectionData, truncation=100,
#                  transect="line",
#                  formula = ~ 1,
#                  cutpoints = c(0,25,50,100),
#                  key = "hn", adjustment = NULL)
#   
#   if(type=="detection"){
#     
#     temp <- summary(hr.model)
#     return(data.frame(Species = myspecies, P = temp$ds$average.p))
#     
#   }
#   
#   if(type =="prediction" | type =="grid"){
#     
#     #make data
#     myObsdata = detectionData %>%
#       mutate(object = 1:nrow(.),
#              Sample.Label = paste(kvadratnr,Year)) %>%
#         filter(kvadratnr %in% environData$kvadratnr)
#     
#     mySegData <- infoS %>%
#       inner_join(.,environData) %>%
#       mutate(Sample.Label = paste(kvadratnr,Year),
#              Effort = 1000)%>%
#       filter(!duplicated(Sample.Label))
#     
#     mySegData<- addCoords(mySegData)
#     
#     # fit a simple smooth of x and y to counts
#     mod1 <- dsm(count ~ buffer_forest + buffer_agri_int +
#                   buffer_urban + buffer_agri_ext + buffer_seawater + 
#                   buffer_wet + buffer_freshwater + yday + sinceSunrise +
#                   skydaekke + regn + vind + s(X,Y), 
#                   hr.model, 
#                   mySegData, myObsdata)
#     
#     #predict to wider area
#     preddata <- dsmPredGrid(myspecies) %>%
#       add_column(Area = 1000 * 1000)
#     
#     if(type=="prediction"){
#       
#    temp <- dsm_var_gam(mod1, preddata, off.set = preddata$Area)
#    return(summary_dsm_var(temp))
#    
#     }else if(type=="grid"){
#       
#     preddata$preds <- predict(mod1, preddata, off.set = preddata$Area)
#     return(preddata)
#       
#     }
#   }
# }


### GBM #####

checkGBM <- function(myspecies, data){
  
  #http://uc-r.github.io/gbm_regression
  
  allData <- getSpeciesData(myspecies, data) 
  
  allDataS <- allData %>%
                dplyr::select(total, Year, yday,sinceSunrise, 
                              X, Y, starts_with("buffer")) %>%
                dplyr::select(-buffer_mapped, -buffer500_mapped, 
                              -buffer500_total,-buffer_total)

  #split into training and test
  train.index <- sample(c(1:nrow(allDataS)), 0.7 * nrow(allDataS), replace=F) 
  train <- allDataS[train.index,]
  test <- allDataS[-train.index,]
  
  #fit model
  gbm1 <- gbm(total ~ . , distribution = "poisson", 
                      data = train, n.trees = 2000, interaction.depth = 5, 
              n.minobsinnode = 10, shrinkage = 0.01, cv.folds = 5) 

  #learning rate/shrinkage: smaller values reduce the chance of overfitting 
  gbm.perf(gbm1, method = "cv")
  
  print(myspecies)
  summary(gbm1)
  
  #check predictions
  train$preds <- predict(gbm1, newdata = train, type="response")
  qplot(total, preds, data = train)#good!!
  test$preds <- predict(gbm1, newdata = test, type="response")
  qplot(total, preds, data = test)
  
  return(data.frame(Species = myspecies,
                    Cor_test = cor(test$total, test$preds),
                    Cor_train = cor(train$total, train$preds),
                    MAD_test = mean(abs(test$preds - test$total)),
                    MAD_train = mean(abs(train$preds - train$total))))

}


checkGBMBoost <- function(myspecies, data){
  
  #http://uc-r.github.io/gbm_regression
  
  allData <- getSpeciesData(myspecies, data) 
  
  allDataS <- allData %>%
    dplyr::select(total, Year, yday,sinceSunrise, 
                  X, Y, starts_with("buffer")) %>%
    dplyr::select(-buffer_mapped, -buffer500_mapped, 
                  -buffer500_total,-buffer_total)
  
  #split into training and test
  #sample along count gradient
  random_split <- rsample::initial_split(allDataS, prop = 0.70, strata = total)
  train <- training(random_split)
  test <- testing(random_split)
  
  #fit model
  gbm1 <- xgboost(data = as.matrix(train[,-1]), label = train$total, 
                  subsample = 0.5, # Lower ratios avoid over-fitting.
                  max.depth = 3, # Lower values avoid over-fitting.
                  min_child_weight = 20, #Larger values avoid over-fitting.
                  eta = 0.15, # Lower values avoid over-fitting.
                  nthread = 3, 
                  nrounds = 500, 
                  gamma = 10, # Larger values avoid over-fitting.
                  colsample_bytree = 0.5, # Lower ratios avoid over-fitting.
                  lambda = 2,
                  alpha = 1,
                  objective = "count:poisson", 
                  booster = 'gbtree')
  
  #check predictions
  train$preds <- predict(gbm1, newdata = as.matrix(train[,-1]), type="response")
  #qplot(total, preds, data = train)#good!!
  test$preds <- predict(gbm1, newdata = as.matrix(test[,-1]), type="response")
  #qplot(total, preds, data = test)
  
  (out <- data.frame(Species = myspecies,
                    Cor_test = cor(test$total, test$preds),
                    Cor_train = cor(train$total, train$preds),
                    MAD_test = mean(abs(test$preds - test$total)),
                    MAD_train = mean(abs(train$preds - train$total))))
  return(out)
  
}

gridsearchGBM <- function(myspecies, data){
  
  #http://uc-r.github.io/gbm_regression
  
  allData <- getSpeciesData(myspecies, data) 
  
  allDataS <- allData %>%
    dplyr::select(total, Year, yday, sinceSunrise, island,
                  X, Y, starts_with("buffer")) %>%
    dplyr::select(-buffer_mapped, -buffer500_mapped, 
                  -buffer500_total,-buffer_total) %>%
    mutate(island = factor(island))
  
  #grid search
  hyper_grid <- expand.grid(
    Species = myspecies,
    shrinkage = c(.005, .01, .05),
    interaction.depth = c(2, 3, 5),
    n.minobsinnode = c(5, 10, 15),
    bag.fraction = c(.5, .75, 1), 
    optimal_trees = 0,               # a place to dump results
    min_RMSE = 0                     # a place to dump results
  )
  
  # total number of combinations
  nrow(hyper_grid)
  
  for(i in 1:nrow(hyper_grid)){
    
  gbm.tune <- gbm(
    formula = total ~ .,
    distribution = "poisson",
    data = allDataS,
    n.trees = 2000,
    interaction.depth = hyper_grid$interaction.depth[i],
    shrinkage = hyper_grid$shrinkage[i],
    n.minobsinnode = hyper_grid$n.minobsinnode[i],
    bag.fraction = hyper_grid$bag.fraction[i],
    train.fraction = .75,
    n.cores = NULL, # will use all cores by default
    verbose = FALSE)

  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
  
  }
  
  return(hyper_grid)

}


#There are in general two ways that you can control overfitting in XGBoost:
#The first way is to directly control model complexity: max_depth, min_child_weight and gamma.
#The second way is to add randomness to make training robust to noise: subsample and colsample_bytree.
#You can also reduce stepsize eta. Remember to increase num_round when you do so.

# lower eta allows deeper interactions
# don't introduce gamma until you see a significant difference in your train and test error.

gridsearchGBMBoost <- function(myspecies, data){
  
  #http://uc-r.github.io/gbm_regression
  
  allData <- getSpeciesData(myspecies, data) 
  
  allDataS <- allData %>%
    dplyr::select(total, Year, yday, sinceSunrise, island,
                  X, Y, starts_with("buffer")) %>%
    dplyr::select(-buffer_mapped, -buffer500_mapped, 
                  -buffer500_total,-buffer_total) %>%
    mutate(island = factor(island))
  
  #grid search
  hyper_grid <- expand.grid(
    Species = myspecies,
    shrinkage = c(.005, .01, .05),
    interaction.depth = c(2, 3, 5),
    n.minobsinnode = c(5, 10, 15),
    bag.fraction = c(.5, .75, 1), 
    optimal_trees = 0,               # a place to dump results
    min_RMSE = 0                     # a place to dump results
  )
  
  # total number of combinations
  nrow(hyper_grid)
  
  for(i in 1:nrow(hyper_grid)){
    
    gbm.tune <- gbm(
      formula = total ~ .,
      distribution = "poisson",
      data = allDataS,
      n.trees = 2000,
      interaction.depth = hyper_grid$interaction.depth[i],
      shrinkage = hyper_grid$shrinkage[i],
      n.minobsinnode = hyper_grid$n.minobsinnode[i],
      bag.fraction = hyper_grid$bag.fraction[i],
      train.fraction = .75,
      n.cores = NULL, # will use all cores by default
      verbose = FALSE)
    
    # add min training error and trees to grid
    hyper_grid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
    hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
    
  }
  
  return(hyper_grid)
  
}

predGBM <- function(myspecies, data){
  
  #http://uc-r.github.io/gbm_regression
  
  allData <- getSpeciesData(myspecies, data) 
  
  allDataS <- allData %>%
    mutate(island = factor(island)) %>%
    dplyr::select(total, Year, yday,sinceSunrise, 
                  X, Y, island, starts_with("buffer")) %>%
    dplyr::select(-buffer_mapped, -buffer500_mapped, 
                  -buffer500_total,-buffer_total)

  #fit model
  #gbm1 <- gbm(total ~ . , distribution = "poisson", 
  #            data = allDataS, n.trees = 2000, interaction.depth = 5, 
  #            n.minobsinnode = 10, shrinkage = 0.01, cv.folds = 5) 
  
  #use gbm step instead
  gbm2 <- gbm.step(data=as.data.frame(allDataS), gbm.x = 2:71, gbm.y = 1, family = "poisson",
                              tree.complexity = 5, learning.rate = 0.005, step.size = 100,
                               bag.fraction = 0.5, n.folds = 5, n.trees = 5000)
  
  #get prediction dataframe for whole country
  bestSeason <- unique(allData$bestSeason[!is.na(allData$bestSeason)])
  gridDataS <- gridData %>%
                  mutate(island = factor(island)) %>%
                  dplyr::select(X, Y, island,starts_with("buffer")) %>%
                  dplyr::select(-buffer_mapped, -buffer500_mapped, 
                  -buffer500_total,-buffer_total) %>%
                  add_column(Year = 2016,
                             yday = median(info$yday[info$type==bestSeason]), 
                             sinceSunrise = median(info$sinceSunrise[info$type==bestSeason]))
  gridDataS <- gridDataS[,names(allDataS)[-1]]
  
  #add predictions
  gridDataS <- gridDataS %>%
                mutate(preds = predict(gbm2, newdata = gridDataS, type="response")) %>%
                add_column(Species = myspecies)
  
  return(gridDataS)
  #ggplot(gridDataS) +
  #  geom_point(aes(x=X, y=Y, colour=(preds)),size=rel(0.9))+
  #  scale_color_viridis_c(direction=-1, trans="log10")
  
}


predGBMBoost <- function(myspecies, data){
  
  #http://uc-r.github.io/gbm_regression
  
  allData <- getSpeciesData(myspecies, data) 
  
  allDataS <- allData %>%
    dplyr::select(total, Year, yday,sinceSunrise, 
                  X, Y, starts_with("buffer")) %>%
    dplyr::select(-buffer_mapped, -buffer500_mapped, 
                  -buffer500_total,-buffer_total)
  
  #model
  xgbm1 <- xgboost(data = as.matrix(allDataS[,-1]), label = allDataS$total, 
                  subsample = 0.5, # Lower ratios avoid over-fitting.
                  max.depth = 3, # Lower values avoid over-fitting.
                  min_child_weight = 20, #Larger values avoid over-fitting.
                  eta = 0.15, # Lower values avoid over-fitting.
                  nthread = 3, 
                  nrounds = 500, 
                  gamma = 10, # Larger values avoid over-fitting.
                  colsample_bytree = 0.5, # Lower ratios avoid over-fitting.
                  lambda = 2,
                  alpha = 1,
                  objective = "count:poisson", 
                  booster = 'gbtree')
  
  #get prediction dataframe for whole country
  bestSeason <- unique(allData$bestSeason[!is.na(allData$bestSeason)])
  
  gridDataS <- gridData %>%
    dplyr::select(X, Y, starts_with("buffer")) %>%
    dplyr::select(-buffer_mapped, -buffer500_mapped, 
                  -buffer500_total,-buffer_total) %>%
    add_column(Year = 2016,
               yday = median(info$yday[info$type==bestSeason]), 
               sinceSunrise = median(info$sinceSunrise[info$type==bestSeason]))
  #reorder to match
  gridDataS <- gridDataS[,names(allDataS)[-1]]
  
  #add predictions
  gridDataS <- gridDataS %>%
    mutate(preds = predict(xgbm1, newdata = as.matrix(gridDataS), type="response")) %>%
    add_column(Species = myspecies)
  
  return(gridDataS)
  #ggplot(gridDataS) +
  #  geom_point(aes(x=X, y=Y, colour=(preds)),size=rel(0.9))+
  #  scale_color_viridis_c(direction=-1, trans="log10")
  
}

### gam ####

identifyCovariates <- function(myspecies, data){
  
  allData <- getSpeciesData(myspecies, data)
  
  if(sum(allData$PA)>20){
    
    require(MuMIn)
    options(na.action = "na.fail")
    
    #full model - habitat vars
    gam1 <- glm(total ~ buffer_forest + buffer_agri_int + buffer_urban +
                  buffer_agri_ext + buffer_freshwater + buffer_open_urban +
                  buffer_wet + buffer_seawater,
                family=poisson, data=allData)#needs to be poisson
    
    dd <- dredge(gam1)
    temp <- subset(dd, delta < 2)
    vars <- names(temp)[colSums(temp,na.rm=T)!=0]
    out <- data.frame(Species = myspecies,
                      Vars = vars[grepl("buffer", vars)])
    
    #other vars
    gam1 <- glm(total ~ skydaekke + regn + vind + 
                  sinceSunrise + yday,
                family=poisson, data=allData)
    
    dd2 <- dredge(gam1)
    temp2 <- subset(dd2, delta < 2)
    vars2 <- names(temp2)[colSums(temp2,na.rm=T)!=0]
    out2 <- data.frame(Species = myspecies,
                      Vars = vars2[vars2 %in% c("skydaekke","regn","vind",
                                                  "sinceSunrise","yday")])
    
    out <- rbind(out,out2)
    
    
  }else{
    
    out <- data.frame(Species = myspecies,
                      Vars = "Not run",
                      Weights = 1)
  }
  return(out)
}


checkGAM <- function(myspecies, data){
  
  allData <-getSpeciesData(myspecies, data)
  
  #get formula
  impVars_sp <- impVars$Vars[impVars$Species==myspecies]
  spCovs <- paste(impVars_sp, collapse = "+")
  formula <- paste("total ~  s(X,Y, bs='ds') + factor(Year)", 
                   spCovs, sep="+")
  
  #add on season if more than one present
  if(length(unique(allData$type))>1){
    formula <- paste0(formula, "+type + type*yday")
  }
  
  #fit gamm  
  
  if(sum(allData$PA)>20){
    
    require(mgcv)
    gam1 <- gam(as.formula(formula),
                 family = quasipoisson,
                 #family=nb,
                 data = allData)
    
  allData$preds <- predict(gam1, type="response")
  
  out <- data.frame(Species = myspecies,
                    cor = cor(allData$total, allData$preds),
                    mad = mean(abs(allData$total-allData$preds)))
  
  }
  
  return(out)

}

dsmGAM <- function(myspecies, data, detDF = detDF){
  
  allData <-getSpeciesData(myspecies, data)
  
  #add on transect area and get effective area
  allData$P <- detDF$P[detDF$Species==myspecies]
  allData$EffectiveArea <- allData$P * 1000 * 1000
  
  #get formula
  impVars_sp <- impVars$Vars[impVars$Species==myspecies]
  spCovs <- paste(impVars_sp, collapse = "+")
  formula <- paste("total ~ offset(log(EffectiveArea)) + s(X,Y)", spCovs, sep="+")
  
  #add on season if more than one present
  if(length(unique(allData$type))>1){
    formula <- paste0(formula, "+type + type*yday")
  }
  
  #fit gamm  
  
  if(sum(allData$PA)>20){
    
    require(mgcv)
    #use zero-inflated
    #if(max(allData$total > 4)){
      
    #gam1 <- gam(list(total ~ factor(Year) + offset(log(EffectiveArea)) + te(X,Y)+buffer_agri_ext+buffer_agri_int+buffer_forest+buffer_freshwater+buffer_urban+skydaekke+vind+buffer_open_urban+buffer_wet+sinceSunrise+yday+regn+buffer_seawater,
    #                 ~ factor(Year) + te(X,Y) + buffer_agri_ext+buffer_agri_int+buffer_forest+buffer_freshwater+buffer_urban+skydaekke+vind+buffer_open_urban+buffer_wet+sinceSunrise+yday+regn+buffer_seawater),
    #            data = allData,
    #            family = ziplss())
    #}else{
      
    gam1 <- gamm(as.formula(formula),
                 random=list(Year=~1),
                 family = quasipoisson,
                 #family=nb,
                 data = allData)
    #}

  }else{
    gam1 <- NA
  }
  
  return(gam1)
  
}

modelCheck <- function(allPreds){
  
  #get mean per grid
  summaryPreds <- allPreds %>%
    group_by(Species, X, Y, kvadratnr) %>%
    summarise(meanPred = mean(preds),
              meantotal = mean(total)) %>%
    ungroup()
  
  #get correlation across species
  speciesCors <- summaryPreds %>%
    group_by(Species) %>%
    summarise(total = mean(meantotal),
              total_sd = sd(meantotal),
              mad = mean(abs(meanPred - meantotal)),
              corr = as.numeric(cor(meanPred, meantotal)))
  
  return(speciesCors)
  
}

plotCheck <- function(allPreds){
  
  #original predictions
  allPredsOriginal <- readRDS("outputs/gamdata_passerines_dsm_constant.rds")
  modelCors <- modelCheck(allPredsOriginal) 
  
  #new predictions
  modelCors2 <- modelCheck(allPreds) %>% rename(corr_new = corr, mad_new = mad)
  modelCors <- subset(modelCors, Species %in% modelCors2$Species)
  
  modelCors2$corr_old <- modelCors$corr  
  q1 <- qplot(corr_new, corr_old, data=modelCors2) + geom_abline(intercept=0, slope=1)
  
  modelCors2$mad_old <- modelCors$mad  
  q2 <- qplot(mad_new, mad_old, data=modelCors2) + geom_abline(intercept=0, slope=1)
  
  cowplot::plot_grid(q1,q2,nrow=1)
  
}

dsmPredGrid <- function(myspecies){
  
  allData <-getSpeciesData(myspecies, data)
  bestSeason <- unique(allData$bestSeason[!is.na(allData$bestSeason)])
  
  #get vars
  newdata = gridData %>%
      dplyr::select(c("buffer_forest","buffer_agri_int",
                      "buffer_urban","buffer_agri_ext",
                      "buffer_freshwater","buffer_wet",
                      "buffer_seawater","buffer_open_urban",
                      "X","Y","island")) %>%
      add_column(skydaekke = 1,
                 regn = 1,
                 vind = 1,
                 type = bestSeason,
                 sinceSunrise = median(allData$sinceSunrise[allData$type==bestSeason]),
                 yday = median(allData$sinceSunrise[allData$type==bestSeason])) %>%
      add_column(EffectiveArea = (1000 * 1000))
    
  #adjust effective area by land cover
  newdata$EffectiveArea <- newdata$EffectiveArea * gridData$LandCover

  #reduce to grids that we want to predict for this species?
  if(myspecies %in% excludeIslands$Species){

    speciesExclusions <- excludeIslands %>%
                            filter(Species == myspecies) %>%
                            pull("island")
    newdata <- newdata %>%
                  filter(!island %in% speciesExclusions)

  }
  
    return(newdata)
    
}

dsmPreds <- function(myspecies, myGAMs, type = "grid"){
  
  allData <-getSpeciesData(myspecies, data)
  
  #get model
  speciesGAM <- myGAMs[[which(myspecies == species)]]
  
  require(mgcv)
  
  if(type=="grid"){
    
    #get prediction grid
    newdata = dsmPredGrid(myspecies)
    
    #get predictions
    preds <- predict(speciesGAM$gam, newdata = newdata, 
                     re_formula = NA, se.fit = TRUE,
                     type="link")
    newdata$preds <- preds$fit
    newdata$preds.se <- preds$se.fit
    
  } else if (type=="data"){
    
    if(any(class(speciesGAM)=="list")){ #if quasi/nb gamm
    preds <- predict(speciesGAM$gam,  
                     re_formula = NA, se.fit = TRUE,
                     type="response")
    }else{ # if zero-inflated gam
    preds <- predict(speciesGAM,se.fit = TRUE,
                       type="response")
    }
    
    allData$preds <- preds$fit
    allData$preds.se <- preds$se.fit
    
    newdata <- allData %>%
      dplyr::select(X,Y,kvadratnr,preds,preds.se,total, Year)
  }
  
  #add species info
  newdata$Species <- myspecies
  
  return(newdata)
  
}

my_dsm_var_gam <- function(myspecies, myGAMs) {
  
  #define prediction grid
  pred.data = dsmPredGrid(myspecies)
  #reduce to grids that we want to predict for this species?
  
  # if we have a gamm, then just pull out the gam object
  dsm.obj <- myGAMs[[which(myspecies==species)]]
  if(any(class(dsm.obj) == "gamm")){
    dsm.obj <- dsm.obj$gam
    is.gamm <- TRUE
  }
  
  # if all the offsets are the same then we can just supply 1 and rep it
  off.set <- pred.data$EffectiveArea

  #get the transformation of the link function i.e, exp
  tmfn <- dsm.obj$family$linkinv
  
  # grab the coefficients
  cft <- coef(dsm.obj)
  preddo <- pred.data
  dpred.db <- matrix(0, 1, length(cft))
  
  ### fancy lp matrix stuff
  
  # set the offset to be zero here so we can use lp
  require(mgcv)
  pred.data$EffectiveArea <- rep(0, nrow(pred.data))
  
  lpmat <- predict(dsm.obj, newdata=pred.data, type='lpmatrix')
  lppred <- lpmat %*% cft
  preddo <-  off.set %*% tmfn(lppred)
  
  # NB previously this was done numerically but things didn't work
  # when se ~=0 and actual 0s were produced. Use analytical expression
  # for the log case here
  dpred.db <- off.set %*% (tmfn(lppred[,1]) * lpmat)
  
  # "'vpred' is the covariance of all the summary-things." - MVB
  # so we want the diagonals if length(pred.data)>1
  # A B A^tr
  vpred <- dpred.db %*% tcrossprod(vcov(dsm.obj), dpred.db)
  
  if(is.matrix(vpred)){
    vpred <- diag(vpred)
  }
  
  result <- list(pred   = preddo,
                 pred.var = vpred)

  return(result)
}


my_summary_dsm_var <- function(myspecies, object, detDF){
  
  # storage
  sinfo <- list()
  # save the alpha value for cis
  sinfo$alpha <- 0.05
  
  
  ### analytical variance estimation (varprop and gam results)
  sinfo$se <- sqrt(object$pred.var)
  sinfo$pred.est <- object$pred
  
  #detecion error
  average.p.se <- detDF$P_se[detDF$Species==myspecies]
  average.p <- detDF$P[detDF$Species==myspecies]
  sinfo$detfct.cv <- c()
  cvp.sq <- 0
  this_cvp.sq <- (average.p.se/average.p)^2
  cvp.sq <- cvp.sq + this_cvp.sq
  sinfo$detfct.cv <- c(sinfo$detfct.cv, sqrt(this_cvp.sq))
  
  #gam error
  sinfo$gam.cv <- sinfo$se/sinfo$pred.est
  
  #combine the two sources of error
  sinfo$cv <- sqrt(cvp.sq+sinfo$gam.cv^2)
  sinfo$se <- sinfo$cv*sinfo$pred.est
  
  #get 95% CI
  # this doesn't transform N, only se(N)
  # this probably should do the following:
  # lower =
  #  qlnorm(alpha/2, log(x$pred.est) - 0.5*log(x$cv^2+1),
  #         sqrt(log(x$cv^2+1)))
  # upper =
  #  qlnorm(1-alpha/2, log(x$pred.est) - 0.5*log(x$cv^2+1),
  #         sqrt(log(x$cv^2+1)))
  unconditional.cv.square <- sinfo$cv^2
  asymp.ci.c.term <- exp(qnorm(1-sinfo$alpha/2) *
                           sqrt(log(1+unconditional.cv.square)))
  asymp.tot <- c(sinfo$pred.est / asymp.ci.c.term,
                 sinfo$pred.est,
                 sinfo$pred.est * asymp.ci.c.term)
  
  names(asymp.tot) <- c(paste0(sinfo$alpha/2*100, "%"),
                        "Mean",
                        paste0((1-sinfo$alpha/2)*100,"%"))
  
  return(data.frame(total = asymp.tot[2],
                    lowerCI = asymp.tot[1],
                    upperCI = asymp.tot[3]))
  
}


dsmPopSize <- function(myspecies, myGAMS, detDF){
  
  var.object1 <- my_dsm_var_gam(myspecies, myGAMs) 
  my_summary_dsm_var(myspecies, var.object1, detDF) %>%
      add_column(Species = myspecies)
  
}

dsm_var_gam <- function(dsm.obj, pred.data, off.set){
  
  # if we have a gamm, then just pull out the gam object
  if(any(class(dsm.obj) == "gamm")){
    dsm.obj <- dsm.obj$gam
    is.gamm <- TRUE
  }
  
  
  # if all the offsets are the same then we can just supply 1 and rep it
  if(length(off.set) == 1){
      off.set <- rep(off.set, nrow(pred.data))
  }
  
  
  # depending on whether we have response or link scale predictions...
  tmfn <- dsm.obj$family$linkinv
 
  
  # grab the coefficients
  cft <- coef(dsm.obj)
  preddo <- list(length(pred.data))
  dpred.db <- matrix(0, length(pred.data), length(cft))
  
  ### fancy lp matrix stuff
  
  # set the offset to be zero here so we can use lp
  pred.data$off.set <- rep(0, nrow(pred.data))
  
  lpmat <- predict(dsm.obj, newdata=pred.data, type='lpmatrix')
  lppred <- lpmat %*% cft
  preddo <-  off.set %*% tmfn(lppred)
  
  # NB previously this was done numerically but things didn't work
  # when se ~=0 and actual 0s were produced. Use analytical expression
  # for the log case here
  dpred.db <- off.set %*% (tmfn(lppred[,1]) * lpmat)
  
  # "'vpred' is the covariance of all the summary-things." - MVB
  # so we want the diagonals if length(pred.data)>1
  # A B A^tr
  vpred <- dpred.db %*% tcrossprod(vcov(dsm.obj), dpred.db)
  
  if(is.matrix(vpred)){
    vpred <- diag(vpred)
  }
  
  result <- list(pred.var       = vpred,
                 pred           = preddo,
                 pred.data      = pred.data,
                 off.set        = off.set,
                 dsm.object     = dsm.obj)
  
  
  return(result)
}

summary_dsm_var <- function(object){
  
  # storage
  sinfo <- list()
  
  # save the alpha value for cis
  sinfo$alpha <- 0.05

  #detection bit
  ddf.summary <- summary(object$dsm.object$ddf)
  sinfo$detfct.cv <- c()
  cvp.sq <- 0
  this_cvp.sq <- (ddf.summary$average.p.se/
                            ddf.summary$average.p)^2
  cvp.sq <- cvp.sq + this_cvp.sq
  sinfo$detfct.cv <- c(sinfo$detfct.cv, sqrt(this_cvp.sq))

  #gam stuff
  sinfo$pred.est <- object$pred
  sinfo$se <- sqrt(object$pred.var)
  sinfo$gam.cv <- sinfo$se/sinfo$pred.est
  sinfo$cv <- sqrt(cvp.sq+sinfo$gam.cv^2)
      
  # total se
  sinfo$se <- sinfo$cv*sinfo$pred.est
  
  #summarise  
  x <- sinfo
  unconditional.cv.square <- x$cv^2
  asymp.ci.c.term <- exp(qnorm(1-x$alpha/2) *
                             sqrt(log(1+unconditional.cv.square)))
  asymp.tot <- c(x$pred.est / asymp.ci.c.term,
                   x$pred.est,
                   x$pred.est * asymp.ci.c.term)
    
  names(asymp.tot) <- c(paste0(x$alpha/2*100, "%"),
                          "Mean",
                          paste0((1-x$alpha/2)*100,"%"))
    
  return(asymp.tot)
    
}

### end ########################