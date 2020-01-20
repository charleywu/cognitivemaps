#Mixed effects models
#All models should be specified in the R notebooks. This code is mainly to facilitate a more convenient workflow
#Charley Wu, 2020
#house keeping
rm(list=ls())

#load packages
packages <- c('dplyr','cowplot', 'Rmisc', 'ggbeeswarm', 'brms', 'sjPlot', 'BayesFactor','scales',  'plyr', 'reshape2', 'ggridges', 'ggplot2', 'jsonlite', 'MASS', 'gridExtra', 'Hmisc', 'lsr', 'pander')
lapply(packages, require, character.only = TRUE)

theme_set(theme_cowplot(font_size=12))
source('dataProcessing.R')
source('statisticalTests.R') 

#Wrapper for brm models such that it saves the full model the first time it is run, otherwise it loads it from disk
run_model <- function(expr, modelName, path='brmsModels', reuse = TRUE) {
  path <- paste0(path,'/', modelName, ".brm")
  if (reuse) {
    fit <- suppressWarnings(try(readRDS(path), silent = TRUE))
  }
  if (is(fit, "try-error")) {
    fit <- eval(expr)
    saveRDS(fit, file = path)
  }
  fit
}
#####################################################################################################################################################################
#Load experiment data
#####################################################################################################################################################################
dataDir <- 'experimentData/full.csv'
df <- dataImport(dataFile =dataDir ,normalize=F)
trajDF <- importTrajData(dataFile =dataDir,normalize=F)
#Trim last rounds
df<- subset(df, round<10)

n_rounds = 9 #without bonus
n_trials = 20

#####################################################################################################################################################################
#Mixed Effects Modeling
#####################################################################################################################################################################
#Previous reward value and distance between selections
distanceRewardMM <- run_model(brm(distance ~ previousReward+context+previousReward*context +(1+previousReward|id), data=subset(df, !is.na(df$distance)), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99)), modelName = 'distanceRewardMM')
bayes_R2(distanceRewardMM)
fixedTerms <- fixef(distanceRewardMM)#Look at fixed terms

#Now generate predictions, removing id as a random effect
xseq <- seq(1,100)
newdat <-data.frame(context = rep(c("Conceptual","Spatial"), each=100), previousReward = rep(xseq,2))
preds <- fitted(distanceRewardMM, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
#create new fixed effects dataframe
fixedDF <- data.frame(context = rep(c("Conceptual","Spatial"), each=100), previousReward = rep(xseq,2),
                      distance = preds[,1], lower = preds[,3], upper = preds[,4] )

ggplot(subset(df, !is.na(df$distance)), aes(previousReward, distance, color = context, fill  = context)) +
  #geom_hline(yintercept = mean(randomDistanceDF$distance, na.rm=T ), size = 1, color = 'black', linetype='dashed')+ 
  geom_line(data = fixedDF,  size = 1)+ #GP is
  geom_ribbon(data = fixedDF, aes(ymin=lower, ymax = upper), color = NA, alpha = 0.4 )+
  stat_summary(fun.y=mean,geom='point', alpha = 0.5)+
  #geom_abline(slope = 1, linetype = 'dashed')+
  #coord_cartesian(xlim = c(0,100))+
  xlim(c(0,100))+
  theme_classic()+
  scale_color_brewer(palette = 'Dark2', name="Task")+
  scale_fill_brewer( palette = 'Dark2', name="Task")+
  #facet_grid(~context, labeller = as_labeller(contextLabels) )+
  xlab("Previous Reward Value")+
  ylab("Distance Between Selections")+
  theme(legend.position=c(0, 0), legend.justification=c(0,0), strip.background=element_blank(), legend.key=element_blank(), legend.background=element_blank())


distanceInitialMM <- run_model(brm(movement ~ previousReward+context+previousReward*context +(1+previousReward|id), data=subset(df, !is.na(df$movement)), cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99)), modelName = 'distanceInitialMM')
bayes_R2(distanceInitialMM)
tab_model(distanceInitialMM)
fixedTerms <- fixef(distanceInitialMM)#Look at fixed terms

#Now generate predictions, removing id as a random effect
xseq <- seq(1,100)
newdat <-data.frame(context = rep(c("Conceptual","Spatial"), each=100), previousReward = rep(xseq,2))
preds <- fitted(distanceInitialMM, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
#create new fixed effects dataframe
fixedDF <- data.frame(context = rep(c("Conceptual","Spatial"), each=100), previousReward = rep(xseq,2),
                      movement = preds[,1], lower = preds[,3], upper = preds[,4] )

ggplot(subset(df, !is.na(df$movement)), aes(previousReward, movement, color = context, fill  = context)) +
  #geom_hline(yintercept = mean(randomDistanceDF$distance, na.rm=T ), size = 1, color = 'black', linetype='dashed')+ 
  geom_line(data = fixedDF,  size = 1)+ #GP is
  geom_ribbon(data = fixedDF, aes(ymin=lower, ymax = upper), color = NA, alpha = 0.4 )+
  stat_summary(fun.y=mean,geom='point', alpha = 0.5)+
  #geom_abline(slope = 1, linetype = 'dashed')+
  #coord_cartesian(xlim = c(0,100))+
  xlim(c(0,100))+
  theme_classic()+
  scale_color_brewer(palette = 'Dark2', name="Task")+
  scale_fill_brewer( palette = 'Dark2', name="Task")+
  #facet_grid(~context, labeller = as_labeller(contextLabels) )+
  xlab("Previous Reward Value")+
  ylab("Distance Between Selections")+
  theme(legend.position=c(0, 0), legend.justification=c(0,0), strip.background=element_blank(), legend.key=element_blank(), legend.background=element_blank())


tab_model(distanceRewardMM, distanceInitialMM)
