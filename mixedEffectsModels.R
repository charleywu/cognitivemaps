#Mixed effects models
#All models are also specified in the R notebooks. This code is mainly to facilitate a more convenient workflow for running multiple models
#Charley Wu, 2020
#house keeping
rm(list=ls())

#load packages
packages <- c('tidyr', 'dplyr','cowplot', 'Rmisc', 'ggbeeswarm', 'brms', 'sjPlot', 'BayesFactor','scales',  'plyr',  'ggplot2', 'jsonlite', 'MASS',  'Hmisc', 'lsr')
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
#Previous reward and distance
#####################################################################################################################################################################
#Previous reward value and distance between selections
prior <- c(set_prior("normal(0,1)", class = "b"),set_prior("normal(0,1)", class = "sd"))
distanceRewardMM <- run_model(brm(distance ~ 0+ intercept+ previousReward+context+previousReward*context +(1+previousReward*context|id), data=subset(df, !is.na(df$distance)), 
                                  prior = prior,cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99)), modelName = 'distanceRewardMM')
#bayes_R2(distanceRewardMM)  R2 = 0.4170613
fixedTerms <- fixef(distanceRewardMM)#Look at fixed terms, reporting medians

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
  stat_summary(fun.y=mean, geom='point', alpha = 0.5)+
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


#####################################################################################################################################################################
#Previous reward and trajectory length
#####################################################################################################################################################################
prior <- c(set_prior("normal(0,1)", class = "b"),set_prior("normal(0,1)", class = "sd"))
distanceInitialMM <- run_model(brm(movement ~0 + intercept +previousReward*context +(1+previousReward*contex|id), data=subset(df, !is.na(df$movement)), 
                                   prior = prior,
                                   cores=4,  iter = 4000, warmup = 1000, control = list(adapt_delta = 0.99)), modelName = 'distanceInitialMM')
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


##########################################################################################################################################
#Bonus round
##########################################################################################################################################
bmtCol <- "#F0E442"
gpCol <- "#D55E00"

bonusDF <- readRDS( "experimentData/bonusRound.Rds")
bonusRoundResults <- readRDS("modelResults/bonusRoundResults.Rds")

predictionDF <- data.frame(id= rep(bonusDF$id, 2), humanPrediction = rep(bonusDF$meanEstimate, 2), 
                           humanConfidence = rep(bonusDF$howCertain, 2), context = rep(bonusDF$context, 2),
                           environment = rep(bonusDF$environment, 2), modelPrediction = c(bonusDF$GPmean, bonusDF$BMTmean), 
                           modelUncertainty = c(bonusDF$GPuncertainty, bonusDF$BMTuncertainty),
                           model = rep(c('GP', 'BMT'),each = nrow(bonusDF) ))

predictionDF$model <- factor(predictionDF$model, level <- c('GP', 'BMT')) 


#Mixed effects models of predictions
prior <- c(set_prior("normal(0,1)", class = "b"),set_prior("normal(0,1)", class = "sd"))
GPjudgments <- run_model(brm(modelPrediction ~ 0+ intercept+ humanPrediction+context+humanPrediction*context +(1+humanPrediction*context|id), data=subset(predictionDF, model=='GP'), 
                             prior = prior, cores=4,  iter = 4000, warmup = 1000,  control = list(adapt_delta = 0.99)), modelName = 'GPbonusJudgments')
bayes_R2(GPjudgments) #R2 = 
tab_model(GPjudgments)
fixedTerms <- fixef(GPjudgments)#Look at fixed terms

#Now generate predictions, removing id as a random effect
xseq <- seq(1,100)
newdat <-data.frame(context = rep(c("Conceptual","Spatial"), each=100), humanPrediction = rep(xseq,2))
preds <- fitted(GPjudgments, re_formula = NA, newdata = newdat, probs = c(0.025, 0.975))
#create new fixed effects dataframe
fixedDF <- data.frame(model='GP', context = rep(c("Conceptual","Spatial"), each=100), humanPrediction = rep(xseq,2),
                      modelPrediction = preds[,1], lower = preds[,3], upper = preds[,4] )
#create BMT df giving the same prediction for all options
bmtDF <- data.frame(model='BMT', context = rep(c("Conceptual","Spatial"), each=100), humanPrediction = rep(xseq,2), modelPrediction = 50, lower = NA, upper = NA)
fixedDF <- rbind(fixedDF, bmtDF)

contextLabels <- contextLabels <- c('Conceptual' = 'Conceptual Task', 'Spatial' = 'Spatial Task')

ggplot(subset(predictionDF, model=='GP'), aes(humanPrediction, modelPrediction, color = model, fill  = model)) +
  geom_point(alpha =.2, color = 'black', fill=NA) +
  #geom_hline(yintercept = 50, size = 1, color = bmtCol)+ #BMT is a flat line
  geom_line(data = fixedDF,  size = 1)+ #GP is
  geom_ribbon(data = fixedDF, aes(ymin=lower, ymax = upper), color = NA, alpha = 0.4 )+
  #geom_abline(slope = 1, linetype = 'dashed')+
  coord_cartesian(xlim = c(0,100), ylim=c(0,100))+
  theme_classic()+
  scale_fill_manual(values = c(gpCol,bmtCol), name='')+
  scale_color_manual(values = c(gpCol,bmtCol), name='')+
  facet_grid(~context, labeller = as_labeller(contextLabels) )+
  xlab("Participant Estimate")+
  ylab("Model Estimate")+
  theme(legend.position=c(0, 1.1), legend.justification=c(0,1), strip.background=element_blank(), legend.key=element_blank(), legend.background=element_blank())



#Model uncertainty to participant confidence
#Mixed effects model
predictionDF$modelUncertainty
prior <- c(set_prior("normal(0,1)", class = "b"),set_prior("normal(0,1)", class = "sd"))
GPconf <- run_model(brm(modelUncertainty~  0 + intercept + humanConfidence*context +(1+humanConfidence*context|id), data=subset(predictionDF, model=='GP'),
                        prior = prior, cores=4,  iter = 4000, warmup = 1000,  control = list(adapt_delta = 0.99)), modelName = 'GPbonusconf')
#bayes_R2(GPconf) #R2 = 
#fixef(GPconf)
#tab_model(GPconf)

#Compute rank-ordered confidence for plot
confidenceDF <- data.frame()
for (pid in unique(bonusDF$id)){
  for (task in c('Conceptual', 'Spatial')){
    for (model in c('GP')){
      dsub <- subset(bonusDF, id == pid & context ==  task)
      modelUncertainty = paste0(model,'uncertainty')
      dummy <- data.frame(model = model,id=pid, rankParticipantConfidence= rank(dsub$howCertain), context = task, rankModelUncertainty = rank(dsub[,modelUncertainty]) )
      confidenceDF <- rbind(dummy,confidenceDF)
    }  
  }
}
confidenceDF$context <- factor(confidenceDF$context , levels = c("Conceptual", "Spatial"))

p5<- ggplot(confidenceDF, aes(x=rankParticipantConfidence, y = rankModelUncertainty,  color= model))+
  #geom_quasirandom(varwidth = T, size = 0.5, alpha = 0.2) +
  #geom_boxplot(width = .25,  outlier.shape = NA, alpha = 0.5, color = 'black') +
  stat_summary(fun.y = mean, geom = "point", color = 'black') + 
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", color = 'black', width = 0)+
  geom_smooth(fill = NA, method = 'lm',formula=y~x,  se=FALSE, fullrange=TRUE)+ 
  #scale_fill_brewer(palette = "Spectral", name="Confidence") +
  #scale_color_viridis(discrete=TRUE, direction = -1)+
  #scale_fill_viridis(discrete=TRUE, direction = -1)+
  #scale_y_continuous(limits = c(0,35))+
  facet_grid(~context, labeller = as_labeller(contextLabels))+
  #scale_x_continuous(limits = c(0.5,11.5), breaks =c(1,3,5,7,9,11))+
  #ylab(expression(paste("GP Uncertainty (", sigma,")")))+
  ylab(expression(paste("GP Uncertainty (rank order)")))+
  xlab('Participant Confidence (rank order)')+
  scale_color_manual(values = c(gpCol,bmtCol), name='')+
  theme_classic() + theme(strip.background = element_blank(), legend.position=c(1,1.1), legend.justification=c(1,1))
p5

tab_model(GPjudgments, GPconf)

