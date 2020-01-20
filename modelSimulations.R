#Charley Wu 2019

#Script to run rational models in comparison to human behavior

#############################################################################################################################
# IMPORT DATA AND ENVIRONMENTS
#############################################################################################################################
#house keeping
rm(list=ls())

#load packages
packages <- c('fgpt', 'plyr', 'jsonlite', 'ggplot2', 'gridExtra', 'reshape2', 'stargazer', 'coefplot', "grid", 'matrixcalc', 'parallel')
lapply(packages, require, character.only = TRUE)==TRUE

source("models.R")
#Participant data
source('dataProcessing.R')
d <- dataImport() #participant data

reps <- 10000 #replications
#cores <- detectCores()-1
cores <- 8

#Environments
setwd("..") #Set into parent folder to read kernel files
setwd("experiment")
#load environments from json, unlist, transform into numeric, convert into matrix, and name dimensions
roughEnvironments <- lapply(fromJSON("roughEnvironment.json"), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=3, byrow=TRUE, dimnames=list(seq(1,64), c('x1', 'x2', 'y'))))
smoothEnvironments <- lapply(fromJSON("smoothEnvironment.json"), FUN=function(x) matrix(as.numeric(unlist(x)), ncol=3, byrow=TRUE, dimnames=list(seq(1,64),  c('x1', 'x2', 'y'))))
setwd("..")#Step back into analysis folder


#parameter estimates
bmtPars <-  read.csv('rationalModels/parameters/BMT.csv')
gpPars <- read.csv('rationalModels/parameters/gpucb.csv')
simKernelPars <- read.csv('rationalModels/parameters/simKernelucb.csv')

#############################################################################################################################
# RANDOM MODEL
# to read the previous simulation from disk:
# randomDF <- read.csv("rationalModels/random.csv")
#############################################################################################################################

randomModel <- function(replications, outputfile){
  #Run for both types of environments
  smoothReward <- mcmapply(1:replications, FUN=function(x) sample(smoothEnvironments[[sample(1:40,1)]][,'y'], 20, replace = TRUE)*100, mc.preschedule=TRUE, mc.cores=cores)
  roughReward <- mcmapply(1:replications, FUN=function(x) sample(roughEnvironments[[sample(1:40,1)]][,'y'], 20, replace = TRUE)*100, mc.preschedule=TRUE, mc.cores=cores)
  #put into dataFrame
  randomDF <- data.frame(trial=rep(seq(1:20), 4),
                         Environment=rep(rep(c("Smooth", "Rough"), each=20),2),
                         Context = c(rep("Spatial", 40), rep("Conceptual", 40)),
                         meanReward = rep(c(rowMeans(smoothReward), rowMeans(roughReward)),2), #reward averaged over replications
                         meanSE = rep(c(apply(smoothReward, 1, FUN = function(x) sd(x)/sqrt(length(x))),
                                        apply(roughReward, 1, FUN = function(x) sd(x)/sqrt(length(x)))),2))
  #add model label
  randomDF$Model <- rep("Random", 80)
  #write output
  if(!is.null(outputfile)){
    write.csv(randomDF, outputfile) 
  }
  return(randomDF)
}

randomDF <- randomModel(reps, "rationalModels/random.csv")

#############################################################################################################################
# BMT-UCB Model
# to read the previous simulation from disk:
# bmtDF <- read.csv("rationalModels/BMTUCB.csv")
#############################################################################################################################

bmtRationalModel <- function(replications, modelName, outputfile, parameters, acq=ucb){
  smoothTotal <- list()
  roughTotal <- list()
  for (task in c('Conceptual', 'Spatial')){
    smoothPars <- subset(parameters, context == task & environment=="Smooth")
    roughPars <- subset(parameters, context == task & environment =="Rough")
    #choice matrix
    choices <- expand.grid(0:7, 0:7) #build choice matrix
    names(choices)<-c("x1","x2")
    #run for smooth environments
    smoothReward <- mcmapply(1:replications, FUN=function(x){
      #sample parameters
      params <- smoothPars[sample(1:nrow(smoothPars), 1),]
      kError <- params$kError
      beta <- params$beta
      tau <- params$tau
      #tau <- median(smoothPars$tau)
      #randomly choose environment
      envNum <- sample(1:40,1) 
      #1st trial is random
      location <- sample(1:64,1)#each location as an integer from 1:25; first value is random, so treat both environments as replications of the same actions
      #Observations
      X1 <- choices[location,'x1']
      X2 <- choices[location,'x2']
      #reward
      reward<- c()
      reward[1] <- Y <- (smoothEnvironments[[envNum]][location,"y"])*100
      #posterior
      prevPost <- NULL  #initialize previous posterior for sequential updating of BMT posterior
      for (j in 2:20){ #after that, loop through remaining trials and make decisions based on BMT-UCB preditions
        #update posterior predictions
        post <- bayesianMeanTracker(x = as.matrix(cbind(X1,X2))[j-1,], y = (reward[j-1]-50)/100, prevPost = prevPost, theta = c(kError))
        prevPost <- post  #save new posterior as prevPost for next round
        #compute acquisition function evaluation
        utilityVec <- acq(post, pars = c(beta))
        #to prevent overflow, subtract max of q(x) vector 
        utilityVec <- utilityVec - max(utilityVec)
        #compute softmax choice probabilities
        p <- exp(utilityVec/tau)
        p <- p/sum(p)
        #Sample next choice
        location <- sample(1:64,1, prob=p, replace=TRUE)
        #update reward, X1, X2, and Y for both smooth and rough
        reward[j] <- (smoothEnvironments[[envNum]][location,"y"]) * 100
        X1 <- c(X1, choices[location, 'x1'])
        X2 <- c(X2, choices[location, 'x2'])
        Y <- c(Y,  reward[j])}
      reward}, mc.preschedule = TRUE, mc.cores=cores)
    #run for smooth environments
    roughReward <- mcmapply(1:replications, FUN=function(x){
      #sample parameters
      params <- roughPars[sample(1:nrow(roughPars), 1),]
      kError <- params$kError
      beta <- params$beta
      tau <- params$tau
      #tau <- median(roughPars$tau)
      #randomly choose environment
      envNum <- sample(1:40,1) 
      #1st trial is random
      location <- sample(1:64,1)#each location as an integer from 1:25; first value is random, so treat both environments as replications of the same actions
      #Observations
      X1 <- choices[location,'x1']
      X2 <- choices[location,'x2']
      #reward
      reward<- c()
      reward[1] <- Y <- (roughEnvironments[[envNum]][location,"y"])*100
      #posterior
      prevPost <- NULL  #initialize previous posterior for sequential updating of BMT posterior
      for (j in 2:20){ #after that, loop through remaining trials and make decisions based on GP preditions
        #update posterior predictions
        post <- bayesianMeanTracker(x = as.matrix(cbind(X1,X2))[j-1,], y = (reward[j-1]-50)/100, prevPost = prevPost, theta = c(kError))
        prevPost <- post  #save new posterior as prevPost for next round
        #compute acquisition function evaluation
        utilityVec <- acq(post, pars = c(beta))
        #prevent overflow by subtracting max
        utilityVec <- utilityVec - max(utilityVec)
        #compute softmax choice probabilities
        p <- exp(utilityVec/tau)
        p <- p/sum(p)
        #Sample next choice
        location <- sample(1:64,1, prob=p, replace=TRUE)
        #update reward, X1, X2, and Y for both smooth and rough
        reward[j] <- (roughEnvironments[[envNum]][location,"y"])* 100
        X1 <- c(X1, choices[location, 'x1'])
        X2 <- c(X2, choices[location, 'x2'])
        Y <- c(Y,  reward[j])}
      reward}, mc.preschedule = TRUE, mc.cores=cores)
    smoothTotal[[length(smoothTotal)+1]] <- smoothReward
    roughTotal[[length(roughTotal)+1]] <-  roughReward
  }
  smoothTotal <- do.call("rbind", smoothTotal)
  roughTotal <- do.call("rbind", roughTotal)
  #put into dataFrame
  bmtDF <-data.frame(trial=rep(seq(1:20), 4),
                     Environment=c(rep("Smooth", 40), rep("Rough", 40)),
                     Context = rep(rep(c("Conceptual", "Spatial"), each=20),2), 
                     meanReward = c(rowMeans(smoothTotal), rowMeans(roughTotal)), #reward averaged over replications
                     meanSE = c(apply(smoothTotal, 1, FUN = function(x) sd(x)/sqrt(length(x))), 
                                apply(roughTotal, 1, FUN = function(x) sd(x)/sqrt(length(x)))))
  #add model label
  bmtDF$Model <- rep(modelName, 80)
  #write to csv
  if(!is.null(outputfile)){
    write.csv(bmtDF, outputfile)
  }
  return(bmtDF)
}

bmtDF <- bmtRationalModel(reps, modelName = 'BMT-UCB', outputfile = "rationalModels/BMTUCB.csv", parameters = bmtPars)
#bmtMean <- bmtRationalModel(reps, modelName = 'BMT-GM', outputfile = "rationalModels/BMTGM.csv", parameters = bmtMeanPars)

#############################################################################################################################
# GP Model
# to read the previous simulation from disk:
# gpDF <- read.csv("rationalModels/GPUCB.csv")
#############################################################################################################################

gpRationalModel <- function(replications, modelName, outputfile, parameters, acq=ucb, kernel=rbf, cores){
  smoothTotal <- list()
  roughTotal <- list()
  for (task in c('Conceptual', 'Spatial')){
    #Extract parameters
    smoothPars <- subset(parameters, context == task & environment=="Smooth")
    roughPars <- subset(parameters, context == task & environment =="Rough")
    #choice matrix
    choices <- expand.grid(0:7, 0:7) #build choice matrix
    names(choices)<-c("x1","x2")
    #run for smooth environments
    smoothReward<- mcmapply(1:replications, FUN=function(x){
      #sample parameters
      params <- smoothPars[sample(1:nrow(smoothPars), 1),]
      lambda <- params$lambda
      beta <- params$beta
      tau <- params$tau
      #tau <- median(smoothPars$tau)
      #randomly choose environment
      envNum <- sample(1:40,1) 
      #1st trial is random
      location <- sample(1:64,1)#each location as an integer from 1:25; first value is random, so treat both environments as replications of the same actions
      #Observations
      X1 <- choices[location,'x1']
      X2 <- choices[location,'x2']
      #reward
      reward<- c()
      reward[1] <- Y <- (smoothEnvironments[[envNum]][location,"y"])*100
      for (j in 2:20){ #after that, loop through remaining trials and make decisions based on GP preditions
        #compute posterior predictions
        post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = (Y-50)/100, k = kernel) #scale observed Y to zero mean, variance of 1
        #compute acquisition function evaluation
        utilityVec <- acq(post, pars = c(beta))
        #scale to max of prevent overflow by subtracting max
        utilityVec <- utilityVec - max(utilityVec)
        #compute softmax choice probabilities
        p <- exp(utilityVec/tau)
        p <- p/sum(p)
        #Sample next choice
        location <- sample(1:64,1, prob=p, replace=TRUE)
        #update reward, X1, X2, and Y for both smooth and rough
        reward[j] <- (smoothEnvironments[[envNum]][location,"y"])* 100
        X1 <- c(X1, choices[location, 'x1'])
        X2 <- c(X2, choices[location, 'x2'])
        Y <- c(Y,  reward[j])}
      reward}, mc.preschedule = TRUE, mc.cores=cores)
    #run for rough environments
    roughReward <- mcmapply(1:replications, FUN=function(x){
      #sample parameters
      params <- roughPars[sample(1:nrow(roughPars), 1),]
      lambda <- params$lambda
      beta <- params$beta
      tau <- params$tau
      #tau <- median(roughPars$tau)
      #randomly choose environment
      envNum <- sample(1:40,1) 
      #1st trial is random
      location <- sample(1:64,1)#each location as an integer from 1:25; first value is random, so treat both environments as replications of the same actions
      #Observations
      X1 <- choices[location,'x1']
      X2 <- choices[location,'x2']
      #reward
      reward<- c()
      reward[1] <- Y <- (roughEnvironments[[envNum]][location,"y"])*100
      for (j in 2:20){ #after that, loop through remaining trials and make decisions based on GP preditions
        #compute posterior predictions
        post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001), X = cbind(X1,X2), Y = (Y-50)/100, k = kernel) #scale observed Y to zero mean, variance of 1
        #compute acquisition function evaluation
        utilityVec <- acq(post, pars = c(beta))
        #scale to max of prevent overflow by subtracting max
        utilityVec <- utilityVec - max(utilityVec)
        #compute softmax choice probabilities
        p <- exp(utilityVec/tau)
        p <- p/sum(p)
        #Sample next choice
        location <- sample(1:64,1, prob=p, replace=TRUE)
        #update reward, X1, X2, and Y for both smooth and rough
        reward[j] <- (roughEnvironments[[envNum]][location,"y"]) * 100
        X1 <- c(X1, choices[location, 'x1'])
        X2 <- c(X2, choices[location, 'x2'])
        Y <- c(Y,  reward[j])}
      reward}, mc.preschedule = TRUE, mc.cores=cores)
    smoothTotal[[length(smoothTotal)+1]] <- smoothReward
    roughTotal[[length(roughTotal)+1]] <-  roughReward
  }
  #put it all together
  smoothTotal <- do.call("rbind", smoothTotal)
  roughTotal <- do.call("rbind", roughTotal)
  gpDF <-data.frame(trial=rep(seq(1:20), 4),
                    Environment=c(rep('Smooth', 40), rep('Rough', 40)),
                    Context =  rep(rep(c("Conceptual", "Spatial"), each=20),2),
                    meanReward = c(rowMeans(smoothTotal), rowMeans(roughTotal)),
                    meanSE = c(apply(smoothTotal, 1, FUN = function(x) sd(x)/sqrt(length(x))), 
                               apply(roughTotal, 1, FUN = function(x) sd(x)/sqrt(length(x)))))
  
  
  #add mmodel label
  gpDF$Model <- rep(modelName, 80)
  #write to csv
  if(!is.null(outputfile)){
    write.csv(gpDF, outputfile)
  }
  return(gpDF)
}


gpDF <- gpRationalModel(reps, modelName = 'GP-UCB', outputfile = "rationalModels/GPUCB.csv", parameters = gpPars, cores = cores)
#gpMeanDF <- gpRationalModel(reps, modelName = 'GP-GM', outputfile = "rationalModels/GPGM.csv", parameters = gpMeanPars)

#############################################################################################################################
#  SimKernel Model
# to read the previous simulation from disk:
# simDF <- read.csv("rationalModels/SimUCB.csv")
#############################################################################################################################

simRationalModel <- function(replications, modelName, outputfile, parameters, acq=ucb, kernel=shepardKernel){
  smoothTotal <- list()
  roughTotal <- list()
  for (task in c('Conceptual', 'Spatial')){
    #Extract parameters
    smoothPars <- subset(parameters, context == task & environment=="Smooth")
    roughPars <- subset(parameters, context == task & environment =="Rough")
    #choice matrix
    choices <- expand.grid(0:7, 0:7) #build choice matrix
    names(choices)<-c("x1","x2")
    #run for smooth environments
    smoothReward<- mcmapply(1:replications, FUN=function(x){
      #sample parameters
      params <- smoothPars[sample(1:nrow(smoothPars), 1),]
      lambda <- params$lambda
      beta <- params$beta
      tau <- params$tau
      p_minkowski <-params$p #minkowski distance parameter
      #tau <- median(smoothPars$tau)
      #randomly choose environment
      envNum <- sample(1:40,1) 
      #1st trial is random
      location <- sample(1:64,1)#each location as an integer from 1:25; first value is random, so treat both environments as replications of the same actions
      #Observations
      X1 <- choices[location,'x1']
      X2 <- choices[location,'x2']
      #reward
      reward<- c()
      reward[1] <- Y <- (smoothEnvironments[[envNum]][location,"y"])*100
      for (j in 2:20){ #after that, loop through remaining trials and make decisions based on GP preditions
        #compute posterior predictions
        post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001, p_minkowski), X = cbind(X1,X2), Y = (Y-50)/100, k = kernel) #scale observed Y to zero mean, variance of 1
        #compute acquisition function evaluation
        utilityVec <- acq(post, pars = c(beta))
        #scale to max of prevent overflow by subtracting max
        utilityVec <- utilityVec - max(utilityVec)
        #compute softmax choice probabilities
        p <- exp(utilityVec/tau)
        p <- p/sum(p)
        #Sample next choice
        location <- sample(1:64,1, prob=p, replace=TRUE)
        #update reward, X1, X2, and Y for both smooth and rough
        reward[j] <- (smoothEnvironments[[envNum]][location,"y"])* 100
        X1 <- c(X1, choices[location, 'x1'])
        X2 <- c(X2, choices[location, 'x2'])
        Y <- c(Y,  reward[j])}
      reward}, mc.preschedule = TRUE, mc.cores=cores)
    #run for rough environments
    roughReward <- mcmapply(1:replications, FUN=function(x){
      #sample parameters
      params <- roughPars[sample(1:nrow(roughPars), 1),]
      lambda <- params$lambda
      beta <- params$beta
      tau <- params$tau
      p_minkowski <- params$p
      #tau <- median(roughPars$tau)
      #randomly choose environment
      envNum <- sample(1:40,1) 
      #1st trial is random
      location <- sample(1:64,1)#each location as an integer from 1:25; first value is random, so treat both environments as replications of the same actions
      #Observations
      X1 <- choices[location,'x1']
      X2 <- choices[location,'x2']
      #reward
      reward<- c()
      reward[1] <- Y <- (roughEnvironments[[envNum]][location,"y"])*100
      for (j in 2:20){ #after that, loop through remaining trials and make decisions based on GP preditions
        #compute posterior predictions
        post <- gpr(X.test = choices, theta = c(lambda, lambda, 1, 0.0001, p_minkowski), X = cbind(X1,X2), Y = (Y-50)/100, k = kernel) #scale observed Y to zero mean, variance of 1
        #compute acquisition function evaluation
        utilityVec <- acq(post, pars = c(beta))
        #scale to max of prevent overflow by subtracting max
        utilityVec <- utilityVec - max(utilityVec)
        #compute softmax choice probabilities
        p <- exp(utilityVec/tau)
        p <- p/sum(p)
        #Sample next choice
        location <- sample(1:64,1, prob=p, replace=TRUE)
        #update reward, X1, X2, and Y for both smooth and rough
        reward[j] <- (roughEnvironments[[envNum]][location,"y"]) * 100
        X1 <- c(X1, choices[location, 'x1'])
        X2 <- c(X2, choices[location, 'x2'])
        Y <- c(Y,  reward[j])}
      reward}, mc.preschedule = TRUE, mc.cores=cores)
    smoothTotal[[length(smoothTotal)+1]] <- smoothReward
    roughTotal[[length(roughTotal)+1]] <-  roughReward
  }
  #put it all together
  smoothTotal <- do.call("rbind", smoothTotal)
  roughTotal <- do.call("rbind", roughTotal)
  simDF <-data.frame(trial=rep(seq(1:20), 4),
                    Environment=c(rep('Smooth', 40), rep('Rough', 40)),
                    Context =  rep(rep(c("Conceptual", "Spatial"), each=20),2),
                    meanReward = c(rowMeans(smoothTotal), rowMeans(roughTotal)),
                    meanSE = c(apply(smoothTotal, 1, FUN = function(x) sd(x)/sqrt(length(x))), 
                               apply(roughTotal, 1, FUN = function(x) sd(x)/sqrt(length(x)))))
  
  
  #add mmodel label
  simDF$Model <- rep(modelName, 80)
  #write to csv
  if(!is.null(outputfile)){
    write.csv(simDF, outputfile)
  }
  return(simDF)
}



simDF <- simRationalModel(reps, modelName = 'SimKernel-UCB', outputfile = "rationalModels/SimUCB.csv", parameters = simKernelPars)

#############################################################################################################################
# # PLOTS
# # Load all models
randomDF <- read.csv("rationalModels/random.csv")
bmtDF <- read.csv("rationalModels/BMTUCB.csv")
gpDF <- read.csv("rationalModels/GPUCB.csv")
simDF <- read.csv("rationalModels/SimUCB.csv")

#############################################################################################################################

#Join all DFs together 
rationalDF <- rbind(randomDF, bmtDF, gpDF, simDF)
rationalDF$Environment <- factor(rationalDF$Environment, levels=c("Rough", "Smooth"))
rationalDF <- rationalDF[ , !(names(rationalDF) %in% c("X"))]

#add human data
d <- dataImport(normalize=FALSE)

dplot<-ddply(d,~trial+environment+context,summarise,meanReward=mean(z), meanSE=se(z))
colnames(dplot)[2]<-"Environment" #rename colname
colnames(dplot)[3] <- "Context"
dplot$Model <- rep("Human", nrow(dplot))
#Join human with rational models
rationalDF <- rbind(rationalDF, dplot)

summary(rationalDF)
#rationalDF$trial <- rationalDF$trial - 1

#reorder factor levels
rationalDF$Model <- factor(rationalDF$Model, levels = c("Random",  "GP-UCB",'SimKernel-UCB', "BMT-UCB" , "Human"))
levels(rationalDF$Model) <- c("Random",  "GP", "SimKernel", "BMT" , "Human")

#Plot of mean Reward
p1<- ggplot(rationalDF, aes(x=trial, y=meanReward, col=Model, shape=Model))+
  stat_summary(fun.y = mean, geom='line', size = 0.8)+
  #stat_summary(fun.y=mean, geom='point', size = 2)+
  #geom_point(size=2) +
  #geom_line(size=.8) +
  #geom_errorbar(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE), width=0.1, size=1) +
  ylab("Average Reward")+xlab("Trial")+
  theme_classic()+
  facet_grid(Environment~ Context) +
  scale_shape_manual(values = c(32, 16,17,15,4,7)) +
  #scale_color_manual(values=c("black",  "#F0E442", "#E69F00", "#009E73", "#56B4E9", "#fb9a99"))+
  #scale_color_manual(values=c("black","#7570B3", "#1B9E77", "#fb9a99"))+
  scale_color_manual(values=c("black","#56B4E9", "#009E73", "#E69F00","#fb9a99"))+
  #scale_linetype_manual(values = c('solid', 'dotted', 'solid', 'dotted', 'solid', 'solid'))+
  #scale_color_brewer(palette="Paired", direction=1)+
  #coord_cartesian(ylim=c(45,85))+
  theme(text = element_text(size=16,  family="sans"), legend.position="top")+
  scale_x_continuous(breaks = round(seq(0,20, by = 5),1))+
  theme(legend.position="right", strip.background=element_blank(), legend.key=element_rect(color=NA))
p1
ggsave(filename = "plots/modelPerformanceAvg.pdf", plot = p1,height = 3.5, width = 9, units = "in", useDingbats=FALSE) 



# #ggsave(filename = "plots/modelPerformanceLegend.pdf", plot = p1 + theme(legend.position="right"),height =2.82, width = 8, units = "in", useDingbats=FALSE) 
# #Plot of mean Reward
# p2<- ggplot(rational5, aes(x=trial5, y=maxReward, col=Model, shape=Model))+
#   geom_point(size=3) +
#   geom_line(size=1) +
#   #geom_errorbar(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE), width=0.1, size=1) +
#   ylab("Maximum Reward")+xlab("Trial")+
#   theme_bw()+
#   facet_wrap(~Environment) +
#   #scale_shape_manual(values = c(32, 16,17,15,4,7)) +
#   scale_color_manual(values=c("black", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#fb9a99"))+
#   #scale_color_brewer(palette="Paired", direction=1)+
#   coord_cartesian(ylim=c(65,100))+
#   theme(text = element_text(size=16,  family="serif"), legend.position="top")+
#   theme(legend.position="none", strip.background=element_blank(), legend.key=element_rect(color=NA))
# p2
# 
# #Barplot of Max Reward
# maxDF <- subset(rational5, trial5==40)
# p3 <- ggplot(maxDF, aes(x=Model, y=maxReward,fill=Model))+
#   stat_summary(fun.y = mean, geom = "bar", position = "dodge", color='black') + 
#   #geom_errorbar(aes(ymin=meanReward-meanSE, ymax=meanReward+meanSE), width=0.1, size=1) +
#   ylab("Max. Reward")+xlab("Trial")+
#   theme_bw()+
#   scale_fill_manual(values=c("black", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#fb9a99"))+
#   coord_cartesian(ylim=c(60,100))+
#   facet_wrap(~Environment)+
#   theme(text = element_text(size=16,  family="serif"))+
#   theme(legend.position="none", strip.background=element_blank(), legend.key=element_rect(color=NA)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
#   theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
# p3
# ggsave(filename = "plots/modelPerformanceMax.pdf", plot = p2, height =2.5, width = 6, units = "in", useDingbats=FALSE) 
