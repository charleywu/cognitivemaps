#Model Comparison 
#Charley Wu Dec 2019

rm(list=ls()) #house keeping

#load packages
#packages <- c('plyr', 'ggplot2', 'jsonlite', 'MASS', 'gridExtra', 'parallel', 'matrixcalc', 'DEoptim', 'fields')
packages <- c('plyr', 'jsonlite', 'DEoptim', "matrixcalc", "fields", "MASS")
lapply(packages, require, character.only = TRUE)

#Source dependencies 

source('dataProcessing.R')
source('models.R') 

participants = 129
horizon = 20 
##############################################################################################################
#Cluster configuration: (1 subject x model combination) per CPU
##############################################################################################################

#IMPORTANT: update batch name
batchName = 'batch5' #saves output csv files in 'modelResults/batchName/*.csv'

#InertiaWeighting of acquisition functions used in Local-GP models
#inertiaTypes <- c('prevChoice', 'initialPoint')
#inertiaTypes <- c('prevChoice')
#inertiaTypes <- c('initialPoint')
inertiaTypes <- c()
inertiaWeight<- FALSE

tasks <- c('Spatial', 'Conceptual')
#tasks <- c('Conceptual')

#create list of all kernel functions
#kernellist<-list(bayesianMeanTracker, rbf, shepardKernel) 
kernellist<-list(rbf, shepardKernel) 
#kernellist<-list(rbf) 
#kernellist<-list(bayesianMeanTracker) 

#names of all kernel functions
kernelnames<-c("RBF", 'SimKernel')
#kernelnames<-c("RBF")
#kernelnames<-c("BMT")

#list of all acquisition functions
acqlist<-list(ucb)
#acqlist<-list(greedyVar, greedyMean, ucb)
#acqlist<-list(exofimp, poi, pmu)

#names of all acquisition functions
acqnames<-c("UCB")
#acqnames<-c("GV", "GM", "UCB")
#acqnames<-c("EXI", "POI", "PMU")


#Cluster id from qsub
#N <-(participants*length(kernellist)*(length(acqlist))*length(tasks) * max(1,length(inertiaTypes)))
# clusterid <- sample(1:N, 1)
clusterid <- as.integer(commandArgs(TRUE)[1]) #Cluster id, corresponds to an integer used to indicate which combination of kernel and acquisition function to simulate

#all combinations of kernels, acquisition functions, and tasks will be needed
combs<-expand.grid(1:length(kernellist), 1:length(acqlist) , 1:length(tasks), 1:max(1,length(inertiaTypes)))

#create a matrix with  combinations of subjectIds and model combinations
subjectComb <- expand.grid(1:participants, 1:(length(kernellist) * length(acqlist) * length(tasks) * max(1,length(inertiaTypes)))) #1:? defines the number of unique models to be analyzed
subjectId <- subjectComb[clusterid,1] #used to identify unique subjects
combId <- subjectComb[clusterid,2] #used to identify a unique model (i.e., combination of GP kernel and acquisition function)

#trim combs if cluster==TRUE
model <- combs[combId,]
task<- tasks[model[[3]]]
if(inertiaWeight==T){
  inertiaType <- inertiaTypes[model[[4]]]  
}else{
  inertiaType <- ''
}


set.seed(clusterid) #set seed as the participant id

#parameter bounds; all parameters expect p are exponentiated inside the optimization function
paramUBound <- list(tau = 5, beta = 5, lambda = 5, errorVariance  = 5, p = 2) 
paramLBound <- list(tau = -5,  beta = -5, lambda = -5, errorVariance  = -5, p = 0.0001)  


params <- paramLabeller(kernellist[[model[[1]]]], acqlist[[model[[2]]]])
allParams <- rep(NA,5)
names(allParams) <- c('tau', 'beta', 'lambda', 'errorVariance', 'p') #all parameters
##############################################################################################################
#Compile Participant Data
##############################################################################################################
data <- dataImport() #sourced from dataMunging.R

uid <- unique(data$id)[subjectId] #convert subjectId to uid

if(inertiaWeight){ #if localized models, initialize inertia weights by reading pre-computed csv files
  blocks<- readBlockDistance(uid, task,inertiaType)
  }
##############################################################################################################
#Model Fitting
##############################################################################################################

#Negative Log Likelihood of model prediction 
#parameters are distinct for each individual model
#subjD = subset of data for a specific subject
#kernel = any of fastKalmanFilter, standard GP (rbf, matern, or oru), or lingp
#acquisition = any of ucb, probofimp, expofimp, or pmu
#rounds list of numbers between 1:4, relative to the order of rounds with the specified horizon length
#inertiaWeight==TRUE weighs the q(x) of the acquisition function by the inverse manhattan distance
#refactorUCB: whether or not to refactor tau and beta parameters in UCB
#returnPredictions: whether or not to also return the dataframe of probabilistic model predictions for each choice. False for MLE, True for out of sample prediction
modelFit<-function(par, subjD, acquisition, k, rounds, inertiaWeight=inertiaWeight, returnPredictions=F){
  #Extract and process parameters
  names(par) <- params
  tau<-exp(par['tau'])  #inverse temperature for softmax; all models have it
  #Which posterior function to use; therefore, which parameters to use
  if (inherits(k, "KalmanFilter")){ #null kernel indicates kalman filter model
    kNoise <- exp(par['errorVariance'])
    parVec <- c(kNoise) #Vector of parameters to send to the KF posterior function
  }else if(inherits(k, "GP")){ #lambda
    lambda <- exp(par['lambda'])
    parVec <- c(lambda, lambda, 1, .0001) # Vector of parameters to send to the GP posterior vector, where sF and sN are fixed
  }else if(inherits(k, "minkowskiGP")){ #lambda
    lambda <- exp(par['lambda'])
    p_minkowski <- par['p'] #Unexponentiated!
    parVec <- c(lambda, lambda, 1, .0001, p_minkowski)
  }
  #Additional acquisition function dependent parameters
  if (inherits(acquisition, "UCB")){ #check if UCB is used
    beta <- exp(par['beta']) #If UCB, beta is always 2nd last
  }
  #Vector to store negative log likelihods
  nLL <- rep(0,length(rounds)) 
  for (r in rounds){ #Begin looping through each round
    #subset of data for round r
    roundD <- subset(subjD, round==r)
    #Observations of subject choice behavior
    chosen <- roundD[2:20,'chosen'] #important! because we don't have a random initially revealed tile, we are only making predictions on 19 of the 20 choices
    y  <- roundD$z[1:horizon] #trim off the last observation, because it was not used to inform a choice (round already over)
    x1 <- roundD[1:horizon,'x']
    x2 <- roundD[1:horizon,'y']
    #create observation matrix
    X<-as.matrix(cbind(x1,x2))
    #initialize Xtest
    Xnew<-as.matrix(expand.grid(0:7,0:7))
    #make sure X is a matrix
    X<-as.matrix(X)
    #Utilties of each choice
    utilities <-NULL
    prevPost <- NULL
    #loop through observations
    for (i in 1:(horizon-1)){ #skip the first observation, since that was completely random
      #new observation
      X1<-matrix(X[1:i,], ncol=2)
      y1<-y[1:i]
      #Which posterior function to use
      if (inherits(k, "KalmanFilter")){# kalman filter model
        out<- bayesianMeanTracker(x = X1[i,], y=(y[i]), prevPost = prevPost, theta = parVec)
        #update prevPost for the next round
        prevPost <- out
      }else if (inherits(k, "linGP")){#linear kernel GP
        out <- lingp(subjD$id[1], i, as.integer(r))
      }else if (inherits(k, "inertia")){#inertia model
        out <- inertia(roundD$chosen[i], distanceMatrix)
        beta <- 0 #set a dummy beta model to make it work with UCB
      }else if (inherits(k, "WSLS")){ #win stay lose go
        out<- WSLS(y1, X1)
      }else{# GP with length-scale parameterized kernel
        out <- gpr(X.test=Xnew, theta=parVec, X=X1, Y=(y1), k=k) #Mu and Sigma predictions for each of the arms; either GP or Kalman filter
      }
      #Slightly different function calls for each acquisition function
      if (inherits(acquisition, "UCB")){ #UCB takes a beta parameter
        utilityVec<-acquisition(out, c(beta))
      }else if(inherits(acquisition, "Imp")){ #ProbImp or ExpectedImp
        y.star <- max(roundD$z[1:i]) #Best revealed solution
        utilityVec <- acquisition(out, y.star)
      }else{ #PMU or any other
        utilityVec <- acquisition(out)
      }
      if (inertiaWeight==TRUE){ #if inertiaWeighting is specified
        #weight the utilityVec by the inverse euclidean distance
        utilityVec <- utilityVec * inertia(blocks[[r]], i)$mu
      }
      utilityVec <- utilityVec - max(utilityVec) #avoid overflow
      utilities <- rbind(utilities, utilityVec) # build horizon_length x 25 matrix, where each row holds the utilities of each choice at each decision time in the search horizon
    }
    #Softmax rule
    p <- exp(utilities/tau)
    p <- p/rowSums(p)
    #avoid underflow by setting a floor and a ceiling
    p <- (pmax(p, 0.00001))
    p <- (pmin(p, 0.99999))
    #Calculate Negative log likelihood
    nLL[which(r==rounds)] <- -sum(log(p[cbind(c(1:(horizon-1)),chosen)]))
  }#end loop through rounds
  if (returnPredictions==FALSE){ #Return only summed log loss for computing MLE
    return(sum(nLL))  #Return negative log likelihoods of all observations 
  }else if (returnPredictions==TRUE){ #Return detailed dataframe of model predictions for outofsample predictions
    detailedPredictions <- list(sum(nLL), p, roundD[2:20,'chosen']) #construct a list of summed log loss, model predictions, and chosen item
    return(detailedPredictions) #return detailed predictions
  }
}

##############################################################################################################
#CROSS VALIDATION FUNCTION
##############################################################################################################
#function to plug in to the optimaztion routine
#selector: scalar, indicates a specific participant
#kernel, function, can be "rbf", "oru", "matern", or NULL (Kalman filter)
#acquisition, function, can be "ucb", "probofimp", "expofimp", or "PMU"
#horizonLength is 20 or 40
#leaveoutindex is 1,2,3, or 4
#inertiaWeight: whether nor not acquisition functions are weighted by inertia
cvfun<-function(selector, kernelFun, acquisition, leaveoutindex, inertiaWeight){
  #subselect participant, horizon and rounds not left out
  d1<-subset(data, id==selector & context==task)
  #training set
  rounds <- seq(1,9) #Trim  last round
  trainingSet <- rounds[! rounds==rounds[leaveoutindex]] #remove round specified by leaveoutindex
  #test set
  testSet <- rounds[leaveoutindex]
  #Compute MLE on the trainingSet
  lbound<-unlist(paramLBound[params])
  ubound<-unlist(paramUBound[params])
  #Begin cross validation routine
  if (length(params)>=2){#if 2 or more parameters
    #TRAINING SET
    fit<-DEoptim(modelFit, lower=lbound, upper=ubound, subjD=d1, k=kernelFun, rounds = trainingSet, acquisition=acquisition,inertiaWeight=inertiaWeight, DEoptim.control(itermax=200))
    paramEstimates <- fit$optim$bestmem #MODEL DEPENDENT PARAMETER ESTIMATES
    allParams[params] <- paramEstimates #index all possible parameters by the parameters of this given model, and insert those values
    #TEST SET
    predict <- modelFit(par=paramEstimates, subjD=d1, acquisition=acquisition, k=kernelFun, rounds=testSet, inertiaWeight=inertiaWeight, returnPredictions=TRUE)
    cvresults <- data.frame(loo = leaveoutindex,nLL =  predict[[1]], tau = exp(allParams['tau']), beta = exp(allParams['beta']), lambda = exp(allParams['lambda']), errorVariance = exp(allParams['errorVariance']), p = allParams['p']) #leaveoutindex, nLL, parameters....
    output <- list(cvresults, predict[[2]], predict[[3]])
  } else{
    #TRAINING SET
    fit<-optimize(modelFit, lower=lbound, upper=ubound, subjD=d1, k=kernelFun, rounds = trainingSet, acquisition=acquisition, inertiaWeight=inertiaWeight)
    paramEstimates <- fit$minimum #MODEL DEPENDENT PARAMETER ESTIMATES
    allParams[params] <- paramEstimates #index all possible parameters by the parameters of this given model, and insert those values
    #TEST SET
    predict <- modelFit(par=paramEstimates, subjD=d1, acquisition=acquisition, k=kernelFun, rounds=testSet, inertiaWeight=inertiaWeight, returnPredictions=TRUE)
    cvresults <- data.frame(loo = leaveoutindex,nLL =  predict[[1]], tau = exp(allParams['tau']), beta = exp(allParams['beta']), lambda = exp(allParams['lambda']), errorVariance = exp(allParams['errorVariance']), p = allParams['p']) #leaveoutindex, nLL, parameters....
    output <- list(cvresults, predict[[2]], predict[[3]])
  }
  return(output) #return optimized value
}

##############################################################################################################
#OPTIMIZATION ROUTINE
##############################################################################################################


#cross-validation routine
crossvalidation <- data.frame() #cross-validation results
modelPredictions <- lapply(1:9, matrix, data= NA, nrow=19, ncol=64)
chosenMatrix <- matrix(NA, nrow=9, ncol=19)


start.time <- Sys.time()
for (loo in 1:9){#leave one out index
  cv <- cvfun(selector = uid, kernelFun=kernellist[[model[[1]]]], acquisition = acqlist[[model[[2]]]], leaveoutindex=loo, inertiaWeight=inertiaWeight)
  crossvalidation <- rbind(crossvalidation, cv[[1]])
  modelPredictions[[loo]] <- cv[[2]]
  chosenMatrix[loo,] <- cv[[3]]
}

output <- list(crossvalidation, modelPredictions, chosenMatrix)
#save the vector with kernel-acquisition-pair as name
name<-paste0("modelResults/", batchName, "/",tasks[model[[3]]],'.', inertiaType, kernelnames[model[[1]]], acqnames[model[[2]]], uid)
save(output, file=paste0(name, '.Rdata'))

print(output)

end.time <- Sys.time()
elapsed <- end.time - start.time
print(elapsed)
##############################################################################################################
#THE END
##############################################################################################################
