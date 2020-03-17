#dataProcessing
#Charley Wu


#load packages
packages <- c('jsonlite')
lapply(packages, require, character.only = TRUE)

num_rounds = 10
num_trials = 20
gridSize_0 = 7 #grid size in base 0
#####################################################################################################################################################################
#Load experment data
#####################################################################################################################################################################
jsonFixer <- function(dataset){
  lapply(dataset, function(j) {
  if(is.list(j)){
    matrix(unlist(j), ncol=2, byrow=T)
  }else{
    j
  }}
)
}

dataImport <- function(dataFile ="experimentData/full.csv", normalize=TRUE){
  #Load data
  rawdata <- read.csv(dataFile, sep=",")
  #remove empty rows
  rawdata <- rawdata[!grepl("NULL",rawdata$grid_assignmentID) & !grepl("NULL",rawdata$gabor_assignmentID),] #Todo: replace this with a more sophisticated way to check that both parts were completed
  #extract sampleSize, records, and dat (deeply nested json data)
  sampleSize <-nrow(rawdata)
  all_opts = expand.grid(0:gridSize_0, 0:gridSize_0)
  #build dataframe
  df<-data.frame()
  
  #Loop through data and build df
  for (i in 1:sampleSize){
    dat = rawdata[i,]
    #general data
    id = dat$id
    age = dat$grid_age #todo: case if the two don't match
    gender = dat$grid_gender
    environment = dat$environment
    contextOrder = dat$scenario
    gridBonus = dat$grid_reward 
    gaborBonus = dat$gabor_reward
    totalBonus = gridBonus + gaborBonus
    
    
    ##Grid data
    gridHistory <- fromJSON(as.character(dat$grid_experimentData))
    gridComprehensionTries <- gridHistory$comprehensionQuestionTries
    grid_x <- as.vector(t(gridHistory$xcollect))
    grid_y <- as.vector(t(gridHistory$ycollect))
    initValues <- do.call(rbind, jsonFixer(gridHistory$initcollect))
    grid_initx<- c(initValues[1:195,1], NA, initValues[196:199,1]) #inserting the null for the bonus round
    grid_inity<- c(initValues[1:195,2], NA, initValues[196:199,2])
    grid_chosen <- apply(cbind(grid_x,grid_y),MARGIN=1, FUN=function(row) which(row[1]==all_opts$Var1 & row[2]==all_opts$Var2))
    grid_ts <- as.vector(t(gridHistory$tscollect))
    grid_z<-as.vector(t(gridHistory$zcollect)) #NOTE ignore bonus round for reaction time data
    if (normalize==TRUE){
      grid_z <- (grid_z-50)/100
    }else{
      grid_z <- grid_z  
    }
    grid_zscaled <- as.vector(t(gridHistory$zscaledcollect))
    grid_previousz <-  as.vector(rbind(rep(NA, num_rounds),matrix(grid_z,nrow = num_trials)[1:(num_trials-1),]))
    gridScale <- as.vector(t(gridHistory$scaleCollect))[1:num_rounds]
    if ('envOrder' %in% names(gridHistory)){#Added late to the second pilot, so it might not exist for that data
      gridEnvOrder <-as.vector(t(gridHistory$envOrder)) 
      gridEnvOrder <- rep(gridEnvOrder, each=num_trials)
    }else {gridEnvOrder=rep(NA, num_rounds*num_trials)}
    #Distance from previous selection
    grid_delta_x <- abs(matrix(grid_x,ncol = num_rounds)[1:(num_trials-1),] - matrix(grid_x,ncol = num_rounds)[2:num_trials,])
    grid_delta_y <- abs(matrix(grid_y,ncol = num_rounds)[1:(num_trials-1),] - matrix(grid_y,ncol = num_rounds)[2:num_trials,])
    gridDistance <- grid_delta_x + grid_delta_y #calculate MANHATTAN distance
    gridDistance <- as.vector(rbind(rep(NA,num_rounds), gridDistance))
    #Distance from initial location
    grid_movement_x <- abs(matrix(grid_x,ncol = num_rounds)- matrix(grid_initx,ncol = num_rounds))
    grid_movement_y <- abs(matrix(grid_y,ncol = num_rounds)- matrix(grid_inity,ncol = num_rounds))
    gridMovement <- as.vector(grid_movement_x + grid_movement_y) #calculate MANHATTAN distance
  
    #Trajectories in main task
    concatedTraj <- lapply(gridHistory$stepcollect, function(j) {
      if(is.list(j)){
        lapply(j, toJSON)
      }else{
        lapply(1:nrow(j), function(k) toJSON(j[k,]))
      }}
    )
    grid_trajectories <- unlist(concatedTraj)
    grid_steps <- sapply(1:length(grid_trajectories), FUN=function(i) length(fromJSON(grid_trajectories[i])))
    #Trajectories in tutorial task (only summary data, since this data frame is designed to be 1 row per trial of the main task)
    gridtraj <-  gridHistory$trajCollect
    gridTrajError <- rowSums((gridtraj$targetcollect - gridtraj$selectioncollect)^2)
    gridTrajCorrect <- sum(gridTrajError==0)/length(gridTrajError)
    gridTrajRMSE <- sqrt(sum(gridTrajError))
    gridTrajAvgSteps <-  mean(unlist(lapply(gridtraj$stepcollect, function(l) length(l))))
    
    #Conceptual data
    gaborHistory <- fromJSON(as.character(dat$gabor_experimentData))
    gaborComprehensionTries <- gaborHistory$comprehensionQuestionTries
    gabor_x <- as.vector(t(gaborHistory$xcollect))
    gabor_y <- as.vector(t(gaborHistory$ycollect))
    initValues <- do.call(rbind, jsonFixer(gaborHistory$initcollect))
    gabor_initx<- c(initValues[1:195,1], NA, initValues[196:199,1]) #inserting the null for the bonus round
    gabor_inity<- c(initValues[1:195,2], NA, initValues[196:199,2]) #inserting the null for the bonus round
    gabor_chosen <- apply(cbind(gabor_x,gabor_y),MARGIN=1, FUN=function(row) which(row[1]==all_opts$Var1 & row[2]==all_opts$Var2))
    gabor_ts <- as.vector(t(gaborHistory$tscollect))
    gabor_z<-as.vector(t(gaborHistory$zcollect)) #NOTE ignore bonus round for reaction time data
    if (normalize==TRUE){
      gabor_z <- (gabor_z-50)/100 #rescale to zero mean
    }else{
      gabor_z <- gabor_z 
    }
    gabor_zscaled <- as.vector(t(gaborHistory$zscaledcollect))
    gabor_previousz <-  as.vector(rbind(rep(NA, num_rounds),matrix(gabor_z,nrow = num_trials)[1:(num_trials-1),]))
    gaborScale <- as.vector(t(gaborHistory$scaleCollect))[1:num_rounds]
    if ('envOrder' %in% names(gaborHistory)){#Added late to the second pilot, so it might not exist for that data
      gaborEnvOrder <-as.vector(t(gaborHistory$envOrder)) 
      gaborEnvOrder <- rep(gaborEnvOrder, each=num_trials)
      }else {gaborEnvOrder=rep(NA, num_rounds*num_trials)}
    #Distance from previous selection
    gabor_delta_x <- abs(matrix(gabor_x,ncol = num_rounds)[1:(num_trials-1),] - matrix(gabor_x,ncol = num_rounds)[2:num_trials,])
    gabor_delta_y <- abs(matrix(gabor_y,ncol = num_rounds)[1:(num_trials-1),] - matrix(gabor_y,ncol = num_rounds)[2:num_trials,])
    gaborDistance <- gabor_delta_x + gabor_delta_y #calculate MANHATTAN distance
    gaborDistance <- as.vector(rbind(rep(NA,num_rounds), gaborDistance))
    #Distance from initial location
    gabor_movement_x <- abs(matrix(gabor_x,ncol = num_rounds)- matrix(gabor_initx,ncol = num_rounds))
    gabor_movement_y <- abs(matrix(gabor_y,ncol = num_rounds)- matrix(gabor_inity,ncol = num_rounds))
    gaborMovement <- gabor_movement_x+ gabor_movement_y #MANHATTAN distance
    #Trajectories
    concatedTraj <- lapply(gaborHistory$stepcollect, function(j) {
      if(is.list(j)){
        lapply(j, toJSON)
      }else{
        lapply(1:nrow(j), function(k) toJSON(j[k,]))
      }}
    )
    gabor_trajectories <- unlist(concatedTraj)
    gabor_steps <- sapply(1:length(gabor_trajectories), FUN=function(i) length(fromJSON(gabor_trajectories[i])))
    #Trajectories in tutorial task (only summary data, since this data frame is designed to be 1 row per trial of the main task)
    gabortraj <-  gaborHistory$trajCollect
    gaborTrajError <- rowSums((gabortraj$targetcollect - gabortraj$selectioncollect)^2)
    gaborTrajCorrect <- sum(gaborTrajError==0)/length(gaborTrajError)
    gaborTrajRMSE <- sqrt(sum(gaborTrajError))
    gaborTrajAvgSteps <-  mean(unlist(lapply(gabortraj$stepcollect, function(l) length(l))))
    
    #Round and trial data
    round<-rep(rep(1:num_rounds, each=num_trials),2)
    trial <- rep(rep(1:num_trials,num_rounds),2)
    n <- length(trial) #number of datapoints
    
    #Start and end times
    grid_start <-  strptime(as.character(dat$grid_task_start),"%Y-%m-%d %H:%M",tz="GMT")
    grid_end <-  strptime(as.character(dat$grid_task_end),"%Y-%m-%d %H:%M",tz="GMT")
    grid_duration <- grid_end - grid_start
    
    gabor_start <-  strptime(as.character(dat$gabor_task_start),"%Y-%m-%d %H:%M",tz="GMT")
    gabor_end <-  strptime(as.character(dat$gabor_task_end),"%Y-%m-%d %H:%M",tz="GMT")
    gabor_duration <- gabor_end - gabor_start
    
    grid_gabor_gap <- ifelse(gabor_start > grid_start, gabor_start - grid_end, grid_start - gabor_end) #hours

    
    dummy<-data.frame(id=rep(id, n), age=rep(age, n), gender=rep(gender, n), environment=rep(environment,n), contextOrder=rep(contextOrder,n), context=rep(c("Spatial", "Conceptual"), each=n/2), 
                      round=round, trial=trial, x=c(grid_x, gabor_x), y=c(grid_y, gabor_y), chosen=c(grid_chosen, gabor_chosen), initx = c(grid_initx, gabor_initx), inity = c(grid_inity, gabor_inity),
                      trajectories = c(grid_trajectories, gabor_trajectories), steps = c(grid_steps, gabor_steps), movement = c(gridMovement, gaborMovement),
                      distance = c(gridDistance, gaborDistance), distance_x =c(as.vector(rbind(rep(NA,num_rounds),grid_delta_x)), as.vector(rbind(rep(NA,num_rounds),gabor_delta_x))), 
                      distance_y =c(as.vector(rbind(rep(NA,num_rounds),grid_delta_y)), as.vector(rbind(rep(NA,num_rounds),gabor_delta_y))),
                      z=c(grid_z, gabor_z), zscaled=c(grid_zscaled, gabor_zscaled), previousReward=c(grid_previousz, gabor_previousz), 
                      ts=c(grid_ts, gabor_ts), scale=rep(c(gridScale, gaborScale), each=num_trials), envOrder=c(gridEnvOrder, gaborEnvOrder), bonus = c(rep(gridBonus, num_rounds*num_trials), rep(gaborBonus, num_rounds*num_trials)), totalBonus = totalBonus,
                      trajCorrect = c(rep(gridTrajCorrect, n/2), rep(gaborTrajCorrect, n/2)), trajRMSE = c(rep(gridTrajRMSE,n/2), rep(gaborTrajRMSE, n/2)), trajAvgSteps = c(rep(gridTrajAvgSteps,n/2), rep(gaborTrajAvgSteps, n/2)),
                      grid_start = rep(grid_start,n), grid_end = rep(grid_end,n), grid_duration = rep(grid_duration,n),
                      gabor_start = rep(gabor_start,n), gabor_end = rep(gabor_end,n), gabor_duration=rep(gabor_duration,n),
                      grid_gabor_gap=rep(grid_gabor_gap,n), comprehensionTries = c(rep(gridComprehensionTries, n/2), rep(gaborComprehensionTries, n/2)))
    df<-rbind(df, dummy)
  }
  
  summary(df)
  #factors
  df$environment<-ifelse(df$environment==0, "Rough", "Smooth")
  df$environment <- factor(df$environment)
  df$contextOrder <- factor(df$contextOrder)
  
  return(df)
}


importTrajData <- function(dataFile ="experimentData/full.csv", normalize=TRUE){
  #Load data
  rawdata <- read.csv(dataFile, sep=",")
  #remove empty rows
  rawdata <- rawdata[!grepl("NULL",rawdata$grid_assignmentID) & !grepl("NULL",rawdata$gabor_assignmentID),] #Todo: replace this with a more sophisticated way to check that both parts were completed
  #extract sampleSize, records, and dat (deeply nested json data)
  sampleSize <-nrow(rawdata)
  all_opts = expand.grid(0:gridSize_0, 0:gridSize_0)
  #build dataframe
  df<-data.frame()
  
  #Loop through data and build df
  for (i in 1:sampleSize){
    dat = rawdata[i,]
    #general data
    id = dat$id
    age = dat$grid_age #todo: case if the two don't match
    gender = dat$grid_gender
    environment = dat$environment
    contextOrder = dat$scenario
    gridBonus = dat$grid_reward 
    gaborBonus = dat$gabor_reward
    totalBonus = gridBonus + gaborBonus
    
    
    ##Grid data
    gridHistory <- fromJSON(as.character(dat$grid_experimentData))
    #Trajectories in tutorial task 
    gridtraj <-  gridHistory$trajCollect
    gridTrajError <- rowSums((gridtraj$targetcollect - gridtraj$selectioncollect)^2)
    gridTrajManhattanDistance <- rowSums(abs(gridtraj$targetcollect - gridtraj$selectioncollect))
    gridTrajManhattanDistance_x <- abs(gridtraj$targetcollect[,1]- gridtraj$selectioncollect[,1]) #Separate X and Y components
    gridTrajManhattanDistance_y <- abs(gridtraj$targetcollect[,2]- gridtraj$selectioncollect[,2])
    gridTrajCorrect <- gridTrajError==0
    gridTrajSteps <-  unlist(lapply(gridtraj$stepcollect, function(l) length(l)))
    gridX <- gridtraj$targetcollect[,1]
    gridY <- gridtraj$targetcollect[,2]
    
    #Conceptual data
    gaborHistory <- fromJSON(as.character(dat$gabor_experimentData))
  
    #Trajectories in tutorial task (only summary data, since this data frame is designed to be 1 row per trial of the main task)
    gabortraj <-  gaborHistory$trajCollect
    gaborTrajError <- rowSums((gabortraj$targetcollect - gabortraj$selectioncollect)^2)
    gaborTrajManhattanDistance <- rowSums(abs(gabortraj$targetcollect - gabortraj$selectioncollect))
    gaborTrajManhattanDistance_x <- abs(gabortraj$targetcollect[,1]- gabortraj$selectioncollect[,1])
    gaborTrajManhattanDistance_y <- abs(gabortraj$targetcollect[,2]- gabortraj$selectioncollect[,2])
    gaborTrajCorrect <- gaborTrajError==0
    gaborTrajSteps <-  unlist(lapply(gabortraj$stepcollect, function(l) length(l)))
    gaborX <- gabortraj$targetcollect[,1]
    gaborY <- gabortraj$targetcollect[,2]
    
    #Round and trial data
    trial <- c(1:length(gridTrajError), 1:length(gaborTrajError))
    n <- length(trial) #number of datapoints
    
    dummy<-data.frame(id=rep(id, n), age=rep(age, n), trial = trial, gender=rep(gender, n), environment=rep(environment,n), 
                      contextOrder=rep(contextOrder,n), context=c(rep("Spatial",length(gridTrajError)), rep("Conceptual", length(gaborTrajError))), 
                      x = c(gridX, gaborX), y = c(gridY, gaborY),trajCorrect = c(gridTrajCorrect, gaborTrajCorrect),
                      trajError = c(gridTrajError, gaborTrajError), manhattanError = c(gridTrajManhattanDistance, gaborTrajManhattanDistance),  
                      manhattanError_x = c(gridTrajManhattanDistance_x, gaborTrajManhattanDistance_x), manhattanError_y = c(gridTrajManhattanDistance_y, gaborTrajManhattanDistance_y),
                      trajSteps = c(gridTrajSteps, gaborTrajSteps))
    df<-rbind(df, dummy)
  }
  
  summary(df)
  #factors
  df$environment<-ifelse(df$environment==0, "Rough", "Smooth")
  df$environment <- factor(df$environment)
  df$contextOrder <- factor(df$contextOrder)
  
  return(df)
}



se <- function(x){sd(x)/sqrt(length(x))}
