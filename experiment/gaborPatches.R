#Gabor patche stimuli

rm(list=ls()) #house keeping

library(grt)
numFeatures <- 8

#define parameter values
#angles <- seq(225, 135, length.out=numFeatures)
angles <- seq(255, 105, length.out=numFeatures)
#angles <- seq(270, 90, length.out=numFeatures+2)
#angles <- angles[2:(numFeatures+1)] #trim first and last , which are visually indistinguishable
cycles <- exp(seq(from = log(1.5), to = log(15), length.out = 8))


#create grid
op <- par(mfrow=c(numFeatures,numFeatures),  
          pty = "m",
          oma = c(0, 0, 0, 0), 
          mai=c(0.01,0.01,0.01,0.01),
          mgp = c(0, 0, 0),
          xpd = NA)
for(i in rev(cycles)){ #Cycles
  for(j in angles){ #rotation  
    gaborPatch(i, theta = j)
  }
}
par(op)


#generate individual outputs

for(i in 1:numFeatures){ #rotation  
  for(j in 1:numFeatures){ #Cycles
    angle <- angles[i]
    cycle <- cycles[j]
    filename = paste0('gabor/gabor.',(i-1),'.',(j-1),'.png')
    png(filename=  filename, width = 240, heigh=240, units = 'px', res = 300)
    par(mar=rep(0,4))
    gaborPatch(cycle, theta = angle)
    dev.off()
  }
}


