#################################################################################################
#This is to compare the Smith Martin model and Simulate data for the CFSE labelling 
#1-We start by but the data in the same shape
#2-Secondly we compute the ARE(no done yet)

##########################################################################################################
#A simple example of how to simulate data from a DDE -- i.e. the Smith Martin model
##########################################################################################################

#model definition
smith.martin.model <- function(t,Y,parms){
  with(as.list(c(parms,Y)), {
    if(t < deltat){
      lag1 <- init.vec[1:(ngen+1)]
    }
    else{
      lag1 <- lagvalue(t-deltat)[1:(ngen+1)]
    }
    temp <- exp(-d_0*deltat)
    dN <- rep(0,length(Y))
    dN[1] <- -(lambda_0 + d_0)*Y[1]
    dN[2] <- 2*lambda_0*temp*lag1[1] - (lambda+d)*Y[2]
    dN[3:(ngen+1)] <- 2*lambda*temp*lag1[2:ngen] - (lambda+d)*Y[3:(ngen+1)]
    
    dN[ngen+2] <- -lambda_0*temp*lag1[1] - d_0*Y[ngen+2]
    dN[ngen+3] <- 2*lambda_0*temp*lag1[1] - lambda*temp*lag1[2] - d*Y[ngen+3]
    dN[(ngen+4):(2*ngen+2)] <- 2*lambda*temp*lag1[2:ngen] - lambda*temp*lag1[3:(ngen+1)] - d*Y[(ngen+4):(2*ngen+2)]
    list(dN)
  })
}

#this function simulates the Smith-Martin model
#and returns the number of cells in each generation
#at each time point found in the vector t_pts.
simulate.smith.martin.model <- function(parms){
  out <- dede(init.vec,t_pts,smith.martin.model,parms)
  return (out[,(ngen+3):(2*ngen+3)])
}

#example run of the model
require(deSolve)
#cells in generations ranging from 0 to ngen will be considered

ngen <- 20
init.vec <- rep(0,2*ngen+2)
#number of cells in A phase in i'th generation is denoted NAi
#total number of cells in i'th generation is denoted Ni
names(init.vec) <- paste(rep(c('NA','N'),each=ngen+1),0:ngen,sep='')
#start with 1e6 cells in 0'th generation
init.vec[c('NA0','N0')] <- 50000
#define parameters of the model
parms <- c(lambda_0=1.876064e-2, lambda=3.425390e-2, d_0=1e-3, d=3.431844e-2, deltat=1e-3)
#set the time points at which data should be collected
t_pts <- seq(0,120,12)
#run the simulation
out <- simulate.smith.martin.model(parms)
out
###################################################

#second data
mydata<- read.delim("~/Dropbox/Data/datasets/agent_based_sim_output_3.txt")
#This command help to delete the column time to the data named my data
mydata <- mydata[,colnames(mydata)!="Time"] 
#This command is to have the total number of cell at time points
ncol(mydata)
nw.data <- data.frame(mydata$N0d + mydata$N0p, mydata$N1d + mydata$N1p, 
                      mydata$N2d + mydata$N2p, mydata$N3d + mydata$N3p, mydata$N4d + mydata$N4p,
                      mydata$N5d + mydata$N5p, mydata$N6d + mydata$N6p, mydata$N7d + mydata$N7p, 
                      mydata$N8d + mydata$N8p, mydata$N9d + mydata$N9p,
                      mydata$N10d + mydata$N10p, mydata$N11d + mydata$N11p,
                      mydata$N12d + mydata$N12p, mydata$N13d + mydata$N13p,
                      mydata$N14d + mydata$N14p, mydata$N15d + mydata$N15p,
                      mydata$N16d + mydata$N16p, mydata$N17d + mydata$N17p,
                      mydata$N18d + mydata$N18p, mydata$N19d + mydata$N19p,
                      mydata$N20d + mydata$N20p) 
names(nw.data)<-c("N0","N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10", 
                  "N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19", "N20")
nw.data <- rbind(c(50000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),nw.data)
View(nw.data,out)
plot(N,)
