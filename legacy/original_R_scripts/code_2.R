#########################################################################################################
#A simple example of how to simulate data from a linear birth death model
##########################################################################################################

#Initialisation of parameters for cells dynamics.
N0<-1e+06
Pa=0.95
Pad=0.90
P1=0.35
p2=0.3
# This is a list which containt the number of cells in different division class at different time point


#model definition
linear.bith.model <- function(t,Y,parms){
  with(as.list(c(parms,Y)), {
    
    
      lag1 <- lagvalue(t)[1:(ngen+1)]
    
    temp <- exp(-d_0*deltat)
    dN <- rep(0,length(Y))
    dN[1] <- -(lambda_0 + d_0)*Y[1]
    dN[2] <- 2*lambda_0*lag1[1] - (lambda+d)*Y[2]
    dN[3:(ngen+1)] <- 2*lambda*lag1[2:ngen] - (lambda+d)*Y[3:(ngen+1)]
    
    dN[(ngen+4):(2*ngen+2)] <- 2*lambda*lag1[2:ngen] - lambda*lag1[3:(ngen+1)] - d*Y[(ngen+4):(2*ngen+2)]
    list(dN)
  })
}




#this function simulates the Smith-Martin model
#and returns the number of cells in each generation
#at each time point found in the vector t_pts.
simulate.linear.bith.model <- function(parms){
  out <- dede(init.vec,t_pts,linear.bith.model,parms)
  return (out[,(ngen+3):(2*ngen+3)])
}

#example run of the model
require(deSolve)
#cells in generations ranging from 0 to ngen will be considered

ngen <- 100
init.vec <- rep(0,2*ngen+2)
#number of cells in A phase in i'th generation is denoted NAi
#total number of cells in i'th generation is denoted Ni
names(init.vec) <- paste(rep(c('N'),each=ngen+1),0:ngen,sep='')
#start with 1e6 cells in 0'th generation
init.vec[c('N0')] <- 1e6
#define parameters of the model
parms <- c(lambda_0=1e-2, lambda=1e-2, d_0=1e-5, d=1e-5)
#set the time points at which data should be collected
t_pts <- seq(0,120,12)
#run the simulation
out <- simulate.linear.bith.model(parms)
out

