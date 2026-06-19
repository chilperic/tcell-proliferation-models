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

#example run of the mod
require(deSolve)
#cells in generations ranging from 0 to ngen will be considered

ngen <- 20
init.vec <- rep(0,2*ngen+2)
#number of cells in A phase in i'th generation is denoted NAi
#total number of cells in i'th generation is denoted Ni
names(init.vec) <- paste(rep(c('NA','N'),each=ngen+1),0:ngen,sep='')
#start with 1e6 cells in 0'th generation
init.vec[c('NA0','N0')] <- 5e4
#define parameters of the model
parms <-c(lambda_0=1.870167e-02, lambda=3.422586e-02, d_0=1.000000e-03, d=3.425349e-02, deltat=1.000000e-3)
#set the time points at which data should be collected
t_pts <- seq(0,120,12)
#run the simulation
out <- simulate.smith.martin.model(parms)
out

View(out)


sum.cols <- function(idx, data) return (rowSums(data[,idx]))


raw.data <- read.table("~/Dropbox/Data/datasets/agent_based_sim_output_2.txt",header = TRUE)
View(raw.data)

View(real.data)
real.data <- matrix(rep(0, (length(t_pts)*(ngen+1))), nrow = length(t_pts))
real.data[1,1] <- init.vec[1]
real.data[-1,]  <- sapply(split(2:43,rep(1:21,each=2)),sum.cols,raw.data)

plot(out[,2],xlab = "time points",ylab = "Generation 0",type="b", col="red")
lines(real.data[,2], type="b", col="blue")
grid(10,10,lwd=1)
legend("topright",legend=c("Real Data", "SM Data", "Cyton Data"),col= c("blue","red","cyan"),yjust=1,lty=1,bty="n",cex=1, inset=0)


#for ( i in seq(1,length(out ),1) ) ts.plot(out[,i],ylab=names(out[1,i]),type="l")
