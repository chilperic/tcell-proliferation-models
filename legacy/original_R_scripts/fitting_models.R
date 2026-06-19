#######################################################################################
## This code fit the cyton model to average of simulate data
#######################################################################################  
ngen <- 20
t_pts <- seq(0,120,12)
#run the simulation

####################CYTON MODELS##########################################################
N <- 50000
parmsc <-  c(f0=1.000000e+00, f1=8.115505e-01, mu0d=4.300041e+02,mud=2.341927e+01,mu0p=5.316612e+01,mup=1.722716e+01)
names(parmsc) <- c('f0', 'f1', 'mu0d','mud','mu0p','mup') #names of the parameters
   
#######################################################################
#   simulate  cyton model -- modified code by Lamin
#######################################################################
#calculations based on explicit analytic formulas
#total number of dividing cell; an explicit expression for r_n^div (defined implicity
#in eqn 13 of Miao et al) was derived and then inserted into eqn 15 of Miao et al to get
#the total number of cells.
sum.rn.div2 <- function(i,a,n,t,theta){
  j <- 0:i
  return ((-1)^(n+i+1)/a^(n-i)*(1/theta^(i+1) - exp(-theta*t)*sum(t^(i-j)/factorial(i-j)/theta^(j+1))))
}
sum.rn.div <- function(n,t,par){
  theta0 <- 1/par['mu0p'] + 1/par['mu0d']
  if(n == 0){
    return (par['f0']*N/par['mu0p']/theta0*(1 - exp(-theta0*t)))
  }
  else{
    theta <- 1/par['mup'] + 1/par['mud']
    a <- theta0 - theta
    temp <- par['f0']/par['mu0p']*(2*par['f1']/par['mup'])^n*N
    temp <- temp*(((-1/a)^n*(1 - exp(-theta0*t))/theta0) +
                    sum(sapply(0:(n-1),sum.rn.div2,a,n,t,theta)))
    return (temp)
  }
}
#calculations based on explicit analytic formulas
#total number of dieing cells; an explicit expression for r_n^die (defined implicity
#in eqn 14 of Miao et al) was derived and then inserted into eqn 15 of Miao et al to get
#the total number of cells.
sum.rn.die2 <- function(i,a,n,t,theta){
  j <- 0:(i+1)
  return ((-1)^(n+i)/a^(n-i-1)*(1/theta^(i+2) - exp(-theta*t)*sum(t^(i-j+1)/factorial(i-j+1)/theta^(j+1))))
}
sum.rn.die <- function(n,t,par){
  theta0 <- 1/par['mu0p'] + 1/par['mu0d']
  if(n == 0){
    return (N/par['mu0d']*((1-par['f0'])*par['mu0d']*(1 - exp(-t/par['mu0d'])) +
                             par['f0']/theta0*(1 - exp(-theta0*t))))
  }
  else{
    theta <- 1/par['mup'] + 1/par['mud']
    a <- theta0 - theta; b <- 1/par['mu0d'] - theta
    temp <- 2*N*par['f0']/par['mu0p']/par['mud']*(2*par['f1']/par['mup'])^(n-1)
    temp2 <- (-1)^(n-1)/a^n*((1 - exp(-theta*t))/theta + (exp(-theta0*t) - 1)/theta0) +
      sum(sapply(0:(n-2),sum.rn.die2,a,n,t,theta))*(n > 1)
    return (temp*temp2)
  }
}

#this function implements eqn (15) of the Miao et al by calling the other functions
#defined above
sum.cells <- function(n,t,par){
  return(ifelse(n == 0, N - sum.rn.div(n,t,par) - sum.rn.die(n,t,par),
                2*sum.rn.div(n-1,t,par) - sum.rn.div(n,t,par) - sum.rn.die(n,t,par)))
}
temp <- rep(0:ngen,each=length(t_pts))
sim.data <- mapply(sum.cells,split(temp,1:length(temp)),split(rep(t_pts,ngen+1),1:length(temp)),rep(list(parmsc),length(temp)))
sim.data <- matrix(sim.data,ncol=ngen+1,byrow=FALSE)


###############################Fitting the models to data#################################



sum.cols <- function(idx, data) return (rowSums(data[,idx]))

###Read the three data set and reshape them to put them in the same shape with the simulated data################################

raw.data1 <- read.table("~/Dropbox/Aims-project/Data/datasets/agent_based_sim_output_1.txt",header = TRUE)

real.data1 <- matrix(rep(0, (length(t_pts)*(ngen+1))), nrow = length(t_pts))
real.data1[1,1] <- init.vec[1]
real.data1[-1,]  <- sapply(split(2:43,rep(1:21,each=2)),sum.cols,raw.data1)

raw.data2 <- read.table("~/Dropbox/Data/datasets/agent_based_sim_output_3.txt",header = TRUE)

real.data2 <- matrix(rep(0, (length(t_pts)*(ngen+1))), nrow = length(t_pts))
real.data2[1,1] <- init.vec[1]
real.data2[-1,]  <- sapply(split(2:43,rep(1:21,each=2)),sum.cols,raw.data2)

raw.data3 <- read.table("~/Dropbox/Data/datasets/agent_based_sim_output_4.txt",header = TRUE)

real.data3 <- matrix(rep(0, (length(t_pts)*(ngen+1))), nrow = length(t_pts))
real.data3[1,1] <- init.vec[1]
real.data3[-1,]  <- sapply(split(2:43,rep(1:21,each=2)),sum.cols,raw.data3)

View(real.data1)
View(real.data2)
View(real.data3)

real.data<-(real.data1+ real.data2+real.data3)/3
View(real.data)






######Plot from genearation 0 up to generation 8 the simulate data VS the real data#######

plot(sim.data[,1],xlab = "time points",ylab = "Number of cell Gen  0",type="o", col="magenta",pch=16)
lines(real.data[,1], type="o", col="blue",pch=16) 
grid(10,10,lwd=1)
#legend("topright",legend=c("Real Data", "SM Data", "Cyton Data"),col= c("blue","red","cyan"),yjust=1,lty=1,bty="n",cex=1, inset=0)

plot(sim.data[,2],xlab = "time points",ylab = "Number of cell Gen 1",type="o", col="magenta",pch=16)
lines(real.data[,2], type="o", col="blue",pch=16)
grid(10,10,lwd=1)
#legend("topright",legend=c("Real Data", "SM Data", "Cyton Data"),col= c("blue","red","cyan"),yjust=1,lty=1,bty="n",cex=1, inset=0)

plot(real.data[,3],xlab = "time points",ylab = "Number of cell Gen 2",type="o", col="blue",pch=16)
lines(sim.data[,3], type="o", col="magenta",pch=16)
grid(10,10,lwd=1)
#legend("topright",legend=c("Real Data", "Cyton Data"),col= c("blue","magenta"),yjust=1,lty=1,bty="n",cex=1, inset=0)

plot(real.data[,4],xlab = "time points",ylab = "Number of cell Gen 3",type="o", col="blue",pch=16)
lines(sim.data[,4], type="o", col="magenta",pch=16)
grid(10,10,lwd=1)
#legend("topright",legend=c("Real Data",  "Cyton Data"),col= c("blue","magenta"),yjust=1,lty=1,bty="n",cex=1, inset=0)

plot(real.data[,5],xlab = "time points",ylab = "Number of cell Gen 4",type="o", col="blue",pch=16)
lines(sim.data[,5], type="o", col="magenta",pch=16)
grid(10,10,lwd=1)
#legend("topright",legend=c("Real Data", "Cyton Data"),col= c("blue","magenta"),yjust=1,lty=1,bty="n",cex=1, inset=0)

plot(real.data[,6],xlab = "time points",ylab = "Number of cell Gen 5",type="o", col="blue",pch=16)
lines(sim.data[,6], type="o", col="magenta",pch=16)
grid(10,10,lwd=1)
#legend("topright",legend=c("Real Data", "Cyton Data"),col= c("blue","magenta"),yjust=1,lty=1,bty="n",cex=1, inset=0)

plot(real.data[,7],xlab = "time points",ylab = "Number of cell Gen 6",type="o", col="blue",pch=16)
lines(sim.data[,7], type="o", col="magenta",pch=16)
grid(10,10,lwd=1)
#legend("topright",legend=c("Real Data", "Cyton Data"),col= c("blue","magenta"),yjust=1,lty=1,bty="n",cex=1, inset=0)

plot(real.data[,8],xlab = "time points",ylab = "Number of cell Gen 7",type="o", col="blue",pch=16)
lines(sim.data[,8], type="o", col="magenta",pch=16)
grid(10,10,lwd=1)
#legend("topright",legend=c("Real Data", "Cyton Data"),col= c("blue","magenta"),yjust=1,lty=1,bty="n",cex=1, inset=0)

plot(real.data[,9],xlab = "time points",ylab = "Number of cell Gen 8",type="o", col="blue",pch=16)
lines(sim.data[,9], type="o", col="magenta",pch=16)
grid(10,10,lwd=1)
#legend("topright",legend=c("Real Data", "Cyton Data"),col= c("blue","magenta"),yjust=1,lty=1,bty="n",cex=1, inset=0)




























