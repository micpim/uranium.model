# uranium.model
model code


CODE:
#install.packages("deSolve")
library(deSolve) #load the deSolve package

## how many iterations ##
nreps <- 500

# length of model run #
duration <-1.8e6 #duration of model run
dt <- 1000 #time step
t.anox.start <- 472000 #time of the onset of perturbation

#set distributions for drawing values for parameters to be estimated
D.other <- runif(n=nreps, min = -0.05, max = 0.05) #fractionation in oxygenated environments
D.anox <- runif(n=nreps, min = 0.2, max = 1.2) #fractionation in anoxic environments
f.anox.b <- runif(n=nreps, min = 0.0, max = 0.01) #anoxic seafloor fraction before/after perturbation 
f.anox.p <- runif(n=nreps, min = 0.0, max = 0.3) #anoxic seafloor fraction during perturbation
dur.p <- runif(n=nreps, min = 0, max = duration-t.anox.start) #duration of perturbation

#values to be used in setting up the model
t.anox.end <- t.anox.start + dur.p #time of the end of the perturbation
steps <- duration/dt #number of time steps in the model

##set up data frame to store output from the model
U.ocean <- data.frame(matrix(nrow=steps+1, ncol=nreps)) #matrix of ocean uranium concentrations
U238.ocean <- data.frame(matrix(nrow=steps+1, ncol=nreps)) #matrix of uranium isotope compositions
params <- data.frame(matrix(nrow=nreps, ncol=5))
names(params) <- c("D.other", "D.anox", "f.anox.b", "f.anox.p", "dur.p")


##set up repetitions of model runs using for loop##
for (i in 1:nreps) {
print(i)

#uranium cycle model#
model <- function (time, y, parms) {
  with(as.list(c(y, parms)), {
    f.anox <- f.an(time)
    dU = J.riv - k.anox*f.anox*y[1] - k.other*(1-f.anox)*y[1]
    dU238 = (J.riv*(deltaU.riv - U238) - k.anox*f.anox*y[1]*D.anox[i] - k.other*(1-f.anox)*y[1]*D.other[i])/y[1]
    list(c(dU, dU238))
  })
}

y <- c(U = 1.96e13, #mol U (Morford and Emerson 1999)
       U238 = -0.165) #Lau et al 2016

parms <- c(J.riv = 0.4e8, #Lau et al 2016
           k.anox = 1.45772594752187e-4, #Lau et al 2016, calculated
           k.other = 1.73834440079269e-6, #Lau et al 2016, to reach steady state
           deltaU.riv=-0.05) #Lau et al 2016, calculated steady state input to produce d238U ~ -0.165

times <- seq(0, duration, dt)

forceanox=c(rep(f.anox.b[i],(t.anox.start/dt+1)), rep(f.anox.p[i], #set up anoxia forcing function
          (t.anox.end[i]-t.anox.start)/dt), rep(f.anox.b[i], 
          (duration-t.anox.end[i])/dt+1))
anox.event=data.frame(times,forceanox)  #put the anoxia forcing into data frame with time
f.an = approxfun(x=anox.event[,1],y=anox.event[,2],method="linear",rule=2) #function for calculating f.anox in model

out <- data.frame(ode(y, times, model, parms))

## need to store output from the model and parameter values in a data frame
U.ocean[,i] <- log10(out[,2]) #put uranium concentrations from model into data frame
U238.ocean[,i] <- out[,3] #put uranium isotopes from model into data frame

params[i,1] <- D.other[i] #store value for D.other, fractionation in oxic environments
params[i,2] <- D.anox[i] #store value for D.anox, fractionation in anoxic environments
params[i,3] <- f.anox.b[i] #store value for f.anox.b, background anoxia level
params[i,4] <- f.anox.p[i] #store value for f.anox.p, perturbation anoxia level
params[i,5] <- dur.p[i] #store value for dur.p, duration of perturbation

}

#add time variable to model output for matching with data
U.ocean$mod.age <- times
U238.ocean$mod.age <- times
