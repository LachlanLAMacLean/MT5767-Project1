install.packages("remotes")
library(remotes)
remotes::install_github("chrissuthy/statsecol")
library(statsecol)
str(wildebeest)
library(tidyverse)

#repurposed general solution from Practical 5 to generate original rainfall-dependent model
rain_rK_nll <- function(pars, years, removals, Nhat, SEhat, rain, type="nll"){
  
  #exponentiate values to match log-likelihood optimisation
  #set up empty parameter placeholders otherwise
  N0 <- exp(pars[1])
  N <- numeric(years)
  r <- numeric(years)
  k <- numeric(years)
  N[1] <- N0
  r[1] <- NA
  k[1] <- NA
  
  #keep a check for model sanity to prevent function errors
  if(length(pars)!=5){stop("par should have 5 values")}
    r[2:years] <- exp(pars[2]+pars[3]*rain[2:years])
    k[2:years] <- exp(pars[4]+pars[5]*rain[2:years])
  
  #generate population dynamics
  for(i in 2:years){
    N[i]=N[i-1] + r[i] * N[i-1] * (1-N[i-1]/k[i]) - removals[i-1]
  }
  
  #get the negative log likelihood
  negloglik <- -sum(dnorm(Nhat,N,SEhat,log=TRUE), na.rm=TRUE)
  
  #returns likelihood by default but can alternatively return Ns
  if(type=="nll"){  return(negloglik)}
  if(type=="proj"){ return(N)}
}

#assign initial values from the actual data
yrs <- nrow(wildebeest)
rmv <- wildebeest$Catch
Nhat <- wildebeest$Nhat
SEhat <- wildebeest$sehat
rain <- wildebeest$rain

#fit the model with suitable initial values
#log values to match log-likelihood
fit_4 <- optim(par = c(log(0.1),log(0.25),0,log(1.5),0), fn = rain_rK_nll, years = yrs, 
               removals = rmv, Nhat = Nhat, SEhat = SEhat, rain = rain)

#Compute the AIC
aic4 <- 2*fit_4$value + 2*length(fit_4$par)

#generate predicted Ns from the same function
proj_4 <- rain_rK_nll(fit_4$par, years = yrs, removals = rmv, Nhat = Nhat, 
                      SEhat = SEhat, rain = rain, type="proj")

#add a column for model identification (after eventual merging)
pred_df <- data.frame(years = wildebeest$year,
                      N = proj_4,
                      Model="Removals Last")

#further modify the provided function to change the order of removals
#now, removals take place before any other step
rain_rK_nll_alt <- function(pars, years, removals, Nhat, SEhat, rain, type="nll"){
  
  #parameter set up: exponentiated values and remaining empty placeholders
  N0 <- exp(pars[1])
  N <- numeric(years)
  r <- numeric(years)
  k <- numeric(years)
  N[1] <- N0 - removals[1]
  r[1] <- NA
  k[1] <- NA
  
  #keep a check for model sanity to prevent function errors
  if(length(pars)!=5){stop("par should have 5 values")}
  r[2:years] <- exp(pars[2]+pars[3]*rain[2:years])
  k[2:years] <- exp(pars[4]+pars[5]*rain[2:years])
  
  #generate population dynamics:
  for(i in 2:years){
    NWithRemovals = N[i-1] - removals[i]
    N[i]= NWithRemovals + r[i] * NWithRemovals * (1-NWithRemovals/k[i]) 
  }
  
  #get the negative log likelihood
  negloglik <- -sum(dnorm(Nhat,N,SEhat,log=TRUE), na.rm=TRUE)
  
  #returns likelihood by default but can alternatively return Ns
  if(type=="nll"){  return(negloglik)}
  if(type=="proj"){ return(N)}
}

#fit the removals first model in the same fashion using previously-defined initial parameters
fit_5 <- optim(par = c(log(0.1),log(0.25),0,log(1.5),0), fn = rain_rK_nll_alt, years = yrs, 
               removals = rmv, Nhat = Nhat, SEhat = SEhat, rain = rain)

#Compute the AIC again to compare
aic5 <- 2*fit_5$value + 2*length(fit_5$par)

#Retrieve the projected N (population sizes)
proj_5 <- rain_rK_nll_alt(fit_5$par, years = yrs, removals = rmv, Nhat = Nhat, 
                      SEhat = SEhat, rain = rain, type="proj")

#add a column for model identification in plots
pred_df2 <- data.frame(years = wildebeest$year,
                      N = proj_5,
                      Model="Removals First")

#combine both model outputs for single plot
pred_df_full <- rbind(pred_df, pred_df2)

#plot the two models for comparison
gpreds <- ggplot(wildebeest, aes(x=year, y=Nhat)) +
  geom_errorbar(aes(ymin=lci,ymax=uci), width=0) +
  geom_point(size=3) +
  geom_line(data=pred_df_full, aes(x=years,y=N,color=Model,group=Model),linewidth=0.8) +
  ylim(0,2.1) + ylab("Abundance (millions)") + xlab("Year") +
  ggtitle("Population models with varying order of removals from illegal harvesting") + 
  theme_bw()

#display plot
gpreds
