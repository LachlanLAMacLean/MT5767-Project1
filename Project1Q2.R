#here is a general solution that can be used to fit all 4 model
rain_rK_nll <- function(pars, years, removals, Nhat, SEhat, rain, type="nll"){
  
  #parameter set up
  N0 <- exp(pars[1])
  N <- numeric(years)
  r <- numeric(years)
  k <- numeric(years)
  N[1] <- N0
  r[1] <- NA
  k[1] <- NA
  
  if(length(pars)!=5){stop("par should have 5 values")}
    r[2:years] <- exp(pars[2]+pars[3]*rain[2:years])
    k[2:years] <- exp(pars[4]+pars[5]*rain[2:years])
  
  #generate population dynamics:
  for(i in 2:years){
    N[i]=N[i-1] + r[i] * N[i-1] * (1-N[i-1]/k[i]) - removals[i-1]
  }
  
  negloglik <- -sum(dnorm(Nhat,N,SEhat,log=TRUE), na.rm=TRUE)
  
  #what should be returned? 
  if(type=="nll"){  return(negloglik)}
  if(type=="proj"){ return(N)}
}


yrs <- nrow(wildebeest)
rmv <- wildebeest$Catch
Nhat <- wildebeest$Nhat
SEhat <- wildebeest$sehat
rain <- wildebeest$rain

#fit the model:
fit_4 <- optim(par = c(log(0.1),log(0.25),0,log(1.5),0), fn = rain_rK_nll, years = yrs, 
               removals = rmv, Nhat = Nhat, SEhat = SEhat, rain = rain)

#Compute the AIC:
aic4 <- 2*fit_4$value + 2*length(fit_4$par)

proj_4 <- rain_rK_nll(fit_4$par, years = yrs, removals = rmv, Nhat = Nhat, 
                      SEhat = SEhat, rain = rain, type="proj")

pred_df <- data.frame(years = wildebeest$year,
                      N = proj_4,
                      Model="Removals Last")

#modified general solution to fit the removals first model
rain_rK_nll_alt <- function(pars, years, removals, Nhat, SEhat, rain, type="nll"){
  
  #parameter set up
  N0 <- exp(pars[1])
  N <- numeric(years)
  r <- numeric(years)
  k <- numeric(years)
  N[1] <- N0 - removals[1]
  r[1] <- NA
  k[1] <- NA
  
  if(length(pars)!=5){stop("par should have 5 values")}
  r[2:years] <- exp(pars[2]+pars[3]*rain[2:years])
  k[2:years] <- exp(pars[4]+pars[5]*rain[2:years])
  
  #generate population dynamics:
  for(i in 2:years){
    NWithRemovals = N[i-1] - removals[i]
    N[i]= NWithRemovals + r[i] * NWithRemovals * (1-NWithRemovals/k[i]) 
  }
  
  negloglik <- -sum(dnorm(Nhat,N,SEhat,log=TRUE), na.rm=TRUE)
  
  #what should be returned? 
  if(type=="nll"){  return(negloglik)}
  if(type=="proj"){ return(N)}
}

#fit the model:
fit_5 <- optim(par = c(log(0.1),log(0.25),0,log(1.5),0), fn = rain_rK_nll_alt, years = yrs, 
               removals = rmv, Nhat = Nhat, SEhat = SEhat, rain = rain)

#Compute the AIC:
aic5 <- 2*fit_5$value + 2*length(fit_5$par)

proj_5 <- rain_rK_nll_alt(fit_4$par, years = yrs, removals = rmv, Nhat = Nhat, 
                      SEhat = SEhat, rain = rain, type="proj")

pred_df2 <- data.frame(years = wildebeest$year,
                      N = proj_5,
                      Model="Removals First")

pred_df_full <- rbind(pred_df, pred_df2)

gpreds <- ggplot(wildebeest, aes(x=year, y=Nhat)) +
  geom_errorbar(aes(ymin=lci,ymax=uci), width=0) +
  geom_point(size=3) +
  geom_line(data=pred_df_full, aes(x=years,y=N,color=Model,group=Model),size=0.8) +
  ylim(0,2.1) + ylab("Abundance (millions)") + xlab("Year") +
  theme_bw()
gpreds
