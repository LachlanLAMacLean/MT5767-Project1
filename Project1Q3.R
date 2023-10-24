install.packages("remotes")
library(remotes)
remotes::install_github("chrissuthy/statsecol")
library(statsecol)
str(wildebeest)
install.packages("gridExtra")
library(gridExtra)


alpha_values <- data.frame(matrix(0, nrow = 2, ncol = 3))
rownames(alpha_values ) = c('r ~ Rt', 'r ~ Rt-1')
colnames(alpha_values )= c('\u03B1 1' , '\u03B1 2', 'AIC')

#---- rt ~ Rt ----
rainr_t <- function(pars, years, removals, Nhat, SEhat, rain, t){
  N0 <- exp(pars[1])
  beta0 <- pars[2]   
  beta1 <- pars[3]   
  k <- exp(pars[4])
  
  N <- numeric(years)
  r <- numeric(years)
  N[1] <- N0
  r[1] <- NA           
  
  #generate population dynamics:
  for(i in 2:years){
    if (t == 't-1'){
      r[i] <- exp(beta0 + beta1*rain[(i-1)]) 
    }
    else{
      r[i] <- exp(beta0 + beta1*rain[i])  
    }
    N[i] = N[i-1] + r[i] * N[i-1] * (1-N[i-1]/k) - removals[i-1]
  }
 
  negloglik <- -sum(dnorm(Nhat,N,SEhat, log=TRUE), na.rm=TRUE)
  
  return(negloglik)   #return the negative log likelihood
}

N0 <- log(0.1)
beta0 <- log(1.5)
beta1 <-  0.2 
k <- log(1.5)
parsr <- c(N0,beta0,beta1,k)

optimised_values <- optim(parsr, 
                   fn = rainr_t,
                   years = nrow(wildebeest), 
                   removals = wildebeest$Catch,
                   Nhat = wildebeest$Nhat,
                   SEhat = wildebeest$sehat,
                   rain = wildebeest$rain,
                   t = 't')

# set up parameters using the optimised values above
N0 <- exp(optimised_values$par[1])
beta0 <- optimised_values$par[2]
beta1 <- optimised_values$par[3]
k <- exp(optimised_values$par[4])
pars <- c(N0,beta0,beta1,k)
N <- numeric(nrow(wildebeest))
r <- numeric(nrow(wildebeest))

alpha_values[1,1] = beta0
alpha_values[1,2] = beta1
alpha_values[1,3] <- 2 * optimised_values$value + 2*(length(optimised_values$par))
BIC1 <- 2 * optimised_values$value + log(12) * (length(optimised_values$par))

#first year
N[1] <- N0
r[1] <- NA #1st K not in the model

#subsequent years
for(i in 2:nrow(wildebeest)){
  r[i] <- exp(beta0 + beta1*wildebeest$rain[i])
  N[i] = N[i-1] + r[i] * N[i-1] * (1-N[i-1]/k) - wildebeest$Catch[i-1]
}

tmp_wilde <- data.frame(Nhat = wildebeest$Nhat,
                        Nproj = N,
                        Year = wildebeest$year,
                        lci = wildebeest$lci,
                        uci = wildebeest$uci)

#plot the projections and the estimates
plot1 <- ggplot(tmp_wilde, aes(x=Year, y=Nproj)) +
  geom_errorbar(aes(ymin=lci,ymax=uci), width=0, color="grey") +
  geom_point(aes(x=Year,y=Nhat), size=3) +
  geom_line(aes(x=Year,y=Nproj),color="blue", size=1) +
  ylim(0,2.1) + ylab("Abundance (millions)") +
  labs(title = expression("Model with r dependent on R"[t])) +
  theme_bw() +
  theme(aspect.ratio = 1)

#---- rt ~ Rt-1 ----
N0 <- log(0.1)
beta0 <- log(1.5)
beta1 <-  0.2 
k <- log(1.5)
parsr <- c(N0,optimised_values$par[2],optimised_values$par[3],optimised_values$par[4])

optimised_values <- optim(parsr, 
                          fn = rainr_t,
                          years = nrow(wildebeest), 
                          removals = wildebeest$Catch,
                          Nhat = wildebeest$Nhat,
                          SEhat = wildebeest$sehat,
                          rain = wildebeest$rain,
                          t = 't-1')

# set up parameters using the optimised values above
N0 <- exp(optimised_values$par[1])
beta0 <- optimised_values$par[2]
beta1 <- optimised_values$par[3]
k <- exp(optimised_values$par[4])
N <- numeric(nrow(wildebeest))
r <- numeric(nrow(wildebeest))

alpha_values[2,1] = beta0
alpha_values[2,2] = beta1
alpha_values[2,3] <- 2 * optimised_values$value + 2*(length(optimised_values$par))
BIC2 <- 2 * optimised_values$value + log(12) * (length(optimised_values$par))


#first year
N[1] <- N0
r[1] <- NA #1st K not in the model

#subsequent years
for(i in 2:nrow(wildebeest)){
  r[i] <- exp(beta0 + beta1*wildebeest$rain[i-1])
  N[i] = N[i-1] + r[i] * N[i-1] * (1-N[i-1]/k) - wildebeest$Catch[i-1]
}

tmp_wilde <- data.frame(Nhat = wildebeest$Nhat,
                        Nproj = N,
                        Year = wildebeest$year,
                        lci = wildebeest$lci,
                        uci = wildebeest$uci)

#plot the projections and the estimates
plot2 <-ggplot(tmp_wilde, aes(x=Year, y=Nproj)) +
  geom_errorbar(aes(ymin=lci,ymax=uci), width=0, color="grey") +
  geom_point(aes(x=Year,y=Nhat), size=3) +
  geom_line(aes(x=Year,y=Nproj),color="blue", size=1) +
  ylim(0,2.1) + ylab("Abundance (millions)") +
  labs(title = "Model with r dependent on R[t-1]") +
  theme_bw() +
  theme(aspect.ratio = 1)

grid.arrange(plot1, plot2, ncol = 2)
#---- k ~ Rt ----

beta_values <- data.frame(matrix(0, ncol = 3, nrow = 2))
rownames(beta_values ) = c('k ~ Rt', 'k ~ Rt-1')
colnames(beta_values )= c('\u03B2_0' , '\u03B2_1', 'AIC')

rainK_t <- function(pars, years, removals, Nhat, SEhat, rain, t){
  
  N0 <- exp(pars[1])
  r <- exp(pars[2])
  beta0 <- pars[3]   #not transformed |
  beta1 <- pars[4]   #not transformed |-> note extra parameter now
  N <- numeric(years)
  k <- numeric(years)
  N[1] <- N0
  k[1] <- NA          #1st K not in the model 
  
  for(i in 2:years){    #generate population dynamics:
    if (t == 't-1'){
      k[i] <- exp(beta0 + beta1*rain[(i-1)]) #link fn of the linear predictor 
    }
    else{
      k[i] <- exp(beta0 + beta1*rain[i]) #link fn of the linear predictor 
    }
    N[i] = N[i-1] + (r * N[i-1] * (1-N[i-1]/k[i])) - removals[i-1]
  }
  negloglik <- -sum(dnorm(Nhat,N,SEhat, log=TRUE), na.rm=TRUE)
  
  return(negloglik) 
}

N0 <- log(0.1)
r <- log(0.25)
beta0 <-  log(0.5)   
beta1 <-  log(0.5)          

parsk <- c(N0,r,beta0,beta1)
fit_rainK <- optim(parsk,
                   fn = rainK_t,
                   years = nrow(wildebeest), 
                   removals = wildebeest$Catch,
                   Nhat = wildebeest$Nhat,
                   SEhat = wildebeest$sehat,
                   rain = wildebeest$rain,
                   t = 't')

# set up parameters
N0 <- exp(fit_rainK$par[1])
r <- exp(fit_rainK$par[2])
beta0 <- fit_rainK$par[3]
beta1 <- fit_rainK$par[4]
pars <- c(N0, r, beta0,beta1)

beta_values[1,1] = beta0
beta_values[1,2] = beta1
N <- numeric(nrow(wildebeest))
k <- numeric(nrow(wildebeest))

beta_values[1,3] <- 2 * fit_rainK$value + 2*length(fit_rainK$par)
BIC3 <- 2 * fit_rainK$value + log(12)*length(fit_rainK$par)

#first year
N[1] <- N0
k[1] <- NA #1st K not in the model

#subsequent years
for(i in 2:nrow(wildebeest)){
  k[i] <- exp(beta0 + beta1*wildebeest$rain[i])
  N[i] = N[i-1] + r * N[i-1] * (1-N[i-1]/k[i]) - wildebeest$Catch[i-1]
}

tmp_wilde <- data.frame(Nhat = wildebeest$Nhat,
                        Nproj = N,
                        Year = wildebeest$year,
                        lci = wildebeest$lci,
                        uci = wildebeest$uci)

#plot the projections and the estimates
plot3 <- ggplot(tmp_wilde, aes(x=Year, y=Nproj)) +
  geom_errorbar(aes(ymin=lci,ymax=uci), width=0, color="grey") +
  geom_point(aes(x=Year,y=Nhat), size=3) +
  geom_line(aes(x=Year,y=Nproj),color="blue", size=1) +
  ylim(0,2.1) + ylab("Abundance (millions)") +
  labs(title = "Model with k dependent on Rt") +
  theme_bw() +
  theme(aspect.ratio = 1)

#---- Kt ~ Rt-1 ----

fit_rainK <- optim(parsk,
                   fn = rainK_t,
                   years = nrow(wildebeest), 
                   removals = wildebeest$Catch,
                   Nhat = wildebeest$Nhat,
                   SEhat = wildebeest$sehat,
                   rain = wildebeest$rain,
                   t = 't-1')

# set up parameters
N0 <- exp(fit_rainK$par[1])
r <- exp(fit_rainK$par[2])
beta0 <- fit_rainK$par[3]
beta1 <- fit_rainK$par[4]
pars <- c(N0, r, beta0,beta1)

beta_values[2,1] = beta0
beta_values[2,2] = beta1
N <- numeric(nrow(wildebeest))
k <- numeric(nrow(wildebeest))
beta_values[2,3] <- 2 * fit_rainK$value + 2*length(fit_rainK$par)
BIC4 <- 2 * fit_rainK$value + log(12)*length(fit_rainK$par)

N[1] <- N0
k[1] <- NA #1st K not in the model

#subsequent years
for(i in 2:nrow(wildebeest)){
  k[i] <- exp(beta0 + beta1*wildebeest$rain[i])
  N[i] = N[i-1] + r * N[i-1] * (1-N[i-1]/k[i]) - wildebeest$Catch[i-1]
}

tmp_wilde <- data.frame(Nhat = wildebeest$Nhat,
                        Nproj = N,
                        Year = wildebeest$year,
                        lci = wildebeest$lci,
                        uci = wildebeest$uci)

#plot the projections and the estimates
plot4 <- ggplot(tmp_wilde, aes(x=Year, y=Nproj)) +
  geom_errorbar(aes(ymin=lci,ymax=uci), width=0, color="grey") +
  geom_point(aes(x=Year,y=Nhat), size=3) +
  geom_line(aes(x=Year,y=Nproj),color="blue", size=1) +
  ylim(0,2.1) + ylab("Abundance (millions)") +
  labs(title = "Model with k dependent on Rt-1") +
  theme_bw()+
  theme(aspect.ratio = 1)


grid.arrange(plot3, plot4, ncol = 2)

print(beta_values)
