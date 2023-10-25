# Load relevant libraries
library(tidyverse)

# Set seed for reproducibility
set.seed(576711)

# a)

# See written report

# b)

# Write an R script to simulate dynamics from this model.
# we do this in a function.

# INPUT : 
# n0 = initial populations
# phi = survival probabilities
# rho = reproduction rates
# p = probability of detection
# nyears = number of years over which population is projected
# PROCESS : 
# calculate the stochastic population dynamics and observations for each year
# OUTPUT :
# y = array of observations for each year, across all age classes
# n = array of abundances for each year, across all age classes

BAS_stoch <- function(n0, phi, rho, p, nyears) {
  
  # initialise
  
  # matrix of all abundances, each column is a year, row is age class
  n <- matrix(data = NA, nrow = length(n0), ncol = (nyears+1))
  
  # matrix of all observations
  y <- n #replicate n
  
  # initial abundance
  # let the first set of abundances (at t = 0) be the initial population size
  n[,1] <- n0
  
  # initial observation
  # ASSUMPTION: binomial distribution as constant probability of detection (see Newman)
  y[1,1] <- rbinom(n = 1, size = n[1,1], prob = p)
  y[2,1] <- rbinom(n = 1, size = n[2,1], prob = p)
  y[3,1] <- rbinom(n = 1, size = n[3,1], prob = p)
  y[4,1] <- rbinom(n = 1, size = n[4,1], prob = p)
  
  # loop over all (other) years
  for (i in 2:(nyears+1)) {
    
    # calculate stochastic sub-processes for BAS model
    
    # Survival process
    u_1s1t <- rbinom(n = 1, size = n[1,i-1], prob = phi[1])
    u_1s2t <- rbinom(n = 1, size = n[2,i-1], prob = phi[2])
    u_1s3t <- rbinom(n = 1, size = n[3,i-1], prob = phi[3])
    u_1s4t <- 0
    
    # Ageing process
    u_2a1t <- 0
    u_2a2t <- u_1s1t
    u_2a3t <- u_1s2t
    u_2a4t <- u_1s3t + u_1s4t
    
    # Reproduction/Birth process
    u_3b1t <- rpois(n = 1, lambda = rho[1] * u_2a2t) + 
              rpois(n = 1, lambda = rho[2] * u_2a3t)
    u_3b2t <- u_2a2t
    u_3b3t <- u_2a3t
    u_3b4t <- u_2a4t
    
    # new abundances
    n[,i] <- c(u_3b1t, u_3b2t, u_3b3t, u_3b4t)
    
    # new observations
    y[1,i] <- rbinom(n = 1, size = n[1,i], prob = p)
    y[2,i] <- rbinom(n = 1, size = n[2,i], prob = p)
    y[3,i] <- rbinom(n = 1, size = n[3,i], prob = p)
    y[4,i] <- rbinom(n = 1, size = n[4,i], prob = p)
  }
  
  # return
  return(list(n=n,y=y))
  
}

# parameters from specification
phi <- c(0.45, 0.7, 0.7)
rho <- c(0.9, 1.9)
p <- 0.5
n0 <- c(150, 70, 50, 30)

# c)

# Simulate 25 years of age-specific population dynamics and observations and 
# produce an informative visualisation of the data.

# specify how many years to project over
nyears <- 25

# run function for 25 years
BAS_proj <- BAS_stoch(n0 = n0, phi = phi, rho = rho, p = p, nyears = nyears)

# visualise the data using a faceted ggplot.

# convert both outputs to dataframes
proj_n_df <- as.data.frame( cbind(0:nyears, t(BAS_proj$n), 
                    rep(2,(nyears+1))) )
proj_y_df <- as.data.frame( cbind(0:nyears, t(BAS_proj$y), 
                    rep(1,(nyears+1))) )

# rename
colnames(proj_n_df) <- c("Year", "First Year Individuals", 
                         "Second Year Individuals", 
                         "Third Year Individuals", 
                         "Fourth Year Individuals", "State")
colnames(proj_y_df) <- c("Year", "First Year Individuals", 
                         "Second Year Individuals", 
                         "Third Year Individuals", 
                         "Fourth Year Individuals", "State")

# pivot longer
proj_n_df_long <- pivot_longer(data = proj_n_df, 
                               cols = c("First Year Individuals", 
                                        "Second Year Individuals", 
                                        "Third Year Individuals", 
                                        "Fourth Year Individuals"))
proj_y_df_long <- pivot_longer(data = proj_y_df, 
                               cols = c("First Year Individuals", 
                                        "Second Year Individuals", 
                                        "Third Year Individuals", 
                                        "Fourth Year Individuals"))
proj_df_long <- rbind(proj_n_df_long, proj_y_df_long)
colnames(proj_df_long) <- c("Year", "State", "Age Class", "Abundance")

# change variable types
proj_df_long$State <- as.factor(proj_df_long$State)
proj_df_long$`Age Class` <- as.factor(proj_df_long$`Age Class`)

# Can now easily use faceted ggplot
plot_BAS <- ggplot(data = proj_df_long, aes(x = Year, y = Abundance)) + 
  geom_line(aes(colour = State), size = 1) +
  scale_colour_manual("State", breaks = c(2, 1),
                      values = c("blue", "grey"), 
                      labels = c("Dynamics","Observations")) +
  facet_wrap(~factor(`Age Class`, levels = c("First Year Individuals", 
                                             "Second Year Individuals", 
                                             "Third Year Individuals", 
                                             "Fourth Year Individuals")), 
             scales = "free_y") +
  ylab("Abundance") +
  ggtitle("Population Dynamics and Observations Projected over Time")
plot_BAS

# Extension to part b) and c)

# Create a new function like in b), but runs poisson instead of Binomial for 
# observations

# INPUT : 
# n0 = initial populations
# phi = survival probabilities
# rho = reproduction rates
# p = probability of detection
# nyears = number of years over which population is projected
# PROCESS : 
# calculate the stochastic population dynamics and observations for each year
# OUTPUT :
# y = array of observations for each year, across all age classes
# n = array of abundances for each year, across all age classes

BAS_extension <- function(n0, phi, rho, p, nyears) {
  
  # initialise
  
  # matrix of all abundances, each column is a year, row is age class
  n <- matrix(data = NA, nrow = length(n0), ncol = (nyears+1))
  
  # matrix of all observations
  y <- n #replicate n
  
  # initial abundance
  # let the first set of abundances (at t = 0) be the initial population size
  n[,1] <- n0
  
  # initial observation
  y[1,1] <- rpois(n = 1, lambda = n[1,1]*p)
  y[2,1] <- rpois(n = 1, lambda = n[2,1]*p)
  y[3,1] <- rpois(n = 1, lambda = n[3,1]*p)
  y[4,1] <- rpois(n = 1, lambda = n[4,1]*p)
  
  # loop over all (other) years
  for (i in 2:(nyears+1)) {
    
    # calculate stochastic sub-processes for BAS model
    
    # Survival process
    u_1s1t <- rbinom(n = 1, size = n[1,i-1], prob = phi[1])
    u_1s2t <- rbinom(n = 1, size = n[2,i-1], prob = phi[2])
    u_1s3t <- rbinom(n = 1, size = n[3,i-1], prob = phi[3])
    u_1s4t <- 0
    
    # Ageing process
    u_2a1t <- 0
    u_2a2t <- u_1s1t
    u_2a3t <- u_1s2t
    u_2a4t <- u_1s3t + u_1s4t
    
    # Reproduction/Birth process
    u_3b1t <- rpois(n = 1, lambda = rho[1] * u_2a2t) + 
      rpois(n = 1, lambda = rho[2] * u_2a3t)
    u_3b2t <- u_2a2t
    u_3b3t <- u_2a3t
    u_3b4t <- u_2a4t
    
    # new abundances
    n[,i] <- c(u_3b1t, u_3b2t, u_3b3t, u_3b4t)
    
    # new observations
    y[1,i] <- rpois(n = 1, lambda = n[1,i]*p)
    y[2,i] <- rpois(n = 1, lambda = n[2,i]*p)
    y[3,i] <- rpois(n = 1, lambda = n[3,i]*p)
    y[4,i] <- rpois(n = 1, lambda = n[4,i]*p)
  }
  
  # return
  return(list(n=n,y=y))
  
}

# Simulate 25 years of age-specific population dynamics and observations and 
# produce an informative visualisation of the data.

# run function for 25 years
# use same parameter values as above
BAS_proj_ext <- BAS_extension(n0 = n0, phi = phi, rho = rho, p = p, 
                              nyears = nyears)

# visualise the data using a faceted ggplot.

# convert both outputs to dataframes
proj_n_df <- as.data.frame( cbind(0:nyears, t(BAS_proj_ext$n), 
                                  rep(2,(nyears+1))) )
proj_y_df <- as.data.frame( cbind(0:nyears, t(BAS_proj_ext$y), 
                                  rep(1,(nyears+1))) )

# rename
colnames(proj_n_df) <- c("Year", "First Year Individuals", 
                         "Second Year Individuals", 
                         "Third Year Individuals", 
                         "Fourth Year Individuals", "State")
colnames(proj_y_df) <- c("Year", "First Year Individuals", 
                         "Second Year Individuals", 
                         "Third Year Individuals", 
                         "Fourth Year Individuals", "State")

# pivot longer
proj_n_df_long <- pivot_longer(data = proj_n_df, 
                               cols = c("First Year Individuals", 
                                        "Second Year Individuals", 
                                        "Third Year Individuals", 
                                        "Fourth Year Individuals"))
proj_y_df_long <- pivot_longer(data = proj_y_df, 
                               cols = c("First Year Individuals", 
                                        "Second Year Individuals", 
                                        "Third Year Individuals", 
                                        "Fourth Year Individuals"))
proj_df_long <- rbind(proj_n_df_long, proj_y_df_long)
colnames(proj_df_long) <- c("Year", "State", "Age Class", "Abundance")

# change variable types
proj_df_long$State <- as.factor(proj_df_long$State)
proj_df_long$`Age Class` <- as.factor(proj_df_long$`Age Class`)

# Can now easily use faceted ggplot
plot_BAS_ext <- ggplot(data = proj_df_long, aes(x = Year, y = Abundance)) + 
  geom_line(aes(colour = State), size = 1) +
  scale_colour_manual("State", breaks = c(2, 1),
                      values = c("blue", "grey"), 
                      labels = c("Dynamics","Observations")) +
  facet_wrap(~factor(`Age Class`, levels = c("First Year Individuals", 
                                             "Second Year Individuals", 
                                             "Third Year Individuals", 
                                             "Fourth Year Individuals")), 
             scales = "free_y") +
  ylab("Abundance") +
  ggtitle("Population Dynamics and Observations Projected over Time")
plot_BAS_ext





# # test
# for (i in 1:4) {
#   plot(0:nyears, BAS_proj$y[i,],ylab="Abundance",xlab="Years",
#        las=1, type="l", pch=16, col = "grey")
#   lines(0:nyears, BAS_proj$n[i,],ylab="Abundance",xlab="Years",
#         las=1, type="l", pch=16, col = "blue")
# }