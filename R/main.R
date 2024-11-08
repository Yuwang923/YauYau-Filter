# Load Rcpp and compile the C++ code
library(Rcpp)
library(ggplot2)
library(gridExtra)

Rcpp::sourceCpp("YauYauFilter/Simulate_State_Obser.cpp")
source("YauYauFilter/plot_State.R")
source("YauYauFilter/plot_Obser.R")

# Define parameters
Dim <- 3
T <- 20
Dt <- 0.001
Dtau <- 5 * Dt
Nt <- as.integer(Dtau / Dt)
Ntau <- as.integer(T / Dtau)
NtNtau <- as.integer(T / Dt)

# Define functions f and h

f <- function(x) {
 return(c(cos(x[1]),
          cos(x[2]),
          cos(x[3])
          ))
}

h <- function(x) {
 return(c(x[2]^3-x[1]^2,
          x[3]^3-x[2]^2,
          x[1]^3-x[3]^2
          ))
}

# Run the simulation
seed_value <- 1234
result <- Simulate_State_Obser(Dt, Ntau, NtNtau, f, h, Dim, seed = seed_value)
 

# Access the results
x <- result$x
y_Dt <- result$y_Dt
y_tau <- result$y_tau

# Generate the plots
plot_State(x)
plot_Obser(y_tau)