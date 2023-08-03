
# Script for generating simulations
#
# - generate matrix based on Ising probabilities
# - insert criteria for choosing the number of iterations
# - create a matrix with a fixed border, as in the Pantanal analysis
# - New simulations for complete 2nd order and incomplete 2nd order
# - Random selection of the evaluated pixel, not linear


#### Run functions to generate samples ####

setwd(directory)
source("pcn_algorithm.R")  


#### Generate matrices of complete order ####


# Auxiliary function
rotate <- function(x) apply(x, 2, rev)

# Define the structure of the dependency tree
# What is the order of dependence
ord <- 2
# How many squares in each order
nviz <- NULL
for (o in 1:ord) {
  config <- 2*(2*o+1) + 2*(2*o-1)
  nviz <- c(nviz, config)
}
# Possible values of s_i - sum of blacks (+1) and whites (-1) - given the order
valor1 <- seq(from = (sum(nviz[1:ord])), to = -(sum(nviz[1:ord])), by = -2)
valor1 # from all blacks to 0 blacks

valor <- (sum(nviz[1:ord]):0) # number of blacks in the neighborhood
valor # from all blacks to 0 blacks

# Note: the order of "valor1" and "valor" must be the same!!! from max number of blacks to min.

# Probability vector (Ising probability). prob (site i being black | neighborhood)
beta <- 0.05
prob <- 1 / (1 + exp(-2 * beta * valor1))  # Ising formula = 1 / (1 + exp(-2 * beta * s))

# Simulate COMPLETE order matrix!
genlattice <- function(valor, prob, ord, n, v1 = 4, v2 = 8, iter = 10,
                       ini = 0, jump = 10, D = trunc(log(n / 2))) {

  # ... [The function code continues here]

}

# ... [The code continues with other functions and simulations]

# Generate matrices of variable order
genmix <- function(valor1 = valor1, prob1 = prob1, ord, n, v1 = 4, v2 = 8, iter = 10,
                   ini = 0, jump = 10, D = trunc(log(n / 2))) {

  # ... [The function code continues here]

}

# ... [The code continues with more simulations and calculations]
