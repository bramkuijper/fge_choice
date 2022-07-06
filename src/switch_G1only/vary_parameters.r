#!/usr/bin/env Rscript

library("simulation.utils")

# start a list of parameters
params <- list()

params$eul <- 0.01
params$c <- c(0.1,0.3,0.4)
params$sigma <- 0.0

params$gamma_G1 <- 5
params$gamma_G2 <- 5

params$psi_G1 <- 10
params$psi_G2 <- seq(0,10,0.2)

params$FS <- 1.0
params$FG1 <- 1.0
params$FG2 <- seq(1.0,10,0.2)
params$FG1G2 <- 0.0

params$dS <- 1.0
params$dG1 <- 1.0
params$dG2 <- 1.0
params$dG1G2 <- 1.0

params$kappa <- 0.001


params$init_vals_S <- 100 
params$init_vals_G1 <- 100 
params$init_vals_G2 <- 100 
params$init_vals_G1G2 <- 0

all.params <- expand.grid(params)

print(nrow(all.params))

make.batch.file(
                parameter_list=all.params
                ,executable_path="./solver.exe"
                ,n_replicates = 1)

