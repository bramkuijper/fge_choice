#!/usr/bin/env Rscript

library("simulation.utils")

# start a list of parameters
params <- list()

params$eul <- 0.001
params$c <- seq(0.01,0.1,0.01)
params$sigma <- 0.0

params$gamma_G1 <- c(0.0,1.0)
params$gamma_G2 <- c(0.0,1.0)

params$psi_G1 <- c(1,5,10)
params$psi_G2 <- c(1,5,10)

params$FS <- 1.0
params$FG1 <- 1.0
params$FG2 <- 1.0
params$FG1G2 <- 1.0

params$dS <- 1.0
params$dG1 <- 1.0
params$dG2 <- c(1,5,10,20,50)
params$dG1G2 <- 1.0

params$kappa <- c(0.001,0.0001,0.00001)


params$init_vals_S <- 100 
params$init_vals_G1 <- 100 
params$init_vals_G2 <- 100 
params$init_vals_G1G2 <- 100 

all.params <- expand.grid(params)


make.batch.file(
                parameter_list=all.params
                ,executable_path="./solver.exe"
                ,n_replicates = 1)

