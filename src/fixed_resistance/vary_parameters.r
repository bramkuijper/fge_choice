#!/usr/bin/env Rscript

library("simulation.utils")

# start a list of parameters
params <- list()

params$gamma_CG1 = 1 
params$gamma_CG2 = 1  
params$gamma_PG1 = 1  
params$gamma_PG2 = 1 

params$psiG1 = 1
params$psiG2 = 1

params$FG1 = 1
params$FG2 = 1

params$dSP = 5
params$dSC = 5 

params$dIPG1 = 5
params$dIPG2 = 1
params$dICG1 = 5
params$dICG2 = 1

params$init_popsize_P <- 100
params$init_popsize_C <- 100
params$init_popsize_PG1 <- 10
params$init_popsize_PG2 <- 1
params$init_popsize_CG1 <- 10
params$init_popsize_CG2 <- 1

params$pi = 1.0
params$c = 0.25
params$kappa = 0.001
params$eul = 0.001
params$demog_feedback = c(1)

all.params <- expand.grid(params)


make.batch.file(
                parameter_list=all.params
                ,executable_path="./solver.exe"
                ,n_replicates = 1)

