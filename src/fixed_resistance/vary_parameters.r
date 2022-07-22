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

params$FG1 = seq(1,10,0.2)
params$FG2 = 10

params$dSP = 1
params$dSC = 1

params$dIPG1 = 1
params$dIPG2 = 1
params$dICG1 = 1
params$dICG2 = 1

params$init_popsize_P <- 100
params$init_popsize_C <- 100
params$init_popsize_PG1 <- 1
params$init_popsize_PG2 <- 1
params$init_popsize_CG1 <- 1 
params$init_popsize_CG2 <- 1

params$pi = c(0.5,1.0)
params$c = seq(0,0.1,0.0005)
params$kappa = 0.001
params$eul = 0.001
params$demog_feedback = c(0,1)

all.params <- expand.grid(params)


make.batch.file(
                parameter_list=all.params
                ,executable_path="./solver.exe"
                ,n_replicates = 1)

