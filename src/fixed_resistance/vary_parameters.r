#!/usr/bin/env Rscript

library("simulation.utils")
library("tibble")

# script to permutate parameter combinations
# that will be fed to the numerical simulation of
# the ODE in C++

# start a list of parameters
params <- list()

params$gamma_CG1 = 1 
params$gamma_CG2 = 1  
params$gamma_PG1 = 1  
params$gamma_PG2 = 1 

params$psiG1 = 10
params$psiG2 = 10

params$FG1 = 2
params$FG2 = 10

d_susceptible <- 5
d_double <- 1

params$dSP = d_susceptible
params$dSC = d_susceptible

params$dIPG1 = d_susceptible
params$dIPG2 = d_double
params$dICG1 = d_susceptible
params$dICG2 = d_double

params$init_popsize_P <- 30
params$init_popsize_C <- 30

params$init_popsize_PG1 <- c(0,1)
params$init_popsize_PG2 <- c(0,1)
params$init_popsize_PpG1G2 <- 0
params$init_popsize_PcG1G2 <- 0

# we just want G1 and G2 to vary initial popsizes
# not all combinations of CG1 and CG2 and PG1 and PG2
# add these columsn later after dataframe expand.grid
#params$init_popsize_CG1 <- 0
#params$init_popsize_CG2 <- 1

params$pi = 1.0
params$c = 0.02
params$kappa = 0.0001
params$sigma = 0.0
params$eul = 0.0001
params$demog_feedback = c(1)

all.params <- as.data.frame(expand.grid(params))

all.params$init_popsize_P <- 60 - all.params$init_popsize_C

# we need to assign the same initial population sizes to CG1 and CG2
all.params <- add_column(all.params, 
        init_popsize_CG1 = 0, 
        init_popsize_CG2 = all.params$init_popsize_PG2, 
        .after="init_popsize_PG2")

# if there is no demographic feedback and you vary infections
# of G1 and G2, then you also need to set psi's to 0
# of those FGEs whose population sizes are set to 0

max.psi <- max(all.params[,c("psiG1","psiG2")])
all.params[all.params$demog_feedback == 0,c("psiG1","psiG2")] <- all.params[all.params$demog_feedback == 0,c("init_popsize_CG1","init_popsize_CG2")]
all.params[all.params$demog_feedback == 0,c("psiG1","psiG2")] <- all.params[all.params$demog_feedback == 0,c("psiG1","psiG2")] * max.psi

make.batch.file(
                parameter_list=all.params
                ,executable_path="./solver.exe"
                ,output_file_prefix="output"
                ,n_replicates = 1)
