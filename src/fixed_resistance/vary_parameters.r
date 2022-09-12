#!/usr/bin/env Rscript

library("simulation.utils")
library("tibble")

# start a list of parameters
params <- list()

params$gamma_CG1 = 1 
params$gamma_CG2 = 1  
params$gamma_PG1 = 1  
params$gamma_PG2 = 1 

params$psiG1 = 10
params$psiG2 = 10

params$FG1 = 6
params$FG2 = seq(6,10,length.out=50)

d_overall <- 5

params$dSP = d_overall
params$dSC = d_overall

params$dIPG1 = d_overall
params$dIPG2 = 1
params$dICG1 = d_overall
params$dICG2 = 1

params$init_popsize_P <- 30
params$init_popsize_C <- 30

params$init_popsize_PG1 <- c(1)
params$init_popsize_PG2 <- c(1)
params$init_popsize_PpG1G2 <- 0
params$init_popsize_PcG1G2 <- 0

# we just want G1 and G2 to vary initial popsizes
# not all combinations of CG1 and CG2 and PG1 and PG2
# add these columsn later after dataframe expand.grid
#params$init_popsize_CG1 <- 0
#params$init_popsize_CG2 <- 1

params$pi = 1.0
params$c = seq(0,0.04,length.out = 50)
params$kappa = 0.001
params$sigma = c(0,1)
params$eul = 0.0001
params$demog_feedback = c(0)

all.params <- as.data.frame(expand.grid(params))

# use tibbles' add_column function
all.params <- add_column(all.params, 
        init_popsize_CG1 = all.params$init_popsize_PG1, 
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
                ,n_replicates = 1)
