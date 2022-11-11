#!/usr/bin/env Rscript

library("simulation.utils")
library("tidyverse")

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

params$FG1 = 6 #seq(0.5,10,length.out=20)
params$FG2 = 10

d_overall <- 1

params$dSP = d_overall
params$dSC = d_overall

params$dIPG1 = d_overall
params$dIPG2 = 1
params$dICG1 = d_overall
params$dICG2 = 1

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
params$c = seq(0,1,length.out=20)
params$kappa = seq(1.0^(-8),1.0^(-4),length.out=20)
params$sigma = 0.0
params$eul = 0.0001
params$demog_feedback = c(1)
params$d_vary <- 1
params$psi_vary <- 20

all.params <- as.data.frame(expand.grid(params))

all.params$dSP <- all.params$d_vary
all.params$dSC <- all.params$d_vary
all.params$dIPG1 <- all.params$d_vary
all.params$dICG1 <- all.params$d_vary

all.params$psiG1 <- all.params$psi_vary
all.params$psiG2 <- all.params$psi_vary

# we need to assign the same initial population sizes to CG1 and CG2
all.params <- add_column(all.params, 
        init_popsize_CG1 = (1.0 - all.params$pi) * all.params$init_popsize_PG1, 
        init_popsize_CG2 = all.params$init_popsize_PG2, 
        .after="init_popsize_PG2")

# remove the d_vary column
all.params <- mutate(all.params
        ,d_vary = NULL
        ,psi_vary = NULL)


write_delim(x=all.params,delim=";",file="all_params.csv")

# if there is no demographic feedback and you vary infections
# of G1 and G2, then you also need to set psi's to 0
# of those FGEs whose population sizes are set to 0

max.psi <- max(all.params[,c("psiG1","psiG2")])
all.params[all.params$demog_feedback == 0,c("psiG1","psiG2")] <- all.params[all.params$demog_feedback == 0,c("init_popsize_CG1","init_popsize_CG2")]
all.params[all.params$demog_feedback == 0,c("psiG1","psiG2")] <- all.params[all.params$demog_feedback == 0,c("psiG1","psiG2")] * max.psi

make.batch.file(
                parameter_list=all.params
                ,executable_path="./solver.exe"
                ,output_file_prefix="competition_output"
                ,n_replicates = 1)
