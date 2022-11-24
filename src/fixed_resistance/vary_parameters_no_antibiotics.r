#!/usr/bin/env Rscript

library("simulation.utils")
library("tibble")
library("utils")

# script to permutate parameter combinations
# that will be fed to the numerical simulation of
# the ODE in C++

# make a 1d list with default parameters in the correct order
# of columns. Later when we specify (subsets of) other parameters
# this will be overwritten for some of those parameters
default_params <- list(
        gamma_CG1=1 # plasmid loss rates
        ,gamma_CG2=1
        ,gamma_PG1=1
        ,gamma_PG2=1
        ,psiG1=10 # force of infection
        ,psiG2=10
        ,FG1=1 # fecundities of the infected hosts
        ,FG2=1
        ,dSP=1 # death rates susceptibe
        ,dSC=1
        ,dIPG1=1 # death rates infecteds
        ,dIPG2=1
        ,dICG1=1
        ,dICG2=1
        ,init_popsize_P=30   # initial population sizes
        ,init_popsize_C=30
        ,init_popsize_PG1=1
        ,init_popsize_PG2=1
        ,init_popsize_PpG1G2=0
        ,init_popsize_PcG1G2=0 
        ,pival = 1.0  # resistance of hosts to G1
        ,c = 0.02 # costs of resistance
        ,kappa = 0.0001     # density dependence
        ,sigma = 0.0   # prob of superinfection
        ,eul = 0.0001  # increment to solve differential equations 
        ,demog_feedback = c(1)
        )

# parameters with antibiotics
FG1_antibiotic <- seq(2,10,0.5)
d_susc <- 5
d_single <- d_susc
d_double <- 1

antibiotics_params <- list(
        FG1=FG1_antibiotic
        ,FG2=10
        ,dSP=d_susc # death rates susceptibe
        ,dSC=d_susc
        ,dIPG1=d_single # death rates infecteds
        ,dIPG2=d_double
        ,dICG1=d_single
        ,dICG2=d_double)

antibiotics_params <- modifyList(default_params, antibiotics_params)

d_susc <- 1
d_single <- 2
d_double <- 2
no_antibiotics_params <- list(
        FG1=1
        ,FG2=1
        ,dSP=d_susc # death rates susceptibe
        ,dSC=d_susc
        ,dIPG1=d_single # death rates infecteds
        ,dIPG2=d_double
        ,dICG1=d_single
        ,dICG2=d_double)

# start the full list of parameters
no_antibiotics_params <- modifyList(default_params, no_antibiotics_params)

# we just want G1 and G2 to vary initial popsizes
# not all combinations of CG1 and CG2 and PG1 and PG2
# add these columsn later after dataframe expand.grid
#params$init_popsize_CG1 <- 0
#params$init_popsize_CG2 <- 1

all.params <- rbind(
        as.data.frame(expand.grid(antibiotics_params))
        ,as.data.frame(expand.grid(no_antibiotics_params))
        )

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
                ,output_file_prefix="output_antibiotic_vs_no_antibiotic"
                ,n_replicates = 1)
