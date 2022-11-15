#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>


class Parameters
{
    public:
        // FGE loss rate
        double gamma[2][2] = {{0.0,0.0},{0.0,0.0}};

        // force of infection
        double psi[2] = {0.0,0.0};

        // fecundity of individuals infected with 
        // G1 or G2
        double F[2] = {0.0,0.0};
        double FGB = 0.0;
        double pi = 1.0;
        double c = 1.0;
        double kappa = 0.001;
        double sigma = 0.0;
        double dS[2] = {0.0,0.0};
        double dI[2][2] = {{0.0,0.0},{0.0,0.0}};
        double dBG[2] = {0.0,0.0};

        std::string base_name = "sim_fixed_resistance";

        // initial population size
        // S_p, S_c, I_pg1, I_pg2, I_cg1, I_cg2
        double init_popsize[2] = {100,100};
        double init_popsize_infected[2][2] = {{1,1},{1,1}};
        double init_popsize_superinfected[2] = {0,0};

        double eul = 0.001;

        long unsigned max_ecol_time = 5e07;

        long unsigned print_interval = 1000;

        double vanish_threshold = 1e-07;

        // whether there should be demographical feedback in the model
        bool demog_feedback = false;
};

#endif
