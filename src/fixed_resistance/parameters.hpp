#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>


class Parameters
{
    public:
        double gamma[2][2] = {{0.0,0.0},{0.0,0.0}};
        double psi[2] = {0.0,0.0};
        double F[2] = {0.0,0.0};
        double pi = 1.0;
        double c = 1.0;
        double kappa = 0.001;
        double dS[2] = {0.0,0.0};
        double dI[2][2] = {{0.0,0.0},{0.0,0.0}};

        std::string base_name = "sim_fixed_resistance";

        // initial population size
        // S_p, S_c, I_pg1, I_pg2, I_cg1, I_cg2
        double init_popsize[2] = {100,100};
        double init_popsize_infected[2][2] = {{1,2},{3,4}};

        double eul = 0.001;

        long unsigned max_ecol_time = 1;

        long unsigned print_interval = 10;

        double vanish_threshold = 1e-07;

        // whether there should be demographical feedback in the model
        bool demog_feedback = false;
};

#endif
