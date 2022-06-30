#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

// Parameters object with default parameters
class Parameters
{
    public:
        // births of individuals and their phages
        // bm, bG1m, bG2m, bG1G2m
        double b[4] = {0.0,0.0,0.0,0.0};

        double d[4] = {0.0,0.0,0.0,0.0};
       
        // density dependent growth parameter
        double kappa = 1.0/1e05;

        // gamma_G1, gamma_G2: plasmid loss rates
        double gamma[4] = {0.0,0.0,0.0,0.0};

        // psi_G1, psi_G2: plasmid infection strengths
        double psi[4] = {0.0,0.0,0.0,0.0};

        // initial value of pi
        double pi_init = 0.0;

        // sigma, susceptibility to co-infection by multiple plasmids
        double sigma = 0.0;

        // initial values of the population level parameters
        double init_vals[4] = {100,1,1,0};

        std::string base_name = {"sim_plasmid_choice"};

        unsigned output_time_interval = 10;

        long unsigned max_time = 100;

        // Euler's constant
        double eul = 0.01;
};


#endif
