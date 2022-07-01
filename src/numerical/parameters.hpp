#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

// Parameters object with default parameters
class Parameters
{
    public:
        // constants of the fecundity function
        // bm, bG1m, bG2m, bG1G2m
        double F[4] = {1.0,0.0,0.0,0.0};

        double bmax = 50.0;
        
        // constants of the fecundity function
        // bm, bG1m, bG2m, bG1G2m
        double d[4] = {1.0,15.0,10.0,7.0};
       
        // density dependent growth parameter
        double kappa = 0.01;

        // gamma_G1, gamma_G2: plasmid loss rates
        double gamma[4] = {0.0,3.0,8.0,0.0};

        // psi_G1, psi_G2: plasmid infection strengths
        double psi[4] = {0.0,10.0,99.0,0.0};

        // initial value of pi
        double pi_init = 0.5;

        // cost of resistance
        double c = 4.0;

        // sigma, susceptibility to co-infection by multiple plasmids
        double sigma = 2.5;

        // initial values of the population level parameters
        double init_vals[4] = {100,1,1,0};

        std::string base_name = {"sim_plasmid_choice"};

        unsigned output_time_interval = 10;

        long unsigned max_time = 1e06;
        long unsigned max_time_ecology = 1e06;

        // boundary at which the difference between 
        // two values can considered to be 0
        double convergence_boundary = 1e-06;

        // Euler's constant
        double eul = 0.01;

};


#endif
