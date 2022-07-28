#ifndef _SOLVER_HPP_
#define _SOLVER_HPP_

#include <iostream>
#include <fstream>
#include "parameters.hpp"

enum HostType
{
    P = 0,
    C = 1
};

enum PhageType
{
    G1 = 0,
    G2 = 1
};

class Solver
{
    public:
        Parameters params;
        std::ofstream data_file;

        double popsize[2] = {0.0,0.0};

        // host type x phage type
        double popsize_infected[2][2] = {
            {0.0,0.0}
            ,{0.0,0.0}
        };

        double N = 0.0;

        unsigned long time_step = 0;

        Solver(Parameters &par);

        // update the population size
        void update_N();

        double dSdt(HostType host_idx) const;
        double dIdt(HostType host_idx, 
                PhageType phage_idx) const;

        double b(HostType host_idx) const;
        double b(HostType host_idx, PhageType phage_idx) const;

        void write_data();
        void write_data_headers();
        void write_parameters();
};

#endif
