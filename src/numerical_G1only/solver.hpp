#include <iostream>
#include <fstream>
#include "parameters.hpp"

enum PopTypes
{
    S_idx = 0,
    G1_idx = 1,
    G2_idx = 2,
    G1G2_idx = 3
};

class Solver
{
    public:
        std::ofstream data_file;

        Parameters par;

        // the evolving resistance value
        double pi = 0.0;
        double pi_tplus1 = 0.0;

        // total population size
        double N = 0.0;

        // population sizes of S, G1, G2 and G1G2
        double popsizes[4] = {0.0,0.0,0.0,0.0};

        // the stable class frequencies, normalized
        double u[4] = {0.0,0.0,0.0,0.0};

        // the reproductive values, normalized so that
        // u^T.v = 1
        double v[4] = {0.0,0.0,0.0,0.0};

        double eigenval = 0.0;

        unsigned long time_step = 0;

        void write_data_headers();
        void write_data();
        void write_parameters(int const state);
        
        Solver(Parameters const &parms);

        double dSdt();
        double dG1dt();
        double dG2dt();
        double dG1G2dt();

        // fecundity function
        double b(PopTypes const idx);
        double dbdpi(PopTypes const idx);

        // solver for the eigenvectors
        void eigenvectors();
        
        void solve_endemic_eq();
        double selgrad_pi();
};


