#include "solver.hpp"

int main(int argc, char **argv)
{
    Parameters pars;

    pars.gamma[C][G1] = std::stod(argv[1]);
    pars.gamma[C][G2] = std::stod(argv[2]);
    pars.gamma[P][G1] = std::stod(argv[3]);
    pars.gamma[P][G2] = std::stod(argv[4]);
    pars.psi[G1] = std::stod(argv[5]);
    pars.psi[G2] = std::stod(argv[6]);
    pars.pi = std::stod(argv[7]);
    pars.c = std::stod(argv[8]);

    Solver sol(pars);
    return 0;
}
