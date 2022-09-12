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

    pars.F[G1] = std::stod(argv[7]);
    pars.F[G2] = std::stod(argv[8]);

    pars.FGB = (pars.F[G1] + pars.F[G2]);

    pars.dS[P] = std::stod(argv[9]);
    pars.dS[C] = std::stod(argv[10]);

    pars.dI[P][G1] = std::stod(argv[11]);
    pars.dI[P][G2] = std::stod(argv[12]);
    pars.dI[C][G1] = std::stod(argv[13]);
    pars.dI[C][G2] = std::stod(argv[14]);

    pars.init_popsize[P] = std::stod(argv[15]);
    pars.init_popsize[C] = std::stod(argv[16]);

    pars.init_popsize_infected[P][G1] = std::stod(argv[17]);
    pars.init_popsize_infected[P][G2] = std::stod(argv[18]);
    pars.init_popsize_infected[C][G1] = std::stod(argv[19]);
    pars.init_popsize_infected[C][G2] = std::stod(argv[20]);

    pars.init_popsize_superinfected[P] = std::stod(argv[21]);
    pars.init_popsize_superinfected[C] = std::stod(argv[22]);

    pars.pi = std::stod(argv[23]);
    pars.c = std::stod(argv[24]);
    pars.kappa = std::stod(argv[25]);
    pars.sigma = std::stod(argv[26]);
    pars.eul = std::stod(argv[27]);
    pars.demog_feedback = std::stod(argv[28]);

    pars.base_name = argv[29];

    Solver sol(pars);
    return 0;
}
