#include <string>
#include "solver.hpp"

int main(int argc, char **argv)
{
    Parameters parms;

    parms.eul = std::stod(argv[1]);
    parms.c = std::stod(argv[2]);
    parms.base_name = argv[3];
    parms.sigma = std::stod(argv[4]);

    parms.gamma[G1_idx] = std::stod(argv[5]);
    parms.gamma[G2_idx] = std::stod(argv[6]);
    
    parms.psi[G1_idx] = std::stod(argv[7]);
    parms.psi[G2_idx] = std::stod(argv[8]);

    parms.F[S_idx]  = std::stod(argv[9]);
    parms.F[G1_idx]  = std::stod(argv[10]);
    parms.F[G2_idx]  = std::stod(argv[11]);
    parms.F[G1G2_idx]  = std::stod(argv[12]);

    parms.d[S_idx]  = std::stod(argv[13]);
    parms.d[G1_idx]  = std::stod(argv[14]);
    parms.d[G2_idx]  = std::stod(argv[15]);
    parms.d[G1G2_idx]  = std::stod(argv[16]);

    parms.kappa  = std::stod(argv[17]);

    parms.init_vals[S_idx]  = std::stod(argv[18]);
    parms.init_vals[G1_idx]  = std::stod(argv[19]);
    parms.init_vals[G2_idx]  = std::stod(argv[20]);
    parms.init_vals[G1G2_idx]  = std::stod(argv[21]);

    Solver sol(parms);
    
    return 0;
}
