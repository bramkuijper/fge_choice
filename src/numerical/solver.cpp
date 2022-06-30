#include <iostream>
#include <fstream>
#include "parameters.hpp"
#include "solver.hpp"

// solver constructur
Solver::Solver(Parameters const &params) :
    data_file{params.base_name.c_str()}
    ,par{params}
    ,popsizes{params.init_vals[0], 
        params.init_vals[1], 
        params.init_vals[2],
        params.init_vals[3]}
{
    for (time_step = 0; 
            time_step < par.max_time; ++time_step)
    {
        popsizes[S_idx] += par.eul * dSdt();
        popsizes[G1_idx] += par.eul * dG1dt();
        popsizes[G2_idx] += par.eul * dG2dt();
        popsizes[G1G2_idx] += par.eul * dG1G2dt();

        if (time_step % par.output_time_interval == 0)
        {
            write_data();
        }
    }
} // end Solver::Solver(Parameters const &params) :


// write the headers to the data file
void Solver::write_data_headers()
{
    data_file << "time_step;S;G1;G2;G1G2;pi;" << std::endl;
} // end write_data_headers

void Solver::write_data()
{
    data_file << time_step << ";" 
        << popsizes[S_idx] << ";"
        << popsizes[G1_idx] << ";"
        << popsizes[G2_idx] << ";"
        << popsizes[G1G2_idx] << ";"
        << pi << ";";
} // write_data()

void Solver::write_parameters()
{
    data_file << std::endl
        << std::endl
        << "bS;" << par.b[S_idx] << std::endl
        << "bG1;" << par.b[G1_idx] << std::endl
        << "bG2;" << par.b[G2_idx] << std::endl
        << "bG1G2;" << par.b[G1G2_idx] << std::endl

        << "init_valsS;" << par.init_vals[S_idx] << std::endl
        << "init_valsG1;" << par.init_vals[G1_idx] << std::endl
        << "init_valsG2;" << par.init_vals[G2_idx] << std::endl
        << "init_valsG1G2;" << par.init_vals[G1G2_idx] << std::endl

        << "dS;" << par.d[S_idx] << std::endl
        << "dG1;" << par.d[G1_idx] << std::endl
        << "dG2;" << par.d[G2_idx] << std::endl
        << "dG1G2;" << par.d[G1G2_idx] << std::endl

        << "kappa;" << par.kappa << std::endl
        << "sigma;" << par.sigma << std::endl
        << "eul;" << par.eul << std::endl
        << "pi_init;" << par.pi_init << std::endl
        << "gammaG1;" << par.gamma[G1_idx] << std::endl
        << "gammaG2;" << par.gamma[G2_idx] << std::endl

        << "psiG1;" << par.psi[G1_idx] << std::endl
        << "psiG2;" << par.psi[G2_idx] << std::endl;
}

double Solver::dSdt()
{
    return(par.b[S_idx] * (1.0 - par.kappa * N) * popsizes[S_idx] 
            + par.gamma[G1_idx] * popsizes[G1_idx] 
            + par.gamma[G2_idx] * popsizes[G2_idx]
            - (par.d[S_idx] + (1.0 - pi) * (
                    par.psi[G1_idx] + par.psi[G2_idx])) * popsizes[S_idx]);
}

double Solver::dG1dt()
{
    return(par.b[S_idx] * (1.0 - par.kappa * N) * popsizes[G1_idx]
            + (1.0 - pi) * par.psi[G1_idx] * popsizes[S_idx]
            - (par.d[G1_idx] + par.gamma[G1_idx] + 
                par.sigma * (1.0 - pi) * par.psi[G2_idx]) * popsizes[G1_idx]
            + par.gamma[G2_idx] * popsizes[G1G2_idx]);

}

double Solver::dG2dt()
{
    return(0.l);
}

double Solver::dG1G2dt()
{
    return(0.0);
}
