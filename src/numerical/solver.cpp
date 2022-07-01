#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "solver.hpp"

// solver constructur
Solver::Solver(Parameters const &params) : // data member initializers
    data_file{params.base_name.c_str()}
    ,par{params}
    ,popsizes{params.init_vals[0], 
        params.init_vals[1], 
        params.init_vals[2],
        params.init_vals[3]}
    ,pi{par.pi_init}
{
    // auxiliary variable that checks 
    // whether population numbers have converged

    write_data_headers();

    eigenvectors();

    exit(1);


    double pi_tplus1 = 0.0;

    for (time_step = 0; 
            time_step <= par.max_time; ++time_step)
    {
        // solve for the ecological equilibrium
        solve_endemic_eq();

        // update the value of pi
        pi_tplus1 = std::clamp(
                pi + par.eul * selgrad_pi()
                ,0.0
                ,1.0);

        if (std::abs(pi_tplus1 - pi) < par.convergence_boundary)
        {
            write_data();
            break;
        }

        // update the population sizes of each category
        // make sure we check for convergence
        
        if (time_step % par.output_time_interval == 0)
        {
            write_data();
        }
    } // end for time_step

    write_parameters();
} // end Solver::Solver(Parameters const &params) :


double Solver::selgrad_pi()
{
    return(0.0);
}

// solve for the endemic equilibrium
void Solver::solve_endemic_eq()
{
    // boolean variable to test for convergence
    bool converged;
    
    double popsizes_tplus1[4] = {0.0,0.0,0.0,0.0};

    for (long unsigned int ecol_time = 0; ecol_time < par.max_time_ecology;
            ++ecol_time)
    {
        popsizes_tplus1[S_idx] = popsizes[S_idx] + 
            par.eul * dSdt();

        popsizes_tplus1[G1_idx] = popsizes[G1_idx] + 
            par.eul * dG1dt();

        popsizes_tplus1[G2_idx] = popsizes[G2_idx] + 
            par.eul * dG2dt();

        popsizes_tplus1[G1G2_idx] = popsizes[G1G2_idx] + 
            par.eul * dG1G2dt();

        for (int idx = 0; idx < 4; ++idx)
        {
            // set limit to the population size;
            if (popsizes_tplus1[idx] < 0.0)
            {
                popsizes_tplus1[idx] = 0.0;
            }


            if (std::abs(popsizes_tplus1[idx] - popsizes[idx]) > 
                    par.convergence_boundary)
            {
                converged = false;
            }

            popsizes[idx] = popsizes_tplus1[idx];
        }

        if (converged)
        {
            break;
        }
    }
} // end  solve_endemic_eq()

// write the headers to the data file
void Solver::write_data_headers()
{
    data_file << "time_step;S;G1;G2;G1G2;N;pi;" << std::endl;
} // end write_data_headers

void Solver::write_data()
{
    data_file << time_step << ";" 
        << popsizes[S_idx] << ";"
        << popsizes[G1_idx] << ";"
        << popsizes[G2_idx] << ";"
        << popsizes[G1G2_idx] << ";"
        << N << ";"
        << pi << ";" << std::endl;
} // write_data()

void Solver::write_parameters()
{
    data_file << std::endl
        << std::endl
        << "FS;" << par.F[S_idx] << std::endl
        << "FG1;" << par.F[G1_idx] << std::endl
        << "FG2;" << par.F[G2_idx] << std::endl
        << "FG1G2;" << par.F[G1G2_idx] << std::endl
        << "bmax;" << par.bmax << std::endl
        << "c;" << par.c << std::endl

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

double Solver::b(PopTypes const idx)
{
    return(par.F[idx] * par.bmax * exp(-par.c * pi));
}

void Solver::eigenvectors()
{
    // dimensions of the matrix as follows: 
    Eigen::MatrixXd m(4,4);

    // row 1:

    m(S_idx,S_idx) = b(S_idx) -
        (par.d[S_idx] + (1.0 - pi) * (par.psi[G1_idx] + par.psi[G2_idx]));

    m(S_idx,G1_idx) = par.gamma[G1_idx];

    m(S_idx,G2_idx) = par.gamma[G2_idx];

    m(S_idx,G1G2_idx) = 0.0;


    // row 2:
    m(G1_idx, S_idx) = (1.0 - pi) * par.psi[G1_idx];

    m(G1_idx, G1_idx) = b(G1_idx) - 
        (par.d[G1_idx] + par.gamma[G1_idx] + par.sigma * (1.0 - pi) * par.psi[G2_idx]);

    m(G1_idx, G2_idx) = 0.0;

    m(G1_idx, G1G2_idx) = par.gamma[G2_idx]; // // los of I_G1G2 due to loss of G2, resulting in IG1


    // row 3:
    m(G2_idx, S_idx) = (1.0 - pi) * par.psi[G2_idx];

    m(G2_idx, G1_idx) = 0.0;

    m(G2_idx, G2_idx) = b(G2_idx) -
        (par.d[G2_idx] + par.gamma[G2_idx] + par.sigma * (1.0 - pi) * par.psi[G1_idx]);

    m(G2_idx, G1G2_idx) = par.gamma[G1_idx]; // los of I_G1G2 due to loss of G1, resulting in IG2

    // row 4:
    m(G1G2_idx, S_idx) = 0.0;

    // infection of a G1 individual with a G2 FGE
    m(G1G2_idx, G1_idx) = (1.0 - pi) * par.sigma * par.psi[G2_idx];
    
    // infection of a G2 individual with a G1 FGE
    m(G1G2_idx, G2_idx) = (1.0 - pi) * par.sigma * par.psi[G1_idx];

    m(G1G2_idx, G1G2_idx) = b(G1G2_idx) - 
        (par.d[G1G2_idx] + par.gamma[G1_idx] + par.gamma[G2_idx]);

    std::cout << m << std::endl;
} // end Solver::eigenvectors()


double Solver::dSdt()
{
    return(b(S_idx) * (1.0 - par.kappa * N) * popsizes[S_idx] 
            + par.gamma[G1_idx] * popsizes[G1_idx] 
            + par.gamma[G2_idx] * popsizes[G2_idx]
            - (par.d[S_idx] + (1.0 - pi) * (
                    par.psi[G1_idx] + par.psi[G2_idx])) * popsizes[S_idx]);
} // end dSdt

double Solver::dG1dt()
{
    return(b(G1_idx) * (1.0 - par.kappa * N) * popsizes[G1_idx]
            + (1.0 - pi) * par.psi[G1_idx] * popsizes[S_idx]
            - (par.d[G1_idx] + par.gamma[G1_idx] + 
                par.sigma * (1.0 - pi) * par.psi[G2_idx]) * popsizes[G1_idx]
            + par.gamma[G2_idx] * popsizes[G1G2_idx]);

} // end dG1dt

double Solver::dG2dt()
{
    return(b(G2_idx) * (1.0 - par.kappa * N) * popsizes[G2_idx]
            + (1.0 - pi) * par.psi[G2_idx] * popsizes[S_idx]
            - (par.d[G2_idx] + par.gamma[G2_idx] + 
                par.sigma * (1.0 - pi) * par.psi[G1_idx]) * popsizes[G2_idx]
            + par.gamma[G1_idx] * popsizes[G1G2_idx]);
} // end dG2dt

double Solver::dG1G2dt()
{
    return(b(G1G2_idx) * (1.0 - par.kappa * N) * popsizes[G1G2_idx]
            + par.sigma * (1.0 - pi) * (
                par.psi[G1_idx] * popsizes[G2_idx] + 
                par.psi[G2_idx] * popsizes[G1_idx])
            - (par.d[G1G2_idx] + par.gamma[G1_idx] + par.gamma[G2_idx]) * popsizes[G1G2_idx]);
} // end dG1G2dt
