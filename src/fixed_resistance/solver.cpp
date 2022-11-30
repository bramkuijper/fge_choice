#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include "solver.hpp"
#include "parameters.hpp"

// constructor of the solver class
Solver::Solver(Parameters &parms) :
    params{parms} // initialize parameters
    ,data_file{parms.base_name.c_str()} // initialize a data file
{
    // write the headers of the data file
    write_data_headers();

    // initialize population sizes
    for (int host_idx = 0; host_idx < 2; ++host_idx)
    {
        // initialize population sizes of the susceptibles
        popsize[host_idx] = params.init_popsize[host_idx];
        
        // initialize population sizes of the superinfected hosts
        popsize_superinfected[host_idx] = 
            params.init_popsize_superinfected[host_idx];

        for (int phage_idx = 0; phage_idx < 2; ++phage_idx)
        {
            // initialize population sizes of the infecteds
            popsize_infected[host_idx][phage_idx] = 
                params.init_popsize_infected[host_idx][phage_idx];
        }
    }

    // variables to store the updated values in the time series
    double popsize_tplus1[2] = {0.0,0.0};
    double popsize_superinfected_tplus1[2] = {0.0,0.0};
    double popsize_infected_tplus1[2][2] = {{0.0,0.0},{0.0,0.0}};

    // aux variables to keep track of indices
    HostType host_type;
    PhageType phage_type;

    bool converged;

    update_N();

    // numerically solve over time 
    for (time_step = 0; 
            time_step <= parms.max_ecol_time; ++time_step)
    {
        // output every n generations, except at the beginning
        if (time_step % parms.print_interval == 0 || time_step < 10000)
        {
            write_data();
        }

        // vary host types
        for (int host_type_idx = 0; 
                host_type_idx < 2; ++host_type_idx)
        {
            // get a host_type variable that is the same as the
            // index but in HostType format (used for accessing functions
            // that accept the HostType)
            host_type = static_cast<HostType>(host_type_idx);

            // perform the actual updating
            popsize_tplus1[host_type_idx] = popsize[host_type_idx] + 
                params.eul * dSdt(host_type);

            // check for boundary
            if (popsize_tplus1[host_type_idx] < 0.0)
            {
                popsize_tplus1[host_type_idx] = 0.0;
            }
            
            // check whether values are sane
            if (!std::isfinite(popsize_tplus1[host_type_idx]))
            {
                std::stringstream msg;

                msg << "isfinite() error at time step "
                        << time_step
                        << " and susceptible type "
                        << (host_type == C ? "C" : "P")
                        << " with value " << popsize_tplus1[host_type_idx];
                        
                throw std::range_error(msg.str());
            }

            // now update the superinfected host numbers
            // perform the actual updating
            popsize_superinfected_tplus1[host_type_idx] = 
                popsize_superinfected[host_type_idx] + 
                    params.eul * dIGBdt(host_type);

            // check for boundary
            if (popsize_superinfected_tplus1[host_type_idx] < 0.0)
            {
                popsize_superinfected_tplus1[host_type_idx] = 0.0;
            }
            
            // check whether values are sane
            if (!std::isfinite(popsize_superinfected_tplus1[host_type_idx]))
            {
                std::stringstream msg;

                msg << "isfinite() error at time step "
                        << time_step
                        << " and superinfected type "
                        << (host_type == C ? "C" : "P")
                        << " with value " << popsize_superinfected_tplus1[host_type_idx];
                        
                throw std::range_error(msg.str());
            }
            
            // now update dynamics for all 
            // Host Type x Phage Type
            // infecteds
            for (int phage_type_idx = 0; 
                    phage_type_idx < 2; ++phage_type_idx)
            {
                phage_type = static_cast<PhageType>(phage_type_idx);

                popsize_infected_tplus1[host_type_idx][phage_type_idx] = 
                    popsize_infected[host_type_idx][phage_type_idx] + 
                    params.eul * dIdt(host_type, phage_type);

                if (popsize_infected_tplus1[host_type_idx][phage_type_idx] < 0.0)
                {
                    popsize_infected_tplus1[host_type_idx][phage_type_idx] = 0.0;
                }

                if (!std::isfinite(
                            popsize_infected_tplus1[host_type_idx][phage_type_idx]))
                {
                    std::stringstream msg;

                    msg << "isfinite() error at time step "
                            << time_step
                            << " and infected type "
                            << (host_type == C ? "C" : "P")
                            << (phage_type == G1 ? "G1" : "G2")
                            << " with value " 
                            << popsize_infected_tplus1[host_type_idx][phage_type_idx];
                            
                    throw std::range_error(msg.str());
                }
            } // end for phage_type_idx
        } // end for host_type_idx

        converged = true;
        
        for (int host_type_idx = 0; 
                host_type_idx < 2; ++host_type_idx)
        {
            if (std::abs(popsize[host_type_idx] - 
                popsize_tplus1[host_type_idx]) > 
                    parms.vanish_threshold)
            {
                converged = false;
            }
            
            popsize[host_type_idx] = 
                popsize_tplus1[host_type_idx];

            if (std::abs(popsize_superinfected[host_type_idx] -
                        popsize_superinfected_tplus1[host_type_idx]) > 
                    parms.vanish_threshold)
            {
                converged = false;
            }
            
            popsize_superinfected[host_type_idx] = 
                popsize_superinfected_tplus1[host_type_idx];

            for (int phage_type_idx = 0; 
                    phage_type_idx < 2; ++phage_type_idx)
            {
                // OK did not converge
                if (std::abs(popsize_infected[host_type_idx][phage_type_idx] - 
                    popsize_infected_tplus1[host_type_idx][phage_type_idx]) > 
                        parms.vanish_threshold)
                {
                    converged = false;
                }

                popsize_infected[host_type_idx][phage_type_idx] = 
                    popsize_infected_tplus1[host_type_idx][phage_type_idx];
            }
        }

        update_N();

        if (converged)
        {
            write_data();
            break;
        }

    } // end for time_step

    write_parameters();
} // end Solver::Solver()

// write a line to the data file
void Solver::write_data()
{
    data_file << time_step << ";";

    // totals for each of the categories c and p
    double totals_pc[] = {0.0,0.0};

    totals_pc[C] = popsize[C] + popsize_superinfected[C] + 
        popsize_infected[C][G1] + popsize_infected[C][G2];

    totals_pc[P] = popsize[P] + popsize_superinfected[P] +
        popsize_infected[P][G1] + popsize_infected[P][G2];

    for (int host_type_idx = 0; 
            host_type_idx < 2; ++host_type_idx)
    {
        // write out numbers of susceptibles that carry either c or b
        data_file << popsize[host_type_idx] << ";";

        // write out conditional frequencies
        // p_S|c and p_S|p
        data_file << popsize[host_type_idx] / totals_pc[host_type_idx] << ";";

        data_file << popsize_superinfected[host_type_idx] << ";";

        for (int phage_type_idx = 0; 
                phage_type_idx < 2; ++phage_type_idx)
        {
            data_file << popsize_infected[host_type_idx][phage_type_idx] << ";";

            // write out conditional frequencies
            // p_G|c, p_B|c, p_G|p, p_B|p
            data_file << popsize_infected[host_type_idx][phage_type_idx] / 
                totals_pc[host_type_idx] << ";";
        }
    }

    data_file << N << ";" << std::endl;
}

void Solver::write_data_headers()
{
    data_file << "time;";

    for (int host_type_idx = 0; 
            host_type_idx < 2; ++host_type_idx)
    {
        data_file << "S" << (host_type_idx == P ? "P" : "C") << ";";
        data_file << "pS" << (host_type_idx == P ? "P" : "C") << ";";
        data_file << "I" << (host_type_idx == P ? "P" : "C") << "G1G2;";
        
        for (int phage_type_idx = 0; 
                phage_type_idx < 2; ++phage_type_idx)
        {
            data_file << "I" << (host_type_idx == P ? "P" : "C") << 
                (phage_type_idx == G1 ? "G1" : "G2") << ";";

            data_file << "p" << (phage_type_idx == G1 ? "G1" : "G2") <<
                (host_type_idx == P ? "P" : "C") << ";";
        }
    }

    data_file << "N;" << std::endl;
} // end write_data_headers()

void Solver::write_parameters()
{
    data_file << std::endl
        << std::endl;

    std::string host_id,phage_id;
    
    for (int host_type_idx = 0; 
            host_type_idx < 2; ++host_type_idx)
    {

        host_id = host_type_idx == C ? "C" : "P";

        data_file << "dS" << host_id << ";" 
            << params.dS[host_type_idx] << std::endl;


        for (int phage_type_idx = 0; 
                phage_type_idx < 2; ++phage_type_idx)
        {
            phage_id = phage_type_idx == G1 ? "G1" : "G2";

            data_file << "gamma" << host_id << phage_id << ";" 
                << params.gamma[host_type_idx][phage_type_idx] << std::endl
                    << "dI" << host_id << phage_id << ";" 
                << params.dI[host_type_idx][phage_type_idx] << std::endl;
        }

    }

    for (int phage_type_idx = 0; 
            phage_type_idx < 2; ++phage_type_idx)
    {
        phage_id = phage_type_idx == G1 ? "G1" : "G2";

        data_file << "psi" << phage_id << ";" << params.psi[phage_type_idx] << std::endl;
        data_file << "F" << phage_id << ";" << params.F[phage_type_idx] << std::endl;
    }

    data_file << "pi;" << params.pi << std::endl
        << "FGB;" << params.FGB << std::endl
        << "c;" << params.c << std::endl
        << "kappa;" << params.kappa << std::endl
        << "sigma;" << params.sigma << std::endl
        << "eul;" << params.eul << std::endl
        << "demog_feedback;" << params.demog_feedback << std::endl
        << "vanish_threshold;" << params.vanish_threshold << std::endl;
} // write_parameters()


void Solver::update_N()
{
    N = 0.0;
    for (int host_type_idx = 0; 
            host_type_idx < 2; ++host_type_idx)
    {
        N += popsize[host_type_idx];
        N += popsize_superinfected[host_type_idx];

        for (int phage_type_idx = 0; 
                phage_type_idx < 2; ++phage_type_idx)
        {
            N += popsize_infected[host_type_idx][phage_type_idx];
        }
    }
} // end Solver::update_N

double Solver::dSdt(HostType host_idx) const
{
    return(b(host_idx) * 
            (1.0 - params.kappa * N) * popsize[host_idx]
            + params.gamma[host_idx][G1] * popsize_infected[host_idx][G1]
            + params.gamma[host_idx][G2] * popsize_infected[host_idx][G2]
            - (params.dS[host_idx] 
                + params.psi[G1] * (host_idx == C ? 1.0 - params.pi : 1.0) 
                    * (params.demog_feedback ? popsize_infected[P][G1] + 
                        popsize_infected[C][G1] : 1.0)
                + params.psi[G2] 
                    * (params.demog_feedback ? popsize_infected[P][G2] + 
                        popsize_infected[C][G2] : 1.0)
                ) * popsize[host_idx]
            );
} // end Solver::dSdt()

double Solver::b(HostType host_idx) const
{
    return(std::exp(-params.c * (
                    host_idx == C ? 
                        params.pi 
                        : 
                        0.0)));

}

double Solver::b(HostType host_idx, PhageType phage_idx) const
{
    return(params.F[phage_idx] * std::exp(-params.c * (
                    host_idx == C ? 
                        params.pi 
                        : 
                        0.0)));
}

double Solver::b(HostType host_idx, PhageType phage1_idx, PhageType phage2_idx) const
{
    return(params.FGB * std::exp(-params.c * (
                    host_idx == C ? 
                        params.pi 
                        : 
                        0.0)));
}

// infected differential
double Solver::dIdt(HostType host_idx, PhageType phage_idx) const
{
    return(b(host_idx, phage_idx) * (1.0 - params.kappa * N) * 
                popsize_infected[host_idx][phage_idx]
            + params.psi[phage_idx] * popsize[host_idx] * 
                (host_idx == C && phage_idx == G1 ? 1.0 - params.pi : 1.0) * 
                (params.demog_feedback ? popsize_infected[P][phage_idx] + 
                                        popsize_infected[C][phage_idx] : 1.0)
                - (params.gamma[host_idx][phage_idx] + params.dI[host_idx][phage_idx]
                        + params.sigma * params.psi[!phage_idx] * 
                        (host_idx == C && phage_idx != G1 ? 1.0 - params.pi : 1.0) *
                            (params.demog_feedback ? 
                                popsize_infected[C][!phage_idx]
                                   + popsize_infected[P][!phage_idx]
                                :
                                1.0)
                    ) * popsize_infected[host_idx][phage_idx]
             + params.gamma[host_idx][!phage_idx] * popsize_superinfected[host_idx]
            );
}

// superinfected differential
double Solver::dIGBdt(HostType host_idx) const
{
    return(b(host_idx, G1, G2) * (1.0 - params.kappa * N) *
            popsize_superinfected[host_idx]
            // infection of a IxG2 host with G1
            + params.sigma * params.psi[G1] * popsize_infected[host_idx][G2] *
                (host_idx == C ? 1.0 - params.pi : 1.0) *
                (params.demog_feedback ? popsize_infected[host_idx][G1] + 
                popsize_infected[!host_idx][G1] : 1.0)
            // infection of a IxG1 host with G2
            + params.sigma * params.psi[G2] * popsize_infected[host_idx][G1] *
                (params.demog_feedback ? popsize_infected[host_idx][G2] +
                 popsize_infected[!host_idx][G2] : 1.0)
            - (params.gamma[host_idx][G1] + params.gamma[host_idx][G2] 
                + params.dBG[host_idx]) * popsize_superinfected[host_idx]);
}
