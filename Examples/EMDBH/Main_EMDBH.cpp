/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "parstream.H" //Gives us pout()
#include <chrono>
#include <iostream>

#include "CallDoAnalysis.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMR.hpp"
#include "GRParmParse.hpp"
#include "SetupFunctions.hpp"
#include "SimulationParameters.hpp"

// builds on GRAMR
#include "BHAMR.hpp"

// Problem specific includes:
#include "EMDBHLevel.hpp"

int runGRChombo(int argc, char *argv[])
{
    // Load the parameter file and construct the SimulationParameter class
    // To add more parameters edit the SimulationParameters file.
    char *in_file = argv[1];
    GRParmParse pp(argc - 2, argv + 2, NULL, in_file);
    SimulationParameters sim_params(pp);

    // The line below selects the problem that is simulated
    // (To simulate a different problem, define a new child of AMRLevel
    // and an associated LevelFactory)
    BHAMR bh_amr;

    // from my old stamr class for tracking stars
    // st_amr.m_star_tracker.initial_setup(
    //     sim_params.star_track_centre, sim_params.do_star_track,
    //     sim_params.number_of_stars, sim_params.initial_star_centres,
    //     sim_params.star_track_resolution, sim_params.star_track_width);

    DefaultLevelFactory<EMDBHLevel> emdbh_level_fact(bh_amr, sim_params);
    setupAMRObject(bh_amr, emdbh_level_fact);


    // Instantiate AMR interpolator for mass/GW extraction
    /*AMRInterpolator<Lagrange<4>> interpolator(
        bh_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);*/
    AMRInterpolator<Lagrange<4>> interpolator(
        bh_amr, sim_params.origin, sim_params.dx, sim_params.boundary_params,
        sim_params.verbosity);
    bh_amr.set_interpolator(&interpolator);



    #ifdef USE_AHFINDER
        if (sim_params.AH_activate)
        {
            AHSphericalGeometry sph(sim_params.horizon_centre_1);
            bh_amr.m_ah_finder.add_ah(sph, sim_params.AH_initial_guess,
                                      sim_params.AH_params);
        }
    #endif



    // Add a scheduler to GRAMR which just calls doAnalysis on every AMRLevel
    // at time 0. It is called later in postTimeStep
    RefCountedPtr<CallDoAnalysis> call_do_analysis_ptr(new CallDoAnalysis);
    RefCountedPtr<Scheduler> scheduler_ptr(new Scheduler);
    scheduler_ptr->schedule(call_do_analysis_ptr, sim_params.max_steps);
    bh_amr.schedule(scheduler_ptr);

    using Clock = std::chrono::steady_clock;
    using Minutes = std::chrono::duration<double, std::ratio<60, 1>>;

    std::chrono::time_point<Clock> start_time = Clock::now();

    // Engage! Run the evolution.
    bh_amr.run(sim_params.stop_time, sim_params.max_steps);

    auto now = Clock::now();
    auto duration = std::chrono::duration_cast<Minutes>(now - start_time);
    pout() << "Total simulation time (mins): " << duration.count() << ".\n";

    bh_amr.conclude();

    // Write Chombo timer report
    CH_TIMER_REPORT();

    return 0;
}

int main(int argc, char *argv[])
{
    mainSetup(argc, argv);

    int status = runGRChombo(argc, argv);

    if (status == 0)
        pout() << "GRChombo finished." << std::endl;
    else
        pout() << "GRChombo failed with return code " << status << std::endl;

    mainFinalize();
    return status;
}
