/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "EMDBHParams.hpp"
#include "EMDCouplingFunction.hpp"

#include "AngMomFluxParams.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // Read the problem specific params
        readParams(pp);
    }

    void readParams(GRParmParse &pp)
    {
        // for regridding
        pp.load("regrid_threshold_A", regrid_threshold_A);
        pp.load("regrid_threshold_phi", regrid_threshold_phi);
        pp.load("regrid_threshold_chi", regrid_threshold_chi);

        // Gravitional constant
        pp.load("G_Newton", G_Newton, 1.0);

        // EMDBH initial data params
        pp.load("emd_data_path", emdbh_params.data_path); // the default should fail
        pp.load("gridpoints", emdbh_params.gridpoints, 40000);
        pp.load("star_centre", emdbh_params.star_centre,
                {0.5 * L, 0.5 * L, 0.5 * L});
        pp.load("binary", emdbh_params.binary, false);
        pp.load("separation", emdbh_params.separation, 0.0);
        pp.load("bh_charge", emdbh_params.bh_charge, 0.0);
        pp.load("bh_mass", emdbh_params.bh_mass, 1.0);
        pp.load("G_Newton", emdbh_params.Newtons_constant, 1.0);

        // Coupling params
        pp.load("emd_alpha", coupling_function_params.alpha, 0.0);
        pp.load("emd_f0", coupling_function_params.f0, 0.0);
        pp.load("emd_f1", coupling_function_params.f1, 0.0);
        pp.load("emd_f2", coupling_function_params.f2, 0.0);

        // Perturbation shell params
        pp.load("Ylm_amplitude", emdbh_params.Ylm_amplitude, 0.);
        pp.load("Ylm_thickness", emdbh_params.Ylm_thickness, 3.);
        pp.load("Ylm_r0", emdbh_params.Ylm_r0 , 0.5 * L);

        // Apparent Horizon stuff
        #ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, emdbh_params.bh_mass*0.5);
        #endif
        pp.load("AH_num_horizons", AH_num_horizons, 0);
        pp.load("AH_expect_merger", AH_expect_merger, 0);
        pp.load("horizon_centre_1", horizon_centre_1,
                {0.5 * L, 0.5 * L, 0.5 * L});
        pp.load("horizon_centre_2", horizon_centre_2,
                {0.5 * L, 0.5 * L, 0.5 * L});

        // Mass extraction
        pp.load("activate_mass_extraction", activate_mass_extraction, 0);
        pp.load("num_mass_extraction_radii",
                mass_extraction_params.num_extraction_radii, 1);
        pp.load("mass_extraction_levels",
                mass_extraction_params.extraction_levels,
                mass_extraction_params.num_extraction_radii, 0);
        pp.load("mass_extraction_radii",
                mass_extraction_params.extraction_radii,
                mass_extraction_params.num_extraction_radii, 0.1);
        pp.load("num_points_phi_mass", mass_extraction_params.num_points_phi,
                2);
        pp.load("num_points_theta_mass",
                mass_extraction_params.num_points_theta, 4);
        pp.load("mass_extraction_center",
                mass_extraction_params.extraction_center,
                {0.5 * L, 0.5 * L, 0.5 * L});

        // Weyl extraction
        pp.load("activate_gw_extraction", activate_weyl_extraction, 0);

        // Em extraction stuff
        pp.load("activate_em_extraction", activate_pheyl_extraction, 0);

        // Em extraction stuff
        pp.load("activate_rs_extraction", activate_realscalar_extraction, 0);

        // Star Tracking
        /*pp.load("do_star_tracking", gaussfit_params.do_star_tracking, 0);
        pp.load("num_points_gaussian_fit", gaussfit_params.num_points, 50);
        // will be replaced
        pp.load("tracked_field_index", gaussfit_params.field_index, 30);
        pp.load("search_width", gaussfit_params.search_width, 16.);
        pp.load("tracking_BH_cutoff", gaussfit_params.BH_cutoff, 0.05);
        pp.load("tracking_AMR_level", gaussfit_params.AMR_level,0);
        pp.load("track_both_centres", gaussfit_params.track_both_centres, true);
        pp.load("track_min_separation",
        gaussfit_params.track_min_separation, 5.); pp.load("tracking_centre",
        gaussfit_params.track_centre, {0.,0.,0.}); pp.load("tracking_centres",
        gaussfit_params.track_centres, {0.,0.,0.,0.,0.,0.});*/ //THIS MIGHT ALL BE TOO OLD

        // Work out the minimum extraction level
        auto min_extraction_level_it =
            mass_extraction_params.min_extraction_level();

        // Do we cant to calculate L2 norms of constraint violations
        pp.load("calculate_constraint_violations",
                calculate_constraint_violations, false);

        // Do we want to calculate and write the Noether Charge to a file
        pp.load("calculate_noether_charge", calculate_noether_charge, false);

        // Variables for outputting to plot files
        // pp.load("num_plot_vars", num_plot_vars, 0);
        // pp.load("plot_vars", plot_vars, num_plot_vars, 0);

        // Variables for outputting inf-norm
        pp.load("num_vars_inf_norm", num_vars_inf_norm, 0);
        pp.load("vars_inf_norm", vars_inf_norm, num_vars_inf_norm, 0);

        pp.load("flux_extraction_level", flux_extraction_level, 0);

        pp.load("star_track_centre", star_track_centre,
                {0.5 * L, 0.5 * L, 0.5 * L});
        pp.load("do_star_track", do_star_track, false);
        pp.load("number_of_stars", number_of_stars, 1);
        pp.load("initial_star_centres", initial_star_centres,
                3 * number_of_stars);
        pp.load("star_track_resolution", star_track_resolution, 21);
        pp.load("star_track_width", star_track_width, 20.);
        pp.load("star_track_level", star_track_level, 0);


        pp.load("robin_manual_origin", robin_manual_origin,
                {0., 0., 0.}); // currently unused

        /*pp.load("flux_number_of_radii", angmomflux_params.number_radii,1);
        pp.load("flux_do", angmomflux_params.do_flux_integration,false);
        pp.load("flux_extraction_level", angmomflux_params.extraction_level,0);
        pp.load("flux_num_theta", angmomflux_params.num_theta,10);
        pp.load("flux_num_phi", angmomflux_params.num_phi,10);
        pp.load("flux_extraction_centre", angmomflux_params.centre,
                                                {0.5 * L, 0.5 * L, 0.5 * L});

        angmomflux_params.radii.resize(angmomflux_params.number_radii);
        pp.load("flux_extraction_radii", angmomflux_params.radii,
                                                angmomflux_params.number_radii);*/
    }

    // Tagging thresholds
    Real regrid_threshold_phi, regrid_threshold_chi, regrid_threshold_A;

    // Initial data for matter and potential
    double G_Newton;
    EMDBH_params_t emdbh_params;
    CouplingFunction::params_t coupling_function_params;
    // GaussFit_params_t gaussfit_params; THIS MIGHT BE TOO OLD!

    // Mass extraction
    int activate_mass_extraction;
    extraction_params_t mass_extraction_params;

    int activate_weyl_extraction;
    int activate_pheyl_extraction;
    int activate_realscalar_extraction;

    // Do we want to write a file with the L2 norms of contraints?
    bool calculate_constraint_violations;

    // Do we want to write the Noether Charge to a file
    bool calculate_noether_charge;

    std::array<double, CH_SPACEDIM> robin_manual_origin; // currently unused

    // Vars for outputting in plot files
    // int num_plot_vars;
    // std::vector<int> plot_vars;

    // Vars for outputting inf-norms
    int num_vars_inf_norm;
    std::vector<int> vars_inf_norm;

    std::array<double, CH_SPACEDIM> star_track_centre;
    bool do_star_track;
    int number_of_stars;
    std::vector<double> initial_star_centres;
    int star_track_resolution;
    double star_track_width;
    int star_track_level;
    int flux_extraction_level;

    #ifdef USE_AHFINDER
    double AH_initial_guess;
    int AH_num_horizons;
    int AH_expect_merger;
    std::array<double, CH_SPACEDIM> horizon_centre_1;
    std::array<double, CH_SPACEDIM> horizon_centre_2;
    #endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
