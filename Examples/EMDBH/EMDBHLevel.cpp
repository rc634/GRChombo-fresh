/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// General includes common to most GR problems
#include "EMDBHLevel.hpp"
#include "BoxLoops.hpp"
#include "GammaCalculator.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "TraceARemoval.hpp"

// For RHS update
#include "MatterCCZ4.hpp"

// For constraints calculation
#include "NewConstraints.hpp"
#include "NewMatterConstraints.hpp"

// For tag cells
#include "EMDExtractionTaggingCriterion.hpp"

// Problem specific includes
#include "EMDBH.hpp"
#include "EMDCouplingFunction.hpp"
#include "EinsteinMaxwellDilatonField.hpp"
#include "ComputePack.hpp"
#include "SetValue.hpp"

// For GW extraction
#include "MatterWeyl4.hpp"
#include "WeylExtraction.hpp"

// FOR EM extraction_level
#include "Pheyl2.hpp"

// For Noether Charge calculation
#include "EMDLorentzScalars.hpp"
#include "SmallDataIO.hpp"

// For EM Tensor
#include "EMTensor.hpp"

// For Star Tracking
// #include "GaussianFitTracking.hpp" // obsolete
// star tracking now included in EMDBHLevel.hpp in STAMR.hpp

// for chombo grid Functions
#include "AMRReductions.hpp"

// Things to do at each advance step, after the RK4 is calculated
void EMDBHLevel::specificAdvance()
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void EMDBHLevel::initialData()
{
    CH_TIME("EMDBHLevel::initialData");
    if (m_verbosity)
        pout() << "EMDBHLevel::initialData " << m_level << endl;

    // First initalise a EMDBH object
    EMDBH emdbh(m_p.emdbh_params, m_p.coupling_function_params,
                         m_p.G_Newton, m_dx, m_verbosity);

   if (m_verbosity)
       pout() << "EMDBHLevel::initialData - Compute_1d_solution " << m_level << endl;
    // the max radius the code might need to calculate out to is L*sqrt(3)
    emdbh.compute_1d_solution(2. * m_p.L);


    if (m_verbosity)
        pout() << "EMDBHLevel::initialData - Interpolate to 3D grid " << m_level << endl;
    // First set everything to zero ... we don't want undefined values in
    // constraints etc, then  initial conditions for EMDBH
    BoxLoops::loop(make_compute_pack(SetValue(0.0), emdbh), m_state_new,
                   m_state_new, INCLUDE_GHOST_CELLS, disable_simd());


    if (m_verbosity)
        pout() << "EMDBHLevel::initialData - GammaCalc " << m_level << endl;
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS, disable_simd());

    fillAllGhosts();
}

// Things to do before outputting a checkpoint file
void EMDBHLevel::preCheckpointLevel()
{
    CH_TIME("EMDBHLevel::preCheckpointLevel");

    fillAllGhosts();
    CouplingFunction coupling_function(m_p.coupling_function_params);
    EinsteinMaxwellDilatonFieldWithCoupling emd_field(coupling_function);
    BoxLoops::loop(
        make_compute_pack(

              MatterWeyl4<EinsteinMaxwellDilatonFieldWithCoupling>(
                              emd_field,
                              m_p.extraction_params.extraction_center, m_dx,
                              m_p.formulation, m_p.G_Newton),

              MatterConstraints<EinsteinMaxwellDilatonFieldWithCoupling>(
                              emd_field, m_dx, m_p.G_Newton, c_Ham,
                              Interval(c_Mom1, c_Mom3)),

              EMDLorentzScalars(m_dx),

              EMTensor<EinsteinMaxwellDilatonFieldWithCoupling>(
                              emd_field, m_dx, c_rho,
                              Interval(c_s1, c_s3), Interval(c_s11, c_s33)),

              Pheyl2(m_p.extraction_params.extraction_center, m_dx)),

        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
}

// Things to do before outputting a plot file
void EMDBHLevel::prePlotLevel()
{
    CH_TIME("EMDBHLevel::prePlotLevel");

    /*fillAllGhosts();
    CouplingFunction coupling_function(m_p.coupling_function_params);
    EinsteinMaxwellDilatonFieldWithCoupling emd_field(coupling_function);
    BoxLoops::loop(make_compute_pack(
                    MatterWeyl4<EinsteinMaxwellDilatonFieldWithCoupling>(
                    emd_field,m_p.extraction_params.extraction_center,
                    m_dx, m_p.formulation, m_p.G_Newton),
                    MatterConstraints<EinsteinMaxwellDilatonFieldWithCoupling>(
                    emd_field, m_dx, m_p.G_Newton, c_Ham,
                    Interval(c_Mom1, c_Mom3)), EMDLorentzScalars(m_dx),
                    EMTensor<EinsteinMaxwellDilatonFieldWithCoupling>(
                    emd_field, m_dx, c_rho, Interval(c_s1,c_s3),
                    Interval(c_s11,c_s33))),
                    m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);*/

    fillAllGhosts();
    CouplingFunction coupling_function(m_p.coupling_function_params);
    EinsteinMaxwellDilatonFieldWithCoupling emd_field(coupling_function);

    BoxLoops::loop(
        make_compute_pack(

            MatterWeyl4<EinsteinMaxwellDilatonFieldWithCoupling>(
                emd_field, m_p.extraction_params.extraction_center,
                m_dx, m_p.formulation, m_p.G_Newton),

            MatterConstraints<EinsteinMaxwellDilatonFieldWithCoupling>(
                emd_field, m_dx, m_p.G_Newton, c_Ham,
                Interval(c_Mom1, c_Mom3)),

            EMDLorentzScalars(m_dx),

            EMTensor<EinsteinMaxwellDilatonFieldWithCoupling>(
                emd_field, m_dx, c_rho,
                Interval(c_s1, c_s3), Interval(c_s11, c_s33)),

            Pheyl2(m_p.extraction_params.extraction_center, m_dx)),

        m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    // older versio with angmomflux GammaCalculator
    //
    //
    // BoxLoops::loop(
    //     make_compute_pack(
    //
    //         MatterWeyl4<EinsteinMaxwellDilatonFieldWithCoupling>(
    //             emd_field, m_p.extraction_params.extraction_center,
    //             m_dx, m_p.formulation, m_p.G_Newton),
    //
    //         MatterConstraints<EinsteinMaxwellDilatonFieldWithCoupling>(
    //             emd_field, m_dx, m_p.G_Newton, c_Ham,
    //             Interval(c_Mom1, c_Mom3)),
    //
    //         EMDLorentzScalars(m_dx),
    //
    //         EMTensor_and_mom_flux<EinsteinMaxwellDilatonFieldWithCoupling>(
    //             emd_field, m_dx, m_p.L, m_p.angmomflux_params.center,
    //             c_Fphi_flux, c_Sphi_source, c_Qphi_density, c_rho,
    //             Interval(c_s1, c_s3), Interval(c_s11, c_s33))),
    //
    //     m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);

    /*BoxLoops::loop(
        SourceIntPreconditioner<EinsteinMaxwellDilatonFieldWithCoupling>(
            emd_field, m_dx, m_p.L, m_p.angmomflux_params.center,
            c_Sphi_source, c_Qphi_density, 10.),
        m_state_diagnostics, m_state_diagnostics, EXCLUDE_GHOST_CELLS);*/
}

// Things to do in RHS update, at each RK4 step
void EMDBHLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                     const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate MatterCCZ4 right hand side with matter_t = EinsteinMaxwellDilatonField
    // We don't want undefined values floating around in the constraints so
    // zero these
    CouplingFunction coupling_function(m_p.coupling_function_params);
    EinsteinMaxwellDilatonFieldWithCoupling emd_field(coupling_function);
    MatterCCZ4RHS<EinsteinMaxwellDilatonFieldWithCoupling> my_ccz4_matter(
        emd_field, m_p.ccz4_params, m_dx, m_p.sigma, m_p.formulation,
        m_p.G_Newton);
    SetValue set_analysis_vars_zero(0.0, Interval(c_Xi + 1, NUM_VARS - 1));
    auto compute_pack =
        make_compute_pack(my_ccz4_matter, set_analysis_vars_zero);
    BoxLoops::loop(compute_pack, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
}

// Things to do at ODE update, after soln + rhs
void EMDBHLevel::specificUpdateODE(GRLevelData &a_soln,
                                       const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce trace free A_ij
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

// Things to do for analysis after each timestep and at the start
void EMDBHLevel::doAnalysis()
{
    CH_TIME("EMDBHLevel::specificPostTimeStep");
    bool first_step = (m_time == 0.0);


    #ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
    {
        CH_TIME("EMDBHLevel::doAnalysis::AH_FINDER");
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
    }
    #endif


    if (m_p.activate_weyl_extraction == 1 &&
        at_level_timestep_multiple(
            m_p.extraction_params.min_extraction_level()))
    {
        CH_TIME("EMDBHLevel::doAnalysis::Weyl4&ADMMass");
        // First compute the ADM Mass integrand values on the grid
        fillAllGhosts();
        CouplingFunction coupling_function(m_p.coupling_function_params);
        EinsteinMaxwellDilatonFieldWithCoupling emd_field(coupling_function);
        auto weyl4_adm_compute_pack = make_compute_pack(
            MatterWeyl4<EinsteinMaxwellDilatonFieldWithCoupling>(
                emd_field, m_p.extraction_params.extraction_center,
                m_dx, m_p.formulation, m_p.G_Newton));
        BoxLoops::loop(weyl4_adm_compute_pack, m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);
        // Do the extraction on the min extraction level
        if (m_level == m_p.extraction_params.min_extraction_level())
        {
            if (m_verbosity)
            {
                pout() << "BinaryBSLevel::specificPostTimeStep:"
                          " Extracting gravitational waves."
                       << endl;
            }

            // Refresh the interpolator and do the interpolation
            m_bh_amr.m_interpolator->refresh();
            WeylExtraction gw_extraction(m_p.extraction_params, m_dt, m_time,
                                         first_step, m_restart_time);
            gw_extraction.execute_query(m_bh_amr.m_interpolator);
        }
    }

    // Do Electromagnetic
    if (m_p.activate_em_extraction == 1)
    {
    }

    // noether charge, max mod phi, min chi, constraint violations
    if (at_level_timestep_multiple(0))
    {
        BoxLoops::loop(EMDLorentzScalars(m_dx), m_state_new, m_state_diagnostics,
                       EXCLUDE_GHOST_CELLS);
    }
    if (m_level == 0)
    {
        AMRReductions<VariableType::diagnostic> amr_reductions(m_bh_amr);

        // EMD scalars should be calculated pre-check and pre plot
        // so automatically here

        // compute integrated volume integrals

        // F_{\mu\nu}F^{\mu\nu}
        // double FF = amr_reductions.sum(c_mod_F);
        // SmallDataIO FF_file("FF", m_dt, m_time,
        //                                 m_restart_time, SmallDataIO::APPEND,
        //                                 first_step);
        // FF_file.remove_duplicate_time_data();
        // if (m_time == 0.)
        // {
        //     FF_file.write_header_line(
        //                             {"Electromagnetic Tensor Modulus Squared"});
        // }
        // FF_file.write_time_data_line({FF});
        //
        // // A_\muA^\mu
        // double AA = amr_reductions.sum(c_mod_A);
        // SmallDataIO AA_file("AA", m_dt, m_time,
        //                                 m_restart_time, SmallDataIO::APPEND,
        //                                 first_step);
        // AA_file.remove_duplicate_time_data();
        // if (m_time == 0.)
        // {
        //     AA_file.write_header_line(
        //                             {"Electromagnetic Vector Modulus Squared"});
        // }
        // AA_file.write_time_data_line({AA});
        //
        // // hamiltonian of scalar field
        // double phiham = amr_reductions.sum(c_phi_ham);
        // SmallDataIO phiham_file("phiham", m_dt, m_time,
        //                                 m_restart_time, SmallDataIO::APPEND,
        //                                 first_step);
        // phiham_file.remove_duplicate_time_data();
        // if (m_time == 0.)
        // {
        //     phiham_file.write_header_line({"Hamiltonian of Scalar Field"});
        // }
        // phiham_file.write_time_data_line({phiham});
        //
        //
        // // Compute the maximum of mod_phi and write it to a file
        // double phi_max = amr_reductions.max(c_phi);
        // SmallDataIO phi_max_file("phi_max", m_dt, m_time,
        //                              m_restart_time, SmallDataIO::APPEND,
        //                              first_step);
        // phi_max_file.remove_duplicate_time_data();
        // if (m_time == 0.)
        // {
        //     phi_max_file.write_header_line({"max phi"});
        // }
        // phi_max_file.write_time_data_line({phi_max});
        //
        // // Compute the min of chi and write it to a file
        // double min_chi = amr_reductions.min(c_chi);
        // SmallDataIO min_chi_file("min_chi", m_dt, m_time, m_restart_time,
        //                          SmallDataIO::APPEND, first_step);
        // min_chi_file.remove_duplicate_time_data();
        // if (m_time == 0.)
        // {
        //     min_chi_file.write_header_line({"min chi"});
        // }
        // min_chi_file.write_time_data_line({min_chi});
        //
        // // constraeints calculated pre check and pre plot so done here already
        //
        // double L2_Ham = amr_reductions.norm(c_Ham);
        // double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3));
        // SmallDataIO constraints_file("constraint_norms", m_dt, m_time,
        //                              m_restart_time, SmallDataIO::APPEND,
        //                              first_step);
        // constraints_file.remove_duplicate_time_data();
        // if (first_step)
        // {
        //     constraints_file.write_header_line({"L^2_Ham", "L^2_Mom"});
        // }
        // constraints_file.write_time_data_line({L2_Ham, L2_Mom});
    }

    // double temp_dx;
    // if (m_p.do_flux_integration &&
    //     at_level_timestep_multiple(m_p.flux_extraction_level))
    // {
    //     CH_TIME("EMDBHLevel::doAnalysis::FphiSphi");
    //     CouplingFunction coupling_function(m_p.coupling_function_params);
    //     EinsteinMaxwellDilatonFieldWithCoupling emd_field(coupling_function);
    //
    //     // execute only on the finest level
    //     if (m_level == m_p.flux_extraction_level)
    //     // if (m_level==m_p.max_level)
    //     {
    //         double S_phi_integral; // integral of angmomsource
    //         double Q_phi_integral; // integral of angmomsource
    //         std::vector<AMRLevel *> all_level_ptrs =
    //             getAMRLevelHierarchy().stdVector();
    //         std::vector<double> S_phi_integrals(
    //             m_p.angmomflux_params
    //                 .num_extraction_radii); // vector storing all integrals
    //         std::vector<double> Q_phi_integrals(
    //             m_p.angmomflux_params
    //                 .num_extraction_radii); // vector storing all integrals
    //
    //         // temp_dx = m_p.coarsest_dx;
    //         // fill grid with angmom variables on every level
    //         for (auto level_ptr : all_level_ptrs)
    //         {
    //             EMDBHLevel *bs_level_ptr =
    //                 dynamic_cast<EMDBHLevel *>(level_ptr);
    //             if (bs_level_ptr == nullptr)
    //             {
    //                 break;
    //             }
    //             temp_dx = bs_level_ptr->m_dx;
    //             BoxLoops::loop(
    //                 EMTensor_and_mom_flux<EinsteinMaxwellDilatonFieldWithCoupling>(
    //                     emd_field, temp_dx, m_p.L,
    //                     m_p.angmomflux_params.center, c_Fphi_flux,
    //                     c_Sphi_source, c_Qphi_density, c_rho,
    //                     Interval(c_s1, c_s3), Interval(c_s11, c_s33)),
    //                 bs_level_ptr->m_state_new,
    //                 bs_level_ptr->m_state_diagnostics, EXCLUDE_GHOST_CELLS);
    //         }
    //         // loop through levels and fill ghosts
    //         for (auto level_ptr : all_level_ptrs)
    //         {
    //             EMDBHLevel *bs_level_ptr =
    //                 dynamic_cast<EMDBHLevel *>(level_ptr);
    //             if (bs_level_ptr == nullptr)
    //             {
    //                 break;
    //             }
    //             bs_level_ptr->fillAllGhosts();
    //         }
    //
    //         for (int i = m_p.angmomflux_params.num_extraction_radii - 1; i >= 0;
    //              i--)
    //         {
    //             temp_dx = m_p.coarsest_dx;
    //             for (auto level_ptr : all_level_ptrs)
    //             {
    //                 EMDBHLevel *bs_level_ptr =
    //                     dynamic_cast<EMDBHLevel *>(level_ptr);
    //                 if (bs_level_ptr == nullptr)
    //                 {
    //                     break;
    //                 }
    //
    //                 // set angmomsource and density to zero outside of
    //                 // extraction radii
    //                 temp_dx = bs_level_ptr->m_dx;
    //                 BoxLoops::loop(
    //                     SourceIntPreconditioner(
    //                         temp_dx, m_p.L, m_p.angmomflux_params.center,
    //                         c_Sphi_source, c_Qphi_density,
    //                         m_p.angmomflux_params.extraction_radii[i]),
    //                     bs_level_ptr->m_state_diagnostics,
    //                     bs_level_ptr->m_state_diagnostics, INCLUDE_GHOST_CELLS);
    //             }
    //             // old code passed this coarse dx, maybe this doesnt work now?
    //             AMRReductions<VariableType::diagnostic> amr_reductions(
    //                 m_bh_amr);
    //             S_phi_integral = amr_reductions.sum(c_Sphi_source);
    //             S_phi_integrals[i] = S_phi_integral;
    //             Q_phi_integral = amr_reductions.sum(c_Qphi_density);
    //             Q_phi_integrals[i] = Q_phi_integral;
    //         }
    //
    //         // save the Source integral to dat file
    //         std::vector<string> title_line(
    //             m_p.angmomflux_params.num_extraction_radii);
    //         string dummy_string;
    //         for (int j = 0; j < m_p.angmomflux_params.num_extraction_radii; j++)
    //         {
    //             dummy_string =
    //                 "r = " +
    //                 to_string(m_p.angmomflux_params.extraction_radii[j]);
    //             title_line[j] = dummy_string;
    //         }
    //
    //         SmallDataIO angmomsource_file("AngMomSource", m_dt, m_time,
    //                                       m_restart_time, SmallDataIO::APPEND,
    //                                       first_step);
    //
    //         if (m_time > 0)
    //             angmomsource_file.remove_duplicate_time_data();
    //
    //         if (m_time == 0.)
    //         {
    //             angmomsource_file.write_header_line(title_line);
    //         }
    //
    //         angmomsource_file.write_time_data_line(S_phi_integrals);
    //
    //         // save the Density integral to dat file
    //         std::vector<string> title_line2(
    //             m_p.angmomflux_params.num_extraction_radii);
    //         string dummy_string2;
    //         for (int j = 0; j < m_p.angmomflux_params.num_extraction_radii; j++)
    //         {
    //             dummy_string2 =
    //                 "r = " +
    //                 to_string(m_p.angmomflux_params.extraction_radii[j]);
    //             title_line2[j] = dummy_string2;
    //         }
    //
    //         SmallDataIO density_file("AngMomDensity", m_dt, m_time,
    //                                  m_restart_time, SmallDataIO::APPEND,
    //                                  first_step);
    //
    //         if (m_time > 0)
    //             density_file.remove_duplicate_time_data();
    //
    //         if (m_time == 0.)
    //         {
    //             density_file.write_header_line(title_line2);
    //         }
    //
    //         density_file.write_time_data_line(Q_phi_integrals);
    //
    //         // Refresh the interpolator and do the interpolation
    //         m_bh_amr.m_interpolator->refresh();
    //         // setup and do angmomflux integral
    //         AngMomFlux ang_mom_flux(m_p.angmomflux_params, m_time, m_dt,
    //                                 m_restart_time, first_step);
    //         ang_mom_flux.run(m_bh_amr.m_interpolator);
    //     }
    // }


    // star tracker not yet compatible with Ah finder,
    // would need to combine st_amr into bh_amr
    // if (m_p.do_star_track && m_level == m_p.star_track_level)
    // {
    //     // if at restart time read data from dat file,
    //     // will default to param file if restart time is 0
    //     if (fabs(m_time - m_restart_time) < m_dt * 1.1)
    //     {
    //         m_st_amr.m_star_tracker.read_old_centre_from_dat(
    //             "StarCentres", m_dt, m_time, m_restart_time, first_step);
    //     }
    //
    //     m_st_amr.m_star_tracker.update_star_centres(c_mod_F);
    //     m_st_amr.m_star_tracker.write_to_dat("StarCentres", m_dt, m_time,
    //                                          m_restart_time, first_step);
    // }
}

void EMDBHLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                             const FArrayBox &current_state)
{
    BoxLoops::loop(EMDExtractionTaggingCriterion(
                       m_dx, m_level, m_p.mass_extraction_params,
                       m_p.regrid_threshold_A, m_p.regrid_threshold_phi,
                       m_p.regrid_threshold_chi), current_state,
                                                  tagging_criterion);
}
