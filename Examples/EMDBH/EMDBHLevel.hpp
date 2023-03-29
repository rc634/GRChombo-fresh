/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMDBHLEVEL_HPP_
#define EMDBHLEVEL_HPP_

#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"
// Problem specific includes
#include "EMDCouplingFunction.hpp"
#include "EinsteinMaxwellDilatonField.hpp"
// includes startracking
#include "STAMR.hpp"

//!  A class for the evolution of a single emdbh.
/*!
     The class takes some initial data for a einstein-maxwell-dilaton black hole
     and evolves it using the CCZ4 equations.
*/
class EMDBHLevel : public GRAMRLevel
{
    friend class DefaultLevelFactory<EMDBHLevel>;
    // Inherit the contructors from GRAMRLevel
    using GRAMRLevel::GRAMRLevel;

    STAMR &m_st_amr = dynamic_cast<STAMR &>(m_gr_amr);

    // Typedef for scalar field
    typedef EinsteinMaxwellDilatonField<CouplingFunction> EinsteinMaxwellDilatonFieldWithCoupling;

    //! Things to do at the end of the advance step, after RK4 calculation
    virtual void specificAdvance() override;

    //! Initialize data for the field and metric variables
    virtual void initialData() override;

    //! routines to do before outputing checkpoint file
    virtual void preCheckpointLevel() override;

    //! routines to do before outputing plot file
    virtual void prePlotLevel() override;

    //! RHS routines used at each RK4 step
    virtual void specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                 const double a_time) override;

    //! Things to do in UpdateODE step, after soln + rhs update
    virtual void specificUpdateODE(GRLevelData &a_soln,
                                   const GRLevelData &a_rhs,
                                   Real a_dt) override;

    //! Specify which variables to write at plot intervals
    // virtual void specificWritePlotHeader(std::vector<int> &plot_states)
    // const; // i think this is obsolete

    //! Tell Chombo how to tag cells for regridding
    virtual void
    computeTaggingCriterion(FArrayBox &tagging_criterion,
                            const FArrayBox &current_state) override;

    //! Things to do for analysis after each timestep and at the start
    virtual void doAnalysis() override;
};

#endif /* EMDBHLEVEL_HPP_ */
