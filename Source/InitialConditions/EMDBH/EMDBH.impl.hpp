/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(EMDBH_HPP_)
#error "This file should only be included through EMDBH.hpp"
#endif

#ifndef EMDBH_IMPL_HPP_
#define EMDBH_IMPL_HPP_

#include "EMDBHSolution.hpp" //for EMDBHSolution class

inline EMDBH::EMDBH(EMDBH_params_t a_params_EMDBH,
                    CouplingFunction::params_t a_params_coupling_function,
                            double a_G_Newton, double a_dx, int a_verbosity)
    : m_dx(a_dx), m_G_Newton(a_G_Newton),
      m_params_EMDBH(a_params_EMDBH),
      m_params_coupling_function(a_params_coupling_function), m_verbosity(a_verbosity)
{
}

void EMDBH::compute_1d_solution(const double max_r)
{
    try
    {
        // set initial parameters and then run the solver
        pout() << "Setting initial conditions" << endl;
        m_1d_sol.set_initialcondition_params(m_params_EMDBH,
                                             m_params_coupling_function, max_r);
        pout() << "run m_1d_sol.main()" << endl;
        m_1d_sol.main();
        pout() << "completed m_1d_sol.main()" << endl;
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

// Compute the value of the initial vars on the grid
template <class data_t> void EMDBH::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<EinsteinMaxwellDilatonField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)
    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_EMDBH.star_centre);


    double separation = m_params_EMDBH.separation;
    bool binary = m_params_EMDBH.binary;

    // boosts and coordinate objects
    double x = (coords.x - separation / 2.);
    double z = coords.z;
    double y = coords.y;
    double r = sqrt(x * x + y * y + z * z);

    // first star physical variables
    double lapse = m_1d_sol.get_lapse_interp(r);
    double psi = m_1d_sol.get_psi_interp(r);
    double At = m_1d_sol.get_At_interp(r);


    if (binary)
    {
        // boosts and coordinate objects
        x = (coords.x + separation / 2.);
        r = sqrt(x * x + y * y + z * z);

        // first star physical variables
        lapse += m_1d_sol.get_lapse_interp(r)-1.;
        psi += m_1d_sol.get_psi_interp(r)-1.;
        At += m_1d_sol.get_At_interp(r);
    }


    vars.chi = pow(psi * psi, -1. / 3.);

    // can pick either lapse maybe
    vars.lapse += sqrt(vars.chi);
    //vars.lapse += lapse;

    vars.At += At;

    double kroneka[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

    FOR2(i, j) vars.h[i][j] = kroneka[i][j];

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* EMDBH_IMPL_HPP_ */
