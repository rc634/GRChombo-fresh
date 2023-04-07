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
    double safe_r = sqrt(x * x + y * y + z * z + 10e-20);


    // first star physical variables
    double psi = m_1d_sol.get_psi_interp(r);
    double dpsi = m_1d_sol.get_dpsi_interp(r);
    double lapse = 1./psi;
    // double At = m_1d_sol.get_At_interp(r);
    double At = m_1d_sol.get_At_interp(r);

    // double Ex = -(At / lapse) * (x / safe_r)
    //           * (1. / safe_r + dpsi / psi + dlapse * safe_inv_lapse);
    // double Ey = -(At / lapse) * (y / safe_r)
    //           * (1. / safe_r + dpsi / psi + dlapse * safe_inv_lapse);
    // double Ez = -(At / lapse) * (z / safe_r)
    //           * (1. / safe_r + dpsi / psi + dlapse * safe_inv_lapse);
    double Er = At / psi / safe_r;
    double Ex = Er * (x / safe_r);
    double Ey = Er * (y / safe_r);
    double Ez = Er * (z / safe_r);


    if (binary)
    {
        // boosts and coordinate objects
        x = (coords.x + separation / 2.);
        r = sqrt(x * x + y * y + z * z);
        safe_r = sqrt(x * x + y * y + z * z + 10e-20);

        // second star physical variables
        double lapse2 = m_1d_sol.get_lapse_interp(r);
        double psi2 = m_1d_sol.get_psi_interp(r);
        double dpsi2 = m_1d_sol.get_dpsi_interp(r);
        double dlapse2 = m_1d_sol.get_dlapse_interp(r);
        double At2 = m_1d_sol.get_At_interp(r);
        //safe_inv_lapse = lapse2 / (sqrt(lapse2*lapse2) + 10e-10);

        double Ex2 = 0.;
        double Ey2 = 0.;
        double Ez2 = 0.;


        At += At2;
        Ex += Ex2;
        Ey += Ey2;
        Ez += Ez2;
        psi = sqrt(psi * psi + psi2 * psi2);
        lapse = sqrt(lapse * lapse + lapse2 * lapse2);
    }


    vars.chi = pow(psi, -2);

    // can pick either lapse maybe
    //vars.lapse += sqrt(vars.chi);
    vars.lapse += lapse;

    vars.At += 0.*At;

    vars.Ex += Ex;

    vars.Ey += Ey;

    vars.Ez += Ez;

    double kroneka[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};

    FOR2(i, j) vars.h[i][j] = kroneka[i][j];

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* EMDBH_IMPL_HPP_ */