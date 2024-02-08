/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(EMDBH_READ_HPP_)
#error "This file should only be included through EMDBH_read.hpp"
#endif

#ifndef EMDBH_READ_IMPL_HPP_
#define EMDBH_READ_IMPL_HPP_

#include "EMDBHSolution_read.hpp" //for EMDBHSolution class

inline EMDBH_read::EMDBH_read(EMDBH_params_t a_params_EMDBH,
                    CouplingFunction::params_t a_params_coupling_function,
                            double a_G_Newton, double a_dx, int a_verbosity)
    : m_dx(a_dx), m_G_Newton(a_G_Newton),
      m_params_EMDBH(a_params_EMDBH),
      m_params_coupling_function(a_params_coupling_function), m_verbosity(a_verbosity)
{
    m_data_path = m_params_EMDBH.data_path;
}

void EMDBH_read::compute_1d_solution(const double max_r)
{
    try
    {
        // set initial parameters and then run the solver
        pout() << "run m_1d_sol.main()" << endl;
        m_1d_sol.main(m_data_path);
        pout() << "completed m_1d_sol.main()" << endl;
    }
    catch (std::exception &exception)
    {
        pout() << exception.what() << "\n";
    }
}

// Compute the value of the initial vars on the grid
template <class data_t> void EMDBH_read::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<EinsteinMaxwellDilatonField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)
    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_EMDBH.star_centre);


    // coordinate objects
    bool binary = m_params_EMDBH.binary;
    double separation = m_params_EMDBH.separation;
    double x = coords.x - 0.5 * separation;
    double z = coords.z;
    double y = coords.y;
    double r = sqrt(x * x + y * y + z * z);
    double safe_r = sqrt(x * x + y * y + z * z + 10e-10);

    // matter conversion from unit conventions
    double root_kappa = 1./sqrt(8.*M_PI);

    //electromagnetism
    double E_r = m_1d_sol.get_value_interp(m_1d_sol.Er,r)*root_kappa; // upstairs
    double X = m_1d_sol.get_value_interp(m_1d_sol.X,r);
    E_r = E_r / ( X * X ) ; // lowered indx with gamma_rr = 1/(X*X)
    // double a_r = m_1d_sol.get_value_interp(m_1d_sol.a_r,r)*root_kappa;
    double E_x = E_r * (x / safe_r);
    double E_y = E_r * (y / safe_r);
    double E_z = E_r * (z / safe_r);


    // check coord transforms laters, for example Ex from Er ..

    vars.chi = X*X;

    vars.lapse = m_1d_sol.get_value_interp(m_1d_sol.lapse,r);

    vars.phi = m_1d_sol.get_value_interp(m_1d_sol.phi,r)*root_kappa;

    vars.Pi = -m_1d_sol.get_value_interp(m_1d_sol.pi,r)*root_kappa;

    vars.Xi = 0.; // magnetic constraint var
    vars.At = 0.; // leccy constraint var

    vars.ax = 0.;
    vars.ay = 0.;
    vars.az = 0.;

    vars.Ex = E_x;
    vars.Ey = E_y;
    vars.Ez = E_z;

    vars.K = m_1d_sol.get_value_interp(m_1d_sol.K,r);

    FOR2(i, j)
    {
        vars.h[i][j] = 0.;
        vars.A[i][j] = 0.;
    }

    // raise indices of shift with inverse metric
    FOR1(i)
    {
        vars.shift[i] = 0.;
        vars.h[i][i] = 1.;
    }


    if (binary)
    {
        x = coords.x + 0.5 * separation;
        z = coords.z;
        y = coords.y;
        r = sqrt(x * x + y * y + z * z);
        safe_r = sqrt(x * x + y * y + z * z + 10e-10);

        // matter conversion from unit conventions
        root_kappa = 1./sqrt(8.*M_PI);

        //electromagnetism
        E_r = m_1d_sol.get_value_interp(m_1d_sol.Er,r)*root_kappa; // upstairs
        double Y = m_1d_sol.get_value_interp(m_1d_sol.X,r);
        E_r = E_r / ( Y * Y ) ; // lowered indx with gamam_rr = 1/(X*X)
        // double a_r = m_1d_sol.get_value_interp(m_1d_sol.a_r,r)*root_kappa;
        E_x = E_r * (x / safe_r);
        E_y = E_r * (y / safe_r);
        E_z = E_r * (z / safe_r);

        // check coord transforms laters, for example Ex from Er ..

        vars.chi = X*X*Y*Y/(X*X+Y*Y);

        vars.lapse = sqrt(vars.chi);

        vars.phi += m_1d_sol.get_value_interp(m_1d_sol.phi,r)*root_kappa;

        vars.Pi += -m_1d_sol.get_value_interp(m_1d_sol.pi,r)*root_kappa;

        vars.Ex += E_x;
        vars.Ey += E_y;
        vars.Ez += E_z;

        vars.K += m_1d_sol.get_value_interp(m_1d_sol.K,r);

    }

    ////////////////////////////
    // flat space spherical harmonic perturbation
    ////////////////////////////

    // reset coords
    x = coords.x;
    z = coords.z;
    y = coords.y;
    r = sqrt(x * x + y * y + z * z);
    safe_r = sqrt(x * x + y * y + z * z + 10e-10);
    // new coords for angles
    double rho = sqrt(x * x + y * y);
    double safe_rho = sqrt(x * x + y * y + 100e-10);
    double sintheta = rho/safe_r;
    double costheta = z/safe_r;
    double sinphi = y/safe_rho;
    double cosphi = x/safe_rho;
    double cos2phi = cosphi * cosphi - sinphi * sinphi;

    //spherial harmonic
    double Y20 = 0.25 * sqrt(5./M_PI) * (3. * costheta * costheta - 1.);
    double Y22 = 0.25 * sqrt(7.5/M_PI) * (cos2phi * sintheta * sintheta);
    // Y22 is actually Y22 + Y2-2 halved

    // inward travelling thin shell params
    // amplitude is A/r in flat space
    double A = m_params_EMDBH.Ylm_amplitude;
    // standard deviation, width of gaussian shell
    double sig = m_params_EMDBH.Ylm_thickness;
    // initial radiause of shell
    double r_0 = m_params_EMDBH.Ylm_r0;

    double psi = (Y22*A/safe_r) * exp(-(r-r_0)*(r-r_0)/(2.*sig*sig));

    vars.phi += psi;
    // boosted inwards with Pi - flat space approx no outgoing wave
    vars.Pi += psi * (r-r_0) / (sig * sig);



    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* EMDBH_READ_IMPL_HPP_ */
