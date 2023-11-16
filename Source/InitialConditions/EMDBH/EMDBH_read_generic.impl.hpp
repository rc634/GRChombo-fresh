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
}

void EMDBH_read::compute_1d_solution(const double max_r)
{
    try
    {
        // set initial parameters and then run the solver
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
template <class data_t> void EMDBH_read::compute(Cell<data_t> current_cell) const
{
    MatterCCZ4<EinsteinMaxwellDilatonField<>>::Vars<data_t> vars;
    // Load variables (should be set to zero if this is a single BS)
    current_cell.load_vars(vars);
    // VarsTools::assign(vars, 0.); // Set only the non-zero components below
    Coordinates<data_t> coords(current_cell, m_dx,
                               m_params_EMDBH.star_centre);


    // boosts and coordinate objects
    double x = coords.x;
    double z = coords.z;
    double y = coords.y;
    double r = sqrt(x * x + y * y + z * z);
    double safe_r = sqrt(x * x + y * y + z * z + 10e-10);
    double rho = sqrt(x * x + y * y);
    double safe_rho = sqrt(x * x + y * y + 10e-10);
    double sintheta = rho/safe_r;
    double costheta = z/safe_r;
    double sinphi = y/safe_rho;
    double cosphi = x/safe_rho;

    // matter conversion from unit conventions
    double root_kappa = 1./sqrt(8.*M_PI);

    // spatial metric and conformal factor
    double X = m_1d_sol.get_value_interp(m_1d_sol.X,r);
    //double dX = m_1d_sol.get_deriv_interp(m_1d_sol.X,r);
    double a = m_1d_sol.get_value_interp(m_1d_sol.a,r);
    double b = m_1d_sol.get_value_interp(m_1d_sol.b,r);
    double chi = X*X/pow(a*b*b,1./3.);
    double h_rr = pow(a/b,2./3.);
    double h_thth = r*r*pow(a/b,-1./3.);
    double h_phph = sintheta*sintheta*h_thth;
    double h_polar[3][3] = {{h_rr, 0., 0.}, {0., h_thth, 0.}, {0., 0., h_phph}};

    // // spatial metric again
    // gamma_polar_LL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    // gamma_polar_LL[0][0] = a/X/X;
    // gamma_polar_LL[1][1] = r*r*b/X/X;
    // gamma_polar_LL[2][2] = r*r*b*sintheta*sintheta/X/X;
    // gamma_polar_UU[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    // gamma_polar_UU[0][0] = 1./gamma_polar_LL[0][0];
    // gamma_polar_UU[1][1] = 1./gamma_polar_LL[1][1];
    // gamma_polar_UU[2][2] = 1./gamma_polar_LL[2][2];
    // gamma_cart_LL[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
    // gamma_cart_UU[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};


    // polar to cart jacobean
    double J[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}; // Jacobean d polar d cart
    J[0][0] = x/safe_r; // dr\dx
    J[0][1] = y/safe_r; // dr\dy
    J[0][2] = z/safe_r; // dr\dz
    J[1][0] = pow(safe_r,-2)*x*z/safe_rho; // dth\dx
    J[1][1] = pow(safe_r,-2)*y*z/safe_rho; // dth\dy
    J[1][2] = -pow(safe_r,-2)*rho; // dth\dz
    J[2][0] = -y/(safe_rho*safe_rho); // dph\dx
    J[2][1] = x/(safe_rho*safe_rho); // dph\dy
    J[2][2] = 0.; // dph\dz




    // gauge
    double lapse = m_1d_sol.get_value_interp(m_1d_sol.lapse,r);
    double shift = m_1d_sol.get_value_interp(m_1d_sol.shift,r);
    double beta_r = shift * chi / h_rr;
    double beta_L_cart[3] = {0.,0.,0.};
    beta_L_cart[0] = x * beta_r / safe_r;
    beta_L_cart[1] = y * beta_r / safe_r;
    beta_L_cart[2] = z * beta_r / safe_r;

    //electromagnetism
    double E_r = m_1d_sol.get_value_interp(m_1d_sol.Er,r)*root_kappa; // upstairs
    E_r = E_r * chi / h_rr;
    // double a_r = m_1d_sol.get_value_interp(m_1d_sol.a_r,r)*root_kappa;
    double E_x = E_r * (x / safe_r);
    double E_y = E_r * (y / safe_r);
    double E_z = E_r * (z / safe_r);
    // assumed no magnetic as B^i = eps^{ijk}(partial_i a_j - partial_j a_i)
    // and only nonzero component of deriv is i=j=r and cancels on eps

    // cruvature
    // components are up and down here, A^i_j
    double AUL_polar[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    AUL_polar[0][0] = pow(a*b*b,1./3.)*m_1d_sol.get_value_interp(m_1d_sol.Aa,r);
    AUL_polar[1][1] = -0.5*AUL_polar[0][0];
    AUL_polar[2][2] = AUL_polar[1][1];
    //put components to downstairs A_ij
    double ALL_polar[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    ALL_polar[0][0] = AUL_polar[0][0] * h_polar[0][0];
    ALL_polar[1][1] = AUL_polar[1][1] * h_polar[1][1];
    ALL_polar[2][2] = AUL_polar[2][2] * h_polar[2][2];



    // check coord transforms laters, for example Ex from Er ..

    vars.chi = chi;

    vars.lapse = lapse;

    //vars.shift[0] = shift * (x / safe_r);
    //vars.shift[1] = shift * (y / safe_r);
    //vars.shift[2] = shift * (z / safe_r);

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

    FOR4(i,j,k,l)
    {
        vars.h[i][j] += h_polar[k][l] * J[k][i] * J[l][j];
        vars.A[i][j] += ALL_polar[k][l] * J[k][i] * J[l][j];
    }

    using namespace TensorAlgebra;
    auto h_UU = compute_inverse_sym(vars.h);

    // raise indices of shift with inverse metric
    FOR1(i) vars.shift[i] = 0.;
    //FOR2(i,j) vars.shift[i] += beta_L_cart[j] * h_UU[i][j] * vars.chi;

    // Store the initial values of the variables
    current_cell.store_vars(vars);
}

#endif /* EMDBH_READ_IMPL_HPP_ */
