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


    // binary parameters (can be NO binary too)
    bool binary = m_params_EMDBH.binary;
    double separation = m_params_EMDBH.separation;

    // matter conversion from unit conventions
    double root_kappa = 1./sqrt(8.*M_PI);

    // the kroneka delta
    double kroneka_delta[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};



    ////////////////////////////
    // read 1st black hole
    ////////////////////////////

    // coord objects
    double x = coords.x - 0.5 * separation;
    double z = coords.z;
    double y = coords.y;
    double cart_coords[3] = {x, y, z};

    // radii and safe (divisible) radii
    double r = sqrt(x * x + y * y + z * z);
    double safe_r = sqrt(x * x + y * y + z * z + 10e-20);
    double rho = sqrt(x * x + y * y);
    double safe_rho = sqrt(x * x + y * y + 10e-20);

    // trig functions
    double sintheta = rho/safe_r;
    double costheta = z/safe_r;
    double sinphi = y/safe_rho;
    double cosphi = x/safe_rho;

    // jacobeans
    double dx_dr = cosphi*sintheta;
    double dy_dr = sinphi*sintheta;
    double dz_dr = costheta;
    double dr_dx = (x / safe_r);
    double dr_dy = (y / safe_r);
    double dr_dz = (z / safe_r);

    // partial cartesian coords (i) by partial polar coords (i)
    // dxc_dxp[i][j]
    // i = {x,y,z}
    // j = {r,th,ph}
    double dxc_dxp[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    dxc_dxp[0][0] = dx_dr;
    dxc_dxp[1][0] = dy_dr;
    dxc_dxp[2][0] = dz_dr;
    dxc_dxp[0][1] = r * cosphi * costheta;
    dxc_dxp[1][1] = r * sinphi * costheta;
    dxc_dxp[2][1] = -r * sintheta;
    dxc_dxp[0][2] = -r * sinphi * sintheta;
    dxc_dxp[1][2] = r * cosphi * sintheta;
    dxc_dxp[2][2] = 0.; // dz/dphi=0

    // partial polar coords (i) by partial cartesian coords (i)
    // dxp_dxc[i][j]
    // i = {r,th,ph}
    // j = {x,y,z}
    double dxp_dxc[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    dxp_dxc[0][0] = dr_dx;
    dxp_dxc[0][1] = dr_dy;
    dxp_dxc[0][2] = dr_dz;
    dxp_dxc[1][0] = x * z / (safe_r * safe_r * safe_rho);
    dxp_dxc[1][1] = y * z / (safe_r * safe_r * safe_rho);
    dxp_dxc[1][2] = - rho / (safe_r * safe_r);
    dxp_dxc[2][0] = - y / (safe_rho * safe_rho);
    dxp_dxc[2][1] = x / (safe_rho * safe_rho);
    dxp_dxc[2][2] = 0.; // dphi/dz=0

    // X^2 = chi
    double X = m_1d_sol.get_value_interp(m_1d_sol.X,r);

    // loading the upstairs E^r then lower with gamma_rr = 1/(X*X)
    double E_r = m_1d_sol.get_value_interp(m_1d_sol.Er,r)*root_kappa;
    E_r = E_r / ( X * X ) ;

    // fabrizio's mixed conformal traceless curvature, then make downstairs
    // THIS ASSUMES DIAGONAL CONFORMALLY FLAT METRIC WITH h_rr=1
    // angular parts must give trace 0
    double AaaUL = m_1d_sol.get_value_interp(m_1d_sol.Aa,r);
    double Arr = AaaUL / ( X * X );
    double Aij_polar[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

    // this is NOT the conformal A, its simply (K_ij-K gamma_ij/3)
    Aij_polar[0][0] = Arr;
    Aij_polar[1][1] = -0.5*Arr*r*r;
    Aij_polar[2][2] = -0.5*Arr*r*r*sintheta*sintheta;

    // upstairs radial shift
    double betaR = m_1d_sol.get_value_interp(m_1d_sol.shift,r);

    // set gauge vars
    vars.lapse = m_1d_sol.get_value_interp(m_1d_sol.lapse,r);
    vars.shift[0] = dx_dr * betaR;
    vars.shift[1] = dy_dr * betaR;
    vars.shift[2] = dz_dr * betaR;

    // set geometry
    vars.chi = X*X;
    vars.K = m_1d_sol.get_value_interp(m_1d_sol.K,r);
    FOR2(i, j)
    {
        vars.h[i][j] = 0.;
        vars.A[i][j] = 0.;
    }
    // this loop is after previous loop to ensure h_ij diag nonzero
    FOR1(i)
    {
        vars.h[i][i] = 1.;
    }
    FOR4(i,j,m,n)
    {
        // conformal decomposition here with chi
        vars.A[i][j] += X*X * dxp_dxc[m][i] * dxp_dxc[n][j] * Aij_polar[m][n];
    }
    // FOR2(i,j)
    // {
    //     vars.A[i][j] = 0.5 * vars.chi * Arr * (
    //           3. * cart_coords[i] * cart_coords[j] / (safe_r * safe_r)
    //           - kroneka_delta[i][j] );
    // }

    // scalar field
    vars.phi = m_1d_sol.get_value_interp(m_1d_sol.phi,r)*root_kappa;
    // minus sign to match Fabrizio's convention
    vars.Pi = -m_1d_sol.get_value_interp(m_1d_sol.pi,r)*root_kappa;

    // electric field
    vars.Ex = dr_dx * E_r;
    vars.Ey = dr_dy * E_r;
    vars.Ez = dr_dz * E_r;


    ////////////////////////////
    // superpose second black hole if required
    ////////////////////////////

    if (binary)
    {
        // coord objects
        x = coords.x + 0.5 * separation;
        z = coords.z;
        y = coords.y;
        cart_coords[0] = x;
        cart_coords[1] = y;
        cart_coords[2] = z;

        // radii and safe (divisible) radii
        r = sqrt(x * x + y * y + z * z);
        safe_r = sqrt(x * x + y * y + z * z + 10e-20);
        rho = sqrt(x * x + y * y);
        safe_rho = sqrt(x * x + y * y + 10e-20);

        // trig functions
        sintheta = rho/safe_r;
        costheta = z/safe_r;
        sinphi = y/safe_rho;
        cosphi = x/safe_rho;

        // jacobeans
        dx_dr = cosphi*sintheta;
        dy_dr = sinphi*sintheta;
        dz_dr = costheta;
        dr_dx = (x / safe_r);
        dr_dy = (y / safe_r);
        dr_dz = (z / safe_r);

        // partial cartesian coords (i) by partial polar coords (i)
        // dxc_dxp[i][j]
        // i = {x,y,z}
        // j = {r,th,ph}
        dxc_dxp[0][0] = dx_dr;
        dxc_dxp[1][0] = dy_dr;
        dxc_dxp[2][0] = dz_dr;
        dxc_dxp[0][1] = r * cosphi * costheta;
        dxc_dxp[1][1] = r * sinphi * costheta;
        dxc_dxp[2][1] = -r * sintheta;
        dxc_dxp[0][2] = -r * sinphi * sintheta;
        dxc_dxp[1][2] = r * cosphi * sintheta;
        dxc_dxp[2][2] = 0.; // dz/dphi=0

        // partial polar coords (i) by partial cartesian coords (i)
        // dxp_dxc[i][j]
        // i = {r,th,ph}
        // j = {x,y,z}
        dxp_dxc[0][0] = dr_dx;
        dxp_dxc[0][1] = dr_dy;
        dxp_dxc[0][2] = dr_dz;
        dxp_dxc[1][0] = x * z / (safe_r * safe_r * safe_rho);
        dxp_dxc[1][1] = y * z / (safe_r * safe_r * safe_rho);
        dxp_dxc[1][2] = - rho / (safe_r * safe_r);
        dxp_dxc[2][0] = - y / (safe_rho * safe_rho);
        dxp_dxc[2][1] = x / (safe_rho * safe_rho);
        dxp_dxc[2][2] = 0.; // dphi/dz=0

        // X^2 = chi
        double Y = m_1d_sol.get_value_interp(m_1d_sol.X,r);

        // loading the upstairs E^r then lower with gamma_rr = 1/(Y*Y)
        E_r = m_1d_sol.get_value_interp(m_1d_sol.Er,r)*root_kappa;
        E_r = E_r / ( Y * Y ) ;

        // fabrizio's mixed conformal traceless curvature, then make downstairs
        // THIS ASSUMES DIAGONAL CONFORMALLY FLAT METRIC WITH h_rr=1
        // angular parts must give trace 0
        AaaUL = m_1d_sol.get_value_interp(m_1d_sol.Aa,r);
        Arr = AaaUL / ( Y * Y );

        // this is NOT the conformal A, its simply (K_ij-K gamma_ij/3)
        Aij_polar[0][0] = Arr;
        Aij_polar[1][1] = -0.5*Arr*r*r;
        Aij_polar[2][2] = -0.5*Arr*r*r*sintheta*sintheta;

        // upstairs radial shift
        betaR = m_1d_sol.get_value_interp(m_1d_sol.shift,r);

        // set gauge vars
        double alpha2 = m_1d_sol.get_value_interp(m_1d_sol.lapse,r);
        // vars.lapse = sqrt(vars.lapse * vars.lapse + alpha2 * alpha2-1.);

        vars.shift[0] += dx_dr * betaR;
        vars.shift[1] += dy_dr * betaR;
        vars.shift[2] += dz_dr * betaR;

        // set geometry
        vars.chi = (X*X*Y*Y)/(X*X + Y*Y - X*X*Y*Y);
        vars.K += m_1d_sol.get_value_interp(m_1d_sol.K,r);

        // lapse after chi as it needs chi
        vars.lapse = sqrt(vars.chi);

        // no need to set h_ij here as its already fine if using kroneka delta

        FOR4(i,j,m,n)
        {
            // conformal decomposition here with chi
            vars.A[i][j] += Y*Y * dxp_dxc[m][i] * dxp_dxc[n][j] * Aij_polar[m][n];
        }
        // hand done algebra check on transformation
        // FOR2(i,j)
        // {
        //     vars.A[i][j] = 0.5 * vars.chi * Arr * (
        //           3. * cart_coords[i] * cart_coords[j] / (safe_r * safe_r)
        //           - kroneka_delta[i][j] );
        // }

        // scalar field
        vars.phi += m_1d_sol.get_value_interp(m_1d_sol.phi,r)*root_kappa;
        // minus sign to match Fabrizio's convention
        vars.Pi += -m_1d_sol.get_value_interp(m_1d_sol.pi,r)*root_kappa;

        // electric field
        vars.Ex += dr_dx * E_r;
        vars.Ey += dr_dy * E_r;
        vars.Ez += dr_dz * E_r;
    }




    ////////////////////////////
    // zeros (magnetics field and electromag constrinats)
    ////////////////////////////

    vars.Xi = 0.; // magnetic constraint var
    vars.At = 0.; // leccy constraint var

    // magnetic fields
    vars.ax = 0.;
    vars.ay = 0.;
    vars.az = 0.;



    ////////////////////////////
    // flat space spherical harmonic perturbation
    ////////////////////////////

    // reset coords
    x = coords.x;
    z = coords.z;
    y = coords.y;
    r = sqrt(x * x + y * y + z * z);
    safe_r = sqrt(x * x + y * y + z * z + 10e-20);

    // new z-plane radius
    rho = sqrt(x * x + y * y);
    safe_rho = sqrt(x * x + y * y + 10e-20);

    // trig
    sintheta = rho/safe_r;
    costheta = z/safe_r;
    sinphi = y/safe_rho;
    cosphi = x/safe_rho;
    double cos2phi = cosphi * cosphi - sinphi * sinphi;

    //spherial harmonic
    // Y22 is actually Y22 + Y2-2 halved
    double Y20 = 0.25 * sqrt(5./M_PI) * (3. * costheta * costheta - 1.);
    double Y22 = 0.25 * sqrt(7.5/M_PI) * (cos2phi * sintheta * sintheta);

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