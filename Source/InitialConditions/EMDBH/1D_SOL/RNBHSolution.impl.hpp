/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(RNBHSOLUTION_HPP_)
#error "This file should only be included through RNBHSolution.hpp"
#endif

#ifndef RNBHSOLUTION_IMPL_HPP_
#define RNBHSOLUTION_IMPL_HPP_

RNBHSolution::RNBHSolution() {}

void RNBHSolution::main()
{
    // create the Reissner Nordstrom solution
    double r=0;
    double eps_0 = 1. ; //permitivity of free space, natural units -> 1
    Q = Qratio / sqrt(8. * M_PI); // note used
    for (int i=0; i<gridsize; i++)
    {
        r = i*dx;
        if (i==0) r = 10e-10;
        omega[i] = (4. * r * r - M * M + Qratio * Qratio)
                      /(4. * r * r + 4. * M * r + M * M - Qratio * Qratio);
        psi[i] = 1. + M / r + (M * M - Qratio * Qratio)/(4. * r * r);
        At[i] =  -Qratio / (sqrt(8. * M_PI) * r ) ;
    }

}

void RNBHSolution::main_polar_areal()
{
    // create the Reissner Nordstrom solution
    double r=0;
    double eps_0 = 1. ; //permitivity of free space, natural units -> 1
    Q = Qratio / sqrt(8. * M_PI);
    for (int i=0; i<gridsize; i++)
    {
        r = i*dx;
        if (i==0) r = 10e-10;
        At[i] = Q / r;
        omega[i] = 1.  -  2. * G * M / r  +  G * Q * Q / (4. * M_PI * r * r );
        psi[i] = 1./omega[i];
    }

}

void RNBHSolution::main_iso_sc()
{
    // create the Reissner Nordstrom solution
    double r=0;
    for (int i=0; i<gridsize; i++)
    {
        r = i*dx;
        if (i==0) r = 10e-10;
        At[i] = 0.;
        omega[i] = (M - 2. * r ) / (M + 2. * r );
        psi[i] = pow( 1. + M / (2. * r)  , 2 );
    }

}


// 4th order error (cubic interpolation) for field. shouts if asked to fetch a
// value outside the ode solution
double RNBHSolution::get_At_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 =
        ((iter == 0)
             ? At[1]
             : At[iter - 1]); // conditionl/ternary imposing zero gradeint at r=0
    f2 = At[iter];
    f3 = At[iter + 1];
    f4 = At[iter + 2];

    if (iter > gridsize - 3)
    {
        std::cout << "Requested Value outside BS initial data domain!"
                  << std::endl;
    }

    // do the cubic spline, from mathematica script written by Robin
    // (rc634@cam.ac.uk)
    interpolated_value =
        (1. / 48.) *
        (f1 * (-3. + 2. * a + 12. * a * a - 8. * a * a * a) +
         (3. + 2. * a) *
             (-(1. + 2. * a) * (-9. * f3 + f4 + 6. * f3 * a - 2 * f4 * a) +
              3. * f2 * (3. - 8. * a + 4. * a * a)));
    return interpolated_value;
}


double RNBHSolution::get_lapse_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter == 0) ? omega[1]
                      : omega[iter - 1]); // conditionl/ternary imposing zero
                                          // gradeint at r=0
    f2 = omega[iter];
    f3 = omega[iter + 1];
    f4 = omega[iter + 2];

    if (iter > gridsize - 3)
    {
        std::cout << "Requested Value outside RN initial data domain!"
                  << std::endl;
    }

    // do the cubic spline, from mathematica script written by Robin
    // (rc634@cam.ac.uk)
    interpolated_value =
        (1. / 48.) *
        (f1 * (-3. + 2. * a + 12. * a * a - 8. * a * a * a) +
         (3. + 2. * a) *
             (-(1. + 2. * a) * (-9. * f3 + f4 + 6. * f3 * a - 2 * f4 * a) +
              3. * f2 * (3. - 8. * a + 4. * a * a)));
    return interpolated_value;
}



double RNBHSolution::get_dlapse_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter == 0) ? omega[1] : omega[iter - 1]);
    f2 = omega[iter];
    f3 = omega[iter + 1];
    f4 = omega[iter + 2];

    if (iter > gridsize - 3)
    {
        std::cout << "Requested Value outside RN initial data domain!"
                  << std::endl;
    }

    // do the cubic spline (for gradient now), from mathematica script written
    // by Robin (rc634@cam.ac.uk)
    interpolated_value =
        (1. / (24. * dx)) *
        ((f1 - 27. * f2 + 27. * f3 - f4) + 12. * a * (f1 - f2 - f3 + f4) -
         12. * a * a * (f1 - 3. * f2 + 3. * f3 - f4));
    return interpolated_value;
}

double RNBHSolution::get_psi_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter == 0) ? psi[1]
                      : psi[iter - 1]); // conditionl/ternary imposing zero
                                          // gradeint at r=0
    f2 = psi[iter];
    f3 = psi[iter + 1];
    f4 = psi[iter + 2];

    if (iter > gridsize - 3)
    {
        std::cout << "Requested Value outside RN initial data domain!"
                  << std::endl;
    }

    // do the cubic spline, from mathematica script written by Robin
    // (rc634@cam.ac.uk)
    interpolated_value =
        (1. / 48.) *
        (f1 * (-3. + 2. * a + 12. * a * a - 8. * a * a * a) +
         (3. + 2. * a) *
             (-(1. + 2. * a) * (-9. * f3 + f4 + 6. * f3 * a - 2 * f4 * a) +
              3. * f2 * (3. - 8. * a + 4. * a * a)));
    return interpolated_value;
}



double RNBHSolution::get_dpsi_interp(const double r) const
{
    int iter = (int)floor(
        r / dx); // index of 2nd (out of 4) gridpoints used for interpolation
    double a =
        (r / dx) - floor(r / dx) - 0.5; // fraction from midpoint of two values,
                                        // a = +- 1/2 is the nearest gridpoints
    double interpolated_value = 0, f1, f2, f3, f4;
    f1 = ((iter == 0) ? psi[1] : psi[iter - 1]);
    f2 = psi[iter];
    f3 = psi[iter + 1];
    f4 = psi[iter + 2];

    if (iter > gridsize - 3)
    {
        std::cout << "Requested Value outside RN initial data domain!"
                  << std::endl;
    }

    // do the cubic spline (for gradient now), from mathematica script written
    // by Robin (rc634@cam.ac.uk)
    interpolated_value =
        (1. / (24. * dx)) *
        ((f1 - 27. * f2 + 27. * f3 - f4) + 12. * a * (f1 - f2 - f3 + f4) -
         12. * a * a * (f1 - 3. * f2 + 3. * f3 - f4));
    return interpolated_value;
}

void RNBHSolution::set_initialcondition_params(
    EMDBH_params_t m_params_EMDBH,
    CouplingFunction::params_t m_params_coupling_function, const double max_r)
{
    gridsize = m_params_EMDBH.gridpoints;
    At.resize(gridsize);            // scalar field modulus
    omega.resize(gridsize);        // lapse
    psi.resize(gridsize);        // conformal factor

    G = m_params_EMDBH.Newtons_constant;
    L = max_r * 1.05; // just to make sure the function domain is slightly
                      // larger than the required cube
    dx = L / (gridsize - 1);
    M = m_params_EMDBH.bh_mass;
    Qratio = M*m_params_EMDBH.bh_charge;
}

#endif /* RNBHSOLUTION_IMPL_HPP_ */