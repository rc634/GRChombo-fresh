/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(EINSTEINMAXWELLDILATONFIELD_HPP_)
#error "This file should only be included through EinsteinMaxwellDilatonField.hpp"
#endif

#ifndef EINSTEINMAXWELLDILATONFIELD_IMPL_HPP_
#define EINSTEINMAXWELLDILATONFIELD_IMPL_HPP_
#include "DimensionDefinitions.hpp"

// Calculate the stress energy tensor elements
template <class coupling_t>
template <class data_t, template <typename> class vars_t>
emtensor_t<data_t> EinsteinMaxwellDilatonField<coupling_t>::compute_emtensor(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<2, data_t> &h_UU, const Tensor<3, data_t> &chris_ULL) const
{
    emtensor_t<data_t> out;

    // Copy the field vars into SFObject
    EMDObject<data_t> vars_emd;
    vars_emd.phi = vars.phi;
    vars_emd.Pi = vars.Pi;
    vars_emd.At = vars.At;
    vars_emd.ax = vars.ax;
    vars_emd.ay = vars.ay;
    vars_emd.az = vars.az;
    vars_emd.Ex = vars.Ex;
    vars_emd.Ey = vars.Ey;
    vars_emd.Ez = vars.Ez;

    // call the function which computes the em tensor excluding the coupling
    emtensor_excl_coupling(out, vars, vars_emd, d1.phi, h_UU,
                            chris_ULL);

    // declare the default coupling values - and some useful auxiliary variables
    data_t f_of_phi = 0.0;
    data_t f_prime_of_phi = 0.0;
    data_t coupling_of_phi = 1.0;

    // some useful variables
    Tensor<1, data_t, 3> a;
    Tensor<1, data_t, 3> Ei;
    Tensor<2, data_t, 3> Bij;
    data_t BB = 0.; // squared variables
    data_t EE = 0.; // squared variables
    data_t PiPi = 0.; // squared variables
    data_t DphiDphi = 0.; // squared variables

    a[0] = vars_emd.ax;
    a[1] = vars_emd.ay;
    a[2] = vars_emd.az;

    Ei[0] = vars_emd.Ex;
    Ei[1] = vars_emd.Ey;
    Ei[2] = vars_emd.Ez;

    Bij[0][0] = 0.;
    Bij[0][1] = d1.ay[0] - d1.ax[1];
    Bij[0][2] = d1.az[0] - d1.ax[2];
    Bij[1][0] = -Bij[0][1];
    Bij[1][1] = 0.;
    Bij[1][2] = d1.az[1] - d1.ay[2];
    Bij[2][0] = -Bij[0][2];
    Bij[1][2] = - Bij[2][1];
    Bij[2][2] = 0.;



    //compute some variables

    PiPi = vars_emd.Pi*vars_emd.Pi;

    FOR2(i, j)
    {
        DphiDphi += d1.Pi[i]*d1.Pi[j]*h_UU[i][j]*vars.chi;
        EE += Ei[i]*Ei[j]*h_UU[i][j]*vars.chi;
    }

    FOR4(i,j,k,l)
    {
        BB += Bij[i][j]*Bij[k][l]*h_UU[i][k]*h_UU[j][l]*vars.chi*vars.chi;
    }



    // compute coupling and add constributions to EM Tensor
    my_coupling.compute_coupling(f_of_phi,f_prime_of_phi,coupling_of_phi,vars);

    out.rho += coupling_of_phi * (BB + 2. * EE);

    FOR2(i, j)
    {
        out.Sij[i][j] += coupling_of_phi
                      * ( -4. * Ei[i] * Ei[j] - vars.h[i][j] * (BB - 2.*EE))
                                /(vars.chi*vars.chi);
    }
    FOR4(i,j,k,l)
    {
        out.Sij[i][j] += 4. * coupling_of_phi
                        * Bij[i][k]*Bij[j][l]*h_UU[k][l]*vars.chi;
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    FOR3(i,j,k)
    {
        out.Si[i] += 4. * coupling_of_phi
                      * Ei[j] * Bij[i][k] * h_UU[k][j] * vars.chi;
    }

    return out;


}

// Calculate the stress energy tensor elements
template <class coupling_t>
template <class data_t, template <typename> class vars_t>
void EinsteinMaxwellDilatonField<coupling_t>::emtensor_excl_coupling(
    emtensor_t<data_t> &out, const vars_t<data_t> &vars,
    const EMDObject<data_t> &vars_emd, const Tensor<1, data_t> &d1_phi,
    const Tensor<2, data_t> &h_UU,
    const Tensor<3, data_t> &chris_ULL)
{

    // compute some useful quantities
    data_t PiPi = vars_emd.Pi*vars_emd.Pi;
    data_t DphiDphi = 0.;

    FOR2(i, j)
    {
        DphiDphi += d1_phi[i]*d1_phi[j]*h_UU[i][j]*vars.chi;
    }

    // Calculate components of EM Tensor
    // S_ij = T_ij
    FOR2(i, j)
    {
        out.Sij[i][j] = 4. * d1_phi[i] * d1_phi[j]
                      - (2./vars.chi)*vars.h[i][j] * (DphiDphi - PiPi);
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // S_i (note lower index) = n^a T_a0
    FOR1(i)
    {
        out.Si[i] = 4.*vars_emd.Pi * d1_phi[i];
    }

    // rho = n^a n^b T_ab
    out.rho = 2. * (PiPi + DphiDphi);
}

// Adds in the RHS for the matter vars
template <class coupling_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void EinsteinMaxwellDilatonField<coupling_t>::add_matter_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const vars_t<data_t> &advec) const
{
    EMDObject<data_t> rhs_emd;
    // advection terms
    EMDObject<data_t> advec_emd;
    advec_emd.phi = advec.phi;
    advec_emd.Pi = advec.Pi;
    advec_emd.At = advec.At;
    advec_emd.ax = advec.ax;
    advec_emd.ay = advec.ay;
    advec_emd.az = advec.az;
    advec_emd.Ex = advec.Ex;
    advec_emd.Ey = advec.Ey;
    advec_emd.Ez = advec.Ez;

    // the vars
    EMDObject<data_t> vars_emd;
    vars_emd.phi = vars.phi;
    vars_emd.Pi = vars.Pi;
    vars_emd.At = vars.At;
    vars_emd.ax = vars.ax;
    vars_emd.ay = vars.ay;
    vars_emd.az = vars.az;
    vars_emd.Ex = vars.Ex;
    vars_emd.Ey = vars.Ey;
    vars_emd.Ez = vars.Ez;

    // call the function for the rhs excluding the coupling
    matter_rhs_excl_coupling(rhs_emd, vars, vars_emd, d1,
                              d1.phi, d1.ax, d1.ay, d1.az,
                              d2.phi, d2.ax, d2.ay, d2.az,
                                          d1.At, advec_emd);

    const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);

    // set the default coupling values
    data_t f_of_phi = 0.;
    data_t f_prime_of_phi = 0.;
    data_t coupling_of_phi = 1.;

    // some useful variables
    Tensor<1, data_t, 3> a;
    Tensor<1, data_t, 3> Ei;
    Tensor<2, data_t, 3> Bij;
    data_t DphiDphi = 0.;
    data_t BB = 0.; // squared variables
    data_t EE = 0.; // squared variables

    a[0] = vars_emd.ax;
    a[1] = vars_emd.ay;
    a[2] = vars_emd.az;

    Ei[0] = vars_emd.Ex;
    Ei[1] = vars_emd.Ey;
    Ei[2] = vars_emd.Ez;

    Bij[0][0] = 0.;
    Bij[0][1] = d1.ay[0] - d1.ax[1];
    Bij[0][2] = d1.az[0] - d1.ax[2];
    Bij[1][0] = -Bij[0][1];
    Bij[1][1] = 0.;
    Bij[1][2] = d1.az[1] - d1.ay[2];
    Bij[2][0] = -Bij[0][2];
    Bij[1][2] = - Bij[2][1];
    Bij[2][2] = 0.;


    FOR2(i, j)
    {
        DphiDphi += d1.Pi[i]*d1.Pi[j]*h_UU[i][j]*vars.chi;
        EE += Ei[i]*Ei[j]*h_UU[i][j]*vars.chi;
    }

    FOR4(i,j,k,l)
    {
        BB += Bij[i][j]*Bij[k][l]*h_UU[i][k]*h_UU[j][l]*vars.chi*vars.chi;
    }

    // compute coupling and add constributions to EOM
    my_coupling.compute_coupling(f_of_phi,f_prime_of_phi,coupling_of_phi,vars);

    // adjust RHS for the coupling term
    total_rhs.phi = rhs_emd.phi;

    total_rhs.Pi = rhs_emd.Pi - 0.5 * f_prime_of_phi * vars.lapse *
                                          coupling_of_phi * (BB - 2. * EE);

    total_rhs.At = rhs_emd.At;

    total_rhs.ax = rhs_emd.ax;

    total_rhs.ay = rhs_emd.ay;

    total_rhs.az = rhs_emd.az;

    total_rhs.Ex = rhs_emd.Ex - 2. * f_prime_of_phi * vars_emd.ax * vars_emd.Pi;

    total_rhs.Ey = rhs_emd.Ey - 2. * f_prime_of_phi * vars_emd.ay * vars_emd.Pi;

    total_rhs.Ez = rhs_emd.Ez - 2. * f_prime_of_phi * vars_emd.az * vars_emd.Pi;

    FOR2(j,k)
    {
        total_rhs.Ex += 2. * vars.lapse * f_prime_of_phi
                           * vars.chi * h_UU[j][k] * Bij[k][0] * d1.phi[j];
        total_rhs.Ey += 2. * vars.lapse * f_prime_of_phi
                           * vars.chi * h_UU[j][k] * Bij[k][1] * d1.phi[j];
        total_rhs.Ez += 2. * vars.lapse * f_prime_of_phi
                           * vars.chi * h_UU[j][k] * Bij[k][2] * d1.phi[j];
    }

}

// the RHS excluding the coupling terms
template <class coupling_t>
template <class data_t, template <typename> class vars_t>
void EinsteinMaxwellDilatonField<coupling_t>::matter_rhs_excl_coupling(
    EMDObject<data_t> &rhs_emd, const vars_t<data_t> &vars,
    const EMDObject<data_t> &vars_emd, const vars_t<Tensor<1, data_t>> &d1,
    const Tensor<1, data_t> &d1_phi,
    const Tensor<1, data_t> &d1_ax,
    const Tensor<1, data_t> &d1_ay,
    const Tensor<1, data_t> &d1_az,
    const Tensor<2, data_t> &d2_phi,
    const Tensor<2, data_t> &d2_ax,
    const Tensor<2, data_t> &d2_ay,
    const Tensor<2, data_t> &d2_az,
    const Tensor<1, data_t> &d1_At,
    const EMDObject<data_t> &advec_emd)
{
    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);
    const auto chris_phys =
                 compute_phys_chris(d1.chi, vars.chi, vars.h, h_UU, chris.ULL);

    // some useful variables
    Tensor<1, data_t, 3> a;
    Tensor<1, data_t, 3> Ei;
    Tensor<2, data_t, 3> Bij;
    Tensor<2, data_t, 3> Da_ij;
    Tensor<3, data_t, 3> DBijk;
    data_t div_a = 0.; // divergence of a

    a[0] = vars_emd.ax;
    a[1] = vars_emd.ay;
    a[2] = vars_emd.az;

    Ei[0] = vars_emd.Ex;
    Ei[1] = vars_emd.Ey;
    Ei[2] = vars_emd.Ez;

    Da_ij[0][0] = d1_ax[0]; // D_i a_j
    Da_ij[1][0] = d1_ax[1];
    Da_ij[2][0] = d1_ax[2];
    Da_ij[0][1] = d1_ay[0];
    Da_ij[1][1] = d1_ay[1];
    Da_ij[2][1] = d1_ay[2];
    Da_ij[0][2] = d1_az[0];
    Da_ij[1][2] = d1_az[1];
    Da_ij[2][2] = d1_az[2];

    FOR3(i,j,k)
    {
        Da_ij[i][j] += - chris_phys[k][i][j]*a[k];

        DBijk[i][j][k] = 0.; // D_i B_jk
    }
    FOR2(i,j)
    {
        div_a += Da_ij[i][j]*h_UU[i][j]*vars.chi;

        Bij[i][j] = Da_ij[i][j] - Da_ij[j][i];

        DBijk[i][j][0] += d2_ax[i][j];
        DBijk[i][0][j] += -d2_ax[i][j];
        DBijk[i][j][1] += d2_ay[i][j];
        DBijk[i][1][j] += -d2_ay[i][j];
        DBijk[i][j][2] += d2_az[i][j];
        DBijk[i][2][j] += -d2_az[i][j];
    }

    FOR4(i,j,k,l)
    {
        DBijk[i][j][k] += - chris_phys[l][i][k] * Bij[j][l]
                          - chris_phys[l][i][j] * Bij[l][k];
    }


    rhs_emd.phi = advec_emd.phi - vars.lapse * vars_emd.Pi;

    rhs_emd.Pi = advec_emd.Pi + vars.lapse * vars.K * vars_emd.Pi;

    rhs_emd.At = advec_emd.At + vars.lapse * (vars_emd.At * vars.K - div_a);

    rhs_emd.ax = advec_emd.ax - vars_emd.At * d1.lapse[0]
                     -vars.lapse * ( vars_emd.Ex + d1_At[0]);

    rhs_emd.ay = advec_emd.ay - vars_emd.At * d1.lapse[1]
                     -vars.lapse * ( vars_emd.Ey + d1_At[1]);

    rhs_emd.az = advec_emd.az - vars_emd.At * d1.lapse[2]
                     -vars.lapse * ( vars_emd.Ez + d1_At[2]);

    rhs_emd.Ex = advec_emd.Ex + vars.lapse * vars_emd.Ex * vars.K;

    rhs_emd.Ey = advec_emd.Ey + vars.lapse * vars_emd.Ey * vars.K;

    rhs_emd.Ez = advec_emd.Ez + vars.lapse * vars_emd.Ez * vars.K;


    FOR2(j,k)
    {
        rhs_emd.Pi += - vars.chi * h_UU[k][j]
                    * (d1_phi[k] * d1.lapse[j] + vars.lapse * d2_phi[k][j]);

        rhs_emd.At += - a[k] * d1.lapse[j] * h_UU[k][j] * vars.chi;

        rhs_emd.Ex += -vars.chi * h_UU[j][k]
                    * ( Bij[k][0] * d1.lapse[j] + vars.lapse * DBijk[j][k][0] );
        rhs_emd.Ey += -vars.chi * h_UU[j][k]
                    * ( Bij[k][1] * d1.lapse[j] + vars.lapse * DBijk[j][k][1] );
        rhs_emd.Ez += -vars.chi * h_UU[j][k]
                    * ( Bij[k][2] * d1.lapse[j] + vars.lapse * DBijk[j][k][2] );
    }


}

#endif /* EINSTEINMAXWELLDILATONFIELD_IMPL_HPP_ */
