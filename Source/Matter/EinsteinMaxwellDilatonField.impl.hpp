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
    vars_emd.Xi = vars.Xi;

    // call the function which computes the em tensor excluding the coupling
    emtensor_excl_coupling(out, vars, vars_emd, d1.phi, h_UU,
                            chris_ULL);

    // declare the default coupling values - and some useful auxiliary variables
    data_t f_of_phi = 0.0;
    data_t f_prime_of_phi = 0.0;
    data_t coupling_of_phi = 1.0;

    // making epsion_i^{\,\,jk}
    auto eps_symbol = TensorAlgebra::epsilon(); //not the tensor, the symbol
    Tensor<3, data_t, 3> epsUUU; //tenosr
    Tensor<3, data_t, 3> epsLUU; //tensor

    // make all upstairs first, divide by root gamma
    FOR3(i,j,k)
    {
      epsUUU[i][j][k] = eps_symbol[i][j][k] / pow(vars.chi,-1.5);
      epsLUU[i][j][k] = 0.;
    }
    // lower first index with physical metric
    FOR4(i,j,k,l) epsLUU[i][j][k] += epsUUU[l][j][k] * vars.h[i][l] / vars.chi;

    // some useful variables
    Tensor<1, data_t, 3> Ei;
    Tensor<1, data_t, 3> Bi;
    data_t BB = 0.; // squared variables
    data_t EE = 0.; // squared variables
    data_t PiPi = 0.; // squared variables
    data_t DphiDphi = 0.; // squared variables

    Bi[0] = vars_emd.ax;
    Bi[1] = vars_emd.ay;
    Bi[2] = vars_emd.az;

    Ei[0] = vars_emd.Ex;
    Ei[1] = vars_emd.Ey;
    Ei[2] = vars_emd.Ez;


    //compute some variables

    PiPi = vars_emd.Pi*vars_emd.Pi;

    FOR2(i, j)
    {
        DphiDphi += d1.phi[i]*d1.phi[j]*h_UU[i][j]*vars.chi;
        EE += Ei[i]*Ei[j]*h_UU[i][j]*vars.chi;
        BB += Bi[i]*Bi[j]*h_UU[i][j]*vars.chi;
    }


    // compute coupling and add constributions to EM Tensor
    my_coupling.compute_coupling(f_of_phi,f_prime_of_phi,coupling_of_phi,vars);

    //out.rho += coupling_of_phi * (BB / 2. + EE);
    out.rho += (BB + EE);

    FOR2(i, j)
    {
        // out.Sij[i][j] += coupling_of_phi * ( -2. * Ei[i] * Ei[j]
        //               -  vars.h[i][j] * (BB / 2. - EE) / vars.chi );
        out.Sij[i][j] +=  - 2. * (Ei[i] * Ei[j] + Bi[i] * Bi[j])
                          + vars.h[i][j] * (BB + EE) / vars.chi ;
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    FOR3(i,j,k)
    {
        out.Si[i] += 2. * epsLUU[i][j][k] * Ei[j] * Bi[k];
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
    data_t PiPi = 0.;
    data_t DphiDphi = 0.;


    PiPi = vars_emd.Pi*vars_emd.Pi;

    FOR2(i, j)
    {
        DphiDphi += d1_phi[i]*d1_phi[j]*h_UU[i][j]*vars.chi;
    }

    // Calculate components of EM Tensor

    // S_ij = T_ij
    FOR2(i, j)
    {
        // out.Sij[i][j] = 2. * d1_phi[i] * d1_phi[j]
        //               - vars.h[i][j] * (DphiDphi - PiPi) / vars.chi;
        out.Sij[i][j] = 0.;
    }

    // S = Tr_S_ij
    // out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);
    out.S = 0.;

    // S_i (note lower index) = n^a T_a0
    FOR1(i)
    {
        // out.Si[i] = 2.*vars_emd.Pi * d1_phi[i];
        out.Si[i] = 0.;
    }

    // rho = n^a n^b T_ab
    // out.rho = PiPi + DphiDphi;
    out.rho = 0.;
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
    advec_emd.Xi = advec.Xi;

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
    vars_emd.Xi = vars.Xi;

    // call the function for the rhs excluding the coupling
    matter_rhs_excl_coupling(rhs_emd, vars, vars_emd, d1,
                              d1.phi, d1.ax, d1.ay, d1.az,
                              d2.phi, d1.Ex, d1.Ey, d1.Ez,
                                  d1.At, d1.Xi, advec_emd);

    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);
    const auto chris_phys =
                 compute_phys_chris(d1.chi, vars.chi, vars.h, h_UU, chris.ULL);

    // set the default coupling values
    data_t f_of_phi = 0.;
    data_t f_prime_of_phi = 0.;
    data_t coupling_of_phi = 1.;

    // some useful variables
    Tensor<1, data_t, 3> Ei;
    Tensor<1, data_t, 3> Bi;
    Tensor<2, data_t, 3> dEij; // partial_i E_j
    data_t DphiDphi = 0.;
    data_t divE = 0.;
    data_t BB = 0.; // squared variables
    data_t EE = 0.; // squared variables

    Bi[0] = vars_emd.ax;
    Bi[1] = vars_emd.ay;
    Bi[2] = vars_emd.az;

    Ei[0] = vars_emd.Ex;
    Ei[1] = vars_emd.Ey;
    Ei[2] = vars_emd.Ez;

    dEij[0][0] = d1.Ex[0];
    dEij[1][0] = d1.Ex[1];
    dEij[2][0] = d1.Ex[2];
    dEij[0][1] = d1.Ey[0];
    dEij[1][1] = d1.Ey[1];
    dEij[2][1] = d1.Ey[2];
    dEij[0][2] = d1.Ez[0];
    dEij[1][2] = d1.Ez[1];
    dEij[2][2] = d1.Ez[2];


    // make DphiDphi and EE and BB
    FOR2(i, j)
    {
        DphiDphi += d1.phi[i] * d1.phi[j] * h_UU[i][j] * vars.chi;
        EE += Ei[i] * Ei[j] * h_UU[i][j] * vars.chi;
        BB += Bi[i] * Bi[j] * h_UU[i][j] * vars.chi;
    }

    // make divE
    FOR2(i, j)
    {
        divE += h_UU[i][j] * vars.chi * dEij[i][j];
        FOR(k)
        {
            divE -= chris_phys[k][i][j] *  Ei[k] * h_UU[i][j] * vars.chi;
        }
    }



    // compute coupling and add constributions to EOM
    my_coupling.compute_coupling(f_of_phi,f_prime_of_phi,coupling_of_phi,vars);

    // adjust RHS for the coupling term
    total_rhs.phi = rhs_emd.phi;

    total_rhs.Pi = rhs_emd.Pi - f_prime_of_phi * vars.lapse *
                                          coupling_of_phi * (BB + EE);

    total_rhs.ax = rhs_emd.ax;

    total_rhs.ay = rhs_emd.ay;

    total_rhs.az = rhs_emd.az;

    total_rhs.Ex = rhs_emd.Ex;

    total_rhs.Ey = rhs_emd.Ey;

    total_rhs.Ez = rhs_emd.Ez;

    total_rhs.Xi = rhs_emd.Xi + vars.lapse * divE;

    total_rhs.At = rhs_emd.At;

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
    const Tensor<1, data_t> &d1_Ex,
    const Tensor<1, data_t> &d1_Ey,
    const Tensor<1, data_t> &d1_Ez,
    const Tensor<1, data_t> &d1_At,
    const Tensor<1, data_t> &d1_Xi,
    const EMDObject<data_t> &advec_emd)
{
    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);
    const auto chris_phys =
                 compute_phys_chris(d1.chi, vars.chi, vars.h, h_UU, chris.ULL);

    // making epsion_i^{\,\,jk}
    auto eps_symbol = epsilon(); //not the tensor, the symbol
    Tensor<3, data_t, 3> epsUUU; //tenosr
    Tensor<3, data_t, 3> epsLUU; //tensor

    // make all upstairs first, divide by root gamma
    FOR3(i,j,k)
    {
      epsUUU[i][j][k] = eps_symbol[i][j][k] / pow(vars.chi,-1.5);
      epsLUU[i][j][k] = 0.;
    }
    // lower first index with physical metric
    FOR4(i,j,k,l) epsLUU[i][j][k] += epsUUU[l][j][k] * vars.h[i][l] / vars.chi;

    // some useful variables
    Tensor<1, data_t, 3> Bi;
    Tensor<1, data_t, 3> Ei;
    Tensor<1, data_t, 3> EKj; // E^\mu K_{\mu\nu}
    Tensor<1, data_t, 3> BKj; // B^\mu K_{\mu\nu}
    Tensor<2, data_t, 3> Kij;
    Tensor<2, data_t, 3> DiEj; // D_i E_j // cov deriv
    Tensor<2, data_t, 3> DiBj; // D_i B_j // cov deriv
    data_t divB = 0.; // divergence of a
    data_t laplace_phi = 0.;

    Bi[0] = vars_emd.ax;
    Bi[1] = vars_emd.ay;
    Bi[2] = vars_emd.az;

    Ei[0] = vars_emd.Ex;
    Ei[1] = vars_emd.Ey;
    Ei[2] = vars_emd.Ez;

    DiBj[0][0] = d1_ax[0]; // partial deriv part
    DiBj[1][0] = d1_ax[1];
    DiBj[2][0] = d1_ax[2];
    DiBj[0][1] = d1_ay[0];
    DiBj[1][1] = d1_ay[1];
    DiBj[2][1] = d1_ay[2];
    DiBj[0][2] = d1_az[0];
    DiBj[1][2] = d1_az[1];
    DiBj[2][2] = d1_az[2];

    DiEj[0][0] = d1_Ex[0]; // partial deriv part
    DiEj[1][0] = d1_Ex[1];
    DiEj[2][0] = d1_Ex[2];
    DiEj[0][1] = d1_Ey[0];
    DiEj[1][1] = d1_Ey[1];
    DiEj[2][1] = d1_Ey[2];
    DiEj[0][2] = d1_Ez[0];
    DiEj[1][2] = d1_Ez[1];
    DiEj[2][2] = d1_Ez[2];

    FOR3(i,j,k) //covariant deriv part
    {
        DiEj[i][j] += - chris_phys[k][i][j] * Ei[k];
        DiBj[i][j] += - chris_phys[k][i][j] * Bi[k];
    }

    FOR2(i,j)
    {
        // spatial laplacian of phi
        laplace_phi += vars.chi * h_UU[i][j] * d2_phi[i][j];
        FOR1(k)
        {
            laplace_phi -= chris_phys[k][i][j] * d1_phi[k]
                         * h_UU[i][j] * vars.chi;
        }

        // div B
        divB += DiBj[i][j]*h_UU[i][j]*vars.chi;

    }

    // Lie derivs with respect to the shift vector
    data_t LieDeriv_Bx = advec_emd.ax;
    data_t LieDeriv_By = advec_emd.ay;
    data_t LieDeriv_Bz = advec_emd.az;
    data_t LieDeriv_Ex = advec_emd.Ex;
    data_t LieDeriv_Ey = advec_emd.Ey;
    data_t LieDeriv_Ez = advec_emd.Ez;

    // partial_j beta^i = d1.shift[i][j]
    FOR1(i)
    {
        LieDeriv_Bx += Bi[i] * d1.shift[i][0];
        LieDeriv_By += Bi[i] * d1.shift[i][1];
        LieDeriv_Bz += Bi[i] * d1.shift[i][2];
        LieDeriv_Ex += Ei[i] * d1.shift[i][0];
        LieDeriv_Ey += Ei[i] * d1.shift[i][1];
        LieDeriv_Ez += Ei[i] * d1.shift[i][2];
    }

    FOR1(i)
    {
        EKj[i] = 0.;
        BKj[i] = 0.;
    }

    FOR2(i,j)
    {
        Kij[i][j] = (vars.A[i][j] + vars.K * vars.h[i][j]/3.)/vars.chi;
    }
    FOR3(i,j,k)
    {
        EKj[i] += Ei[j] * Kij[k][i] * h_UU[j][k] * vars.chi;
        BKj[i] += Bi[j] * Kij[k][i] * h_UU[j][k] * vars.chi;
    }

    // the right hans side updates, no coupling function to dilaton here

    rhs_emd.phi = advec_emd.phi - vars.lapse * vars_emd.Pi;

    rhs_emd.Pi = advec_emd.Pi + vars.lapse * vars.K * vars_emd.Pi
                              - vars.lapse * laplace_phi;

    rhs_emd.ax = LieDeriv_Bx + vars.lapse * vars_emd.ax * vars.K
                             + vars.lapse * d1_At[0]
                             - 2. * vars.lapse * BKj[0];

    rhs_emd.ay = LieDeriv_By + vars.lapse * vars_emd.ay * vars.K
                             + vars.lapse * d1_At[1]
                             - 2. * vars.lapse * BKj[1];

    rhs_emd.az = LieDeriv_Bz + vars.lapse * vars_emd.az * vars.K
                             + vars.lapse * d1_At[2]
                             - 2. * vars.lapse * BKj[2];

    rhs_emd.Ex = LieDeriv_Ex + vars.lapse * vars_emd.Ex * vars.K
                             + vars.lapse * d1_Xi[0]
                             - 2. * vars.lapse * EKj[0];

    rhs_emd.Ey = LieDeriv_Ey + vars.lapse * vars_emd.Ey * vars.K
                             + vars.lapse * d1_Xi[1]
                             - 2. * vars.lapse * EKj[1];

    rhs_emd.Ez = LieDeriv_Ez + vars.lapse * vars_emd.Ez * vars.K
                             + vars.lapse * d1_Xi[2]
                             - 2. *vars.lapse *  EKj[2];

    data_t kappa_E = 1. ; // maxwell E damping param, should be positive

    data_t kappa_B = 1. ; // maxwell B damping param, should be positive

    rhs_emd.Xi = advec_emd.Xi - vars.lapse * kappa_E * vars_emd.Xi;

    rhs_emd.At = advec_emd.At + vars.lapse * (divB - kappa_B * vars_emd.At);


    FOR2(j,k)
    {
        rhs_emd.Pi += - vars.chi * h_UU[k][j] * d1_phi[k] * d1.lapse[j];

        rhs_emd.Ex += epsLUU[0][j][k]
                       * (Bi[k] * d1.lapse[j] + vars.lapse * DiBj[j][k]);

        rhs_emd.Ey += epsLUU[1][j][k]
                       * (Bi[k] * d1.lapse[j] + vars.lapse * DiBj[j][k]);

        rhs_emd.Ez += epsLUU[2][j][k]
                       * (Bi[k] * d1.lapse[j] + vars.lapse * DiBj[j][k]);

        rhs_emd.ax += - epsLUU[0][j][k]
                       * (Ei[k] * d1.lapse[j] + vars.lapse * DiEj[j][k]);

        rhs_emd.ay += - epsLUU[1][j][k]
                       * (Ei[k] * d1.lapse[j] + vars.lapse * DiEj[j][k]);

        rhs_emd.az += - epsLUU[2][j][k]
                       * (Ei[k] * d1.lapse[j] + vars.lapse * DiEj[j][k]);
    }


}

#endif /* EINSTEINMAXWELLDILATONFIELD_IMPL_HPP_ */
