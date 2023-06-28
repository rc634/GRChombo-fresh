/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMDLORENTZSCALARS_HPP_
#define EMDLORENTZSCALARS_HPP_

#include "ADMConformalVars.hpp" // needed for CCz4 and matter variables
#include "Cell.hpp"
#include "EinsteinMaxwellDilatonField.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp"
#include "simd.hpp"

//! Calculates the Noether Charge integrand values and the modulus of the
//! complex scalar field on the grid
class EMDLorentzScalars
{
    // Need matter variables and chi
    template <class data_t>
    using ADMVars = ADMConformalVars::VarsWithGauge<data_t>;
    template <class data_t>
    using MatterVars = EinsteinMaxwellDilatonField<>::Vars<data_t>;

    double m_dx;

  public:
    EMDLorentzScalars(double a_dx) : m_dx(a_dx){}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // load vars locally
        const FourthOrderDerivatives m_deriv(m_dx);
        const auto adm_vars = current_cell.template load_vars<ADMVars>();
        const auto matter_vars = current_cell.template load_vars<MatterVars>();
        const auto adm_d1 = m_deriv.template diff1<ADMVars>(current_cell);
        const auto matter_d1 = m_deriv.template diff1<MatterVars>(current_cell);
        const auto matter_d2 = m_deriv.template diff2<MatterVars>(current_cell);

        // root minus g
        data_t root_minus_g = pow(adm_vars.chi, -1.5)*adm_vars.lapse;

        // inverse conformal metric
        const auto h_UU = TensorAlgebra::compute_inverse_sym(adm_vars.h);
        auto gamma_UU = h_UU;
        FOR2(i,j) gamma_UU[i][j] = h_UU[i][j]*adm_vars.chi;
        const auto chris = TensorAlgebra::compute_christoffel(adm_d1.h, h_UU);
        const auto chris_phys =
                           TensorAlgebra::compute_phys_chris(adm_d1.chi, adm_vars.chi,
                                                       adm_vars.h, h_UU, chris.ULL);


        ////////////////////////////////////////////
        // calculate A_\mu A^\mu = -At^2 + a_i a^i
        // (At = A^\mu n_\mu = -\alpha A^t)
        ////////////////////////////////////////////
        data_t H2norm_damping = matter_vars.At * matter_vars.At
                              + matter_vars.Xi * matter_vars.Xi;


        ////////////////////////////////////////////
        // calculate F_{\mu\nu} F^{\mu\nu} = 2 B_i B^i - 2 E_i E^i
        ////////////////////////////////////////////

        data_t FF = 0.;
        Tensor<1, data_t, 3> Ei;
        Tensor<1, data_t, 3> Bi;
        data_t BB = 0.; // squared variables
        data_t EE = 0.; // squared variables

        Ei[0] = matter_vars.Ex;
        Ei[1] = matter_vars.Ey;
        Ei[2] = matter_vars.Ez;

        Bi[0] = matter_vars.ax;
        Bi[1] = matter_vars.ay;
        Bi[2] = matter_vars.az;


        FOR2(i, j)
        {
            EE += Ei[i] * Ei[j] * h_UU[i][j] * adm_vars.chi;
            BB += Bi[i] * Bi[j] * h_UU[i][j] * adm_vars.chi;
        }

        FF = 2. * BB - 2. * EE;


        ////////////////////////////////////////////
        // calculate hamiltonian of scalar field
        ////////////////////////////////////////////

        data_t phi_hamiltonian = matter_vars.Pi * matter_vars.Pi;

        FOR1(i)
        {
            phi_hamiltonian += 2. * matter_vars.Pi
                                  * adm_vars.shift[i] * matter_d1.phi[i];
        }

        FOR2(i,j)
        {
            phi_hamiltonian += matter_d1.phi[i] * matter_d1.phi[j]
                             * h_UU[i][j] * adm_vars.chi;
        }

        phi_hamiltonian *= root_minus_g;





        ////////////////////////////////////////////
        // maxwell constriants
        ////////////////////////////////////////////

        data_t H2norm_maxwell_constraints = 0.;
        Tensor<2, data_t, 3> DiBj;
        Tensor<2, data_t, 3> DiEj;
        data_t magnetic_constraint=0.;
        data_t electric_constraint=0.;

        FOR1(i)
        {
            DiBj[i][0] = matter_d1.ax[i];
            DiBj[i][1] = matter_d1.ay[i];
            DiBj[i][2] = matter_d1.az[i];

            DiEj[i][0] = matter_d1.Ex[i];
            DiEj[i][1] = matter_d1.Ey[i];
            DiEj[i][2] = matter_d1.Ez[i];
        }

        // partial derivs
        FOR2(i,j)
        {
            electric_constraint += DiEj[i][j] * h_UU[i][j] * adm_vars.chi;
            magnetic_constraint += DiBj[i][j] * h_UU[i][j] * adm_vars.chi;
        }

        // covariant corrections
        FOR3(i,j,k)
        {
            magnetic_constraint += - chris_phys[k][i][j] * Bi[k]
                                     * h_UU[i][j] * adm_vars.chi;
            electric_constraint += - chris_phys[k][i][j] * Ei[k]
                                     * h_UU[i][j] * adm_vars.chi;
        }

        H2norm_maxwell_constraints = electric_constraint * electric_constraint
                                   + magnetic_constraint * magnetic_constraint;



        ////////////////////////////////////////////
        // maxwell constriants
        ////////////////////////////////////////////


        // store variables
        //current_cell.store_vars(phi_hamiltonian, c_phi_ham);
        current_cell.store_vars(sqrt(H2norm_damping), c_phi_ham);
        current_cell.store_vars(FF, c_mod_F);
        current_cell.store_vars(sqrt(H2norm_maxwell_constraints), c_mod_A);
    }
};

#endif /* EMDLORENTZSCALARS_HPP_ */
