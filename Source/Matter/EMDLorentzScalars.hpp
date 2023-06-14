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
        data_t AA = -matter_vars.At*matter_vars.At;

        Tensor<1, data_t, 3> a;

        a[0] = matter_vars.ax;
        a[1] = matter_vars.ay;
        a[2] = matter_vars.az;

        FOR2(i,j)
        {
            AA += a[i] * a[j] * h_UU[i][j] * adm_vars.chi;
        }

        ////////////////////////////////////////////
        // calculate F_{\mu\nu} F^{\mu\nu} = B_{ij}B^{ij} - 2 E_i E^i
        ////////////////////////////////////////////
        data_t FF = 0.;
        Tensor<1, data_t, 3> Ei;
        Tensor<2, data_t, 3> Bij;
        data_t BB = 0.; // squared variables
        data_t EE = 0.; // squared variables

        Ei[0] = matter_vars.Ex;
        Ei[1] = matter_vars.Ey;
        Ei[2] = matter_vars.Ez;

        Bij[0][0] = 0.;
        Bij[0][1] = matter_d1.ay[0] - matter_d1.ax[1];
        Bij[0][2] = matter_d1.az[0] - matter_d1.ax[2];
        Bij[1][0] = -Bij[0][1];
        Bij[1][1] = 0.;
        Bij[1][2] = matter_d1.az[1] - matter_d1.ay[2];
        Bij[2][0] = -Bij[0][2];
        Bij[2][1] = -Bij[1][2];
        Bij[2][2] = 0.;

        FOR2(i, j)
        {
            EE += Ei[i] * Ei[j] * h_UU[i][j] * adm_vars.chi;
        }

        FOR4(i,j,k,l)
        {
            BB += Bij[i][j] * Bij[k][l]
                * h_UU[i][k] * h_UU[j][l]
                * adm_vars.chi * adm_vars.chi;
        }

        FF = BB - 2. * EE;


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

        auto epsilon = TensorAlgebra::epsilon();
        Tensor<3, data_t, 3> Dijak;
        Tensor<2, data_t, 3> Diaj;
        Tensor<2, data_t, 3> DiEj;
        data_t magnetic_constraint=0.;
        data_t electric_constraint=0.;

        FOR1(i)
        {
            Diaj[i][0] = matter_d1.ax[i];
            Diaj[i][1] = matter_d1.ay[i];
            Diaj[i][2] = matter_d1.az[i];

            DiEj[i][0] = matter_d1.Ex[i];
            DiEj[i][1] = matter_d1.Ey[i];
            DiEj[i][2] = matter_d1.Ez[i];
        }
        FOR2(i,j)
        {
            Dijak[i][j][0] = matter_d2.ax[i][j];
            Dijak[i][j][1] = matter_d2.ay[i][j];
            Dijak[i][j][2] = matter_d2.az[i][j];

            electric_constraint += DiEj[i][j] * h_UU[i][j] * adm_vars.chi;
        }

        FOR4(i,j,k,l)
        {
            Dijak[i][j][k] += - chris_phys[l][i][j] * Diaj[l][k]
                              - chris_phys[l][i][k] * Diaj[j][l];
        }

        FOR3(i,j,k)
        {
            magnetic_constraint += epsilon[i][j][k]*Dijak[i][j][k];
            electric_constraint += - chris_phys[k][i][j] * Ei[k]
                                     * h_UU[i][j] * adm_vars.chi;
        }
        current_cell.store_vars(electric_constraint, c_phi_ham);


        // store variables
        //current_cell.store_vars(phi_hamiltonian, c_phi_ham);
        current_cell.store_vars(FF, c_mod_F);
        current_cell.store_vars(AA, c_mod_A);
    }
};

#endif /* EMDLORENTZSCALARS_HPP_ */
