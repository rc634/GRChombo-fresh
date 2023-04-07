/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EINSTEINMAXWELLDILATONFIELD_HPP_
#define EINSTEINMAXWELLDILATONFIELD_HPP_

#include "CCZ4Geometry.hpp"
#include "DefaultEMDCouplingFunction.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

//!  Calculates the matter type specific elements such as the EMTensor and
//   matter evolution

template <class coupling_t = DefaultEMDCouplingFunction>
                                              class EinsteinMaxwellDilatonField
{
  protected:
    //! The local copy of the coupling
    coupling_t my_coupling;

  public:
    //!  Constructor of class EinsteinMaxwellDilaton field, inputs are the matter parameters.
    EinsteinMaxwellDilatonField(const coupling_t a_coupling)
        : my_coupling(a_coupling)
    {
    }

    //! Structure containing the variables for the matter fields
    template <class data_t> struct EMDObject
    {
        data_t phi;
        data_t Pi;
        data_t At;
        data_t ax;
        data_t ay;
        data_t az;
        data_t Ex;
        data_t Ey;
        data_t Ez;
        data_t Xi;
    };

    //! Structure containing the rhs variables for the matter fields
    template <class data_t> struct Vars
    {
      data_t phi;
      data_t Pi;
      data_t At;
      data_t ax;
      data_t ay;
      data_t az;
      data_t Ex;
      data_t Ey;
      data_t Ez;
      data_t Xi;

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi, phi);
            VarsTools::define_enum_mapping(mapping_function, c_Pi, Pi);
            VarsTools::define_enum_mapping(mapping_function, c_At, At);
            VarsTools::define_enum_mapping(mapping_function, c_ax, ax);
            VarsTools::define_enum_mapping(mapping_function, c_ay, ay);
            VarsTools::define_enum_mapping(mapping_function, c_az, az);
            VarsTools::define_enum_mapping(mapping_function, c_Ex, Ex);
            VarsTools::define_enum_mapping(mapping_function, c_Ey, Ey);
            VarsTools::define_enum_mapping(mapping_function, c_Ez, Ez);
            VarsTools::define_enum_mapping(mapping_function, c_Xi, Xi);
        }
    };

    //! Structure containing the rhs variables for the matter fields requiring
    //!  2nd derivs
    template <class data_t> struct Diff2Vars
    {
        data_t phi;
        data_t ax;
        data_t ay;
        data_t az;

        /// Defines the mapping between members of Vars and Chombo grid
        ///  variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi, phi);
            VarsTools::define_enum_mapping(mapping_function, c_ax, ax);
            VarsTools::define_enum_mapping(mapping_function, c_ay, ay);
            VarsTools::define_enum_mapping(mapping_function, c_az, az);
        }
    };

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, including the coupling
    template <class data_t, template <typename> class vars_t>
    emtensor_t<data_t> compute_emtensor(
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices)
        const Tensor<3, data_t> &chris_ULL)
        const; //!< the conformal christoffel symbol

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, excluding the coupling
    template <class data_t, template <typename> class vars_t>
    static void emtensor_excl_coupling(
        emtensor_t<data_t> &out,           //!< the em tensor output
        const vars_t<data_t> &vars,        //!< the value of the variables
        const EMDObject<data_t> &vars_emd, //!< the value of the emd variables
        const Tensor<1, data_t>
            &d1_phi, //!< the value of the first deriv of phi
        const Tensor<2, data_t> &h_UU, //!< the inverse metric (raised indices).
        const Tensor<3, data_t>
            &chris_ULL); //!< the conformal christoffel symbol

    //! The function which adds in the RHS for the matter field vars,
    //! including the coupling
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void add_matter_rhs(
        rhs_vars_t<data_t> &total_rhs,       //!< value of the RHS for all vars
        const vars_t<data_t> &vars,          //!< value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< value of the 2nd derivs
        const vars_t<data_t> &advec) //!< the value of the advection terms
        const;

    //! The function which calculates the RHS for the matter field vars
    //! excluding the coupling
    template <class data_t, template <typename> class vars_t>
    static void matter_rhs_excl_coupling(
        EMDObject<data_t>
            &rhs_emd, //!< the value of the RHS terms for the sf vars
        const vars_t<data_t> &vars,        //!< the values of all the variables
        const EMDObject<data_t> &vars_emd, //!< the value of the sf variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const Tensor<1, data_t>
            &d1_phi, //!< the value of the 1st derivs of phi
        const Tensor<1, data_t>
            &d1_ax,
        const Tensor<1, data_t>
            &d1_ay,
        const Tensor<1, data_t>
            &d1_az,
        const Tensor<2, data_t>
            &d2_phi, //!< the value of the 2nd derivs of phi
        const Tensor<2, data_t>
            &d2_ax,
        const Tensor<2, data_t>
            &d2_ay,
        const Tensor<2, data_t>
            &d2_az,
        const Tensor<1, data_t>
            &d1_At,
        const Tensor<1, data_t>
            &d1_Xi,
        const EMDObject<data_t>
            &advec_emd); //!< advection terms for the emd vars
};

#include "EinsteinMaxwellDilatonField.impl.hpp"

#endif /* EINSTEINMAXWELLDILATONFIELD_HPP_ */
