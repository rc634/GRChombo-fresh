/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

//need this to stop double defining tetrad if we use weyl 4 aswell
#ifndef WEYL4_HPP_

//! Struct for the null tetrad
template <class data_t> struct Tetrad_t
{
    Tensor<1, data_t> u; //!< the vector u^i
    Tensor<1, data_t> v; //!< the vector v^i
    Tensor<1, data_t> w; //!< the vector w^i
};

#endif

//! Struct for the EM Radiation scalar
template <class data_t> struct EMRScalar_t
{
    data_t Real; // Real component
    data_t Im;   // Imaginary component
};




#ifndef PHEYL2_HPP_
#define PHEYL2_HPP_

#include "CCZ4.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs c_NUM - total number of components
#include "simd.hpp"
#include <array>




//!  Calculates the Pheyl2 scalar for spacetimes with electromagnetism
/*!
   This class calculates the Phi2 (called Pheyl2) scalar real and im parts using
   definitions from "https://arxiv.org/pdf/1205.1063.pdf"

   Tetrads defined in "gr-qc/0104063" to calculate Pheyl2.
*/
class Pheyl2
{
  public:
    // Use the variable definitions containing the needed quantities
    template <class data_t> using Vars = CCZ4Vars::VarsWithGauge<data_t>;
    template <class data_t> using EMDVars =
                               EinsteinMaxwellDilatonField<>::Vars<data_t>;

    //! Constructor of class Pheyl2
    /*!
        Takes in the centre for the calculation of the tetrads and grid spacing.
    */
    Pheyl2(const std::array<double, CH_SPACEDIM> a_center, const double a_dx)
        : m_center(a_center), m_dx(a_dx)
    {
    }

    //! The compute member which calculates the wave quantities at each point on
    //! the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    const std::array<double, CH_SPACEDIM> m_center; //!< The grid center
    const double m_dx;                              //!< the grid spacing

    //! Compute spatial volume element
    template <class data_t>
    Tensor<3, data_t> compute_epsilon3_LLU(const Vars<data_t> &vars,
                                           const Tensor<2, data_t> &h_UU) const;

    //! Calculation of Pheyl_2 scalar
    template <class data_t>
    EMRScalar_t<data_t> compute_Pheyl2(const Vars<data_t> &vars,
                                     const EMDVars<data_t> &emdvars,
                                     const Tensor<3, data_t> &epsilon3_LUU,
                                     const Tensor<2, data_t> &h_UU,
                                     const Coordinates<data_t> &coords) const;

    //! Calculation of the tetrads
    template <class data_t>
    Tetrad_t<data_t>
    compute_null_tetrad(const Vars<data_t> &vars,
                        const Tensor<2, data_t> &h_UU,
                        const Coordinates<data_t> &coords) const;


};

#include "Pheyl2.impl.hpp"

#endif /* PHEYL2_HPP_ */
