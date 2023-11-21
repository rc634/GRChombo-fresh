/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef EMDBH_READ_HPP_
#define EMDBH_READ_HPP_


#include "EMDBHSolution_read.hpp"
#include "Cell.hpp"
#include "EMDCouplingFunction.hpp"
#include "EinsteinMaxwellDilatonField.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "parstream.H" //gives pout
#include "simd.hpp"
#include <vector>

//! Class which solves for the initial data for a spherically symmetric emdbh

class EMDBH_read
{

  public:
    //! The constructor
    EMDBH_read(EMDBH_params_t a_params_EMDBH,
              CouplingFunction::params_t a_params_coupling_function,
              double a_G_Newton, double a_dx, int a_verbosity);

    //! Computes the 1d solution and stores in m_1d_sol
    void compute_1d_solution(const double max_r);

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    EMDBHSolution_read m_1d_sol; /*<
    The object that stores the solution found by reading data file */

  protected:
    double m_dx;
    double m_G_Newton;
    EMDBH_params_t m_params_EMDBH;  //!< The complex scalar field params
    CouplingFunction::params_t m_params_coupling_function; //!< The potential params
    int m_verbosity;
};

#include "EMDBH_read.impl.hpp"

#endif /* EMDBH_READ_HPP_ */
