/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include <array>

#ifndef EMDBHPARAMS_HPP_
#define EMDBHPARAMS_HPP_

//! A structure for the input params for the emdbh
struct EMDBH_params_t
{
    int gridpoints;
    double separation;
    double bh_charge;
    double bh_mass;
    bool binary;
    double Newtons_constant;
    double Ylm_amplitude;
    double Ylm_thickness;
    double Ylm_r0;
    std::string data_path;
    std::array<double, CH_SPACEDIM>
        star_centre; //!< coordinates of the centre of the emdbh
};

#endif /* EMDBHPARAMS_HPP_ */
