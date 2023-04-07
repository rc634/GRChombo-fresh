/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef USERVARIABLES_HPP
#define USERVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4UserVariables.hpp"
#include "DiagnosticVariables.hpp"
#include <array>
#include <string>

/// This enum gives the index of every variable stored in the grid
enum
{

    c_phi = NUM_CCZ4_VARS, // scalar field
    c_Pi,                    // scalar field momentum
    c_At,                 // time component of A
    c_ax,                  // xcomponent of A
    c_ay,                  // ycomponent of A
    c_az,                  // zcomponent of A
    c_Ex,                  // x electric field
    c_Ey,                  // y electric field
    c_Ez,                  // z electric field
    c_Xi,                  // maxwell constriant


    // Note that it is important that the first enum value is set to 1 more than
    // the last CCZ4 var enum
    NUM_VARS
};

namespace UserVariables
{
static const std::array<std::string, NUM_VARS - NUM_CCZ4_VARS>
    matterfield_variable_names =
    {"phi", "Pi", "At",
               "ax", "ay", "az", "Ex", "Ey", "Ez", "Xi"};

static const std::array<std::string, NUM_VARS> variable_names =
    ArrayTools::concatenate(ccz4_variable_names, matterfield_variable_names);
} // namespace UserVariables

#include "UserVariables.inc.hpp"

#endif /* USERVARIABLES_HPP */
