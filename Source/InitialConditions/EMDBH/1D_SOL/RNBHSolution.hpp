
#ifndef RNBHSOLUTION_HPP_
#define RNBHSOLUTION_HPP_

// a lightweight and basic polar areal reissner nordstrom black hole
class RNBHSolution
{

  private:

    double G, L, dx, M=1., Qratio, Q; // Qratio = Q/Qextremal
    int gridsize;

    std::vector<double> At;            // electric potential
    std::vector<double> omega;          // lapse
    std::vector<double> psi;            // conformal factor

  public:

    RNBHSolution();
    void set_initialcondition_params(EMDBH_params_t m_params_EMDBH,
                        CouplingFunction::params_t m_params_coupling_function,
                                                          const double max_r);
    double get_At_interp(const double r) const;
    double get_lapse_interp(const double r) const;
    double get_dlapse_interp(const double r) const;
    double get_psi_interp(const double r) const;
    double get_dpsi_interp(const double r) const;

    void main();
    void main_polar_areal();
};

#include "RNBHSolution.impl.hpp"

#endif /* RNBHSOLUTION_HPP_ */
