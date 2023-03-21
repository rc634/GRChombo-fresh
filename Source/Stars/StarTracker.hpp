/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef STARTRACKER_HPP_
#define STARTRACKER_HPP_

#include "AMRInterpolator.hpp"
#include "AlwaysInline.hpp"
#include "Lagrange.hpp"

//!  The class tracks the puncture locations by integrating the shift at
//!  The puncture position
class StarTracker
{
  private:
    //! Params for puncture tracking
    int m_num_stars;
    std::vector<double> m_star_coords;
    std::array<double, CH_SPACEDIM> m_centre;
    int m_tracking_level; // level (i.e. times) to execute tracking
    int m_resolution;     // number of points n, (2n + 1 points to integrate)
    double m_width;
    std::string m_punctures_filename;

    // saved pointer to external interpolator
    AMRInterpolator<Lagrange<4>> *m_interpolator;

  public:
    //! The constructor
    StarTracker() : m_interpolator(nullptr) {}

    //! set puncture locations on start (or restart)
    //! this needs to be done before 'setupAMRObject'
    //! if the puncture locations are required for Tagging Criteria
    void initial_setup(std::array<double, CH_SPACEDIM> a_star_track_centre,
                       bool a_do_star_track, int a_number_of_stars,
                       std::vector<double> a_initial_star_centres,
                       int a_star_track_resolution, double a_star_track_width)
    {
        m_num_stars = a_number_of_stars;
        m_star_coords = a_initial_star_centres;
        m_resolution = a_star_track_resolution;
        m_width = a_star_track_width;
        m_centre = a_star_track_centre;

        for (int n = 0; n < m_num_stars; n++)
        {
            for (int i = 0; i < CH_SPACEDIM; i++)
            {
                m_star_coords[n * CH_SPACEDIM + i] += m_centre[n];
            }
        }
    }

    void test();

    ALWAYS_INLINE void
    set_interpolator(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        m_interpolator = a_interpolator;
    }

    void update_star_centres(int a_field_index);

    void get_star_centres(std::vector<double> &a_centre);

    void write_to_dat(std::string a_filename, double a_dt, double a_time,
                      double a_restart_time, bool a_first_step);

    void read_old_centre_from_dat(std::string a_filename, double a_dt,
                                  double a_time, double a_restart_time,
                                  bool a_first_step);

    void
    get_field_value_at_centres(int a_field_index,
                               std::vector<double> &a_out_data,
                               AMRInterpolator<Lagrange<4>> *a_interpolator);

  private:
};

#endif /* STARTRACKER_HPP_ */
