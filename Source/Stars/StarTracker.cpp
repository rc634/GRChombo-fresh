/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "StarTracker.hpp"
#include "ChomboParameters.hpp" // for writing data
#include "DimensionDefinitions.hpp"
#include "InterpolationQuery.hpp"
#include "SmallDataIO.hpp"   // for writing data
#include "UserVariables.hpp" // for writing data

//! Set punctures post restart
void StarTracker::test()
{
    std::cout << "test shout for STAMR ! " << std::endl;
}

void StarTracker::update_star_centres(int a_field_index)
{
    int i_tot = m_num_stars * m_resolution * 3;
    int i_max = m_num_stars * 3;
    int n1, n2, n3;
    double delta;
    std::vector<double> x_coords(i_tot);
    std::vector<double> y_coords(i_tot);
    std::vector<double> z_coords(i_tot);
    std::vector<double> vals(i_tot);

    // setup positions in a 3d cross about old centre
    for (int n = 0; n < m_num_stars; n++)
    {
        n1 = 3 * n;     // index of x gaussian
        n2 = 3 * n + 1; // index of y gaussian
        n3 = 3 * n + 2; // index of z gaussian

        for (int i = 0; i < m_resolution; i++)
        {

            delta = m_width * (double(i) / double(m_resolution - 1) - 0.5);

            x_coords[n1 * m_resolution + i] = m_star_coords[3 * n] + delta;
            y_coords[n1 * m_resolution + i] = m_star_coords[3 * n + 1];
            z_coords[n1 * m_resolution + i] = m_star_coords[3 * n + 2];

            x_coords[n2 * m_resolution + i] = m_star_coords[3 * n];
            y_coords[n2 * m_resolution + i] = m_star_coords[3 * n + 1] + delta;
            z_coords[n2 * m_resolution + i] = m_star_coords[3 * n + 2];

            x_coords[n3 * m_resolution + i] = m_star_coords[3 * n];
            y_coords[n3 * m_resolution + i] = m_star_coords[3 * n + 1];
            z_coords[n3 * m_resolution + i] = m_star_coords[3 * n + 2] + delta;
        }
    }

    // interpolate field values on 3d cross
    m_interpolator->refresh();
    InterpolationQuery query(i_tot);
    query.setCoords(0, x_coords.data());
    query.setCoords(1, y_coords.data());
    query.setCoords(2, z_coords.data());
    query.addComp(a_field_index, vals.data(), Derivative::LOCAL,
                  VariableType::diagnostic);
    m_interpolator->interp(query);

    // calculate expectations of
    double x_int, x_weighted_int, y_int, y_weighted_int, z_int, z_weighted_int;
    for (int n = 0; n < m_num_stars; n++)
    {
        x_int = 0;
        x_weighted_int = 0;
        y_int = 0;
        y_weighted_int = 0;
        z_int = 0;
        z_weighted_int = 0;
        n1 = 3 * n;     // index of x gaussian
        n2 = 3 * n + 1; // index of y gaussian
        n3 = 3 * n + 2; // index of z gaussian

        for (int i = 0; i < m_resolution; i++)
        {
            x_int += vals[n1 * m_resolution + i];
            y_int += vals[n2 * m_resolution + i];
            z_int += vals[n3 * m_resolution + i];
            x_weighted_int +=
                x_coords[n1 * m_resolution + i] * vals[n1 * m_resolution + i];
            y_weighted_int +=
                y_coords[n2 * m_resolution + i] * vals[n2 * m_resolution + i];
            z_weighted_int +=
                z_coords[n3 * m_resolution + i] * vals[n3 * m_resolution + i];
        }

        m_star_coords[n1] = x_weighted_int / x_int;
        m_star_coords[n2] = y_weighted_int / y_int;
        m_star_coords[n3] = z_weighted_int / z_int;
    }
}

void StarTracker::write_to_dat(std::string a_filename, double a_dt,
                               double a_time, double a_restart_time,
                               bool a_first_step)
{
    double eps = 10e-8;
    SmallDataIO star_centre_file(a_filename, a_dt, a_time, a_restart_time,
                                 SmallDataIO::APPEND, a_first_step);

    if (a_time > a_restart_time + eps)
        star_centre_file.remove_duplicate_time_data();

    std::vector<string> header_line(3. * m_num_stars);
    int n1, n2, n3;
    for (int n = 0; n < m_num_stars; n++)
    {
        n1 = 3 * n;
        n2 = 3 * n + 1;
        n3 = 3 * n + 2;
        header_line[n1] = "Star " + to_string(n) + " x";
        header_line[n2] = "Star " + to_string(n) + " y";
        header_line[n3] = "Star " + to_string(n) + " z";
    }

    if (a_time == 0.)
    {
        star_centre_file.write_header_line(header_line);
    }
    star_centre_file.write_time_data_line(m_star_coords);
}

void StarTracker::get_star_centres(std::vector<double> &a_centre)
{
    int i_max = m_num_stars * 3;
    a_centre.resize(i_max);
    for (int i = 0; i < i_max; i++)
    {
        a_centre[i] = m_star_coords[i];
    }
}

// read a data line from the previous timestep
void StarTracker::read_old_centre_from_dat(std::string a_filename, double a_dt,
                                           double a_time, double a_restart_time,
                                           bool a_first_step)
{
    std::vector<double> data_line;
    if (a_time > a_dt / 3.)
    {
        SmallDataIO star_centre_file(a_filename, a_dt, a_time, a_restart_time,
                                     SmallDataIO::READ, a_first_step);
        star_centre_file.get_specific_data_line(data_line, a_time - a_dt);

        bool length_match = data_line.size() == m_num_stars * CH_SPACEDIM;

        if (length_match)
        {
            for (int i = 0; i < data_line.size(); i++)
            {
                m_star_coords[i] = data_line[i];
            }
        }
        else
        {
            for (int i = 0; i < m_star_coords.size(); i++)
            {
                m_star_coords[i] = NAN;
            }
            std::cout
                << "Array Size Mismatch While Loading Star Centre From File ! "
                << std::endl;
        }
    }
}

// pass this an empty std::vector and it will resize and fill with centre values
// of field
void StarTracker::get_field_value_at_centres(
    int a_field_index, std::vector<double> &a_out_data,
    AMRInterpolator<Lagrange<4>> *a_interpolator)
{
    bool gone_NAN = false;
    a_out_data.resize(m_num_stars);

    // detect if gone stars centres have gone nan
    for (int i = 0; i < m_num_stars * CH_SPACEDIM; i++)
    {
        if (std::isnan(m_star_coords[i]))
        {
            gone_NAN = true;
        }
    }

    if (gone_NAN == false)
    {

        std::vector<double> x(m_num_stars);
        std::vector<double> y(m_num_stars);
        std::vector<double> z(m_num_stars);
        std::vector<double> f(m_num_stars);

        for (int n = 0; n < m_num_stars; n++)
        {
            x[n] = m_star_coords[n * CH_SPACEDIM];
            y[n] = m_star_coords[n * CH_SPACEDIM + 1];
            z[n] = m_star_coords[n * CH_SPACEDIM + 2];
        }

        a_interpolator->refresh();
        InterpolationQuery query(m_num_stars);
        query.setCoords(0, x.data());
        query.setCoords(1, y.data());
        query.setCoords(2, z.data());
        query.addComp(a_field_index, f.data());
        a_interpolator->interp(query);

        for (int n = 0; n < m_num_stars; n++)
        {
            a_out_data[n] = f[n];
        }
    }
    else
    {
        for (int n = 0; n < m_num_stars; n++)
        {
            a_out_data[n] = NAN;
        }
    }
}
