/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

 This file is part of ASPECT.

 ASPECT is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 ASPECT is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ASPECT; see the file LICENSE.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/particle/interpolator/quadratic_least_squares.h>
#include <aspect/postprocess/particles.h>
#include <aspect/simulator.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/qr.h>

#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      template <int dim>
      std::vector<std::vector<double>>
                                    QuadraticLeastSquares<dim>::properties_at_points(const ParticleHandler<dim> &particle_handler,
                                                                                     const std::vector<Point<dim>> &positions,
                                                                                     const ComponentMask &selected_properties,
                                                                                     const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        const unsigned int n_particle_properties = particle_handler.n_properties_per_particle();

        const unsigned int property_index = selected_properties.first_selected_component(selected_properties.size());

        AssertThrow(property_index != numbers::invalid_unsigned_int,
                    ExcMessage("Internal error: the particle property interpolator was "
                               "called without a specified component to interpolate."));

        const Point<dim> approximated_cell_midpoint = std::accumulate (positions.begin(), positions.end(), Point<dim>())
                                                      / static_cast<double> (positions.size());

        typename parallel::distributed::Triangulation<dim>::active_cell_iterator found_cell;

        if (cell == typename parallel::distributed::Triangulation<dim>::active_cell_iterator())
          {
            // We can not simply use one of the points as input for find_active_cell_around_point
            // because for vertices of mesh cells we might end up getting ghost_cells as return value
            // instead of the local active cell. So make sure we are well in the inside of a cell.
            Assert(positions.size() > 0,
                   ExcMessage("The particle property interpolator was not given any "
                              "positions to evaluate the particle properties at."));


            found_cell =
              (GridTools::find_active_cell_around_point<> (this->get_mapping(),
                                                           this->get_triangulation(),
                                                           approximated_cell_midpoint)).first;
          }
        else
          found_cell = cell;

        const typename ParticleHandler<dim>::particle_iterator_range particle_range =
          particle_handler.particles_in_cell(found_cell);


        std::vector<std::vector<double>> cell_properties(positions.size(),
                                                         std::vector<double>(n_particle_properties,
                                                                             numbers::signaling_nan<double>()));

        const unsigned int n_particles = std::distance(particle_range.begin(), particle_range.end());

        const unsigned int n_matrix_columns = (dim == 2) ? 6 : 10;
        AssertThrow(n_particles >= n_matrix_columns,
                    ExcMessage("At least one cell didn't contain enough particles to interpolate a unique solution for the cell. "
                               "The 'quadratic' interpolation scheme does not support this case."));


        // Notice that the size of matrix A is n_particles x n_matrix_columns
        // which usually is not a square matrix. Therefore, we find the
        // least squares solution of Ac=r by solving the reduced QR factorization
        // Ac = QRc = b -> Q^TQRc = Rc =Q^Tb
        // A is a std::vector of Vectors(which are it's columns) so that we create what the ImplicitQR
        // class needs.
        std::vector<Vector<double>> A(n_matrix_columns, Vector<double>(n_particles));
        std::vector<Vector<double>> b(n_particle_properties, Vector<double>(n_particles));

        unsigned int particle_index = 0;
        const double cell_diameter = found_cell->diameter();
        for (typename ParticleHandler<dim>::particle_iterator particle = particle_range.begin();
             particle != particle_range.end(); ++particle, ++particle_index)
          {
            const auto &particle_property_value = particle->get_properties();
            for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
              if (selected_properties[property_index])
                b[property_index][particle_index] = particle_property_value[property_index];
            const Tensor<1, dim, double> relative_particle_position = (particle->get_location() - approximated_cell_midpoint) / cell_diameter;
            // A is accessed by A[column][row] here since we need to append
            // columns into the qr matrix.

            // There is a potential that in the future we may change to
            // interpolate to $Q_2$ instead of the current $P_2$. This would
            // involve adding more terms to the interpolation for some problems.
            // We for now leave those terms out of the interpolation because they
            // would require more particles per cell and we are not yet aware of
            // a benchmark where those terms have affected the convegence.
            A[0][particle_index] = 1;
            A[1][particle_index] = relative_particle_position[0];
            A[2][particle_index] = relative_particle_position[1];
            if (dim == 2)
              {
                A[3][particle_index] = relative_particle_position[0] * relative_particle_position[1];
                A[4][particle_index] = relative_particle_position[0] * relative_particle_position[0];
                A[5][particle_index] = relative_particle_position[1] * relative_particle_position[1];
              }
            else
              {
                A[3][particle_index] = relative_particle_position[2];
                A[4][particle_index] = relative_particle_position[0] * relative_particle_position[1];
                A[5][particle_index] = relative_particle_position[0] * relative_particle_position[2];
                A[6][particle_index] = relative_particle_position[1] * relative_particle_position[2];
                A[7][particle_index] = relative_particle_position[0] * relative_particle_position[0];
                A[8][particle_index] = relative_particle_position[1] * relative_particle_position[1];
                A[9][particle_index] = relative_particle_position[2] * relative_particle_position[2];
              }
          }

        ImplicitQR<Vector<double>> qr;
        for (const auto &column : A)
          qr.append_column(column);
        // If A is rank deficent then qr.append_column will not append
        // the first column that can be written as a linear combination of
        // other columns. We check that all columns were added through
        // this one assertion
        AssertThrow(qr.size() == n_matrix_columns,
                    ExcMessage("The matrix A was rank deficent during quadratic least squares interpolation."));

        std::vector<Vector<double>> QTb(n_particle_properties, Vector<double>(n_matrix_columns));
        std::vector<Vector<double>> c(n_particle_properties, Vector<double>(n_matrix_columns));
        for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
          {
            if (selected_properties[property_index])
              {
                qr.multiply_with_QT(QTb[property_index], b[property_index]);
                qr.solve(c[property_index], QTb[property_index]);
              }
          }
        unsigned int index_positions = 0;
        for (typename std::vector<Point<dim>>::const_iterator itr = positions.begin(); itr != positions.end(); ++itr, ++index_positions)
          {
            const Tensor<1, dim, double> relative_support_point_location = (*itr - approximated_cell_midpoint) / cell_diameter;
            for (unsigned int property_index = 0; property_index < n_particle_properties; ++property_index)
              {
                double interpolated_value = c[property_index][0] +
                                            c[property_index][1] * relative_support_point_location[0] +
                                            c[property_index][2] * relative_support_point_location[1];
                if (dim == 2)
                  {
                    interpolated_value += c[property_index][3] * relative_support_point_location[0] * relative_support_point_location[1] +
                                          c[property_index][4] * relative_support_point_location[0] * relative_support_point_location[0] +
                                          c[property_index][5] * relative_support_point_location[1] * relative_support_point_location[1];
                  }
                else
                  {
                    interpolated_value += c[property_index][3] * relative_support_point_location[2] +
                                          c[property_index][4] * relative_support_point_location[0] * relative_support_point_location[1] +
                                          c[property_index][5] * relative_support_point_location[0] * relative_support_point_location[2] +
                                          c[property_index][6] * relative_support_point_location[1] * relative_support_point_location[2] +
                                          c[property_index][7] * relative_support_point_location[0] * relative_support_point_location[0] +
                                          c[property_index][8] * relative_support_point_location[1] * relative_support_point_location[1] +
                                          c[property_index][9] * relative_support_point_location[2] * relative_support_point_location[2];
                  }

                // Overshoot and undershoot correction of interpolated particle property.
                if (use_global_values_limiter)
                  {
                    interpolated_value = std::min(interpolated_value, global_maximum_particle_properties[property_index]);
                    interpolated_value = std::max(interpolated_value, global_minimum_particle_properties[property_index]);
                  }

                cell_properties[index_positions][property_index] = interpolated_value;
              }
          }
        return cell_properties;
      }



      template <int dim>
      void
      QuadraticLeastSquares<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Interpolator");
            {
              prm.enter_subsection("Quadratic least squares");
              {
                prm.declare_entry ("Global particle property maximum",
                                   boost::lexical_cast<std::string>(std::numeric_limits<double>::max()),
                                   Patterns::List(Patterns::Double ()),
                                   "The maximum global particle property values that will be used as a "
                                   "limiter for the quadratic least squares interpolation. The number of the input "
                                   "'Global particle property maximum' values separated by ',' has to be "
                                   "the same as the number of particle properties.");
                prm.declare_entry ("Global particle property minimum",
                                   boost::lexical_cast<std::string>(std::numeric_limits<double>::lowest()),
                                   Patterns::List(Patterns::Double ()),
                                   "The minimum global particle property that will be used as a "
                                   "limiter for the quadratic least squares interpolation. The number of the input "
                                   "'Global particle property minimum' values separated by ',' has to be "
                                   "the same as the number of particle properties.");
                prm.declare_entry("Use limiter", "false",
                                  Patterns::Bool (),
                                  "Whether to apply a global particle property limiting scheme to the interpolated "
                                  "particle properties.");

              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      QuadraticLeastSquares<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Interpolator");
            {
              prm.enter_subsection("Quadratic least squares");
              {
                use_global_values_limiter = prm.get_bool("Use limiter");
                if (use_global_values_limiter)
                  {
                    global_maximum_particle_properties = Utilities::string_to_double(Utilities::split_string_list(prm.get("Global particle property maximum")));
                    global_minimum_particle_properties = Utilities::string_to_double(Utilities::split_string_list(prm.get("Global particle property minimum")));

                    const Postprocess::Particles<dim> &particle_postprocessor =
                      this->get_postprocess_manager().template get_matching_postprocessor<const Postprocess::Particles<dim>>();
                    const auto &particle_property_information = particle_postprocessor.get_particle_world().get_property_manager().get_data_info();
                    const unsigned int n_property_components = particle_property_information.n_components();
                    const unsigned int n_internal_components = particle_property_information.get_components_by_field_name("internal: integrator properties");

                    // Check that if a global limiter is used, we were given the minimum and maximum value for each
                    // particle property that is not an internal property.
                    AssertThrow(global_minimum_particle_properties.size() == n_property_components - n_internal_components,
                                ExcMessage("Make sure that the size of list 'Global minimum particle property' "
                                           "is equivalent to the number of particle properties."));

                    // Check that if a global limiter is used, we were given the minimum and maximum value for each
                    // particle property that is not an internal property.
                    AssertThrow(global_maximum_particle_properties.size() == n_property_components - n_internal_components,
                                ExcMessage("Make sure that the size of list 'Global maximum particle property' "
                                           "is equivalent to the number of particle properties."));
                  }
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(QuadraticLeastSquares,
                                            "quadratic least squares",
                                            "Interpolates particle properties onto a vector of points using a "
                                            "quadratic least squares method. Note that deal.II must be configured "
                                            "with BLAS/LAPACK.")
    }
  }
}
