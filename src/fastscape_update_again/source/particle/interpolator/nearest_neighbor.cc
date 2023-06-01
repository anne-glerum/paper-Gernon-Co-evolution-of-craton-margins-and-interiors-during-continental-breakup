/*
 Copyright (C) 2017 - 2019 by the authors of the ASPECT code.

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

#include <aspect/particle/interpolator/nearest_neighbor.h>
#include <aspect/postprocess/particles.h>
#include <aspect/simulator.h>

#include <deal.II/grid/grid_tools.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      template <int dim>
      std::vector<std::vector<double>>
                                    NearestNeighbor<dim>::properties_at_points(const ParticleHandler<dim> &particle_handler,
                                                                               const std::vector<Point<dim>> &positions,
                                                                               const ComponentMask &selected_properties,
                                                                               const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        typename parallel::distributed::Triangulation<dim>::active_cell_iterator found_cell;

        if (cell->state() == IteratorState::invalid)
          {
            // We can not simply use one of the points as input for find_active_cell_around_point
            // because for vertices of mesh cells we might end up getting ghost_cells as return value
            // instead of the local active cell. So make sure we are well in the inside of a cell.
            Assert(positions.size() > 0,
                   ExcMessage("The particle property interpolator was not given any "
                              "positions to evaluate the particle properties at."));

            const Point<dim> approximated_cell_midpoint = std::accumulate (positions.begin(), positions.end(), Point<dim>())
                                                          / static_cast<double> (positions.size());

            found_cell =
              (GridTools::find_active_cell_around_point<> (this->get_mapping(),
                                                           this->get_triangulation(),
                                                           approximated_cell_midpoint)).first;
          }
        else
          found_cell = cell;

        const typename ParticleHandler<dim>::particle_iterator_range particle_range = particle_handler.particles_in_cell(found_cell);

        const unsigned int n_particles = std::distance(particle_range.begin(),particle_range.end());
        const unsigned int n_particle_properties = particle_handler.n_properties_per_particle();

        std::vector<double> temp(n_particle_properties, 0.0);
        std::vector<std::vector<double>> point_properties(positions.size(), temp);

        for (unsigned int pos_idx=0; pos_idx < positions.size(); ++pos_idx)
          {
            double minimum_distance = std::numeric_limits<double>::max();
            if (n_particles > 0)
              {
                typename ParticleHandler<dim>::particle_iterator nearest_neighbor;
                for (typename ParticleHandler<dim>::particle_iterator particle = particle_range.begin();
                     particle != particle_range.end(); ++particle)
                  {
                    const double dist = (positions[pos_idx] - particle->get_location()).norm_square();
                    if (dist < minimum_distance)
                      {
                        minimum_distance = dist;
                        nearest_neighbor = particle;
                      }
                  }
                const dealii::ArrayView<const double> neighbor_props = nearest_neighbor->get_properties();
                for (unsigned int i = 0; i < n_particle_properties; ++i)
                  if (selected_properties[i])
                    point_properties[pos_idx][i] = neighbor_props[i];
              }
            else
              {
                std::vector<typename parallel::distributed::Triangulation<dim>::active_cell_iterator> neighbors;
                GridTools::get_active_neighbors<parallel::distributed::Triangulation<dim>>(found_cell,neighbors);

                unsigned int nearest_neighbor_cell = numbers::invalid_unsigned_int;
                for (unsigned int i=0; i<neighbors.size(); ++i)
                  {
                    // Only recursively call this function if the neighbor cell contains
                    // particles (else we end up in an endless recursion)
                    if (particle_handler.n_particles_in_cell(neighbors[i]) == 0)
                      continue;

                    const double dist = (positions[pos_idx] - neighbors[i]->center()).norm_square();
                    if (dist < minimum_distance)
                      {
                        minimum_distance = dist;
                        nearest_neighbor_cell = i;
                      }
                  }

                if (!allow_cells_without_particles)
                  {
                    AssertThrow(nearest_neighbor_cell != numbers::invalid_unsigned_int,
                                ExcMessage("A cell and all of its neighbors do not contain any particles. "
                                           "The `nearest neighbor' interpolation scheme does not support this case unless specified "
                                           "in Allow cells without particles."));
                  }

                if (nearest_neighbor_cell != numbers::invalid_unsigned_int)
                  {
                    point_properties[pos_idx] = properties_at_points(particle_handler,
                                                                     std::vector<Point<dim>> (1,positions[pos_idx]),
                                                                     selected_properties,
                                                                     neighbors[nearest_neighbor_cell])[0];
                  }
                else if (allow_cells_without_particles && nearest_neighbor_cell == numbers::invalid_unsigned_int)
                  {
                    for (unsigned int i = 0; i < n_particle_properties; ++i)
                      if (selected_properties[i])
                        point_properties[pos_idx][i] = 0;
                  }

              }

          }

        return point_properties;
      }



      template <int dim>
      void
      NearestNeighbor<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.declare_entry ("Allow cells without particles", "false",
                               Patterns::Bool (),
                               "By default, every cell needs to contain particles to use this interpolator "
                               "plugin. If this parameter is set to true, cells are allowed to have no particles. "
                               "In case both the current cell and its neighbors are empty, "
                               "the interpolator will return 0 for the current cell's properties.");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }



      template <int dim>
      void
      NearestNeighbor<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            allow_cells_without_particles = prm.get_bool("Allow cells without particles");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
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
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(NearestNeighbor,
                                            "nearest neighbor",
                                            "Return the properties of the nearest neighboring particle "
                                            "in the current cell, or nearest particle in nearest neighboring "
                                            "cell if current cell is empty. "
                                            "In case the neighboring cells are also empty, and 'Allow cells "
                                            "without particles' is set to true, the interpolator returns 0. "
                                            "Otherwise, an exception is thrown. ")
    }
  }
}
