/*
  Copyright (C) 2015 - 2018 by the authors of the ASPECT code.

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

#include <aspect/particle/integrator/rk_4.h>
#include <aspect/particle/property/interface.h>
#include <aspect/particle/world.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      template <int dim>
      RK4<dim>::RK4()
        :
        integrator_substep(0)
      {}



      template <int dim>
      void
      RK4<dim>::initialize ()
      {
        const auto &property_information = this->get_particle_world().get_property_manager().get_data_info();

        property_index_k[0] = property_information.get_position_by_field_name("internal: integrator properties");
        property_index_k[1] = property_index_k[0] + dim;
        property_index_k[2] = property_index_k[1] + dim;
        property_index_k[3] = property_index_k[2] + dim;
      }



      template <int dim>
      void
      RK4<dim>::local_integrate_step(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                     const typename ParticleHandler<dim>::particle_iterator &end_particle,
                                     const std::vector<Tensor<1,dim>> &old_velocities,
                                     const std::vector<Tensor<1,dim>> &velocities,
                                     const double dt)
      {
        Assert(static_cast<unsigned int> (std::distance(begin_particle, end_particle)) == old_velocities.size(),
               ExcMessage("The particle integrator expects the old velocity vector to be of equal size "
                          "to the number of particles to advect. For some unknown reason they are different, "
                          "most likely something went wrong in the calling function."));

        Assert(old_velocities.size() == velocities.size(),
               ExcMessage("The particle integrator expects the velocity vector to be of equal size "
                          "to the number of particles to advect. For some unknown reason they are different, "
                          "most likely something went wrong in the calling function."));

        const bool geometry_has_periodic_boundary = (this->get_geometry_model().get_periodic_boundary_pairs().size() != 0);

        typename std::vector<Tensor<1,dim>>::const_iterator old_velocity = old_velocities.begin();
        typename std::vector<Tensor<1,dim>>::const_iterator velocity = velocities.begin();

        std::array<Point<dim>,5> k;
        for (typename ParticleHandler<dim>::particle_iterator it = begin_particle;
             it != end_particle; ++it, ++velocity, ++old_velocity)
          {
            ArrayView<double> properties = it->get_properties();

            if (integrator_substep == 0)
              {
                k[0] = it->get_location();
                k[1] = dt * (*old_velocity);

                Point<dim> new_location = k[0] + 0.5 * k[1];

                // Check if we crossed a periodic boundary and if necessary adjust positions
                if (geometry_has_periodic_boundary)
                  this->get_geometry_model().adjust_positions_for_periodicity(new_location,
                                                                              ArrayView<Point<dim>>(&k[0],2));

                for (unsigned int i=0; i<dim; ++i)
                  {
                    properties[property_index_k[0] + i] = k[0][i];
                    properties[property_index_k[1] + i] = k[1][i];
                  }

                it->set_location(new_location);
              }
            else if (integrator_substep == 1)
              {
                k[2] = dt * ((*old_velocity) + (*velocity)) / 2.0;
                for (unsigned int i=0; i<dim; ++i)
                  k[0][i] = properties[property_index_k[0] + i];

                Point<dim> new_location = k[0] + 0.5 * k[2];

                if (geometry_has_periodic_boundary)
                  {
                    for (unsigned int i=0; i<dim; ++i)
                      k[1][i] = properties[property_index_k[1] + i];

                    this->get_geometry_model().adjust_positions_for_periodicity(new_location,
                                                                                ArrayView<Point<dim>>(&k[0],3));

                    for (unsigned int i=0; i<dim; ++i)
                      {
                        properties[property_index_k[0] + i] = k[0][i];
                        properties[property_index_k[1] + i] = k[1][i];
                      }
                  }

                for (unsigned int i=0; i<dim; ++i)
                  properties[property_index_k[2] + i] = k[2][i];

                it->set_location(new_location);
              }
            else if (integrator_substep == 2)
              {
                k[3] = dt * (*old_velocity + *velocity) / 2.0;
                for (unsigned int i=0; i<dim; ++i)
                  k[0][i] = properties[property_index_k[0] + i];

                Point<dim> new_location = k[0] + k[3];

                if (geometry_has_periodic_boundary)
                  {
                    for (unsigned int i=0; i<dim; ++i)
                      {
                        k[1][i] = properties[property_index_k[1] + i];
                        k[2][i] = properties[property_index_k[2] + i];
                      }

                    this->get_geometry_model().adjust_positions_for_periodicity(new_location,
                                                                                ArrayView<Point<dim>>(&k[0],4));

                    for (unsigned int i=0; i<dim; ++i)
                      {
                        properties[property_index_k[0] + i] = k[0][i];
                        properties[property_index_k[1] + i] = k[1][i];
                        properties[property_index_k[2] + i] = k[2][i];
                      }
                  }

                for (unsigned int i=0; i<dim; ++i)
                  {
                    properties[property_index_k[3] + i] = k[3][i];
                  }

                it->set_location(new_location);
              }
            else if (integrator_substep == 3)
              {
                k[4] = dt * (*velocity);

                for (unsigned int i=0; i<dim; ++i)
                  {
                    k[0][i] = properties[property_index_k[0] + i];
                    k[1][i] = properties[property_index_k[1] + i];
                    k[2][i] = properties[property_index_k[2] + i];
                    k[3][i] = properties[property_index_k[3] + i];
                  }

                Point<dim> new_location = k[0] + (k[1] + 2.0*k[2] + 2.0*k[3] + k[4])/6.0;

                // No need to fix intermediate values, this is the last integrator step
                if (geometry_has_periodic_boundary)
                  this->get_geometry_model().adjust_positions_for_periodicity(new_location);

                it->set_location(new_location);
              }
            else
              {
                Assert(false,
                       ExcMessage("The RK4 integrator should never continue after four integration steps."));
              }
          }
      }



      template <int dim>
      bool
      RK4<dim>::new_integration_step()
      {
        integrator_substep = (integrator_substep+1)%4;

        // Continue until we're at the last step
        return (integrator_substep != 0);
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      ASPECT_REGISTER_PARTICLE_INTEGRATOR(RK4,
                                          "rk4",
                                          "Runge Kutta fourth order integrator, where "
                                          "$y_{n+1} = y_n + \\frac{1}{6} k_1 + \\frac{1}{3} k_2 "
                                          "+ \\frac{1}{3} k_3 + \\frac{1}{6} k_4$ "
                                          "and $k_1$, $k_2$, $k_3$, $k_4$ are defined as usual.")
    }
  }
}
