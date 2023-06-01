/*
  Copyright (C) 2015 - 2021 by the authors of the ASPECT code.

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

#include <aspect/particle/property/reference_position.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      ReferencePosition<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                               std::vector<double> &data) const
      {
        for (unsigned int i = 0; i < dim; ++i)
          data.push_back(position[i]);
      }

      template <int dim>
      void
      ReferencePosition<dim>::update_particle_property(const unsigned int data_position,
                                                       const Vector<double> &/*solution*/,
                                                       const std::vector<Tensor<1,dim>> &/*gradients*/,
                                                       typename ParticleHandler<dim>::particle_iterator &particle) const
      {
        for (unsigned int i = 0; i < dim; ++i)
          particle->get_properties()[data_position+i] = particle->get_reference_location()[i];
      }

      template <int dim>
      UpdateTimeFlags
      ReferencePosition<dim>::need_update() const
      {
        return update_output_step;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
                                                     ReferencePosition<dim>::get_property_information() const
      {
        const std::vector<std::pair<std::string,unsigned int>> property_information (1,std::make_pair("reference position",dim));
        return property_information;
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(ReferencePosition,
                                        "reference position",
                                        "Implementation of a plugin in which the particle "
                                        "property is defined as the current reference position.")
    }
  }
}
