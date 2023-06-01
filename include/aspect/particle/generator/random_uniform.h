/*
 Copyright (C) 2015 - 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_generator_random_uniform_h
#define _aspect_particle_generator_random_uniform_h

#include <aspect/particle/generator/probability_density_function.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Generate random uniform distribution of particles over entire simulation domain.
       *
       * @ingroup ParticleGenerators
       */
      template <int dim>
      class RandomUniform : public ProbabilityDensityFunction<dim>
      {
        public:

        private:
          /**
           * Returns the weight of one cell, which is interpreted as the probability
           * to generate particles in this cell.
           */
          double
          get_cell_weight (const typename DoFHandler<dim>::active_cell_iterator &cell) const override;
      };

    }
  }
}

#endif
