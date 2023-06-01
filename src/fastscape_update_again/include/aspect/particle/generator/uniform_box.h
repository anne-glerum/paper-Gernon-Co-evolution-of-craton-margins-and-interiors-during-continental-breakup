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

#ifndef _aspect_particle_generator_uniform_box_h
#define _aspect_particle_generator_uniform_box_h

#include <aspect/particle/generator/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Generate a uniform distribution of particles in a box region in the
       * model domain. Uniform here means the particles will be generated with
       * an equal spacing in each spatial dimension. Note that in order
       * to produce a regular distribution the number of generated
       * particles might not exactly match the one specified in the
       * input file.
       *
       * @ingroup ParticleGenerators
       */
      template <int dim>
      class UniformBox : public Interface<dim>
      {
        public:
          /**
           * Generate a uniformly randomly distributed set of particles in a
           * box-like subdomain of the global domain.
           *
           * @param [in,out] particles A multimap between cells and their
           * particles. This map will be filled in this function.
           */
          void
          generate_particles(std::multimap<Particles::internal::LevelInd, Particle<dim>> &particles) override;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * Number of initial particles to create.
           */
          types::particle_index n_particles;

          /**
           * The minimum coordinates of the particle region, i.e. one corner of
           * the n-dimensional box region in which particles are generated.
           */
          Point<dim> P_min;

          /**
           * The maximum coordinates of the particle region, i.e. the opposite
           * corner of the n-dimensional box region in which particles are
           * generated.
           */
          Point<dim> P_max;
      };

    }
  }
}

#endif
