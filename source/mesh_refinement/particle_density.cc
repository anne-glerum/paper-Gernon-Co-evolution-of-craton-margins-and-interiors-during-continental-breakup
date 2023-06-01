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


#include <aspect/mesh_refinement/particle_density.h>

#include <aspect/postprocess/particles.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    ParticleDensity<dim>::execute(Vector<float> &indicators) const
    {
      AssertThrow(this->get_postprocess_manager().template has_matching_postprocessor<Postprocess::Particles<dim>>(),
                  ExcMessage("The mesh refinement plugin `particle density' requires the "
                             "postprocessor plugin `particles' to be selected. Please activate the "
                             "particles or deactivate this mesh refinement plugin."));

      const Postprocess::Particles<dim> &particle_postprocessor =
        this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Particles<dim>>();

      const Particle::ParticleHandler<dim> &particle_handler = particle_postprocessor.get_particle_world().get_particle_handler();

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            // Note  that this refinement indicator will level out the number
            // of particles per cell, therefore creating fine cells in regions
            // of high particle density and coarse cells in low particle
            // density regions.
            indicators(cell->active_cell_index()) = static_cast<float>(particle_handler.n_particles_in_cell(cell));
          }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(ParticleDensity,
                                              "particle density",
                                              "A mesh refinement criterion that computes "
                                              "refinement indicators based on the density "
                                              "of particles. In practice this plugin "
                                              "equilibrates the number of particles per cell, "
                                              "leading to fine cells in high particle density regions "
                                              "and coarse cells in low particle density regions. "
                                              "This plugin is mostly useful for models with inhomogeneous "
                                              "particle density, e.g. when tracking an initial interface "
                                              "with a high particle density, or when the spatial particle "
                                              "density denotes the region of interest. Additionally, this "
                                              "plugin tends to balance the computational load between "
                                              "processes in parallel computations, because the particle "
                                              "and mesh density is more aligned.")
  }
}
