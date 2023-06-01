/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#include <aspect/mesh_refinement/velocity.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/error_estimator.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    Velocity<dim>::execute(Vector<float> &indicators) const
    {
      indicators = 0;

      KellyErrorEstimator<dim>::estimate (this->get_mapping(),
                                          this->get_dof_handler(),
                                          QGauss<dim-1>(this->introspection().polynomial_degree.velocities+1),
                                          std::map<types::boundary_id,const Function<dim>*>(),
                                          this->get_solution(),
                                          indicators,
                                          this->introspection().component_masks.velocities,
                                          nullptr,
                                          0,
                                          this->get_triangulation().locally_owned_subdomain());
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Velocity,
                                              "velocity",
                                              "A mesh refinement criterion that computes "
                                              "refinement indicators from the velocity field."
                                              "\n\n"
                                              "The way these indicators are computed is by "
                                              "evaluating the `Kelly error indicator' on the "
                                              "velocity field. This error indicator takes the "
                                              "finite element approximation of the velocity "
                                              "field and uses it to compute an approximation "
                                              "of the second derivatives of the velocity for "
                                              "each cell. This approximation is then multiplied "
                                              "by an appropriate power of the cell's diameter "
                                              "to yield an indicator for how large the error "
                                              "is likely going to be on this cell. This "
                                              "construction rests on the observation that for "
                                              "many partial differential equations, the error "
                                              "on each cell is proportional to some power of "
                                              "the cell's diameter times the second derivatives "
                                              "of the solution on that cell."
                                              "\n\n"
                                              "For complex equations such as those we solve "
                                              "here, this observation may not be strictly "
                                              "true in the mathematical sense, but it often "
                                              "yields meshes that are surprisingly good.")
  }
}
