/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#ifndef _aspect_mesh_deformation_free_surface_h
#define _aspect_mesh_deformation_free_surface_h

#include <aspect/mesh_deformation/interface.h>

#include <aspect/simulator_access.h>
#include <aspect/simulator/assemblers/interface.h>


namespace aspect
{
  using namespace dealii;

  namespace MeshDeformation
  {
    /**
     * A plugin that computes the deformation of surface
     * vertices according to the solution of the flow problem.
     * In particular this means if the surface of the domain is
     * left open to flow, this flow will carry the mesh with it.
     */
    template<int dim>
    class FreeSurface : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize function, which connects the set_assemblers function
         * to the appropriate Simulator signal.
         */
        void initialize() override;

        /**
         * A function that creates constraints for the velocity of certain mesh
         * vertices (e.g. the surface vertices) for a specific boundary.
         * The calling class will respect
         * these constraints when computing the new vertex positions.
         */
        void
        compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                 AffineConstraints<double> &mesh_velocity_constraints,
                                                 const std::set<types::boundary_id> &boundary_id) const override;

        /**
         * Declare parameters for the free surface handling.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters for the free surface handling.
         */
        void parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Project the Stokes velocity solution onto the
         * free surface. Called by make_constraints()
         */
        void project_velocity_onto_boundary (const DoFHandler<dim> &free_surface_dof_handler,
                                             const IndexSet &mesh_locally_owned,
                                             const IndexSet &mesh_locally_relevant,
                                             LinearAlgebra::Vector &output) const;

        /**
         * A struct for holding information about how to advect the free surface.
         */
        struct SurfaceAdvection
        {
          enum Direction { normal, vertical };
        };

        /**
         * Stores whether to advect the free surface in the normal direction
         * or the direction of the local vertical.
         */
        typename SurfaceAdvection::Direction advection_direction;
    };
  }
}


#endif
