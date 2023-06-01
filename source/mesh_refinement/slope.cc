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


#include <aspect/mesh_refinement/slope.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    Slope<dim>::execute(Vector<float> &indicators) const
    {
      const types::boundary_id top_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      indicators = 0;

      QMidpoint<dim-1> quadrature;
      UpdateFlags update_flags = UpdateFlags(update_normal_vectors | update_quadrature_points );
      FEFaceValues<dim> fe_face_values (this->get_mapping(), this->get_fe(), quadrature, update_flags);

      // iterate over all of the cells, get a face normal, then
      // dot it with the local gravity

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          {
            const unsigned int idx = cell->active_cell_index();
            for (const unsigned int face_no : cell->face_indices())
              if (cell->face(face_no)->at_boundary())
                {
                  const types::boundary_id boundary_indicator
                    = cell->face(face_no)->boundary_id();

                  // Use cases for this plugin include a deforming mesh,
                  // or a fixed mesh with initial topography
                  if ( (this->get_parameters().mesh_deformation_enabled &&
                        this->get_mesh_deformation_boundary_indicators().find(boundary_indicator) !=
                        this->get_mesh_deformation_boundary_indicators().end()) ||
                       (Plugins::plugin_type_matches<const InitialTopographyModel::ZeroTopography<dim>>(this->get_initial_topography_model()) &&
                        boundary_indicator == top_boundary_id)  )
                    {
                      fe_face_values.reinit(cell, face_no);

                      const Tensor<1,dim> normal( fe_face_values.normal_vector(0) ); // Only one q point
                      const Point<dim> midpoint = fe_face_values.quadrature_point(0);
                      const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(midpoint);

                      indicators(idx) = std::acos( std::abs ( normal * gravity / gravity.norm() ) ) // Don't care whether gravity is in the opposite direction
                                        * std::pow( cell->diameter(), double(dim-1)); // scale with approximate surface area of the cell
                      break;  // no need to loop over the rest of the faces
                    }
                }
          }

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Slope,
                                              "slope",
                                              "A class that implements a mesh refinement criterion intended for "
                                              "use with deforming mesh boundaries, like the free surface. "
                                              "It calculates a local slope based on "
                                              "the angle between the surface normal and the local gravity vector. "
                                              "Cells with larger angles are marked for refinement."
                                              "\n\n"
                                              "To use this refinement criterion, you may want to combine "
                                              "it with other refinement criteria, setting the 'Normalize "
                                              "individual refinement criteria' flag and using the `max' "
                                              "setting for 'Refinement criteria merge operation'.")
  }
}
