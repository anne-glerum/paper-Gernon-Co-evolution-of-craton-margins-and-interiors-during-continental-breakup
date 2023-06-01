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


#ifndef _aspect_mesh_refinement_boundary_h
#define _aspect_mesh_refinement_boundary_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
namespace MeshRefinement
{

  /**
   * A class that implements a mesh refinement criterion that refines the
   * mesh on selected boundaries. This is useful for cases where one wants
   * to accurately model processes at a boundary.  Frequently, there are
   * also strong temperature or compositional gradients at boundaries, so
   * this refinement criterion can be redundant.
   *
   * @ingroup MeshRefinement
   */
  template <int dim>
  class Boundary : public Interface<dim>,
    public SimulatorAccess<dim>
  {
    public:

      /**
       * Execute this mesh refinement criterion.
       *
       * @param[out] error_indicators A vector that for every active cell of
       * the current mesh (which may be a partition of a distributed mesh)
       * provides an error indicator. This vector will already have the
       * correct size when the function is called.
       */
      void
      execute (Vector<float> &error_indicators) const override;

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
       * A set of boundary indicators marked for refinement
       */
      std::set<types::boundary_id> boundary_refinement_indicators;
  };
}
}

#endif
