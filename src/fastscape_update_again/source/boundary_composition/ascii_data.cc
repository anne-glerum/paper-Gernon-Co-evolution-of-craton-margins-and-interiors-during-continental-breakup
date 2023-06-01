/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/boundary_composition/ascii_data.h>

#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace BoundaryComposition
  {
    template <int dim>
    AsciiData<dim>::AsciiData ()
    {}

    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      const std::set<types::boundary_id> boundary_ids = this->get_fixed_composition_boundary_indicators();

      Utilities::AsciiDataBoundary<dim>::initialize(boundary_ids,
                                                    this->n_compositional_fields());
    }


    template <int dim>
    void
    AsciiData<dim>::update ()
    {
      Interface<dim>::update ();

      Utilities::AsciiDataBoundary<dim>::update();
    }


    template <int dim>
    double
    AsciiData<dim>::
    boundary_composition (const types::boundary_id boundary_indicator,
                          const Point<dim> &position,
                          const unsigned int compositional_field) const
    {
      return Utilities::AsciiDataBoundary<dim>::get_data_component(boundary_indicator,
                                                                   position,
                                                                   compositional_field);
    }

    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                              "$ASPECT_SOURCE_DIR/data/boundary-composition/ascii-data/test/",
                                                              "box_2d_%s.%d.txt");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary composition model");
      {
        Utilities::AsciiDataBoundary<dim>::parse_parameters(prm);
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryComposition
  {
    ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(AsciiData,
                                               "ascii data",
                                               "Implementation of a model in which the boundary "
                                               "composition is derived from files containing data "
                                               "in ascii format. Note the required format of the "
                                               "input data: The first lines may contain any number of comments "
                                               "if they begin with `#', but one of these lines needs to "
                                               "contain the number of grid points in each dimension as "
                                               "for example `# POINTS: 3 3'. "
                                               "The order of the data columns "
                                               "has to be `x', `composition1', `composition2', "
                                               "etc. in a 2d model and "
                                               "`x', `y', `composition1', `composition2', "
                                               "etc., in a 3d model, according "
                                               "to the number of compositional fields, which means that "
                                               "there has to be a single column "
                                               "for every composition in the model. "
                                               "Note that the data in the input "
                                               "files need to be sorted in a specific order: "
                                               "the first coordinate needs to ascend first, "
                                               "followed by the second in order to "
                                               "assign the correct data to the prescribed coordinates."
                                               "If you use a spherical model, "
                                               "then the assumed grid changes. `x' will be replaced by "
                                               "the radial distance of the point to the bottom of the model, "
                                               "`y' by the azimuth angle and `z' by the polar angle measured "
                                               "positive from the north pole. The grid will be assumed to be "
                                               "a latitude-longitude grid. Note that the order "
                                               "of spherical coordinates is `r', `phi', `theta' "
                                               "and not `r', `theta', `phi', since this allows "
                                               "for dimension independent expressions.")
  }
}
