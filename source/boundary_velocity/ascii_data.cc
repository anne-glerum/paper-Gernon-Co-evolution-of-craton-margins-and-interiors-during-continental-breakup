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


#include <aspect/global.h>
#include <aspect/boundary_velocity/ascii_data.h>

#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace BoundaryVelocity
  {
    template <int dim>
    AsciiData<dim>::AsciiData ()
    {}


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      for (const auto &p : this->get_boundary_velocity_manager().get_active_boundary_velocity_conditions())
        {
          for (const auto &plugin : p.second)
            if (plugin.get() == this)
              boundary_ids.insert(p.first);
        }
      AssertThrow(*(boundary_ids.begin()) != numbers::invalid_boundary_id,
                  ExcMessage("Did not find the boundary indicator for the prescribed data plugin."));

      Utilities::AsciiDataBoundary<dim>::initialize(boundary_ids,
                                                    dim);
    }



    template <int dim>
    void
    AsciiData<dim>::update ()
    {
      Interface<dim>::update ();

      Utilities::AsciiDataBoundary<dim>::update();
    }


    template <int dim>
    Tensor<1,dim>
    AsciiData<dim>::
    boundary_velocity (const types::boundary_id ,
                       const Point<dim> &position) const
    {
      Tensor<1,dim> velocity;
      for (unsigned int i = 0; i < dim; i++)
        velocity[i] = Utilities::AsciiDataBoundary<dim>::get_data_component(*(boundary_ids.begin()),
                                                                            position,
                                                                            i);
      if (use_spherical_unit_vectors)
        velocity = Utilities::Coordinates::spherical_to_cartesian_vector(velocity, position);

      return velocity;
    }


    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary velocity model");
      {
        Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                              "$ASPECT_SOURCE_DIR/data/boundary-velocity/ascii-data/test/",
                                                              "box_2d_%s.%d.txt");
        prm.enter_subsection("Ascii data model");
        {
          prm.declare_entry ("Use spherical unit vectors", "false",
                             Patterns::Bool (),
                             "Specify velocity as r, phi, and theta components "
                             "instead of x, y, and z. Positive velocities point up, east, "
                             "and north (in 3D) or out and clockwise (in 2D). "
                             "This setting only makes sense for spherical geometries."
                            );
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary velocity model");
      {
        Utilities::AsciiDataBoundary<dim>::parse_parameters(prm);

        if (this->convert_output_to_years() == true)
          {
            this->scale_factor               /= year_in_seconds;
          }
        prm.enter_subsection("Ascii data model");
        {
          use_spherical_unit_vectors = prm.get_bool("Use spherical unit vectors");
          if (use_spherical_unit_vectors)
            AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical,
                         ExcMessage ("Spherical unit vectors should not be used "
                                     "when geometry model is not spherical."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryVelocity
  {
    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(AsciiData,
                                            "ascii data",
                                            "Implementation of a model in which the boundary "
                                            "velocity is derived from files containing data "
                                            "in ascii format. Note the required format of the "
                                            "input data: The first lines may contain any number of comments "
                                            "if they begin with `#', but one of these lines needs to "
                                            "contain the number of grid points in each dimension as "
                                            "for example `# POINTS: 3 3'. "
                                            "The order of the data columns "
                                            "has to be `x', `velocity${}_x$', `velocity${}_y$' in a 2d model "
                                            "or `x', `y', `velocity${}_x$', `velocity${}_y$', "
                                            "`velocity${}_z$' in a 3d model. "
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
                                            "for dimension independent expressions. "
                                            "Velocities can be specified using cartesian "
                                            "(by default) or spherical unit vectors. "
                                            "No matter which geometry model is chosen, "
                                            "the unit of the velocities is assumed to be "
                                            "m/s or m/yr depending on the `Use years in output "
                                            "instead of seconds' flag. If you provide velocities "
                                            "in cm/yr, set the `Scale factor' option to 0.01.")
  }
}
