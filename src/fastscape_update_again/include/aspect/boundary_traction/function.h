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


#ifndef _aspect_boundary_traction_function_h
#define _aspect_boundary_traction_function_h

#include <aspect/boundary_traction/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace BoundaryTraction
  {
    using namespace dealii;

    /**
     * A class that implements traction boundary conditions based on a
     * functional description provided in the input file.
     *
     * @ingroup BoundaryTractions
     */
    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Function ();

        /**
         * Return the boundary traction as a function of position. The
         * (outward) normal vector to the domain is also provided as
         * a second argument.
         */
        Tensor<1,dim>
        boundary_traction (const types::boundary_id boundary_indicator,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const override;

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the boundary values will
         * next be evaluated. For the current class, the function passes to
         * the parsed function what the current time is.
         */
        void
        update () override;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * A function object representing the components of the traction.
         */
        Functions::ParsedFunction<dim> boundary_traction_function;

        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are depth, cartesian and spherical.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}


#endif
