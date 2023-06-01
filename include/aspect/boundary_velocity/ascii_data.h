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


#ifndef _aspect_boundary_velocity_ascii_data_h
#define _aspect_boundary_velocity_ascii_data_h

#include <aspect/boundary_velocity/interface.h>

#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace BoundaryVelocity
  {
    using namespace dealii;

    /**
     * A class that implements prescribed velocity boundary conditions
     * determined from a AsciiData input file.
     *
     * @ingroup BoundaryVelocities
     */
    template <int dim>
    class AsciiData : public Utilities::AsciiDataBoundary<dim>, public Interface<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        AsciiData ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataBoundary<dim>::initialize;

        /**
         * A function that is called at the beginning of each time step. For
         * the current plugin, this function loads the next data files if
         * necessary and outputs a warning if the end of the set of data files
         * is reached.
         */
        void
        update () override;

        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function returns value from the text files.
         */
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
                           const Point<dim> &position) const override;

        // avoid -Woverloaded-virtual warning until the deprecated function
        // is removed from the interface:
        using Interface<dim>::boundary_velocity;

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
        std::set<types::boundary_id> boundary_ids;

        /**
         * Whether to specify velocity in x, y, z components, or
         * r, phi, theta components.
         */
        bool use_spherical_unit_vectors;
    };
  }
}


#endif
