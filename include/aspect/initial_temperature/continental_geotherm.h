/*
  Copyright (C) 2014 - 2019 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/



#ifndef _aspect_initial_temperature_continental_geotherm_h
#define _aspect_initial_temperature_continental_geotherm_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace InitialTemperature
  {

    /**
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class ContinentalGeotherm : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        ContinentalGeotherm ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        /**
         * Return the initial temperature as a function of depth,
         * based on the solution of the steady-state heat conduction
         * differential equation.
         */
        double initial_temperature (const Point<dim> &position) const override;

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
         * Surface temperature
         */
        double T0;

        /**
         * Value of the isotherm defining the
         * Lithosphere-Asthenosphere boundary.
         */
        double LAB_isotherm;

        /**
         * Vector for the thicknesses of the compositional fields
         * representing the layers 'upper_crust', 'lower_crust' and 'lithospheric_mantle'.
         */
        std::vector<double> thicknesses;

        /**
         * Vector for the heat production rates of the different
         * compositional fields, read from parameter file.
         */
        std::vector<double> heat_productivities;

        /**
         * Vector for the thermal conductivities of the different
         * compositional fields, read from parameter file.
         */
        std::vector<double> conductivities;

        /**
         * Vector for the densities of the different compositional
         * fields, read from parameter file.
         */
        std::vector<double> densities;
    };
  }
}


#endif
