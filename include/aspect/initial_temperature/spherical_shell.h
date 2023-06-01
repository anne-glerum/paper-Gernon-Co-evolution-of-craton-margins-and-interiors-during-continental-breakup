/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_temperature_spherical_shell_h
#define _aspect_initial_temperature_spherical_shell_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

    /**
     * A class that describes a perturbed initial temperature field for the
     * spherical shell.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class SphericalHexagonalPerturbation : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
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
         * The angular mode is the number of perturbations to apply to the
         * spherical shell. Historically, this was permanently set to 6 (hence
         * the class name SphericalHexagonalPerturbation) The default is 6 in
         * order to provide backwards compatibility.
         */
        unsigned int angular_mode;

        /**
         * The rotation offset describes the number of degrees to rotate the
         * perturbation counterclockwise. Setting the rotation offset to 0
         * will cause one of the perturbations to point north/up. Rotation
         * offset is set to -45 degrees by default in order to provide
         * backwards compatibility.
         */
        double rotation_offset;

        /**
         * Outer radius.
         */
        double R1;

    };


    /**
     * A class that describes a perturbed initial temperature field for the
     * spherical shell.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class SphericalGaussianPerturbation : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        /**
         * Constructor.
         */
        SphericalGaussianPerturbation<dim>();

        /**
         * Return the initial temperature as a function of position.
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
        double angle;
        double depth;
        double amplitude;
        double sigma;
        double sign;
        unsigned int npoint;

        std::vector<double> radial_position;
        std::vector<double> geotherm;

        /**
         * Inner radius.
         */
        double R0;

        /**
         * Outer radius.
         */
        double R1;
    };
  }
}

#endif
