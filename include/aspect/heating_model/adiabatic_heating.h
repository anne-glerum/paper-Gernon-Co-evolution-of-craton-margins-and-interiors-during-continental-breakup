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


#ifndef _aspect_heating_model_adiabatic_heating_h
#define _aspect_heating_model_adiabatic_heating_h

#include <aspect/heating_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    /**
     * A class that implements a standard adiabatic heating rate.
     *
     * This adds the term from adiabatic compression heating
     *    $ \alpha T (\mathbf u \cdot \nabla p) $
     * where we use the definition of
     *    $ \alpha = - \frac{1}{\rho} \frac{\partial \rho}{\partial T} $
     * Note: this term is often simplified using the relationship
     *    $ \rho \mathbf g = - \nabla p $
     * to yield
     *    $ - \alpha \rho T (\mathbf u \cdot \mathbf g) $
     *
     * The user can specify if the simplification should be used
     * by setting the corresponding flag in the input file.
     * Also see the Equations section in the manual.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class AdiabaticHeating : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Returns whether we compute a simplified adiabatic heating term by
         * neglecting dynamic pressure effects. I.e. we use
         *     $ \alpha T (\mathbf u \cdot \nabla p) $
         * as adiabatic heating term if this function returns false, and
         *    $ - \alpha \rho T (\mathbf u \cdot \mathbf g) $
         * if it returns true.
         */
        bool
        use_simplified_adiabatic_heating() const;

        /**
         * Compute the heating model outputs for this class.
         */
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */

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

        /**
         * @}
         */

      private:
        bool simplified_adiabatic_heating;
    };
  }
}


#endif
