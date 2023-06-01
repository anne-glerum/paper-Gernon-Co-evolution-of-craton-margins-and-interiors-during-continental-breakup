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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/



#ifndef _aspect_heating_model_radioactive_decay_h
#define _aspect_heating_model_radioactive_decay_h

#include <aspect/simulator_access.h>
#include <aspect/heating_model/interface.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    /**
     * A class that implements a heating model based on radioactive decay.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class RadioactiveDecay : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        RadioactiveDecay ();

        /**
         * Return the specific heating rate as calculated by radioactive
         * decay.
         */
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;

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
         * Number of radio active heating elements.
         */
        unsigned int                   n_radio_heating_elements;

        /**
         * Store the half life of different elements.
         */
        std::vector<double>            half_decay_times;

        /**
         * Store the unit heating rate of different elements.
         */
        std::vector<double>            radioactive_heating_rates;

        /**
         * Store the initial concentration in the crust.
         */
        std::vector<double>            radioactive_initial_concentrations_crust;

        /**
         * Store the initial concentration in the mantle.
         */
        std::vector<double>            radioactive_initial_concentrations_mantle;

        /**
         * Whether crust defined by composition or depth
         */
        bool                           is_crust_defined_by_composition;

        /**
         * Depth of the crust.
         */
        double                         crust_depth;

        /**
         * Composition number of crust.
         */
        unsigned int                   crust_composition_num;
    };
  }
}


#endif
