/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/heating_model/constant_heating.h>


namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    ConstantHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &/*material_model_inputs*/,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          // return a constant value
          heating_model_outputs.heating_source_terms[q] = radiogenic_heating_rate
                                                          * material_model_outputs.densities[q];
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }



    template <int dim>
    void
    ConstantHeating<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Constant heating");
        {
          prm.declare_entry ("Radiogenic heating rate", "0.",
                             Patterns::Double (0.),
                             "The specific rate of heating due to radioactive decay (or other bulk sources "
                             "you may want to describe). This parameter corresponds to the variable "
                             "$H$ in the temperature equation stated in the manual, and the heating "
                             "term is $\rho H$. Units: W/kg.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ConstantHeating<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Constant heating");
        {
          radiogenic_heating_rate    = prm.get_double ("Radiogenic heating rate");
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
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ConstantHeating,
                                  "constant heating",
                                  "Implementation of a model in which the heating "
                                  "rate is constant.")
  }
}
