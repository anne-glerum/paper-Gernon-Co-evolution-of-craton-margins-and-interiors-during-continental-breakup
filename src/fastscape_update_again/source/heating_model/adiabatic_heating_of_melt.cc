/*
  Copyright (C) 2015 - 2020 by the authors of the ASPECT code.

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


#include <aspect/heating_model/adiabatic_heating_of_melt.h>
#include <aspect/melt.h>
#include <aspect/simulator.h>
#include <deal.II/numerics/fe_field_function.h>

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    AdiabaticHeatingMelt<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      AssertThrow(this->include_melt_transport(),
                  ExcMessage ("Heating model Adiabatic heating with melt only works if melt transport is enabled."));

      // get the melt velocity
      const MaterialModel::MeltInputs<dim> *melt_in = material_model_inputs.template get_additional_input<MaterialModel::MeltInputs<dim>>();
      AssertThrow(melt_in != nullptr,
                  ExcMessage ("Need MeltInputs from the material model for adiabatic heating with melt!"));

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          const double porosity = material_model_inputs.composition[q][this->introspection().compositional_index_for_name("porosity")];

          if (!simplified_adiabatic_heating)
            heating_model_outputs.heating_source_terms[q] = ((-porosity * material_model_inputs.velocity[q] + porosity * melt_in->fluid_velocities[q])
                                                             * material_model_inputs.pressure_gradient[q])
                                                            * material_model_outputs.thermal_expansion_coefficients[q]
                                                            * material_model_inputs.temperature[q];
          else
            {
              const MaterialModel::MeltOutputs<dim> *melt_outputs = material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim>>();
              Assert(melt_outputs != nullptr, ExcMessage("Need MeltOutputs from the material model for adiabatic heating with melt."));

              heating_model_outputs.heating_source_terms[q] = (-porosity * material_model_inputs.velocity[q]
                                                               * this->get_gravity_model().gravity_vector(material_model_inputs.position[q]))
                                                              * material_model_outputs.thermal_expansion_coefficients[q]
                                                              * material_model_inputs.temperature[q]
                                                              * material_model_outputs.densities[q]
                                                              +
                                                              ((porosity * melt_in->fluid_velocities[q])
                                                               * this->get_gravity_model().gravity_vector(material_model_inputs.position[q]))
                                                              * material_model_outputs.thermal_expansion_coefficients[q]
                                                              * material_model_inputs.temperature[q]
                                                              * melt_outputs->fluid_densities[q];
            }

          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }

    template <int dim>
    void
    AdiabaticHeatingMelt<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Adiabatic heating of melt");
        {
          prm.declare_entry ("Use simplified adiabatic heating", "false",
                             Patterns::Bool (),
                             "A flag indicating whether the adiabatic heating should be simplified "
                             "from $\\alpha T (\\mathbf u \\cdot \\nabla p)$ to "
                             "$ \\alpha \\rho T (\\mathbf u \\cdot \\mathbf g) $.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    AdiabaticHeatingMelt<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Adiabatic heating of melt");
        {
          simplified_adiabatic_heating = prm.get_bool ("Use simplified adiabatic heating");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    AdiabaticHeatingMelt<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &output) const
    {
      MeltHandler<dim>::create_material_model_outputs(output);
    }



    template <int dim>
    void
    AdiabaticHeatingMelt<dim>::
    create_additional_material_model_inputs(MaterialModel::MaterialModelInputs<dim> &inputs) const
    {
      // we need the melt inputs for this adiabatic heating of melt
      if (inputs.template get_additional_input<MaterialModel::MeltInputs<dim>>() != nullptr)
        return;

      inputs.additional_inputs.push_back(
        std::make_unique<MaterialModel::MeltInputs<dim>> (inputs.position.size()));
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(AdiabaticHeatingMelt,
                                  "adiabatic heating of melt",
                                  "Implementation of a standard and a simplified model of "
                                  "adiabatic heating of melt. The full model implements the "
                                  "heating term \n"
                                  "$\\alpha T (-\\phi \\mathbf u_s \\cdot \\nabla p) "
                                  "+ \\alpha T (\\phi \\mathbf u_f \\cdot \\nabla p)$.\n"
                                  "For full adiabatic heating, "
                                  "this has to be used in combination with the heating model "
                                  "`adiabatic heating' to also include adiabatic heating for "
                                  "the solid part, and the full heating term is then "
                                  "$\\alpha T ((1-\\phi) \\mathbf u_s \\cdot \\nabla p) "
                                  "+ \\alpha T (\\phi \\mathbf u_f \\cdot \\nabla p)$.")
  }
}
