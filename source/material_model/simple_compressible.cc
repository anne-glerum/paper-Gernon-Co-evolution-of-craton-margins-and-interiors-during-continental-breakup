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

#include <aspect/material_model/simple_compressible.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    SimpleCompressible<dim>::
    evaluate(const MaterialModelInputs<dim> &in,
             MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          const Point<dim> position = in.position[i];
          const double temperature = in.temperature[i];
          const double pressure = in.pressure[i];

          out.viscosities[i] = constant_rheology.compute_viscosity();
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = k_value;
          out.thermal_expansion_coefficients[i] = thermal_alpha;

          double rho = reference_rho * std::exp(reference_compressibility * (pressure - this->get_surface_pressure()));
          rho *= (1 - thermal_alpha * (temperature - this->get_adiabatic_conditions().temperature(position)));

          out.densities[i] = rho;
          out.compressibilities[i] = reference_compressibility; // 1/rho drho/dp
          out.entropy_derivative_pressure[i] = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;
          // Change in composition due to chemical reactions at the
          // given positions. The term reaction_terms[i][c] is the
          // change in compositional field c at point i.
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }

    template <int dim>
    double
    SimpleCompressible<dim>::
    reference_viscosity () const
    {
      return constant_rheology.compute_viscosity();
    }



    template <int dim>
    bool
    SimpleCompressible<dim>::
    is_compressible () const
    {
      return (reference_compressibility != 0);
    }



    template <int dim>
    void
    SimpleCompressible<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple compressible model");
        {
          prm.declare_entry ("Reference density", "3300.",
                             Patterns::Double (0.),
                             "Reference density $\\rho_0$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Reference specific heat", "1250.",
                             Patterns::Double (0.),
                             "The value of the specific heat $C_p$. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0.),
                             "The value of the thermal expansion coefficient $\\alpha$. "
                             "Units: \\si{\\per\\kelvin}.");
          prm.declare_entry ("Reference compressibility", "4e-12",
                             Patterns::Double (0.),
                             "The value of the reference compressibility. "
                             "Units: \\si{\\per\\pascal}.");

          Rheology::ConstantViscosity::declare_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SimpleCompressible<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple compressible model");
        {
          reference_rho              = prm.get_double ("Reference density");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          reference_compressibility  = prm.get_double ("Reference compressibility");

          constant_rheology.parse_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::none;

      if (thermal_alpha != 0)
        this->model_dependence.density |= NonlinearDependence::temperature;
      if (reference_compressibility != 0)
        this->model_dependence.density |= NonlinearDependence::pressure;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SimpleCompressible,
                                   "simple compressible",
                                   "A material model that has constant values "
                                   "for all coefficients but the density. The defaults for all "
                                   "coefficients are chosen to be similar to what is believed to be correct "
                                   "for Earth's mantle. All of the values that define this model are read "
                                   "from a section ``Material model/Simple compressible model'' in the input file, see "
                                   "Section~\\ref{parameters:Material_20model/Simple_20compressible_20model}."
                                   "\n\n"
                                   "This model uses the following equations for the density: "
                                   "\\begin{align}"
                                   "  \\rho(p,T) = \\rho_0"
                                   "              \\left(1-\\alpha (T-T_a)\\right) "
                                   "              \\exp{\\beta (P-P_0))}"
                                   "\\end{align}"
                                   "This formulation for the density assumes that the compressibility "
                                   "provided by the user is the adiabatic compressibility ($\\beta_S$). "
                                   "The thermal expansivity and isentropic compressibility implied by "
                                   "the pressure and temperature dependence are equal to the "
                                   "user-defined constant values only along the reference isentrope, and "
                                   "there is also an implicit pressure dependence to the heat capacity "
                                   "$C_p$ via Maxwell's relations.")
  }
}
