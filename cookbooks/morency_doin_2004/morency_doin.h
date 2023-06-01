/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_morency_doin_h
#define _aspect_material_model_morency_doin_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model based on the rheology described in (Morency and Doin,
     * 2004): Brittle-ductile rheology with a viscosity strongly depending on
     * temperature and composition. Using pseudo-brittle rheology limits the
     * strength of the lithosphere at low temperature.
     *
     * The effective viscosity is defined as the harmonic mean of two stress-
     * dependent viscosity functions: a simple temperature/pressure-dependent,
     * non-Newtonian viscosity, and a more strongly stress-dependent,
     * "plastic" viscosity.
     *
     * @f[ v_{eff}^v = B \left(\frac{\dot{\varepsilon}}{\dot{\varepsilon}_{ref}}\right)^{-1+1/n_v}
     * exp\left(\frac{E_a +V_a \rho_m g z}{n_v R T}\right) @f]
     * @f[ v_{eff}^p = (\tau_0 + \gamma \rho_m g z) \left( \frac{\dot{\varepsilon}^{-1+1/n_p}}
     * {\dot{\varepsilon}_{ref}^{1/n_p}} \right) @f]
     * @f[ v_{eff} = \left(\frac{1}{v_{eff}^v}+\frac{1}{v_{eff}^p}\right)^{-1} @f]
     *
     * Where $v_{eff}$ is the effective viscosity, $B$ is a scaling constant,
     * $\dot{\varepsilon}$ is related to the second invariant of the strain
     * rate tensor, $\dot{\varepsilon}_{ref}$ is a reference strain rate,
     * $n_v$ and $n_p$ are stress exponents, $E_a$ is the activation energy,
     * $V_a$ is the activation volume, $\rho_m$ is the mantle density, $R$ is
     * the gas constant, $T$ is temperature, $\tau_0$ is the cohesive
     * strength of rocks at the surface, $\gamma$ is a coefficient of yield
     * stress increase with depth, and $z$ is depth.
     *
     * Several model parameters (reference densities, activation energies,
     * thermal expansivities, and stress exponents-- both viscous and plastic)
     * can be defined per-compositional field. If a list of values is given
     * for any of these parameters, the weighted sum of the values based on
     * volume fractions of the compositional fields is used in their place. If
     * only one value is given for any of these parameters, all compositions
     * are assigned the same value. The first value in the list is the value
     * assigned to "background mantle" (regions where the sum of the
     * compositional fields is < 1.0).
     *
     * For more on the material model and its applications, see: Morency, C.,
     * and M‐P. Doin. "Numerical simulations of the mantle lithosphere
     * delamination." Journal of Geophysical Research: Solid Earth
     * (1978–2012) 109.B3 (2004)
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MorencyDoin : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Evaluate material properties.
         */
        virtual void evaluate(const MaterialModelInputs<dim> &in,
                              MaterialModelOutputs<dim> &out) const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;

        virtual double reference_viscosity () const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        std::vector<double> densities;
        std::vector<double> activation_energies;
        std::vector<double> thermal_expansivities;
        std::vector<double> nvs; // Stress exponent, viscous rheology
        std::vector<double> nps;//Stress exponent, plastic rheology
        double thermal_diffusivity;
        double gamma; // Coefficient of yield stress increase with depth
        double heat_capacity;
        double activation_volume;
        double ref_strain_rate;
        double B; // Preexponential constant in the viscous rheology law B
        double tau_0; // cohesive strength of rocks at the surface
        double reference_T;
        double min_strain_rate;

        double ref_visc;
    };

  }
}

#endif
