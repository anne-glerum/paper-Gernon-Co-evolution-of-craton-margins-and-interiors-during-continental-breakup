/*
  Copyright (C) 2015 - 2019 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_drucker_prager_h
#define _aspect_material_model_drucker_prager_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/rheology/drucker_prager.h>
#include <aspect/material_model/equation_of_state/linearized_incompressible.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that consists of globally constant values for all
     * material parameters except density and viscosity.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible.
     *
     * The viscosity is computed according to the Drucker Prager frictional
     * plasticity criterion based on a user-defined internal angle of friction $\phi$
     * and cohesion $C$. In 3D:
     * $\sigma_y = \frac{6 C \cos(\phi)}{\sqrt{3} (3+\sin(\phi))} +
     * \frac{6 P \sin(\phi)}{\sqrt{3} (3+\sin(\phi))}$,
     * where $P$ is the pressure.
     * See for example Zienkiewicz, O. C., Humpheson, C. and Lewis, R. W. (1975),
     * G&eacute;otechnique 25, No. 4, 671-689.
     * With this formulation we circumscribe instead of inscribe the Mohr Coulomb
     * yield surface.
     * In 2D the Drucker Prager yield surface is the same
     * as the Mohr Coulomb surface:
     * $\sigma_y = P \sin(\phi) + C \cos(\phi)$.
     * Note that in 2D for $\phi=0$, these criteria
     * revert to the von Mises criterion (no pressure dependence).
     * See for example Thieulot, C. (2011), PEPI 188, 47-68.
     *
     * Note that we enforce the pressure to be positive in the computation of
     * the yield strength by replacing it with
     * a zero value whenever it is negative to prevent negative
     * yield strengths and viscosities.
     * We then use the computed yield strength to scale back the viscosity on
     * to the yield surface using the Viscosity Rescaling Method described in
     * Kachanov, L. M. (2004), Fundamentals of the Theory of Plasticity,
     * Dover Publications, Inc.
     *
     * To avoid numerically unfavourably large (or even negative) viscosity ranges,
     * we cut off the viscosity with a user-defined minimum and maximum viscosity:
     * $\eta_eff = \frac{1}{\frac{1}{\eta_min + \eta}+\\
     * \frac{1}{\eta_max}}$.
     *
     * Note that this model uses the formulation that assumes an incompressible
     * medium despite the fact that the density follows the law
     * $\rho(T)=\rho_0(1-\beta(T-T_{\text{ref}}))$.
     *
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class DruckerPrager : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        bool is_compressible () const override;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        double reference_viscosity () const override;
        /**
         * @}
         */

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

        double reference_T;
        double reference_eta;
        double thermal_conductivities;

        EquationOfState::LinearizedIncompressible<dim> equation_of_state;

        /*
         * Objects for computing plastic stresses, viscosities, and additional outputs
         */
        Rheology::DruckerPrager<dim> drucker_prager_plasticity;

        /**
         * The angle of internal friction
         */
        double angle_of_internal_friction;

        /**
         * The cohesion
         */
        double cohesion;

        /**
         * The applied viscosity bounds
         */
        double minimum_viscosity;
        double maximum_viscosity;

        /**
         * The reference strain rate used as a first estimate
         */
        double reference_strain_rate;
    };

  }
}

#endif
