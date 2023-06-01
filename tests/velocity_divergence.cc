#include <aspect/material_model/interface.h>
#include <aspect/global.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  template <int dim>
  class TestMaterial:
    public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_viscosity () const
      {
        return 1e19;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            out.viscosities[i] = 1e19;
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1e-4;
            out.compressibilities[i] = 0.0;
            out.densities[i] = 3000;
            for (unsigned int c=0; c<in.composition[i].size(); ++c)
              {
                if (in.requests_property(MaterialModel::MaterialProperties::reaction_terms))
                  out.reaction_terms[i][c] = trace(in.strain_rate[i]) * this->get_timestep();
                else
                  out.reaction_terms[i][c] = 0.0;
              }

          }
      }
  };
}

// explicit instantiations
namespace aspect
{

  ASPECT_REGISTER_MATERIAL_MODEL(TestMaterial,
                                 "test material",
                                 "")

}
