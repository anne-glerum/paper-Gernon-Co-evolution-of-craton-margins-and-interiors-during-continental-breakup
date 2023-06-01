#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/material_model/simpler.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    class CellMaterial :
      public aspect::MaterialModel::Simpler<dim>, aspect::SimulatorAccess<dim>
    {
      public:
        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                 MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          Simpler<dim>::evaluate(in,out);

          if (in.current_cell.state() == IteratorState::valid)
            {
              std::cout << "Level: " << in.current_cell->level() << " Index: " << in.current_cell->index() << std::endl;
            }
        }
    };
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(CellMaterial,
                                   "cell",
                                   "")
  }
}
