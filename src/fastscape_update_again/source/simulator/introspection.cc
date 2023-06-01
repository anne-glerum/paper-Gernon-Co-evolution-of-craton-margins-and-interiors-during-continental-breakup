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


#include <aspect/introspection.h>
#include <aspect/global.h>


#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <tuple>

namespace aspect
{
  namespace internal
  {
    /**
     * Return a filled ComponentIndices structure.
     */
    template <int dim>
    typename Introspection<dim>::ComponentIndices
    setup_component_indices (FEVariableCollection<dim> &fevs)
    {
      typename Introspection<dim>::ComponentIndices ci;
      for (unsigned int i=0; i<dim; ++i)
        ci.velocities[i] = fevs.variable("velocity").first_component_index+i;

      ci.pressure = fevs.variable("pressure").first_component_index;
      ci.temperature = fevs.variable("temperature").first_component_index;
      unsigned int n_compositional_fields = fevs.variable("compositions").n_components();
      for (unsigned int i=0; i<n_compositional_fields; ++i)
        ci.compositional_fields.push_back(fevs.variable("compositions").first_component_index+i);

      return ci;
    }

    /**
     * return filled BlockIndices structure.
     */
    template <int dim>
    typename Introspection<dim>::BlockIndices
    setup_blocks (FEVariableCollection<dim> &fevs)
    {
      typename Introspection<dim>::BlockIndices b;
      b.velocities = fevs.variable("velocity").block_index;
      b.pressure = fevs.variable("pressure").block_index;
      b.temperature = fevs.variable("temperature").block_index;

      unsigned int n_compositional_fields = fevs.variable("compositions").n_components();
      for (unsigned int i=0; i<n_compositional_fields; ++i)
        b.compositional_fields.push_back(fevs.variable("compositions").block_index+i);

      return b;
    }



    template <int dim>
    typename Introspection<dim>::BaseElements
    setup_base_elements (FEVariableCollection<dim> &fevs)
    {
      typename Introspection<dim>::BaseElements base_elements;

      base_elements.velocities = fevs.variable("velocity").base_index;
      base_elements.pressure = fevs.variable("pressure").base_index;
      base_elements.temperature = fevs.variable("temperature").base_index;
      base_elements.compositional_fields = fevs.variable("compositions").base_index;

      return base_elements;
    }



    template <int dim>
    typename Introspection<dim>::PolynomialDegree
    setup_polynomial_degree (const Parameters<dim> &parameters)
    {
      typename Introspection<dim>::PolynomialDegree polynomial_degree;

      polynomial_degree.velocities = parameters.stokes_velocity_degree;
      polynomial_degree.temperature = parameters.temperature_degree;
      polynomial_degree.compositional_fields = parameters.composition_degree;

      return polynomial_degree;
    }



    template <int dim>
    std::shared_ptr<FiniteElement<dim>>
                                     new_FE_Q_or_DGP(const bool discontinuous,
                                                     const unsigned int degree)
    {
      if (discontinuous)
        return std::make_shared<FE_DGP<dim>>(degree);
      else
        return std::make_shared<FE_Q<dim>>(degree);
    }



    template <int dim>
    std::shared_ptr<FiniteElement<dim>>
                                     new_FE_Q_or_DGQ(const bool discontinuous,
                                                     const unsigned int degree)
    {
      if (discontinuous)
        return std::make_shared<FE_DGQ<dim>>(degree);
      else
        return std::make_shared<FE_Q<dim>>(degree);
    }
  }




  template <int dim>
  std::vector<VariableDeclaration<dim>>
                                     construct_default_variables (const Parameters<dim> &parameters)
  {
    std::vector<VariableDeclaration<dim>> variables;

    const unsigned int n_velocity_blocks = parameters.use_direct_stokes_solver ? 0 : 1;
    variables.push_back(
      VariableDeclaration<dim>("velocity",
                               std::make_shared<FE_Q<dim>>(parameters.stokes_velocity_degree),
                               dim,
                               n_velocity_blocks));

    if (parameters.use_equal_order_interpolation_for_stokes == false)
      variables.push_back(
        VariableDeclaration<dim>(
          "pressure",
          internal::new_FE_Q_or_DGP<dim>(parameters.use_locally_conservative_discretization,
                                         parameters.stokes_velocity_degree-1),
          1,
          1));
    else
      variables.push_back(
        VariableDeclaration<dim>(
          "pressure",
          std::shared_ptr<FiniteElement<dim>>(new FE_Q<dim>(parameters.stokes_velocity_degree)),
          1,
          1));


    variables.push_back(
      VariableDeclaration<dim>(
        "temperature",
        internal::new_FE_Q_or_DGQ<dim>(parameters.use_discontinuous_temperature_discretization,
                                       parameters.temperature_degree),
        1,
        1));

    variables.push_back(
      VariableDeclaration<dim>(
        "compositions",
        internal::new_FE_Q_or_DGQ<dim>(parameters.use_discontinuous_composition_discretization,
                                       parameters.composition_degree),
        parameters.n_compositional_fields,
        parameters.n_compositional_fields));

    return variables;
  }


  template <int dim>
  Introspection<dim>::Introspection(const std::vector<VariableDeclaration<dim>> &variable_definition,
                                    const Parameters<dim> &parameters
                                   )
    :
    FEVariableCollection<dim>(variable_definition),
    n_components (FEVariableCollection<dim>::n_components()),
    n_compositional_fields (parameters.n_compositional_fields),
    use_discontinuous_temperature_discretization (parameters.use_discontinuous_temperature_discretization),
    use_discontinuous_composition_discretization (parameters.use_discontinuous_composition_discretization),
    component_indices (internal::setup_component_indices<dim>(*this)),
    n_blocks(FEVariableCollection<dim>::n_blocks()),
    block_indices (internal::setup_blocks<dim>(*this)),
    extractors (component_indices),
    base_elements (internal::setup_base_elements<dim>(*this)),
    polynomial_degree (internal::setup_polynomial_degree<dim>(parameters)),
    component_masks (*this),
    system_dofs_per_block (n_blocks),
    compositional_field_methods(parameters.compositional_field_methods),
    composition_names(parameters.names_of_compositional_fields),
    composition_descriptions(parameters.composition_descriptions)
  {}



  template <int dim>
  Introspection<dim>::~Introspection ()
  {}



  namespace
  {
    std::vector<FEValuesExtractors::Scalar>
    make_extractor_sequence (const std::vector<unsigned int> &compositional_fields)
    {
      std::vector<FEValuesExtractors::Scalar> x;
      for (unsigned int i=0; i<compositional_fields.size(); ++i)
        x.emplace_back(compositional_fields[i]);
      return x;
    }
  }



  template <int dim>
  Introspection<dim>::Extractors::Extractors (const Introspection<dim>::ComponentIndices &component_indices)
    :
    velocities (component_indices.velocities[0]),
    pressure (component_indices.pressure),
    temperature (component_indices.temperature),
    compositional_fields (make_extractor_sequence (component_indices.compositional_fields))
  {}



  namespace
  {
    template <int dim>
    std::vector<ComponentMask>
    make_component_mask_sequence(const FEVariable<dim> &variable)
    {
      std::vector<ComponentMask> result;
      for (unsigned int i=0; i<variable.multiplicity; ++i)
        {
          result.push_back(ComponentMask(variable.component_mask.size(), false));
          result.back().set(variable.first_component_index+i, true);
        }
      return result;
    }
  }



  template <int dim>
  Introspection<dim>::ComponentMasks::ComponentMasks (FEVariableCollection<dim> &fevs)
    :
    velocities (fevs.variable("velocity").component_mask),
    pressure (fevs.variable("pressure").component_mask),
    temperature (fevs.variable("temperature").component_mask),
    compositional_fields (make_component_mask_sequence (fevs.variable("compositions")))
  {}



  template <int dim>
  unsigned int
  Introspection<dim>::compositional_index_for_name (const std::string &name) const
  {
    const std::vector<std::string>::const_iterator
    it = std::find(composition_names.begin(), composition_names.end(), name);
    AssertThrow (it != composition_names.end(),
                 ExcMessage ("The compositional field " + name +
                             " you asked for is not used in the simulation."));

    return (it - composition_names.begin());
  }



  template <int dim>
  std::string
  Introspection<dim>::name_for_compositional_index (const unsigned int index) const
  {
    // make sure that what we get here is really an index of one of the compositional fields
    AssertIndexRange(index,composition_names.size());
    return composition_names[index];
  }



  template <int dim>
  const std::vector<std::string> &
  Introspection<dim>::get_composition_names () const
  {
    // Simply return the full list of composition names
    return composition_names;
  }



  template <int dim>
  const std::vector<typename Parameters<dim>::CompositionalFieldDescription> &
  Introspection<dim>::get_composition_descriptions () const
  {
    return composition_descriptions;
  }



  template <int dim>
  bool
  Introspection<dim>::composition_type_exists (const typename Parameters<dim>::CompositionalFieldDescription::Type &type) const
  {
    for (unsigned int c=0; c<composition_descriptions.size(); ++c)
      if (composition_descriptions[c].type == type)
        return true;
    return false;
  }



  template <int dim>
  unsigned int
  Introspection<dim>::find_composition_type (const typename Parameters<dim>::CompositionalFieldDescription::Type &type) const
  {
    for (unsigned int c=0; c<composition_descriptions.size(); ++c)
      if (composition_descriptions[c].type == type)
        return c;
    return composition_descriptions.size();
  }



  template <int dim>
  bool
  Introspection<dim>::compositional_name_exists (const std::string &name) const
  {
    return (std::find(composition_names.begin(), composition_names.end(), name) != composition_names.end()
            ?
            true
            :
            false);
  }



  template <int dim>
  bool
  Introspection<dim>::is_stokes_component (const unsigned int component_index) const
  {
    if (component_index == component_indices.pressure)
      return true;

    for (unsigned int i=0; i<dim; ++i)
      if (component_index == component_indices.velocities[i])
        return true;

    return false;
  }


}


// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template struct Introspection<dim>; \
  template \
  std::vector<VariableDeclaration<dim>> \
  construct_default_variables (const Parameters<dim> &parameters);

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
