/*
  Copyright (C) 2014 - 2021 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>
#include <aspect/geometry_model/box.h>
#include <aspect/simulator.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1_eulerian.h>
#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include <aspect/simulator_signals.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/simulator/assemblers/interface.h>
#include <aspect/melt.h>

#ifdef ASPECT_USE_PETSC
#  include <deal.II/lac/sparsity_tools.h>
#endif



namespace aspect
{

namespace Assemblers
  {
    template <int dim>
    ApplyStabilization<dim>::ApplyStabilization(const double stabilization_theta)
      :
      surface_theta(stabilization_theta)
    {}

    template <int dim>
    void
    ApplyStabilization<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>       &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim>      &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      AssertThrow(!this->get_mesh_deformation_handler().get_active_mesh_deformation_boundary_indicators().empty(),
                  ExcMessage("Applying free surface stabilization, even though no mesh deformation is active. "));

      if (this->get_parameters().include_melt_transport)
        {
          this->get_melt_handler().apply_free_surface_stabilization_with_melt (surface_theta,
                                                                               scratch.cell,
                                                                               scratch,
                                                                               data);
          return;
        }

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const typename DoFHandler<dim>::active_cell_iterator cell (&this->get_triangulation(),
                                                                 scratch.finite_element_values.get_cell()->level(),
                                                                 scratch.finite_element_values.get_cell()->index(),
                                                                 &this->get_dof_handler());

      const unsigned int n_face_q_points = scratch.face_finite_element_values.n_quadrature_points;
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();

      // Get the boundary indicators of those boundaries with a free surface
      const std::set<types::boundary_id> tmp_free_surface_boundary_indicators = this->get_mesh_deformation_handler().get_active_mesh_deformation_boundary_indicators();


      // only apply on free surface faces
      if (cell->at_boundary() && cell->is_locally_owned())
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
          if (cell->face(face_no)->at_boundary())
            {
              const types::boundary_id boundary_indicator
                = cell->face(face_no)->boundary_id();

              if (tmp_free_surface_boundary_indicators.find(boundary_indicator)
                  == tmp_free_surface_boundary_indicators.end())
                continue;

              scratch.face_finite_element_values.reinit(cell, face_no);

              this->compute_material_model_input_values (this->get_solution(),
                                                         scratch.face_finite_element_values,
                                                         cell,
                                                         false,
                                                         scratch.face_material_model_inputs);

              this->get_material_model().evaluate(scratch.face_material_model_inputs, scratch.face_material_model_outputs);

              for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point)
                {
                  for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
                    {
                      if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                        {
                          scratch.phi_u[i_stokes] = scratch.face_finite_element_values[introspection.extractors.velocities].value(i, q_point);
                          ++i_stokes;
                        }
                      ++i;
                    }

                  const Tensor<1,dim>
                  gravity = this->get_gravity_model().gravity_vector(scratch.face_finite_element_values.quadrature_point(q_point));
                  const double g_norm = gravity.norm();

                  // construct the relevant vectors
                  const Tensor<1,dim> n_hat = scratch.face_finite_element_values.normal_vector(q_point);
                  const Tensor<1,dim> g_hat = (g_norm == 0.0 ? Tensor<1,dim>() : gravity/g_norm);

                  const double pressure_perturbation = scratch.face_material_model_outputs.densities[q_point] *
                                                       this->get_timestep() *
                                                       surface_theta *
                                                       g_norm;

                  // see Kaus et al 2010 for details of the stabilization term
                  for (unsigned int i=0; i< stokes_dofs_per_cell; ++i)
                    for (unsigned int j=0; j< stokes_dofs_per_cell; ++j)
                      {
                        // The fictive stabilization stress is (phi_u[i].g)*(phi_u[j].n)
                        const double stress_value = -pressure_perturbation*
                                                    (scratch.phi_u[i]*g_hat) * (scratch.phi_u[j]*n_hat)
                                                    *scratch.face_finite_element_values.JxW(q_point);

                        data.local_matrix(i,j) += stress_value;
                      }
                }
            }
    }
  }

  namespace MeshDeformation
  {
    template <int dim>
    void
    Interface<dim>::initialize ()
    {}



    template <int dim>
    void
    Interface<dim>::update ()
    {}



    template <int dim>
    Tensor<1,dim>
    Interface<dim>::
    compute_initial_deformation_on_boundary(const types::boundary_id /*boundary_indicator*/,
                                            const Point<dim> &/*position*/) const
    {
      return Tensor<1,dim>();
    }



    template <int dim>
    void
    Interface<dim>::
    compute_velocity_constraints_on_boundary(const DoFHandler<dim> &/*mesh_deformation_dof_handler*/,
                                             AffineConstraints<double> &/*mesh_velocity_constraints*/,
                                             const std::set<types::boundary_id> &/*boundary_id*/) const
    {}



    template <int dim>
    void
    Interface<dim>::
    declare_parameters (ParameterHandler &)
    {}



    template <int dim>
    void
    Interface<dim>::parse_parameters (ParameterHandler &)
    {}



    template <int dim>
    MeshDeformationHandler<dim>::MeshDeformationHandler (Simulator<dim> &simulator)
      : sim(simulator),  // reference to the simulator that owns the MeshDeformationHandler
        mesh_deformation_fe (FE_Q<dim>(1),dim), // Q1 elements which describe the mesh geometry
        mesh_deformation_dof_handler (sim.triangulation),
        include_initial_topography(false)
    {
      // Now reset the mapping of the simulator to be something that captures mesh deformation in time.
      sim.mapping.reset (new MappingQ1Eulerian<dim, LinearAlgebra::Vector> (mesh_deformation_dof_handler,
                                                                            mesh_displacements));
    }



    template <int dim>
    MeshDeformationHandler<dim>::~MeshDeformationHandler ()
    {
      // Free the Simulator's mapping object, otherwise
      // when the MeshDeformationHandler gets destroyed,
      // the mapping's reference to the mesh displacement
      // vector will be invalid.
      sim.mapping.reset();
    }



    namespace
    {
      std::tuple
      <void *,
      void *,
      aspect::internal::Plugins::PluginList<Interface<2>>,
      aspect::internal::Plugins::PluginList<Interface<3>>> registered_plugins;
    }



    template <int dim>
    void
    MeshDeformationHandler<dim>::initialize ()
    {
      // In case we prescribed initial topography, we should take this into
      // account. However, it is not included in the mesh displacements,
      // so we need to fetch it separately.
      if (!Plugins::plugin_type_matches<InitialTopographyModel::ZeroTopography<dim>>(this->get_initial_topography_model()))
        include_initial_topography = true;


      this->get_signals().set_assemblers.connect(std::bind(&MeshDeformationHandler<dim>::set_assemblers,
                                                           std::cref(*this),
                                                           std::placeholders::_1,
                                                           std::placeholders::_2));

    }

    template <int dim>
    void MeshDeformationHandler<dim>::set_assemblers(const SimulatorAccess<dim> &,
                                          aspect::Assemblers::Manager<dim> &assemblers) const
    {
      aspect::Assemblers::ApplyStabilization<dim> *surface_stabilization
        = new aspect::Assemblers::ApplyStabilization<dim>(surface_theta);

      assemblers.stokes_system.push_back(
        std::unique_ptr<aspect::Assemblers::ApplyStabilization<dim> > (surface_stabilization));

      // Note that we do not want face_material_model_data, because we do not
      // connect to a face assembler. We instead connect to a normal assembler,
      // and compute our own material_model_inputs in apply_stabilization
      // (because we want to use the solution instead of the current_linearization_point
      // to compute the material properties).
      assemblers.stokes_system_assembler_on_boundary_face_properties.needed_update_flags |= (update_values  |
          update_gradients |
          update_quadrature_points |
          update_normal_vectors |
          update_JxW_values);
    }



    template <int dim>
    void
    MeshDeformationHandler<dim>::register_mesh_deformation (const std::string &name,
                                                            const std::string &description,
                                                            void (*declare_parameters_function) (ParameterHandler &),
                                                            Interface<dim> *(*factory_function) ())
    {
      std::get<dim>(registered_plugins).register_plugin (name,
                                                         description,
                                                         declare_parameters_function,
                                                         factory_function);
    }



    template <int dim>
    void
    MeshDeformationHandler<dim>::update ()
    {
      AssertThrow(sim.parameters.mesh_deformation_enabled, ExcInternalError());

      for (const auto &boundary_and_deformation_objects : mesh_deformation_objects)
        {
          for (const auto &model : boundary_and_deformation_objects.second)
            model->update();
        }
    }


    template <int dim>
    void MeshDeformationHandler<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry ("Additional tangential mesh velocity boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "where there the mesh is allowed to move tangential to the "
                           "boundary. All tangential mesh movements along "
                           "those boundaries that have tangential material velocity "
                           "boundary conditions are allowed by default, this parameters "
                           "allows to generate mesh movements along other boundaries that are "
                           "open, or have prescribed material velocities or tractions."
                           "\n\n"
                           "The names of the boundaries listed here can either be "
                           "numbers (in which case they correspond to the numerical "
                           "boundary indicators assigned by the geometry object), or they "
                           "can correspond to any of the symbolic names the geometry object "
                           "may have provided for each part of the boundary. You may want "
                           "to compare this with the documentation of the geometry model you "
                           "use in your model.");
        prm.declare_entry ("Mesh deformation boundary indicators", "",
                           Patterns::List (Patterns::Anything()),
                           "A comma separated list of names denoting those boundaries "
                           "where there the mesh is allowed to move according to the "
                           "specified mesh deformation objects. "
                           "\n\n"
                           "The names of the boundaries listed here can either be "
                           "numbers (in which case they correspond to the numerical "
                           "boundary indicators assigned by the geometry object), or they "
                           "can correspond to any of the symbolic names the geometry object "
                           "may have provided for each part of the boundary. You may want "
                           "to compare this with the documentation of the geometry model you "
                           "use in your model. "
                           "\n\n"
                           "The format is id1: object1 \\& object2, id2: object3 \\& object2, where "
                           "objects are one of " + std::get<dim>(registered_plugins).get_description_string());
          prm.declare_entry("Surface stabilization theta", "0.5",
                            Patterns::Double(0., 1.),
                            "Theta parameter described in \\cite{KMM2010}. "
                            "An unstabilized free surface can overshoot its "
                            "equilibrium position quite easily and generate "
                            "unphysical results.  One solution is to use a "
                            "quasi-implicit correction term to the forces near the "
                            "free surface.  This parameter describes how much "
                            "the free surface is stabilized with this term, "
                            "where zero is no stabilization, and one is fully "
                            "implicit.");

      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }

    template <int dim>
    void MeshDeformationHandler<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        surface_theta = prm.get_double("Surface stabilization theta");

        // Create the map of prescribed mesh movement boundary indicators
        // Each boundary indicator can carry a number of mesh deformation plugin names.
        const std::vector<std::string> x_mesh_deformation_boundary_indicators
          = Utilities::split_string_list(prm.get("Mesh deformation boundary indicators"),",");

        for (const auto &entry : x_mesh_deformation_boundary_indicators)
          {
            // each entry has the format (white space is optional):
            // <boundary_id> : <object_name & object_name, ...>
            const std::vector<std::string> split_parts = Utilities::split_string_list (entry, ':');
            AssertThrow (split_parts.size() == 2,
                         ExcMessage ("The format for mesh deformation indicators "
                                     "requires that each entry has the form `"
                                     "<id> : <value & value & ...>', but there does not "
                                     "appear to be a colon in the entry <"
                                     + entry
                                     + ">."));

            // Get the values, i.e. the mesh deformation plugin names
            const std::vector<std::string> object_names = Utilities::split_string_list(split_parts[1],"&");

            // Try to translate the id into a boundary_id.
            types::boundary_id boundary_id;
            try
              {
                boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id(split_parts[0]);
              }
            catch (const std::string &error)
              {
                AssertThrow (false, ExcMessage ("While parsing the entry <Mesh deformation/"
                                                "Mesh deformation boundary indicators>, there was an error. Specifically, "
                                                "the conversion function complained as follows:\n\n"
                                                + error));
              }

            // Store the boundary indicator. If the entry exists this does nothing.
            prescribed_mesh_deformation_boundary_indicators.insert(boundary_id);

            for (const auto &object_name : object_names)
              {
                // Make sure there are no duplicated entries. If this boundary is not
                // already in the map the first call to map[key] will create an empty entry.
                AssertThrow(std::find(mesh_deformation_object_names[boundary_id].begin(),
                                      mesh_deformation_object_names[boundary_id].end(), object_name)
                            == mesh_deformation_object_names[boundary_id].end(),
                            ExcMessage("The current mesh deformation object is listed twice for boundary indicator "
                                       + dealii::Utilities::int_to_string(boundary_id)));

                mesh_deformation_object_names[boundary_id].push_back(object_name);

                if (object_name == "free surface")
                  free_surface_boundary_indicators.insert(boundary_id);
              }
          }

        // Create the list of tangential mesh movement boundary indicators
        try
          {
            const std::vector<types::boundary_id> x_additional_tangential_mesh_boundary_indicators
              = this->get_geometry_model().translate_symbolic_boundary_names_to_ids(Utilities::split_string_list
                                                                                    (prm.get ("Additional tangential mesh velocity boundary indicators")));

            tangential_mesh_deformation_boundary_indicators.insert(x_additional_tangential_mesh_boundary_indicators.begin(),
                                                                   x_additional_tangential_mesh_boundary_indicators.end());
          }
        catch (const std::string &error)
          {
            AssertThrow (false, ExcMessage ("While parsing the entry <Mesh deformation/Additional tangential "
                                            "mesh velocity boundary indicators>, there was an error. Specifically, "
                                            "the conversion function complained as follows:\n\n"
                                            + error));
          }

        // Boundaries with tangential Stokes velocity boundary conditions are implicitly
        // treated as tangential mesh boundaries, but only if they do not have
        // assigned mesh deformation objects.
        for (const auto &boundary_id : this->get_boundary_velocity_manager().get_tangential_boundary_velocity_indicators())
          tangential_mesh_deformation_boundary_indicators.insert(boundary_id);

        // The tangential mesh boundaries can accidentally contain prescribed mesh
        // boundaries if they were in the list of tangential Stokes boundaries.
        // If so remove them.
        for (const auto &boundary_id : prescribed_mesh_deformation_boundary_indicators)
          tangential_mesh_deformation_boundary_indicators.erase(boundary_id);

        // All periodic boundaries are implicitly treated as tangential mesh deformation boundaries.
        using periodic_boundary_pair = std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int>;
        for (const periodic_boundary_pair &p : this->get_geometry_model().get_periodic_boundary_pairs())
          {
            tangential_mesh_deformation_boundary_indicators.insert(p.first.first);
            tangential_mesh_deformation_boundary_indicators.insert(p.first.second);
          }

        // Create the list of zero mesh movement (fixed) boundary indicators, these are
        // all boundaries, which are not prescribed or tangential.
        zero_mesh_deformation_boundary_indicators = this->get_geometry_model().get_used_boundary_indicators();

        for (const auto &boundary_id : prescribed_mesh_deformation_boundary_indicators)
          zero_mesh_deformation_boundary_indicators.erase(boundary_id);

        for (const auto &boundary_id : tangential_mesh_deformation_boundary_indicators)
          zero_mesh_deformation_boundary_indicators.erase(boundary_id);
      }
      prm.leave_subsection ();

      // go through the list of object names, create objects and let them parse
      // their own parameters
      for (const auto &boundary_and_object_names : mesh_deformation_object_names)
        {
          for (const auto &object_name : boundary_and_object_names.second)
            {
              mesh_deformation_objects[boundary_and_object_names.first].push_back(
                std::unique_ptr<Interface<dim>> (std::get<dim>(registered_plugins)
                                                 .create_plugin (object_name,
                                                                 "Mesh deformation::Model names")));

              if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(mesh_deformation_objects[boundary_and_object_names.first].back().get()))
                sim->initialize_simulator (this->get_simulator());

              mesh_deformation_objects[boundary_and_object_names.first].back()->parse_parameters (prm);
              mesh_deformation_objects[boundary_and_object_names.first].back()->initialize ();
            }
        }
    }



    template <int dim>
    void MeshDeformationHandler<dim>::execute()
    {
      AssertThrow(sim.parameters.mesh_deformation_enabled, ExcInternalError());

      TimerOutput::Scope timer (sim.computing_timer, "Mesh deformation");

      old_mesh_displacements = mesh_displacements;

      // Make the constraints for the elliptic problem.
      make_constraints();

      // Assemble and solve the vector Laplace problem which determines
      // the mesh displacements in the interior of the domain
      compute_mesh_displacements();

      // Interpolate the mesh velocity into the same
      // finite element space as used in the Stokes solve, which
      // is needed for the ALE corrections.
      interpolate_mesh_velocity();

      // After changing the mesh we need to rebuild things
      sim.rebuild_stokes_matrix = sim.rebuild_stokes_preconditioner = true;
    }



    template <int dim>
    const Mapping<dim> &
    MeshDeformationHandler<dim>::get_level_mapping(const unsigned int level) const
    {
      return *level_mappings[level].get();
    }


    template <int dim>
    double
    MeshDeformationHandler<dim>::get_surface_theta()const
    {
      return surface_theta;
    }




    template <int dim>
    void MeshDeformationHandler<dim>::make_constraints()
    {
      AssertThrow(sim.parameters.mesh_deformation_enabled, ExcInternalError());

      // Now construct the mesh displacement constraints
      mesh_velocity_constraints.clear();
      mesh_velocity_constraints.reinit(mesh_locally_relevant);

      // mesh_velocity_constraints can use the same hanging node
      // information that was used for mesh_vertex constraints.
      mesh_velocity_constraints.merge(mesh_vertex_constraints);

      // Add the vanilla periodic boundary constraints
      using periodic_boundary_pairs = std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int>>;
      const periodic_boundary_pairs pbp = this->get_geometry_model().get_periodic_boundary_pairs();
      for (const auto &p : pbp)
        DoFTools::make_periodicity_constraints(mesh_deformation_dof_handler,
                                               p.first.first,
                                               p.first.second,
                                               p.second,
                                               mesh_velocity_constraints);

      // Zero out the displacement for the zero-velocity boundaries
      // if the boundary is not in the set of tangential mesh boundaries and not in the set of mesh deformation boundary indicators
      for (const auto &boundary_id : zero_mesh_deformation_boundary_indicators)
        {
          VectorTools::interpolate_boundary_values (this->get_mapping(),
                                                    mesh_deformation_dof_handler,
                                                    boundary_id,
                                                    Functions::ZeroFunction<dim>(dim),
                                                    mesh_velocity_constraints);
        }

      this->get_signals().pre_compute_no_normal_flux_constraints(sim.triangulation);
      // Make the no flux boundary constraints
      VectorTools::compute_no_normal_flux_constraints (mesh_deformation_dof_handler,
                                                       /* first_vector_component= */
                                                       0,
                                                       tangential_mesh_deformation_boundary_indicators,
                                                       mesh_velocity_constraints,
                                                       this->get_mapping());

      this->get_signals().post_compute_no_normal_flux_constraints(sim.triangulation);

      // Ask all plugins to add their constraints.
      // For the moment add constraints from all plugins into one matrix, then
      // merge that matrix with the existing constraints (respecting the existing
      // constraints as more important)
      AffineConstraints<double> plugin_constraints(mesh_vertex_constraints.get_local_lines());

      for (const auto &boundary_id : mesh_deformation_objects)
        {
          std::set<types::boundary_id> boundary_id_set;
          boundary_id_set.insert(boundary_id.first);

          for (const auto &model : boundary_id.second)
            {
              AffineConstraints<double> current_plugin_constraints(mesh_vertex_constraints.get_local_lines());

              model->compute_velocity_constraints_on_boundary(mesh_deformation_dof_handler,
                                                              current_plugin_constraints,
                                                              boundary_id_set);

              const IndexSet local_lines = current_plugin_constraints.get_local_lines();
              for (dealii::IndexSet::size_type local_line : local_lines)
                {
                  if (current_plugin_constraints.is_constrained(local_line))
                    {
                      if (plugin_constraints.is_constrained(local_line) == false)
                        {
                          plugin_constraints.add_line(local_line);
                          plugin_constraints.set_inhomogeneity(local_line, current_plugin_constraints.get_inhomogeneity(local_line));
                        }
                      else
                        {
                          // Add the inhomogeneity of the current plugin to the existing constraints
                          const double inhomogeneity = plugin_constraints.get_inhomogeneity(local_line);
                          plugin_constraints.set_inhomogeneity(local_line, current_plugin_constraints.get_inhomogeneity(local_line) + inhomogeneity);
                        }
                    }
                }
            }
        }

      mesh_velocity_constraints.merge(plugin_constraints,AffineConstraints<double>::left_object_wins);
      mesh_velocity_constraints.close();
    }



    template <int dim>
    void MeshDeformationHandler<dim>::make_initial_constraints()
    {
      AssertThrow(this->get_parameters().mesh_deformation_enabled, ExcInternalError());

      // This might look incorrect at first glance, but it is okay to
      // overwrite the velocity constraints with our displacements
      // because this object is used for updating the displacement in
      // compute_mesh_displacements().
      mesh_velocity_constraints.clear();
      mesh_velocity_constraints.reinit(mesh_locally_relevant);

      // mesh_velocity_constraints can use the same hanging node
      // information that was used for mesh_vertex constraints.
      mesh_velocity_constraints.merge(mesh_vertex_constraints);

      // Add the vanilla periodic boundary constraints
      std::set< types::boundary_id > periodic_boundaries;
      using periodic_boundary_pairs = std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int>>;
      const periodic_boundary_pairs pbp = this->get_geometry_model().get_periodic_boundary_pairs();
      for (const auto &p : pbp)
        {
          periodic_boundaries.insert(p.first.first);
          periodic_boundaries.insert(p.first.second);

          DoFTools::make_periodicity_constraints(mesh_deformation_dof_handler,
                                                 p.first.first,
                                                 p.first.second,
                                                 p.second,
                                                 mesh_velocity_constraints);
        }

      // Zero out the displacement for the fixed boundaries
      for (const types::boundary_id &boundary_id : zero_mesh_deformation_boundary_indicators)
        {
          VectorTools::interpolate_boundary_values (this->get_mapping(),
                                                    mesh_deformation_dof_handler,
                                                    boundary_id,
                                                    Functions::ZeroFunction<dim>(dim),
                                                    mesh_velocity_constraints);
        }

      // Make tangential deformation constraints for tangential boundaries
      this->get_signals().pre_compute_no_normal_flux_constraints(sim.triangulation);
      VectorTools::compute_no_normal_flux_constraints (mesh_deformation_dof_handler,
                                                       /* first_vector_component= */
                                                       0,
                                                       tangential_mesh_deformation_boundary_indicators,
                                                       mesh_velocity_constraints,
                                                       this->get_mapping());
      this->get_signals().post_compute_no_normal_flux_constraints(sim.triangulation);

      // Ask all plugins to add their constraints.
      // For the moment add constraints from all plugins into one matrix, then
      // merge that matrix with the existing constraints (respecting the existing
      // constraints as more important)
      AffineConstraints<double> plugin_constraints(mesh_vertex_constraints.get_local_lines());

      for (const auto &boundary_id_and_deformation_objects: mesh_deformation_objects)
        {
          for (const auto &deformation_object : boundary_id_and_deformation_objects.second)
            {
              AffineConstraints<double> current_plugin_constraints(mesh_vertex_constraints.get_local_lines());

              Utilities::VectorFunctionFromVelocityFunctionObject<dim> vel
              (dim,
               [&] (const Point<dim> &x) -> Tensor<1,dim>
              {
                return deformation_object->compute_initial_deformation_on_boundary(boundary_id_and_deformation_objects.first, x);
              });

              VectorTools::interpolate_boundary_values (this->get_mapping(),
                                                        mesh_deformation_dof_handler,
                                                        boundary_id_and_deformation_objects.first,
                                                        vel,
                                                        current_plugin_constraints);

              const IndexSet local_lines = current_plugin_constraints.get_local_lines();
              for (dealii::IndexSet::size_type local_line : local_lines)
                {
                  if (current_plugin_constraints.is_constrained(local_line))
                    {
                      if (plugin_constraints.is_constrained(local_line) == false)
                        {
                          plugin_constraints.add_line(local_line);
                          plugin_constraints.set_inhomogeneity(local_line, current_plugin_constraints.get_inhomogeneity(local_line));
                        }
                      else
                        {
                          // Add the current plugin constraints to the existing inhomogeneity
                          const double inhomogeneity = plugin_constraints.get_inhomogeneity(local_line);
                          plugin_constraints.set_inhomogeneity(local_line, current_plugin_constraints.get_inhomogeneity(local_line) + inhomogeneity);
                        }
                    }
                }
            }
        }

      mesh_velocity_constraints.merge(plugin_constraints,
                                      AffineConstraints<double>::left_object_wins);
      mesh_velocity_constraints.close();
    }



    template <int dim>
    void MeshDeformationHandler<dim>::compute_mesh_displacements()
    {
      // This functions updates the mesh displacement of the whole
      // domain (stored in the vector mesh_displacements) based on
      // information on the boundary.
      //
      // Each step, we get the velocity specified on the free surface
      // boundary (stored in mesh_velocity_constraints) and solve for
      // the velocity in the interior by solving a vector Laplace
      // problem. This velocity is then used to update the
      // displacement vector.
      //
      // This is different in timestep 0. Here, the information on the
      // boundary is actually a displacement (given initial
      // topography), which is used to set the initial
      // displacement. The process in this function is otherwise
      // identical.

      QGauss<dim> quadrature(mesh_deformation_fe.degree + 1);
      UpdateFlags update_flags = UpdateFlags(update_values | update_JxW_values | update_gradients);
      FEValues<dim> fe_values (*sim.mapping, mesh_deformation_fe, quadrature, update_flags);

      const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                         dofs_per_face = sim.finite_element.dofs_per_face,
                         n_q_points    = fe_values.n_quadrature_points;

      std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
      std::vector<unsigned int> face_dof_indices (dofs_per_face);
      Vector<double> cell_vector (dofs_per_cell);
      FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

      // We are just solving a Laplacian in each spatial direction, so
      // the degrees of freedom for different dimensions do not couple.
      Table<2,DoFTools::Coupling> coupling (dim, dim);
      coupling.fill(DoFTools::none);

      for (unsigned int c=0; c<dim; ++c)
        coupling[c][c] = DoFTools::always;

      LinearAlgebra::SparseMatrix mesh_matrix;
      TrilinosWrappers::SparsityPattern sp (mesh_locally_owned,
                                            mesh_locally_owned,
                                            mesh_locally_relevant,
                                            sim.mpi_communicator);
      DoFTools::make_sparsity_pattern (mesh_deformation_dof_handler,
                                       coupling, sp,
                                       mesh_velocity_constraints, false,
                                       Utilities::MPI::
                                       this_mpi_process(sim.mpi_communicator));
      sp.compress();
      mesh_matrix.reinit (sp);

      // carry out the solution
      FEValuesExtractors::Vector extract_vel(0);

      LinearAlgebra::Vector rhs, solution;
      rhs.reinit(mesh_locally_owned, sim.mpi_communicator);
      solution.reinit(mesh_locally_owned, sim.mpi_communicator);

      typename DoFHandler<dim>::active_cell_iterator cell = mesh_deformation_dof_handler.begin_active(),
                                                     endc= mesh_deformation_dof_handler.end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            cell->get_dof_indices (cell_dof_indices);
            fe_values.reinit (cell);

            cell_vector = 0;
            cell_matrix = 0;
            for (unsigned int point=0; point<n_q_points; ++point)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    cell_matrix(i,j) += scalar_product( fe_values[extract_vel].gradient(i,point),
                                                        fe_values[extract_vel].gradient(j,point) ) *
                                        fe_values.JxW(point);
                }

            mesh_velocity_constraints.distribute_local_to_global (cell_matrix, cell_vector,
                                                                  cell_dof_indices, mesh_matrix, rhs, false);
          }

      rhs.compress (VectorOperation::add);
      mesh_matrix.compress (VectorOperation::add);

      // Make the AMG preconditioner
      std::vector<std::vector<bool>> constant_modes;
      DoFTools::extract_constant_modes (mesh_deformation_dof_handler,
                                        ComponentMask(dim, true),
                                        constant_modes);
      // TODO: think about keeping object between time steps
      LinearAlgebra::PreconditionAMG preconditioner_stiffness;
      LinearAlgebra::PreconditionAMG::AdditionalData Amg_data;
      Amg_data.constant_modes = constant_modes;
      Amg_data.elliptic = true;
      Amg_data.higher_order_elements = false;
      Amg_data.smoother_sweeps = 2;
      Amg_data.aggregation_threshold = 0.02;
      preconditioner_stiffness.initialize(mesh_matrix);

      // we solve with higher accuracy in the initial timestep:
      const double tolerance
        = sim.parameters.linear_stokes_solver_tolerance
          * ((this->simulator_is_past_initialization()) ? 1.0 : 1e-5);

      SolverControl solver_control(5*rhs.size(), tolerance * rhs.l2_norm());
      SolverCG<LinearAlgebra::Vector> cg(solver_control);

      cg.solve (mesh_matrix, solution, rhs, preconditioner_stiffness);
      this->get_pcout() << "   Solving mesh displacement system... " << solver_control.last_step() <<" iterations."<< std::endl;

      mesh_velocity_constraints.distribute (solution);

      // Update the mesh velocity vector
      fs_mesh_velocity = solution;

      // Update the mesh displacement vector
      if (this->simulator_is_past_initialization())
        {
          // during the simulation, we add dt*solution
          LinearAlgebra::Vector distributed_mesh_displacements(mesh_locally_owned, sim.mpi_communicator);
          distributed_mesh_displacements = mesh_displacements;
          distributed_mesh_displacements.add(this->get_timestep(), solution);
          mesh_displacements = distributed_mesh_displacements;
        }
      else
        {
          // In the initial step we apply 100% of the initial displacement
          mesh_displacements = solution;
        }

      if (this->is_stokes_matrix_free())
        update_multilevel_deformation();
    }



    template <int dim>
    void MeshDeformationHandler<dim>::set_initial_topography()
    {
      LinearAlgebra::Vector distributed_initial_topography;
      distributed_initial_topography.reinit(mesh_locally_owned, sim.mpi_communicator);

      if (!include_initial_topography)
        distributed_initial_topography = 0.;
      else
        {
          const std::vector<Point<dim>> support_points
                                     = mesh_deformation_fe.base_element(0).get_unit_support_points();

          const Quadrature<dim> quad(support_points);
          const UpdateFlags update_flags = UpdateFlags(update_quadrature_points);
          FEValues<dim> fs_fe_values (*sim.mapping, mesh_deformation_fe, quad, update_flags);

          const unsigned int n_q_points = fs_fe_values.n_quadrature_points,
                             dofs_per_cell = fs_fe_values.dofs_per_cell;

          std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);

          typename DoFHandler<dim>::active_cell_iterator cell = mesh_deformation_dof_handler.begin_active(),
                                                         endc = mesh_deformation_dof_handler.end();

          for (; cell!=endc; ++cell)
            if (cell->is_locally_owned())
              {
                cell->get_dof_indices (cell_dof_indices);

                fs_fe_values.reinit (cell);
                for (unsigned int j=0; j<n_q_points; ++j)
                  {
                    Point<dim-1> surface_point;
                    std::array<double, dim> natural_coord = this->get_geometry_model().cartesian_to_natural_coordinates(fs_fe_values.quadrature_point(j));
                    if (const GeometryModel::Box<dim> *geometry = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model()))
                      {
                        for (unsigned int d=0; d<dim-1; ++d)
                          surface_point[d] = natural_coord[d];
                      }
                    else
                      {
                        for (unsigned int d=1; d<dim; ++d)
                          surface_point[d-1] = natural_coord[d];
                      }
                    // Get the topography at this point.
                    const double topo = this->get_initial_topography_model().value(surface_point);


                    // TODO adapt to radial topography
                    const unsigned int support_point_index
                      = mesh_deformation_fe.component_to_system_index(dim-1,/*dof index within component=*/ j);
                    distributed_initial_topography[cell_dof_indices[support_point_index]] = topo;
                  }
              }
        }

      distributed_initial_topography.compress(VectorOperation::insert);
      initial_topography = distributed_initial_topography;
    }


    template <int dim>
    void MeshDeformationHandler<dim>::interpolate_mesh_velocity()
    {
      // Interpolate the mesh vertex velocity onto the Stokes velocity system for use in ALE corrections
      LinearAlgebra::BlockVector distributed_mesh_velocity;
      distributed_mesh_velocity.reinit(sim.introspection.index_sets.system_partitioning, sim.mpi_communicator);

      const std::vector<Point<dim>> support_points
                                 = sim.finite_element.base_element(sim.introspection.component_indices.velocities[0]).get_unit_support_points();

      const Quadrature<dim> quad(support_points);
      const UpdateFlags update_flags = UpdateFlags(update_values | update_JxW_values);
      FEValues<dim> fs_fe_values (*sim.mapping, mesh_deformation_fe, quad, update_flags);
      FEValues<dim> fe_values (*sim.mapping, sim.finite_element, quad, update_flags);
      const unsigned int n_q_points = fe_values.n_quadrature_points,
                         dofs_per_cell = fe_values.dofs_per_cell;

      std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
      FEValuesExtractors::Vector extract_vel(0);
      std::vector<Tensor<1,dim>> velocity_values(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator
      cell = sim.dof_handler.begin_active(), endc= sim.dof_handler.end();
      typename DoFHandler<dim>::active_cell_iterator
      fscell = mesh_deformation_dof_handler.begin_active();

      for (; cell!=endc; ++cell, ++fscell)
        if (cell->is_locally_owned())
          {
            cell->get_dof_indices (cell_dof_indices);

            fe_values.reinit (cell);
            fs_fe_values.reinit (fscell);
            fs_fe_values[extract_vel].get_function_values(fs_mesh_velocity, velocity_values);
            for (unsigned int j=0; j<n_q_points; ++j)
              for (unsigned int dir=0; dir<dim; ++dir)
                {
                  const unsigned int support_point_index
                    = sim.finite_element.component_to_system_index(/*velocity component=*/ sim.introspection.component_indices.velocities[dir],
                                                                                           /*dof index within component=*/ j);
                  distributed_mesh_velocity[cell_dof_indices[support_point_index]] = velocity_values[j][dir];
                }
          }

      distributed_mesh_velocity.compress(VectorOperation::insert);
      mesh_velocity = distributed_mesh_velocity;
    }


    template <int dim>
    void MeshDeformationHandler<dim>::setup_dofs()
    {
      AssertThrow(sim.parameters.mesh_deformation_enabled, ExcInternalError());

      // these live in the same FE as the velocity variable:
      mesh_velocity.reinit(sim.introspection.index_sets.system_partitioning,
                           sim.introspection.index_sets.system_relevant_partitioning,
                           sim.mpi_communicator);

      mesh_deformation_dof_handler.distribute_dofs(mesh_deformation_fe);

      if (this->is_stokes_matrix_free())
        {
          mesh_deformation_dof_handler.distribute_mg_dofs();

          const unsigned int n_levels = this->get_triangulation().n_global_levels();

          level_displacements.resize(0, n_levels-1);
          // Important! Preallocate level vectors with all needed ghost
          // entries. While interpolate_to_mg can create these vectors
          // automatically, they will not contain all ghost values that we
          // need to evaluate the mapping later.
          for (unsigned int level = 0; level < n_levels; ++level)
            {
              IndexSet relevant_mg_dofs;
              DoFTools::extract_locally_relevant_level_dofs(mesh_deformation_dof_handler,
                                                            level,
                                                            relevant_mg_dofs);
              level_displacements[level].reinit(mesh_deformation_dof_handler.locally_owned_mg_dofs(level),
                                                relevant_mg_dofs,
                                                sim.mpi_communicator);
              level_displacements[level].update_ghost_values();
            }

          // create the mappings on each level:
          level_mappings.resize(0, n_levels-1);
          level_mappings.apply([&](const unsigned int level, std::unique_ptr<Mapping<dim>> &object)
          {
            object = std::make_unique<MappingQEulerian<dim,
            dealii::LinearAlgebra::distributed::Vector<double>>>(
              /* degree = */ 1,
              mesh_deformation_dof_handler,
              level_displacements[level],
              level);
          });

          mg_transfer.build(mesh_deformation_dof_handler);

        }

      this->get_pcout() << "Number of mesh deformation degrees of freedom: "
                        << mesh_deformation_dof_handler.n_dofs()
                        << std::endl;

      // Renumber the DoFs hierarchical so that we get the
      // same numbering if we resume the computation. This
      // is because the numbering depends on the order the
      // cells are created.
      DoFRenumbering::hierarchical (mesh_deformation_dof_handler);

      mesh_locally_owned = mesh_deformation_dof_handler.locally_owned_dofs();
      DoFTools::extract_locally_relevant_dofs (mesh_deformation_dof_handler,
                                               mesh_locally_relevant);

      // This will initialize the mesh displacement and free surface
      // mesh velocity vectors with zero-valued entries.
      mesh_displacements.reinit(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);
      old_mesh_displacements.reinit(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);
      initial_topography.reinit(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);
      fs_mesh_velocity.reinit(mesh_locally_owned, mesh_locally_relevant, sim.mpi_communicator);

      // if we are just starting, we need to set the initial topography.
      if (this->simulator_is_past_initialization() == false ||
          this->get_timestep_number() == 0)
        set_initial_topography();

      // We would like to make sure that the mesh stays conforming upon
      // redistribution, so we construct mesh_vertex_constraints, which
      // keeps track of hanging node constraints.
      // Note: this would be a more natural fit in make_constraints(),
      // but we would like to be able to apply vertex constraints directly
      // after setup_dofs(), as is done, for instance, during mesh
      // refinement.
      mesh_vertex_constraints.clear();
      mesh_vertex_constraints.reinit(mesh_locally_relevant);

      DoFTools::make_hanging_node_constraints(mesh_deformation_dof_handler, mesh_vertex_constraints);

      // We can safely close this now
      mesh_vertex_constraints.close();

      // if we are just starting, we need to prescribe the initial deformation
      if (this->simulator_is_past_initialization() == false ||
          this->get_timestep_number() == 0)
        {
          TimerOutput::Scope timer (sim.computing_timer, "Mesh deformation initialize");

          make_initial_constraints();
          compute_mesh_displacements();
        }

      if (this->is_stokes_matrix_free())
        update_multilevel_deformation();
    }



    template <int dim>
    void MeshDeformationHandler<dim>::update_multilevel_deformation ()
    {
      Assert(this->is_stokes_matrix_free(), ExcInternalError());

      // Convert the mesh_displacements to a d:Vector that we can use
      // to transfer to the MG levels below. The conversion is done by
      // going through a ReadWriteVector.
      dealii::LinearAlgebra::distributed::Vector<double> displacements(mesh_deformation_dof_handler.locally_owned_dofs(),
                                                                       this->get_triangulation().get_communicator());
      dealii::LinearAlgebra::ReadWriteVector<double> rwv;
      rwv.reinit(mesh_displacements);
      displacements.import(rwv, VectorOperation::insert);

      mg_transfer.interpolate_to_mg(mesh_deformation_dof_handler,
                                    level_displacements,
                                    displacements);
    }



    template <int dim>
    const std::map<types::boundary_id, std::vector<std::string>> &
                                                              MeshDeformationHandler<dim>::get_active_mesh_deformation_names () const
    {
      return mesh_deformation_object_names;
    }



    template <int dim>
    const std::map<types::boundary_id,std::vector<std::unique_ptr<Interface<dim>>>> &
    MeshDeformationHandler<dim>::get_active_mesh_deformation_models () const
    {
      return mesh_deformation_objects;
    }



    template <int dim>
    const std::set<types::boundary_id> &
    MeshDeformationHandler<dim>::get_active_mesh_deformation_boundary_indicators () const
    {
      return prescribed_mesh_deformation_boundary_indicators;
    }



    template <int dim>
    const std::set<types::boundary_id> &
    MeshDeformationHandler<dim>::get_free_surface_boundary_indicators () const
    {
      return free_surface_boundary_indicators;
    }



    template <int dim>
    const LinearAlgebra::Vector &
    MeshDeformationHandler<dim>::get_mesh_displacements () const
    {
      return mesh_displacements;
    }



    template <int dim>
    const LinearAlgebra::Vector &
    MeshDeformationHandler<dim>::get_initial_topography () const
    {
      return initial_topography;
    }



    template <int dim>
    std::string
    get_valid_model_names_pattern ()
    {
      return std::get<dim>(registered_plugins).get_pattern_of_names ();
    }



    template <int dim>
    void
    MeshDeformationHandler<dim>::write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Mesh deformation interface",
                                                            out);
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<MeshDeformation::Interface<2>>::PluginInfo> *
                                                                          internal::Plugins::PluginList<MeshDeformation::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<MeshDeformation::Interface<3>>::PluginInfo> *
                                                                          internal::Plugins::PluginList<MeshDeformation::Interface<3>>::plugins = nullptr;
    }
  }

  namespace MeshDeformation
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class MeshDeformationHandler<dim>; \
  \
  template \
  std::string \
  get_valid_model_names_pattern<dim> ();

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
