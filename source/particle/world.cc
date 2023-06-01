/*
  Copyright (C) 2015 - 2021 by the authors of the ASPECT code.

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

#include <aspect/particle/world.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/citation_info.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/matrix_free/fe_point_evaluation.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <boost/serialization/map.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace aspect
{
  namespace Particle
  {
    template <int dim>
    World<dim>::World()
    {}

    template <int dim>
    World<dim>::~World()
    {}

    template <int dim>
    void
    World<dim>::initialize()
    {
      CitationInfo::add("particles");
      if (particle_load_balancing & ParticleLoadBalancing::repartition)
        this->get_triangulation().signals.cell_weight.connect(
          [&] (const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
               const typename parallel::distributed::Triangulation<dim>::CellStatus status)
          -> unsigned int
        {
          return this->cell_weight(cell, status);
        });

      // Create a particle handler that stores the future particles.
      // If we restarted from a checkpoint we will fill this particle handler
      // later with its serialized variables and stored particles
      particle_handler = std::make_unique<ParticleHandler<dim>>(this->get_triangulation(),
                                                                this->get_mapping(),
                                                                property_manager->get_n_property_components());

      particle_handler_backup.initialize(this->get_triangulation(),
                                         this->get_mapping(),
                                         property_manager->get_n_property_components());

      connect_to_signals(this->get_signals());
    }

    template <int dim>
    const Property::Manager<dim> &
    World<dim>::get_property_manager() const
    {
      return *property_manager;
    }



    template <int dim>
    const Particles::ParticleHandler<dim> &
    World<dim>::get_particle_handler() const
    {
      return *particle_handler.get();
    }



    template <int dim>
    Particles::ParticleHandler<dim> &
    World<dim>::get_particle_handler()
    {
      return *particle_handler.get();
    }



    template <int dim>
    void
    World<dim>::copy_particle_handler (const Particles::ParticleHandler<dim> &from_particle_handler,
                                       Particles::ParticleHandler<dim> &to_particle_handler) const
    {
      {
        TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Copy");

        to_particle_handler.copy_from(from_particle_handler);
      }
    }



    template <int dim>
    void
    World<dim>::backup_particles ()
    {
      copy_particle_handler (*particle_handler.get(), particle_handler_backup);
    }



    template <int dim>
    void
    World<dim>::restore_particles ()
    {
      copy_particle_handler (particle_handler_backup, *particle_handler.get());
    }



    template <int dim>
    const Interpolator::Interface<dim> &
    World<dim>::get_interpolator() const
    {
      return *interpolator;
    }



    template <int dim>
    types::particle_index
    World<dim>::n_global_particles() const
    {
      return particle_handler->n_global_particles();
    }



    template <int dim>
    void
    World<dim>::connect_to_signals(aspect::SimulatorSignals<dim> &signals)
    {
      signals.post_set_initial_state.connect(
        [&] (const SimulatorAccess<dim> &)
      {
        this->setup_initial_state();
      });

      connect_particle_handler_signals(signals,*particle_handler);
      // Particle handler backup will not be stored for checkpointing
      connect_particle_handler_signals(signals, particle_handler_backup, false);

      signals.post_refinement_load_user_data.connect(
        [&] (typename parallel::distributed::Triangulation<dim> &)
      {
        this->apply_particle_per_cell_bounds();
      });

      signals.post_resume_load_user_data.connect(
        [&] (typename parallel::distributed::Triangulation<dim> &)
      {
        this->apply_particle_per_cell_bounds();
      });
    }



    template <int dim>
    void
    World<dim>::connect_particle_handler_signals(aspect::SimulatorSignals<dim> &signals,
                                                 ParticleHandler<dim> &particle_handler_,
                                                 const bool connect_to_checkpoint_signals) const
    {
      signals.pre_refinement_store_user_data.connect(
        [&] (typename parallel::distributed::Triangulation<dim> &)
      {
        particle_handler_.register_store_callback_function();
      });

      signals.post_refinement_load_user_data.connect(
        [&] (typename parallel::distributed::Triangulation<dim> &)
      {
        particle_handler_.register_load_callback_function(false);
      });

      // Only connect to checkpoint signals if requested
      if (connect_to_checkpoint_signals)
        {
          signals.pre_checkpoint_store_user_data.connect(
            [&] (typename parallel::distributed::Triangulation<dim> &)
          {
            particle_handler_.register_store_callback_function();
          });

          signals.post_resume_load_user_data.connect(
            [&] (typename parallel::distributed::Triangulation<dim> &)
          {
            particle_handler_.register_load_callback_function(true);
          });
        }

      if (update_ghost_particles &&
          dealii::Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) > 1)
        {
          auto do_ghost_exchange = [&] (typename parallel::distributed::Triangulation<dim> &)
          {
            particle_handler_.exchange_ghost_particles();
          };
          signals.post_refinement_load_user_data.connect(do_ghost_exchange);
          signals.post_resume_load_user_data.connect(do_ghost_exchange);
        }

      signals.post_mesh_deformation.connect(
        [&] (const SimulatorAccess<dim> &)
      {
        particle_handler->sort_particles_into_subdomains_and_cells();
      },
      boost::signals2::at_front);
    }



    template <int dim>
    void
    World<dim>::apply_particle_per_cell_bounds()
    {
      // If any load balancing technique is selected that creates/destroys particles
      if (particle_load_balancing & ParticleLoadBalancing::remove_and_add_particles)
        {
          // First do some preparation for particle generation in poorly
          // populated areas. For this we need to know which particle ids to
          // generate so that they are globally unique.
          // Ensure this by communicating the number of particles that every
          // process is going to generate.
          particle_handler->update_cached_numbers();
          types::particle_index local_next_particle_index = particle_handler->get_next_free_particle_index();
          if (particle_load_balancing & ParticleLoadBalancing::add_particles)
            {
              types::particle_index particles_to_add_locally = 0;

              // Loop over all cells and determine the number of particles to generate
              for (const auto &cell : this->get_dof_handler().active_cell_iterators())
                if (cell->is_locally_owned())
                  {
                    const unsigned int particles_in_cell = particle_handler->n_particles_in_cell(cell);

                    if (particles_in_cell < min_particles_per_cell)
                      particles_to_add_locally += static_cast<types::particle_index> (min_particles_per_cell - particles_in_cell);
                  }

              // Determine the starting particle index of this process, which
              // is the highest currently existing particle index plus the sum
              // of the number of newly generated particles of all
              // processes with a lower rank.

              types::particle_index local_start_index = 0.0;

              const int ierr = MPI_Scan(&particles_to_add_locally, &local_start_index, 1, DEAL_II_PARTICLE_INDEX_MPI_TYPE, MPI_SUM, this->get_mpi_communicator());
              AssertThrowMPI(ierr);

              local_start_index -= particles_to_add_locally;
              local_next_particle_index += local_start_index;

              const types::particle_index globally_generated_particles =
                dealii::Utilities::MPI::sum(particles_to_add_locally,this->get_mpi_communicator());

              AssertThrow (particle_handler->get_next_free_particle_index()
                           <= std::numeric_limits<types::particle_index>::max() - globally_generated_particles,
                           ExcMessage("There is no free particle index left to generate a new particle id. Please check if your "
                                      "model generates unusually many new particles (by repeatedly deleting and regenerating particles), or "
                                      "recompile deal.II with the DEAL_II_WITH_64BIT_INDICES option enabled, to use 64-bit integers for "
                                      "particle ids."));
            }

          std::mt19937 random_number_generator;

          // Loop over all cells and generate or remove the particles cell-wise
          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                const unsigned int n_particles_in_cell = particle_handler->n_particles_in_cell(cell);

                // Add particles if necessary
                if ((particle_load_balancing & ParticleLoadBalancing::add_particles) &&
                    (n_particles_in_cell < min_particles_per_cell))
                  {
                    for (unsigned int i = n_particles_in_cell; i < min_particles_per_cell; ++i,++local_next_particle_index)
                      {
                        std::pair<Particles::internal::LevelInd,Particles::Particle<dim>> new_particle = generator->generate_particle(cell,local_next_particle_index);

                        const std::vector<double> particle_properties =
                          property_manager->initialize_late_particle(new_particle.second.get_location(),
                                                                     *particle_handler,
                                                                     *interpolator,
                                                                     cell);

                        typename ParticleHandler<dim>::particle_iterator particle = particle_handler->insert_particle(new_particle.second,
                                                                                    typename parallel::distributed::Triangulation<dim>::cell_iterator (&this->get_triangulation(),
                                                                                        new_particle.first.first,
                                                                                        new_particle.first.second));
                        particle->set_properties(particle_properties);
                      }
                  }

                // Remove particles if necessary
                else if ((particle_load_balancing & ParticleLoadBalancing::remove_particles) &&
                         (n_particles_in_cell > max_particles_per_cell))
                  {
                    const unsigned int n_particles_to_remove = n_particles_in_cell - max_particles_per_cell;

#if DEAL_II_VERSION_GTE(10,0,0)
                    for (unsigned int i=0; i < n_particles_to_remove; ++i)
                      {
                        const unsigned int current_n_particles_in_cell = particle_handler->n_particles_in_cell(cell);
                        const unsigned int index_to_remove = std::uniform_int_distribution<unsigned int>
                                                             (0,current_n_particles_in_cell-1)(random_number_generator);

                        auto particle_to_remove = particle_handler->particles_in_cell(cell).begin();
                        std::advance(particle_to_remove, index_to_remove);
                        particle_handler->remove_particle(particle_to_remove);
                      }
#else
                    const boost::iterator_range<typename ParticleHandler<dim>::particle_iterator> particles_in_cell
                      = particle_handler->particles_in_cell(cell);

                    std::set<unsigned int> particle_ids_to_remove;
                    while (particle_ids_to_remove.size() < n_particles_to_remove)
                      particle_ids_to_remove.insert(random_number_generator() % n_particles_in_cell);

                    std::vector<typename ParticleHandler<dim>::particle_iterator> particles_to_remove;
                    particles_to_remove.reserve(n_particles_to_remove);

                    for (const auto id : particle_ids_to_remove)
                      {
                        typename ParticleHandler<dim>::particle_iterator particle_to_remove = particles_in_cell.begin();
                        std::advance(particle_to_remove, id);

                        particles_to_remove.push_back(particle_to_remove);
                      }

                    for (const auto &particle : particles_to_remove)
                      {
                        particle_handler->remove_particle(particle);
                      }
#endif
                  }
              }

          particle_handler->update_cached_numbers();
        }
    }

    template <int dim>
    unsigned int
    World<dim>::cell_weight(const typename parallel::distributed::Triangulation<dim>::cell_iterator &cell,
                            const typename parallel::distributed::Triangulation<dim>::CellStatus status)
    {
      if (cell->is_active() && !cell->is_locally_owned())
        return 0;

      if (status == parallel::distributed::Triangulation<dim>::CELL_PERSIST
          || status == parallel::distributed::Triangulation<dim>::CELL_REFINE)
        {
          const unsigned int n_particles_in_cell = particle_handler->n_particles_in_cell(cell);
          return n_particles_in_cell * particle_weight;
        }
      else if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN)
        {
          unsigned int n_particles_in_cell = 0;

          for (unsigned int child_index = 0; child_index < GeometryInfo<dim>::max_children_per_cell; ++child_index)
            n_particles_in_cell += particle_handler->n_particles_in_cell(cell->child(child_index));

          return n_particles_in_cell * particle_weight;
        }

      Assert (false, ExcInternalError());
      return 0;
    }


    template <int dim>
    std::map<types::subdomain_id, unsigned int>
    World<dim>::get_subdomain_id_to_neighbor_map() const
    {
      std::map<types::subdomain_id, unsigned int> subdomain_id_to_neighbor_map;
      const std::set<types::subdomain_id> ghost_owners = this->get_triangulation().ghost_owners();
      std::set<types::subdomain_id>::const_iterator ghost_owner = ghost_owners.begin();

      for (unsigned int neighbor_id=0; neighbor_id<ghost_owners.size(); ++neighbor_id,++ghost_owner)
        {
          subdomain_id_to_neighbor_map.insert(std::make_pair(*ghost_owner,neighbor_id));
        }
      return subdomain_id_to_neighbor_map;
    }



    template <int dim>
    void
    World<dim>::local_initialize_particles(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                           const typename ParticleHandler<dim>::particle_iterator &end_particle)
    {
      for (typename ParticleHandler<dim>::particle_iterator it = begin_particle; it!=end_particle; ++it)
        property_manager->initialize_one_particle(it);
    }



    template <int dim>
    void
    World<dim>::local_update_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                       const typename ParticleHandler<dim>::particle_iterator &end_particle,
                                       internal::SolutionEvaluators<dim> &evaluators)
    {
#if DEAL_II_VERSION_GTE(10,0,0)
      const unsigned int n_particles_in_cell = particle_handler->n_particles_in_cell(cell);
#else
      const unsigned int n_particles_in_cell = std::distance(begin_particle,end_particle);
#endif

      std::vector<Point<dim> > positions;
      positions.reserve(n_particles_in_cell);

      for (auto particle = begin_particle; particle!=end_particle; ++particle)
        positions.push_back(particle->get_reference_location());

      const UpdateFlags update_flags = property_manager->get_needed_update_flags();

      boost::container::small_vector<double, 100> solution_values(this->get_fe().dofs_per_cell);

      cell->get_dof_values(this->get_solution(),
                           solution_values.begin(),
                           solution_values.end());

      if (update_flags & (update_values | update_gradients))
        evaluators.reinit(cell, positions, {solution_values.data(), solution_values.size()}, update_flags);

      Vector<double> solution;
      if (update_flags & update_values)
        solution.reinit(this->introspection().n_components);

      std::vector<Tensor<1,dim>> gradients;
      if (update_flags & update_gradients)
        gradients.resize(this->introspection().n_components);

      auto particle = begin_particle;
      for (unsigned int i = 0; particle!=end_particle; ++particle,++i)
        {
          // Evaluate the solution, but only if it is requested in the update_flags
          if (update_flags & update_values)
            evaluators.get_solution(i, solution);

          // Evaluate the gradients, but only if they are requested in the update_flags
          if (update_flags & update_gradients)
            evaluators.get_gradients(i, gradients);

          property_manager->update_one_particle(particle,
                                                solution,
                                                gradients);
        }
    }



    template <int dim>
    void
    World<dim>::local_update_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                       const typename ParticleHandler<dim>::particle_iterator &end_particle)
    {
#if DEAL_II_VERSION_GTE(10,0,0)
      const unsigned int n_particles_in_cell = particle_handler->n_particles_in_cell(cell);
#else
      const unsigned int n_particles_in_cell = std::distance(begin_particle,end_particle);
#endif
      const unsigned int solution_components = this->introspection().n_components;

      Vector<double>              value (solution_components);
      std::vector<Tensor<1,dim>> gradient (solution_components,Tensor<1,dim>());

      std::vector<Vector<double>>              values(n_particles_in_cell,value);
      std::vector<std::vector<Tensor<1,dim>>> gradients(n_particles_in_cell,gradient);
      std::vector<Point<dim>>                  positions(n_particles_in_cell);

      typename ParticleHandler<dim>::particle_iterator it = begin_particle;
      for (unsigned int i = 0; it!=end_particle; ++it,++i)
        {
          positions[i] = it->get_reference_location();
        }

      const Quadrature<dim> quadrature_formula(positions);
      const UpdateFlags update_flags = property_manager->get_needed_update_flags();
      FEValues<dim> fe_value (this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_flags);

      fe_value.reinit (cell);
      if (update_flags & update_values)
        fe_value.get_function_values (this->get_solution(),
                                      values);
      if (update_flags & update_gradients)
        fe_value.get_function_gradients (this->get_solution(),
                                         gradients);

      it = begin_particle;
      for (unsigned int i = 0; it!=end_particle; ++it,++i)
        {
          property_manager->update_one_particle(it,
                                                values[i],
                                                gradients[i]);
        }
    }




    template <int dim>
    void
    World<dim>::local_advect_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                       const typename ParticleHandler<dim>::particle_iterator &end_particle,
                                       internal::SolutionEvaluators<dim> &evaluators)
    {
#if DEAL_II_VERSION_GTE(10,0,0)
      const unsigned int n_particles_in_cell = particle_handler->n_particles_in_cell(cell);
#else
      const unsigned int n_particles_in_cell = std::distance(begin_particle,end_particle);
#endif

      boost::container::small_vector<Point<dim>, 100>   positions;
      positions.reserve(n_particles_in_cell);
      for (auto particle = begin_particle; particle!=end_particle; ++particle)
        positions.push_back(particle->get_reference_location());

      boost::container::small_vector<double, 100> solution_values(this->get_fe().dofs_per_cell);
      boost::container::small_vector<double, 100> old_solution_values(this->get_fe().dofs_per_cell);

      cell->get_dof_values(this->get_current_linearization_point(),
                           solution_values.begin(),
                           solution_values.end());

      cell->get_dof_values(this->get_old_solution(),
                           old_solution_values.begin(),
                           old_solution_values.end());

      const bool use_fluid_velocity = this->include_melt_transport() &&
                                      property_manager->get_data_info().fieldname_exists("melt_presence");
      auto &evaluator = evaluators.get_velocity_or_fluid_velocity_evaluator(use_fluid_velocity);

      evaluator.reinit (cell, {positions.data(),positions.size()});

      evaluator.evaluate({solution_values.data(),solution_values.size()},
                         EvaluationFlags::values);

      std::vector<Tensor<1,dim>> velocities;
      velocities.reserve(n_particles_in_cell);
      for (unsigned int i=0; i<n_particles_in_cell; ++i)
        velocities.push_back(evaluator.get_value(i));

      evaluator.evaluate({old_solution_values.data(),old_solution_values.size()},
                         EvaluationFlags::values);

      std::vector<Tensor<1,dim>> old_velocities;
      old_velocities.reserve(n_particles_in_cell);
      for (unsigned int i=0; i<n_particles_in_cell; ++i)
        old_velocities.push_back(evaluator.get_value(i));

      integrator->local_integrate_step(begin_particle,
                                       end_particle,
                                       old_velocities,
                                       velocities,
                                       this->get_timestep());
    }



    template <int dim>
    void
    World<dim>::local_advect_particles(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                       const typename ParticleHandler<dim>::particle_iterator &end_particle)
    {
#if DEAL_II_VERSION_GTE(10,0,0)
      const unsigned int n_particles_in_cell = particle_handler->n_particles_in_cell(cell);
#else
      const unsigned int n_particles_in_cell = std::distance(begin_particle,end_particle);
#endif

      std::vector<Tensor<1,dim>>  velocity(n_particles_in_cell);
      std::vector<Tensor<1,dim>>  old_velocity(n_particles_in_cell);

      // Below we manually evaluate the solution at all support points of the
      // current cell, and then use the shape functions to interpolate the
      // solution to the particle points. All of this can be done with less
      // code using an FEValues object, but since this object initializes a lot
      // of memory for other purposes and we can not reuse the FEValues object
      // for other cells, it is much faster to do the work manually. Also this
      // function is quite performance critical.

      std::vector<types::global_dof_index> cell_dof_indices (this->get_fe().dofs_per_cell);
      cell->get_dof_indices (cell_dof_indices);

      const FiniteElement<dim> &velocity_fe = this->get_fe().base_element(this->introspection()
                                                                          .base_elements.velocities);

      const bool compute_fluid_velocity = this->include_melt_transport() &&
                                          property_manager->get_data_info().fieldname_exists("melt_presence");

      const unsigned int fluid_component_index = (compute_fluid_velocity ?
                                                  this->introspection().variable("fluid velocity").first_component_index
                                                  :
                                                  numbers::invalid_unsigned_int);

      // In regions without melt, the fluid velocity equals the solid velocity, so we can use it for all particles.
      std::vector<bool> use_fluid_velocity((compute_fluid_velocity ?
                                            n_particles_in_cell
                                            :
                                            0), compute_fluid_velocity);

      for (unsigned int j=0; j<velocity_fe.dofs_per_cell; ++j)
        {
          Tensor<1,dim> velocity_at_support_point;
          Tensor<1,dim> old_velocity_at_support_point;

          for (unsigned int dir=0; dir<dim; ++dir)
            {
              const unsigned int support_point_index
                = this->get_fe().component_to_system_index(this->introspection()
                                                           .component_indices.velocities[dir],j);

              velocity_at_support_point[dir] = this->get_current_linearization_point()[cell_dof_indices[support_point_index]];
              old_velocity_at_support_point[dir] = this->get_old_solution()[cell_dof_indices[support_point_index]];
            }

          Tensor<1,dim> fluid_velocity_at_support_point;
          Tensor<1,dim> old_fluid_velocity_at_support_point;

          if (compute_fluid_velocity)
            for (unsigned int dir=0; dir<dim; ++dir)
              {
                const unsigned int support_point_index
                  = this->get_fe().component_to_system_index(fluid_component_index + dir,j);

                fluid_velocity_at_support_point[dir] = this->get_solution()[cell_dof_indices[support_point_index]];
                old_fluid_velocity_at_support_point[dir] = this->get_old_solution()[cell_dof_indices[support_point_index]];
              }

          typename ParticleHandler<dim>::particle_iterator it = begin_particle;
          for (unsigned int particle_index = 0; it!=end_particle; ++it,++particle_index)
            {
              // melt FE uses the same FE so the shape value is the same
              const double shape_value = velocity_fe.shape_value(j,it->get_reference_location());

              if (compute_fluid_velocity && use_fluid_velocity[particle_index])
                {
                  velocity[particle_index] += fluid_velocity_at_support_point * shape_value;
                  old_velocity[particle_index] += old_fluid_velocity_at_support_point * shape_value;
                }
              else
                {
                  velocity[particle_index] += velocity_at_support_point * shape_value;
                  old_velocity[particle_index] += old_velocity_at_support_point * shape_value;
                }
            }
        }

      integrator->local_integrate_step(begin_particle,
                                       end_particle,
                                       old_velocity,
                                       velocity,
                                       this->get_timestep());
    }



    template <int dim>
    void
    World<dim>::setup_initial_state ()
    {
      // If we are in the first adaptive refinement cycle generate particles
      if (this->get_pre_refinement_step() == 0)
        generate_particles();

      // And initialize the particle properties according to the initial
      // conditions on the current mesh
      initialize_particles();
    }



    template <int dim>
    void
    World<dim>::generate_particles()
    {
      TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Generate");

      std::multimap<Particles::internal::LevelInd, Particles::Particle<dim>> particles;
      generator->generate_particles(particles);

      std::multimap<typename Triangulation<dim>::active_cell_iterator, Particles::Particle<dim>> new_particles;

      for (const auto &particle : particles)
        new_particles.insert(new_particles.end(),
                             std::make_pair(typename Triangulation<dim>::active_cell_iterator(&this->get_triangulation(),
                                            particle.first.first, particle.first.second),
                                            particle.second));

      particle_handler->insert_particles(new_particles);
    }



    template <int dim>
    void
    World<dim>::initialize_particles()
    {
#if !DEAL_II_VERSION_GTE(10,0,0)
      // Initialize the particle's access to the property_pool. This is necessary
      // even if the Particle do not carry properties, because they need a
      // way to determine the number of properties they carry.
      for (ParticleIterator<dim> particle = particle_handler->begin(); particle!=particle_handler->end(); ++particle)
        particle->set_property_pool(particle_handler->get_property_pool());
#endif

      // TODO: Change this loop over all cells to use the WorkStream interface
      if (property_manager->get_n_property_components() > 0)
        {
          TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Initialize properties");

          particle_handler->get_property_pool().reserve(2 * particle_handler->n_locally_owned_particles());


          if (particle_handler->n_locally_owned_particles() > 0)
            local_initialize_particles(particle_handler->begin(),
                                       particle_handler->end());

          if (update_ghost_particles &&
              dealii::Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) > 1)
            {
              TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Exchange ghosts");
              particle_handler->exchange_ghost_particles();
            }
        }
    }



    namespace internal
    {
      // This class evaluates the solution vector at arbitrary positions inside a cell.
      // This base class only provides the interface for SolutionEvaluatorsImplementation.
      // See there for more details.
      template <int dim>
      class SolutionEvaluators
      {
        public:
          // virtual Destructor.
          virtual ~SolutionEvaluators() = default;

          // Reinitialize all variables to evaluate the given solution for the given cell
          // and the given positions. The update flags control if only the solution or
          // also the gradients should be evaluated.
          // If other flags are set an assertion is triggered.
          virtual
          void
          reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
                 const ArrayView<Point<dim>> &positions,
                 const ArrayView<double> &solution_values,
                 const UpdateFlags update_flags) = 0;

          // Fill @p solution with all solution components at the given @p evaluation_point. Note
          // that this function only works after a successful call to reinit(),
          // because this function only returns the results of the computation that
          // happened in reinit().
          virtual
          void get_solution(const unsigned int evaluation_point,
                            Vector<double> &solution) = 0;

          // Fill @p gradients with all solution gradients at the given @p evaluation_point. Note
          // that this function only works after a successful call to reinit(),
          // because this function only returns the results of the computation that
          // happened in reinit().
          virtual
          void get_gradients(const unsigned int evaluation_point,
                             std::vector<Tensor<1,dim>> &gradients) = 0;

          // Return the evaluator for velocity or fluid velocity. This is the only
          // information necessary for advecting particles.
          virtual
          FEPointEvaluation<dim, dim> &
          get_velocity_or_fluid_velocity_evaluator(const bool use_fluid_velocity) = 0;
      };

      // This class evaluates the solution vector at arbitrary positions inside a cell.
      // It uses the deal.II class FEPointEvaluation to do this efficiently. Because
      // FEPointEvaluation only supports a single finite element, but ASPECT uses a FESystem with
      // many components, this class creates several FEPointEvaluation objects that are used for
      // the individual finite elements of our solution (pressure, velocity, temperature, and
      // all other optional variables). Because FEPointEvaluation is templated based on the
      // number of components, but ASPECT only knows the number of components at runtime
      // we create this derived class with an additional template. This makes it possible
      // to access the functionality through the base class, but create an object of this
      // derived class with the correct number of components at runtime.
      template <int dim, int n_compositional_fields>
      class SolutionEvaluatorsImplementation: public SolutionEvaluators<dim>
      {
        public:
          // Constructor. Create the member variables given a simulator and a set of
          // update flags. The update flags control if only the solution or also the gradients
          // should be evaluated.
          SolutionEvaluatorsImplementation(const SimulatorAccess<dim> &simulator,
                                           const UpdateFlags update_flags);

          // Reinitialize all variables to evaluate the given solution for the given cell
          // and the given positions. The update flags control if only the solution or
          // also the gradients should be evaluated.
          // If other flags are set an assertion is triggered.
          void
          reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
                 const ArrayView<Point<dim>> &positions,
                 const ArrayView<double> &solution_values,
                 const UpdateFlags update_flags) override;

          // Return the value of all solution components at the given evaluation point. Note
          // that this function only works after a successful call to reinit(),
          // because this function only returns the results of the computation that
          // happened in reinit().
          void get_solution(const unsigned int evaluation_point,
                            Vector<double> &solution) override;

          // Return the value of all solution gradients at the given evaluation point. Note
          // that this function only works after a successful call to reinit(),
          // because this function only returns the results of the computation that
          // happened in reinit().
          void get_gradients(const unsigned int evaluation_point,
                             std::vector<Tensor<1,dim>> &gradients) override;

          // Return the evaluator for velocity or fluid velocity. This is the only
          // information necessary for advecting particles.
          FEPointEvaluation<dim, dim> &
          get_velocity_or_fluid_velocity_evaluator(const bool use_fluid_velocity) override;

        private:
          // FEPointEvaluation objects for all common
          // components of ASPECT's finite element solution.
          // These objects are used inside of the member functions of this class.
          FEPointEvaluation<dim, dim> velocity;
          FEPointEvaluation<1, dim> pressure;
          FEPointEvaluation<1, dim> temperature;

          // If instantiated evaluate multiple compositions at once, if
          // not fall back to evaluating them individually.
          FEPointEvaluation<n_compositional_fields, dim> compositions;
          std::vector<FEPointEvaluation<1, dim>> additional_compositions;

          // Pointers to FEPointEvaluation objects for all melt
          // components of ASPECT's finite element solution, which only
          // point to valid objects in case we use melt transport. Other
          // documentation like for the objects directly above.
          std::unique_ptr<FEPointEvaluation<dim, dim>> fluid_velocity;
          std::unique_ptr<FEPointEvaluation<1, dim>> compaction_pressure;
          std::unique_ptr<FEPointEvaluation<1, dim>> fluid_pressure;

          // The component indices for the three melt formulation
          // variables fluid velocity, compaction pressure, and
          // fluid pressure (in this order). They are cached
          // to avoid repeated expensive lookups.
          std::array<unsigned int, 3> melt_component_indices;

          // Reference to the active simulator access object. Provides
          // access to the general simulation variables.
          const SimulatorAccess<dim> &simulator_access;
      };



      template <int dim, int n_compositional_fields>
      SolutionEvaluatorsImplementation<dim, n_compositional_fields>::SolutionEvaluatorsImplementation(const SimulatorAccess<dim> &simulator,
          const UpdateFlags update_flags)
        :
        velocity(simulator.get_mapping(),
                 simulator.get_fe(),
                 update_flags,
                 simulator.introspection().component_indices.velocities[0]),
        pressure(simulator.get_mapping(),
                 simulator.get_fe(),
                 update_flags,
                 simulator.introspection().component_indices.pressure),
        temperature(simulator.get_mapping(),
                    simulator.get_fe(),
                    update_flags,
                    simulator.introspection().component_indices.temperature),
        compositions(simulator.get_mapping(),
                     simulator.get_fe(),
                     update_flags,
                     simulator.n_compositional_fields() > 0 ? simulator.introspection().component_indices.compositional_fields[0] : simulator.introspection().component_indices.temperature),
        melt_component_indices(),
        simulator_access(simulator)
      {
        // Create the evaluators for all compositional fields beyond the ones this class was
        // instantiated for
        const unsigned int n_total_compositional_fields = simulator_access.n_compositional_fields();
        const auto &component_indices = simulator_access.introspection().component_indices.compositional_fields;
        for (unsigned int composition = n_compositional_fields; composition < n_total_compositional_fields; ++composition)
          additional_compositions.emplace_back(FEPointEvaluation<1, dim>(simulator_access.get_mapping(),
                                                                         simulator_access.get_fe(),
                                                                         update_flags,
                                                                         component_indices[composition]));

        // Create the melt evaluators, but only if we use melt transport in the model
        if (simulator_access.include_melt_transport())
          {
            // Store the melt component indices to avoid repeated string lookups later on
            melt_component_indices[0] = simulator_access.introspection().variable("fluid velocity").first_component_index;
            melt_component_indices[1] = simulator_access.introspection().variable("fluid pressure").first_component_index;
            melt_component_indices[2] = simulator_access.introspection().variable("compaction pressure").first_component_index;

            fluid_velocity = std::make_unique<FEPointEvaluation<dim, dim>>(simulator_access.get_mapping(),
                                                                           simulator_access.get_fe(),
                                                                           update_flags,
                                                                           melt_component_indices[0]);
            fluid_pressure = std::make_unique<FEPointEvaluation<1, dim>>(simulator_access.get_mapping(),
                                                                         simulator_access.get_fe(),
                                                                         update_flags,
                                                                         melt_component_indices[1]);
            compaction_pressure = std::make_unique<FEPointEvaluation<1, dim>>(simulator_access.get_mapping(),
                                                                              simulator_access.get_fe(),
                                                                              update_flags,
                                                                              melt_component_indices[2]);

          }
      }



      template <int dim, int n_compositional_fields>
      void
      SolutionEvaluatorsImplementation<dim, n_compositional_fields>::reinit(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                                            const ArrayView<Point<dim>> &positions,
                                                                            const ArrayView<double> &solution_values,
                                                                            const UpdateFlags update_flags)
      {
        // FEPointEvaluation uses different evaluation flags than the common UpdateFlags.
        // Translate between the two.
        EvaluationFlags::EvaluationFlags evaluation_flags = EvaluationFlags::nothing;

        if (update_flags & update_values)
          evaluation_flags = evaluation_flags | EvaluationFlags::values;

        if (update_flags & update_gradients)
          evaluation_flags = evaluation_flags | EvaluationFlags::gradients;

        // Make sure only the flags are set that we can deal with at the moment
        Assert ((update_flags & ~(update_gradients | update_values)) == false,
                ExcNotImplemented());

        // Reinitialize and evaluate all evaluators.
        // TODO: It would be nice to be able to hand over a ComponentMask
        // to specify which evaluators to use. Currently, this is only
        // possible by manually accessing the public members of this class.
        velocity.reinit (cell, positions);

#if DEAL_II_VERSION_GTE(10,0,0)
        // Only compute the mapping data once for velocity,
        // and reuse it for the other components.
        const auto &mapping_data = velocity.get_mapping_data();

        pressure.reinit (cell, positions, mapping_data);
        temperature.reinit (cell, positions, mapping_data);
        compositions.reinit (cell, positions, mapping_data);

        for (auto &evaluator_composition: additional_compositions)
          evaluator_composition.reinit (cell, positions, mapping_data);

        if (simulator_access.include_melt_transport())
          {
            fluid_velocity->reinit (cell, positions, mapping_data);
            fluid_pressure->reinit (cell, positions, mapping_data);
            compaction_pressure->reinit (cell, positions, mapping_data);
          }
#else
        pressure.reinit (cell, positions);
        temperature.reinit (cell, positions);
        compositions.reinit (cell, positions);

        for (auto &evaluator_composition: additional_compositions)
          evaluator_composition.reinit (cell, positions);

        if (simulator_access.include_melt_transport())
          {
            fluid_velocity->reinit (cell, positions);
            fluid_pressure->reinit (cell, positions);
            compaction_pressure->reinit (cell, positions);
          }
#endif

        velocity.evaluate (solution_values, evaluation_flags);
        pressure.evaluate (solution_values, evaluation_flags);
        temperature.evaluate (solution_values, evaluation_flags);
        compositions.evaluate (solution_values, evaluation_flags);

        for (auto &evaluator_composition: additional_compositions)
          evaluator_composition.evaluate (solution_values, evaluation_flags);

        if (simulator_access.include_melt_transport())
          {
            fluid_velocity->evaluate (solution_values, evaluation_flags);
            fluid_pressure->evaluate (solution_values, evaluation_flags);
            compaction_pressure->evaluate (solution_values, evaluation_flags);
          }
      }



      template <int dim, int n_compositional_fields>
      void
      SolutionEvaluatorsImplementation<dim, n_compositional_fields>::get_solution(const unsigned int evaluation_point,
                                                                                  Vector<double> &solution)
      {
        Assert(solution.size() == simulator_access.introspection().n_components,
               ExcDimensionMismatch(solution.size(), simulator_access.introspection().n_components));

        const auto &component_indices = simulator_access.introspection().component_indices;

        const Tensor<1,dim> velocity_value = velocity.get_value(evaluation_point);
        for (unsigned int j=0; j<dim; ++j)
          solution[component_indices.velocities[j]] = velocity_value[j];

        solution[component_indices.pressure] = pressure.get_value(evaluation_point);
        solution[component_indices.temperature] = temperature.get_value(evaluation_point);

        const typename FEPointEvaluation<n_compositional_fields, dim>::value_type composition_values = compositions.get_value(evaluation_point);
        for (unsigned int j=0; j<n_compositional_fields; ++j)
          solution[component_indices.compositional_fields[j]] = dealii::internal::FEPointEvaluation::EvaluatorTypeTraits<dim, n_compositional_fields, double>::access(composition_values,j);

        const unsigned int n_additional_compositions = additional_compositions.size();
        for (unsigned int j=0; j<n_additional_compositions; ++j)
          solution[component_indices.compositional_fields[n_compositional_fields+j]] = additional_compositions[j].get_value(evaluation_point);

        if (simulator_access.include_melt_transport())
          {
            const Tensor<1,dim> fluid_velocity_value = velocity.get_value(evaluation_point);
            for (unsigned int j=0; j<dim; ++j)
              solution[melt_component_indices[0]+j] = fluid_velocity_value[j];

            solution[melt_component_indices[1]] = fluid_pressure->get_value(evaluation_point);
            solution[melt_component_indices[2]] = compaction_pressure->get_value(evaluation_point);
          }
      }



      template <int dim, int n_compositional_fields>
      void
      SolutionEvaluatorsImplementation<dim, n_compositional_fields>::get_gradients(const unsigned int evaluation_point,
                                                                                   std::vector<Tensor<1,dim>> &gradients)
      {
        Assert(gradients.size() == simulator_access.introspection().n_components,
               ExcDimensionMismatch(gradients.size(), simulator_access.introspection().n_components));

        const auto &component_indices = simulator_access.introspection().component_indices;

        const Tensor<2,dim> velocity_gradient = velocity.get_gradient(evaluation_point);
        for (unsigned int j=0; j<dim; ++j)
          gradients[component_indices.velocities[j]] = velocity_gradient[j];

        gradients[component_indices.pressure] = pressure.get_gradient(evaluation_point);
        gradients[component_indices.temperature] = temperature.get_gradient(evaluation_point);

        const typename FEPointEvaluation<n_compositional_fields, dim>::gradient_type composition_gradients = compositions.get_gradient(evaluation_point);
        for (unsigned int j=0; j<n_compositional_fields; ++j)
          gradients[component_indices.compositional_fields[j]] = dealii::internal::FEPointEvaluation::EvaluatorTypeTraits<dim, n_compositional_fields, double>::access(composition_gradients,j);

        const unsigned int n_additional_compositions = additional_compositions.size();
        for (unsigned int j=0; j<n_additional_compositions; ++j)
          gradients[component_indices.compositional_fields[n_compositional_fields+j]] = additional_compositions[j].get_gradient(evaluation_point);

        if (simulator_access.include_melt_transport())
          {
            const Tensor<2,dim> fluid_velocity_gradient = velocity.get_gradient(evaluation_point);
            for (unsigned int j=0; j<dim; ++j)
              gradients[melt_component_indices[0]+j] = fluid_velocity_gradient[j];

            gradients[melt_component_indices[1]] = fluid_pressure->get_gradient(evaluation_point);
            gradients[melt_component_indices[2]] = compaction_pressure->get_gradient(evaluation_point);
          }
      }


      template <int dim, int n_compositional_fields>
      FEPointEvaluation<dim, dim> &
      SolutionEvaluatorsImplementation<dim, n_compositional_fields>::get_velocity_or_fluid_velocity_evaluator(const bool use_fluid_velocity)
      {
        if (use_fluid_velocity)
          return *fluid_velocity;
        else
          return velocity;

        return velocity;
      }

      // A function to create a pointer to a SolutionEvaluators object.
      template <int dim>
      std::unique_ptr<internal::SolutionEvaluators<dim>>
                                                      construct_solution_evaluators (const SimulatorAccess<dim> &simulator_access,
                                                                                     const UpdateFlags update_flags)
      {
        switch (simulator_access.n_compositional_fields())
          {
            case 0:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,0>>(simulator_access, update_flags);
            case 1:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,1>>(simulator_access, update_flags);
            case 2:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,2>>(simulator_access, update_flags);
            case 3:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,3>>(simulator_access, update_flags);
            case 4:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,4>>(simulator_access, update_flags);
            case 5:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,5>>(simulator_access, update_flags);
            case 6:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,6>>(simulator_access, update_flags);
            case 7:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,7>>(simulator_access, update_flags);
            case 8:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,8>>(simulator_access, update_flags);
            case 9:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,9>>(simulator_access, update_flags);
            case 10:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,10>>(simulator_access, update_flags);
            case 11:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,11>>(simulator_access, update_flags);
            case 12:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,12>>(simulator_access, update_flags);
            case 13:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,13>>(simulator_access, update_flags);
            case 14:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,14>>(simulator_access, update_flags);
            case 15:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,15>>(simulator_access, update_flags);
            case 16:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,16>>(simulator_access, update_flags);
            case 17:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,17>>(simulator_access, update_flags);
            case 18:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,18>>(simulator_access, update_flags);
            case 19:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,19>>(simulator_access, update_flags);
            // Return the maximally instantiated object. The class will handle additional compositional fields
            // by dynamically allocating additional scalar evaluators.
            default:
              return std::make_unique<SolutionEvaluatorsImplementation<dim,20>>(simulator_access, update_flags);
          }
      }
    }



    template <int dim>
    void
    World<dim>::update_particles()
    {
      // TODO: Change this loop over all cells to use the WorkStream interface

      if (property_manager->get_n_property_components() > 0)
        {
          TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Update properties");

          const UpdateFlags update_flags = property_manager->get_needed_update_flags();

          // Only use deal.II FEPointEvaluation if its fast path is used. Prior to deal.II 10.0
          // FEPointEvaluation did not support MappingCartesian for box geometries, and there was
          // a bug for dynamically allocating scalar evaluators for individual components of a
          // base element with multiplicity (see https://github.com/dealii/dealii/pull/12786).
          bool use_fast_path = false;
#if DEAL_II_VERSION_GTE(10,0,0)
          if (dynamic_cast<const MappingQGeneric<dim> *>(&this->get_mapping()) != nullptr ||
              dynamic_cast<const MappingCartesian<dim> *>(&this->get_mapping()) != nullptr)
            use_fast_path = true;
#else
          if (dynamic_cast<const MappingQGeneric<dim> *>(&this->get_mapping()) != nullptr &&
              this->n_compositional_fields() <= 20)
            use_fast_path = true;
#endif
          std::unique_ptr<internal::SolutionEvaluators<dim>> evaluators;

          if (use_fast_path == true)
            evaluators = internal::construct_solution_evaluators(*this,
                                                                 update_flags);


          // Loop over all cells and update the particles cell-wise
          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                typename ParticleHandler<dim>::particle_iterator_range
                particles_in_cell = particle_handler->particles_in_cell(cell);

                // Only update particles, if there are any in this cell
                if (particles_in_cell.begin() != particles_in_cell.end())
                  {

                    if (use_fast_path)
                      local_update_particles(cell,
                                             particles_in_cell.begin(),
                                             particles_in_cell.end(),
                                             *evaluators);
                    else
                      local_update_particles(cell,
                                             particles_in_cell.begin(),
                                             particles_in_cell.end());
                  }

              }
        }
    }



    template <int dim>
    void
    World<dim>::advect_particles()
    {
      {
        // TODO: Change this loop over all cells to use the WorkStream interface
        TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Advect");

        std::unique_ptr<internal::SolutionEvaluators<dim>> evaluators =
                                                          std::make_unique<internal::SolutionEvaluatorsImplementation<dim, 0>>(*this,
                                                              update_values);

        // Loop over all cells and advect the particles cell-wise
        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned())
            {
              const typename ParticleHandler<dim>::particle_iterator_range
              particles_in_cell = particle_handler->particles_in_cell(cell);

              // Only advect particles, if there are any in this cell
              if (particles_in_cell.begin() != particles_in_cell.end())
                {
                  // Only use deal.II FEPointEvaluation if it's fast path is used
                  bool use_fast_path = false;
#if DEAL_II_VERSION_GTE(10,0,0)
                  if (dynamic_cast<const MappingQGeneric<dim> *>(&this->get_mapping()) != nullptr ||
                      dynamic_cast<const MappingCartesian<dim> *>(&this->get_mapping()) != nullptr)
                    use_fast_path = true;
#else
                  if (dynamic_cast<const MappingQGeneric<dim> *>(&this->get_mapping()) != nullptr)
                    use_fast_path = true;
#endif

                  if (use_fast_path)
                    local_advect_particles(cell,
                                           particles_in_cell.begin(),
                                           particles_in_cell.end(),
                                           *evaluators);
                  else
                    local_advect_particles(cell,
                                           particles_in_cell.begin(),
                                           particles_in_cell.end());


                }
            }
      }

      {
        TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Sort");
        // Find the cells that the particles moved to
        particle_handler->sort_particles_into_subdomains_and_cells();
      }
    }



    template <int dim>
    void
    World<dim>::advance_timestep()
    {
      do
        {
          advect_particles();
        }
      // Keep calling the integrator until it indicates it is finished
      while (integrator->new_integration_step());

      apply_particle_per_cell_bounds();

      // Update particle properties
      if (property_manager->need_update() == Property::update_time_step)
        update_particles();

      // Now that all particle information was updated, exchange the new
      // ghost particles.
      if (update_ghost_particles &&
          dealii::Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) > 1)
        {
          TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Exchange ghosts");
          particle_handler->exchange_ghost_particles();
        }
    }



    template <int dim>
    void
    World<dim>::save (std::ostringstream &os) const
    {
      aspect::oarchive oa (os);
      oa << (*this);
    }



    template <int dim>
    void
    World<dim>::load (std::istringstream &is)
    {
      aspect::iarchive ia (is);
      ia >> (*this);
    }



    template <int dim>
    void
    World<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particles");
        {
          prm.declare_entry ("Load balancing strategy", "repartition",
                             Patterns::MultipleSelection ("none|remove particles|add particles|"
                                                          "remove and add particles|repartition"),
                             "Strategy that is used to balance the computational "
                             "load across processors for adaptive meshes.");
          prm.declare_entry ("Minimum particles per cell", "0",
                             Patterns::Integer (0),
                             "Lower limit for particle number per cell. This limit is "
                             "useful for adaptive meshes to prevent fine cells from being empty "
                             "of particles. It will be checked and enforced after mesh "
                             "refinement and after particle movement. "
                             "If there are "
                             "\\texttt{n\\_number\\_of\\_particles} $<$ \\texttt{min\\_particles\\_per\\_cell} "
                             "particles in one cell then "
                             "\\texttt{min\\_particles\\_per\\_cell} - \\texttt{n\\_number\\_of\\_particles} "
                             "particles are generated and randomly placed in "
                             "this cell. If the particles carry properties the "
                             "individual property plugins control how the "
                             "properties of the new particles are initialized.");
          prm.declare_entry ("Maximum particles per cell", "100",
                             Patterns::Integer (0),
                             "Upper limit for particle number per cell. This limit is "
                             "useful for adaptive meshes to prevent coarse cells from slowing down "
                             "the whole model. It will be checked and enforced after mesh "
                             "refinement, after MPI transfer of particles and after particle "
                             "movement. If there are "
                             "\\texttt{n\\_number\\_of\\_particles} $>$ \\texttt{max\\_particles\\_per\\_cell} "
                             "particles in one cell then "
                             "\\texttt{n\\_number\\_of\\_particles} - \\texttt{max\\_particles\\_per\\_cell} "
                             "particles in this cell are randomly chosen and destroyed.");
          prm.declare_entry ("Particle weight", "10",
                             Patterns::Integer (0),
                             "Weight that is associated with the computational load of "
                             "a single particle. The sum of particle weights will be added "
                             "to the sum of cell weights to determine the partitioning of "
                             "the mesh if the `repartition' particle load balancing strategy "
                             "is selected. The optimal weight depends on the used "
                             "integrator and particle properties. In general for a more "
                             "expensive integrator and more expensive properties a larger "
                             "particle weight is recommended. Before adding the weights "
                             "of particles, each cell already carries a weight of 1000 to "
                             "account for the cost of field-based computations.");
          prm.declare_entry ("Update ghost particles", "false",
                             Patterns::Bool (),
                             "Some particle interpolation algorithms require knowledge "
                             "about particles in neighboring cells. To allow this, "
                             "particles in ghost cells need to be exchanged between the "
                             "processes neighboring this cell. This parameter determines "
                             "whether this transport is happening.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      Generator::declare_parameters<dim>(prm);
      Integrator::declare_parameters<dim>(prm);
      Interpolator::declare_parameters<dim>(prm);
      Property::Manager<dim>::declare_parameters(prm);
    }



    template <int dim>
    void
    World<dim>::parse_parameters (ParameterHandler &prm)
    {
      // First do some error checking. The current algorithm does not find
      // the cells around particles, if the particles moved more than one
      // cell in one timestep and we are running in parallel, because they
      // skip the layer of ghost cells around our local domain. Assert this
      // is not possible.
      const double CFL_number = prm.get_double ("CFL number");
      const unsigned int n_processes = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());

      AssertThrow((n_processes == 1) || (CFL_number <= 1.0),
                  ExcMessage("The current particle algorithm does not work in "
                             "parallel if the CFL number is larger than 1.0, because "
                             "in this case particles can move more than one cell "
                             "diameter in one time step and therefore skip the layer "
                             "of ghost cells around the local subdomain."));

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particles");
        {
          min_particles_per_cell = prm.get_integer("Minimum particles per cell");
          max_particles_per_cell = prm.get_integer("Maximum particles per cell");

          AssertThrow(min_particles_per_cell <= max_particles_per_cell,
                      ExcMessage("Please select a 'Minimum particles per cell' parameter "
                                 "that is smaller than or equal to the 'Maximum particles per cell' parameter."));

          particle_weight = prm.get_integer("Particle weight");

          update_ghost_particles = prm.get_bool("Update ghost particles");

          const std::vector<std::string> strategies = Utilities::split_string_list(prm.get ("Load balancing strategy"));
          AssertThrow(Utilities::has_unique_entries(strategies),
                      ExcMessage("The list of strings for the parameter "
                                 "'Postprocess/Particles/Load balancing strategy' contains entries more than once. "
                                 "This is not allowed. Please check your parameter file."));

          particle_load_balancing = ParticleLoadBalancing::no_balancing;

          for (std::vector<std::string>::const_iterator strategy = strategies.begin(); strategy != strategies.end(); ++strategy)
            {
              if (*strategy == "remove particles")
                particle_load_balancing = typename ParticleLoadBalancing::Kind(particle_load_balancing | ParticleLoadBalancing::remove_particles);
              else if (*strategy == "add particles")
                particle_load_balancing = typename ParticleLoadBalancing::Kind(particle_load_balancing | ParticleLoadBalancing::add_particles);
              else if (*strategy == "remove and add particles")
                particle_load_balancing = typename ParticleLoadBalancing::Kind(particle_load_balancing | ParticleLoadBalancing::remove_and_add_particles);
              else if (*strategy == "repartition")
                particle_load_balancing = typename ParticleLoadBalancing::Kind(particle_load_balancing | ParticleLoadBalancing::repartition);
              else if (*strategy == "none")
                {
                  particle_load_balancing = ParticleLoadBalancing::no_balancing;
                  AssertThrow(strategies.size() == 1,
                              ExcMessage("The particle load balancing strategy `none' is not compatible "
                                         "with any other strategy, yet it seems another is selected as well. "
                                         "Please check the parameter file."));
                }
              else
                AssertThrow(false,
                            ExcMessage("The 'Load balancing strategy' parameter contains an unknown value: <" + *strategy
                                       + ">. This value does not correspond to any known load balancing strategy. Possible values "
                                       "are listed in the corresponding manual subsection."));
            }

        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      TimerOutput::Scope timer_section(this->get_computing_timer(), "Particles: Initialization");

      // Create a generator object depending on what the parameters specify
      generator.reset(Generator::create_particle_generator<dim> (prm));
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(generator.get()))
        sim->initialize_simulator (this->get_simulator());
      generator->parse_parameters(prm);
      generator->initialize();

      // Create a property_manager object and initialize its properties
      property_manager = std::make_unique<Property::Manager<dim>> ();
      SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(property_manager.get());
      sim->initialize_simulator (this->get_simulator());
      property_manager->parse_parameters(prm);
      property_manager->initialize();

      // Create an integrator object depending on the specified parameter
      integrator.reset(Integrator::create_particle_integrator<dim> (prm));
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(integrator.get()))
        sim->initialize_simulator (this->get_simulator());
      integrator->parse_parameters(prm);
      integrator->initialize();

      // Create an interpolator object depending on the specified parameter
      interpolator.reset(Interpolator::create_particle_interpolator<dim> (prm));
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(interpolator.get()))
        sim->initialize_simulator (this->get_simulator());
      interpolator->parse_parameters(prm);
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Particle
  {
#define INSTANTIATE(dim) \
  template class World<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
