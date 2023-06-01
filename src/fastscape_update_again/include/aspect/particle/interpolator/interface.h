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

#ifndef _aspect_particle_interpolator_interface_h
#define _aspect_particle_interpolator_interface_h

#include <aspect/plugins.h>
#include <aspect/global.h>

#include <deal.II/particles/particle.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/component_mask.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      using namespace dealii;
      using namespace dealii::Particles;

      /**
       * An abstract class defining virtual methods for performing
       * interpolation of particle properties to arbitrary points.
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      class Interface
      {
        public:
          /**
           * Destructor. Made virtual so that derived classes can be created
           * and destroyed through pointers to the base class.
           */
          virtual ~Interface () = default;

          /**
           * Perform an interpolation of the properties of the particles in
           * this cell onto a vector of positions in this cell.
           * Implementations of this function must return a vector of a vector
           * of doubles which contains a somehow computed
           * value of all particle properties at all given positions.
           *
           * @param [in] particle_handler Reference to the particle handler
           * that allows accessing the particles in the domain.
           * @param [in] positions The vector of positions where the properties
           * should be evaluated.
           * @param [in] cell An optional iterator to the cell containing the
           * particles. Not all callers will know the cell of the particles,
           * but providing the cell when known speeds up the interpolation
           * significantly.
           * @return A vector with as many entries as @p positions. Every entry
           * is a vector of interpolated particle properties at this position.
           */
          DEAL_II_DEPRECATED
          virtual
          std::vector<std::vector<double>>
                                        properties_at_points(const ParticleHandler<dim> &particle_handler,
                                                             const std::vector<Point<dim>> &positions,
                                                             const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell = typename parallel::distributed::Triangulation<dim>::active_cell_iterator()) const;

          /**
           * Perform an interpolation of the properties of the particles in
           * this cell onto a vector of positions in this cell.
           * Implementations of this function must return a vector of a vector
           * of doubles with as many entries as positions in @p positions.
           * Each entry is a vector with as many entries as there are  particle
           * properties in this computation.
           * All in @p selected_properties selected components
           * will be filled with computed properties, all other components
           * are not filled (or filled with invalid values).
           *
           * @param [in] particle_handler Reference to the particle handler
           * that allows accessing the particles in the domain.
           * @param [in] positions The vector of positions where the properties
           * should be evaluated.
           * @param [in] selected_properties A component mask that determines
           * which particle properties are interpolated in this function.
           * @param [in] cell An optional iterator to the cell containing the
           * particles. Not all callers will know the cell of the particles,
           * but providing the cell when known speeds up the interpolation
           * significantly.
           * @return A vector with as many entries as @p positions. Every entry
           * is a vector of interpolated particle properties at this position.
           */
          virtual
          std::vector<std::vector<double>>
                                        properties_at_points(const ParticleHandler<dim> &particle_handler,
                                                             const std::vector<Point<dim>> &positions,
                                                             const ComponentMask &selected_properties,
                                                             const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell = typename parallel::distributed::Triangulation<dim>::active_cell_iterator()) const = 0;

          /**
           * Declare the parameters this class takes through input files. The
           * default implementation of this function does not describe any
           * parameters. Consequently, derived classes do not have to overload
           * this function if they do not take any runtime parameters.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           * The default implementation of this function does not read any
           * parameters. Consequently, derived classes do not have to overload
           * this function if they do not take any runtime parameters.
           */
          virtual
          void
          parse_parameters (ParameterHandler &prm);
      };


      /**
       * Return a list of names (separated by '|') of possible interpolator
       * classes for particles.
       */
      std::string
      interpolator_object_names ();


      /**
       * Register a particle interpolator so that it can be selected from
       * the parameter file.
       *
       * @param name A string that identifies the particle interpolator
       * @param description A text description of what this interpolator does and that
       * will be listed in the documentation of the parameter file.
       * @param declare_parameters_function A pointer to a function that can be
       * used to declare the parameters that this particle interpolator wants to read
       * from input files.
       * @param factory_function A pointer to a function that can create an
       * object of this particle interpolator.
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      void
      register_particle_interpolator (const std::string &name,
                                      const std::string &description,
                                      void (*declare_parameters_function) (ParameterHandler &),
                                      Interface<dim> *(*factory_function) ());

      /**
       * A function that given the name of a model returns a pointer to an
       * object that describes it. Ownership of the pointer is transferred to
       * the caller.
       *
       * The model object returned is not yet initialized and has not
       * read its runtime parameters yet.
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      Interface<dim> *
      create_particle_interpolator (ParameterHandler &prm);


      /**
       * Declare the runtime parameters of the registered particle interpolators.
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      void
      declare_parameters (ParameterHandler &prm);


      /**
       * For the current plugin subsystem, write a connection graph of all of the
       * plugins we know about, in the format that the
       * programs dot and neato understand. This allows for a visualization of
       * how all of the plugins that ASPECT knows about are interconnected, and
       * connect to other parts of the ASPECT code.
       *
       * @param output_stream The stream to write the output to.
       */
      template <int dim>
      void
      write_plugin_graph (std::ostream &output_stream);

      /**
       * Given a class name, a name, and a description for the parameter file
       * for a particle interpolator, register it with the functions that
       * can declare their parameters and create these objects.
       *
       * @ingroup ParticleInterpolators
       */
#define ASPECT_REGISTER_PARTICLE_INTERPOLATOR(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_PARTICLE_INTERPOLATOR_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Interpolator::Interface<2>,classname<2>> \
        dummy_ ## classname ## _2d (&aspect::Particle::Interpolator::register_particle_interpolator<2>, \
                                    name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Particle::Interpolator::Interface<3>,classname<3>> \
        dummy_ ## classname ## _3d (&aspect::Particle::Interpolator::register_particle_interpolator<3>, \
                                    name, description); \
  }
    }
  }
}


#endif
