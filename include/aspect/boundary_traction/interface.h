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


#ifndef _aspect_boundary_traction_interface_h
#define _aspect_boundary_traction_interface_h

#include <aspect/plugins.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  /**
   * A namespace in which we define everything that has to do with defining
   * traction boundary conditions for the Stokes equations.
   *
   * @ingroup BoundaryTractions
   */
  namespace BoundaryTraction
  {
    using namespace dealii;

    /**
     * A base class for parameterizations of traction boundary conditions.
     *
     * @ingroup BoundaryTractions
     */
    template <int dim>
    class Interface
    {
      public:
        /**
         * Destructor. Made virtual to enforce that derived classes also have
         * virtual destructors.
         */
        virtual ~Interface() = default;

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after the
         * SimulatorAccess (if applicable) is initialized.
         */
        virtual
        void
        initialize ();

        /**
         * A function that is called at the beginning of each time step.
         * The default implementation of the function does nothing, but
         * derived classes that need more elaborate setups for a given time
         * step may overload the function.
         *
         * The point of this function is to allow complex boundary traction
         * models to do an initialization step once at the beginning of each
         * time step. An example would be a model that needs to call an
         * external program to compute stresses for a set of plates.
         */
        virtual
        void
        update ();

        /**
         * Return the boundary traction as a function of position.
         *
         * @deprecated Use boundary_traction(const types::boundary_id boundary_indicator,
         * const Point<dim> &position, const Tensor<1,dim> &normal_vector) const instead.
         */
        DEAL_II_DEPRECATED
        virtual
        Tensor<1,dim>
        traction (const Point<dim> &position,
                  const Tensor<1,dim> &normal_vector) const;

        /**
         * Return the boundary traction as a function of position. The
         * (outward) normal vector to the domain is also provided as
         * a second argument.
         */
        virtual
        Tensor<1,dim>
        boundary_traction (const types::boundary_id boundary_indicator,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const;

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

      protected:
        /**
         * Pointer to the geometry object in use.
         */
        const GeometryModel::Interface<dim> *geometry_model;
    };




    /**
     * Register a traction boundary conditions model so that it can be
     * selected from the parameter file.
     *
     * @param name A string that identifies the traction boundary conditions
     * model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this traction boundary conditions
     * model wants to read from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this traction boundary conditions model.
     *
     * @ingroup BoundaryTractions
     */
    template <int dim>
    void
    register_boundary_traction (const std::string &name,
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
     * @ingroup BoundaryTractions
     */
    template <int dim>
    Interface<dim> *
    create_boundary_traction (const std::string &name);

    /**
     * Return a list of names of all implemented boundary traction models,
     * separated by '|' so that it can be used in an object of type
     * Patterns::Selection.
     */
    template <int dim>
    std::string
    get_names ();

    /**
     * Declare the runtime parameters of the registered traction boundary
     * conditions models.
     *
     * @ingroup BoundaryTractions
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
     * for a traction boundary conditions model, register it with the
     * functions that can declare their parameters and create these objects.
     *
     * @ingroup BoundaryTractions
     */
#define ASPECT_REGISTER_BOUNDARY_TRACTION_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_BOUNDARY_TRACTION_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryTraction::Interface<2>,classname<2>> \
        dummy_ ## classname ## _2d (&aspect::BoundaryTraction::register_boundary_traction<2>, \
                                    name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryTraction::Interface<3>,classname<3>> \
        dummy_ ## classname ## _3d (&aspect::BoundaryTraction::register_boundary_traction<3>, \
                                    name, description); \
  }
  }
}


#endif
