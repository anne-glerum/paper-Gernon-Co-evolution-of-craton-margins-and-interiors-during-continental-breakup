/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_boundary_temperature_interface_h
#define _aspect_boundary_temperature_interface_h

#include <aspect/plugins.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/distributed/tria.h>

#include <boost/core/demangle.hpp>
#include <typeinfo>


namespace aspect
{
  template <int dim> class SimulatorAccess;

  /**
   * A namespace for the definition of things that have to do with describing
   * the boundary values for the temperature.
   *
   * @ingroup BoundaryTemperatures
   */
  namespace BoundaryTemperature
  {
    using namespace dealii;

    /**
     * Base class for classes that describe temperature boundary values.
     *
     * @ingroup BoundaryTemperatures
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
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        virtual void initialize ();

        /**
         * Return the temperature that is to hold at a particular position on
         * the boundary of the domain.
         *
         * @param boundary_indicator The boundary indicator of the part of the
         * boundary of the domain on which the point is located at which we
         * are requesting the temperature.
         * @param position The position of the point at which we ask for the
         * temperature.
         * @return Boundary temperature at position @p position.
         */
        virtual
        double boundary_temperature (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position) const;

        /**
         * Return the minimal temperature on that part of the boundary on
         * which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        virtual
        double minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids
                                    = std::set<types::boundary_id>()) const = 0;

        /**
         * Return the maximal temperature on that part of the boundary on
         * which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        virtual
        double maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids
                                    = std::set<types::boundary_id>()) const = 0;

        /**
         * A function that is called at the beginning of each time step. The
         * default implementation of the function does nothing, but derived
         * classes that need more elaborate setups for a given time step may
         * overload the function.
         *
         * The point of this function is to allow complex boundary temperature
         * models to do an initialization step once at the beginning of each
         * time step. An example would be a model that needs to call an
         * external program to compute temperature change at bottom.
         */
        virtual
        void
        update ();

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
     * A class that manages all boundary temperature objects.
     *
     * @ingroup BoundaryTemperatures
     */
    template <int dim>
    class Manager : public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Destructor. Made virtual since this class has virtual member
         * functions.
         */
        ~Manager () override;

        /**
         * A function that is called at the beginning of each time step and
         * calls the corresponding functions of all created plugins.
         *
         * The point of this function is to allow complex boundary temperature
         * models to do an initialization step once at the beginning of each
         * time step. An example would be a model that needs to call an
         * external program to compute the temperature change at a boundary.
         */
        virtual
        void
        update ();

        /**
         * Declare the parameters of all known boundary temperature plugins, as
         * well as the ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which boundary temperature objects will be created;
         * then let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * A function that calls the boundary_temperature functions of all the
         * individual boundary temperature objects and uses the stored operators
         * to combine them.
         */
        double
        boundary_temperature (const types::boundary_id boundary_indicator,
                              const Point<dim> &position) const;

        /**
         * Return the minimal temperature of all selected plugins on that
         * part of the boundary on which Dirichlet conditions are posed.
         */
        double minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids
                                    = std::set<types::boundary_id>()) const;

        /**
         * Return the maximal temperature of all selected plugins on that
         * part of the boundary on which Dirichlet conditions are posed.
         */
        double maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids
                                    = std::set<types::boundary_id>()) const;

        /**
         * A function that is used to register boundary temperature objects in such
         * a way that the Manager can deal with all of them without having to
         * know them by name. This allows the files in which individual
         * plugins are implemented to register these plugins, rather than also
         * having to modify the Manager class by adding the new boundary
         * temperature plugin class.
         *
         * @param name A string that identifies the boundary temperature model
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this boundary temperature model
         * wants to read from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this boundary temperature model.
         */
        static
        void
        register_boundary_temperature (const std::string &name,
                                       const std::string &description,
                                       void (*declare_parameters_function) (ParameterHandler &),
                                       Interface<dim> *(*factory_function) ());


        /**
         * Return a list of names of all boundary temperature models currently
         * used in the computation, as specified in the input file.
         */
        const std::vector<std::string> &
        get_active_boundary_temperature_names () const;

        /**
         * Return a list of pointers to all boundary temperature models
         * currently used in the computation, as specified in the input file.
         */
        const std::vector<std::unique_ptr<Interface<dim>>> &
        get_active_boundary_temperature_conditions () const;

        /**
         * Go through the list of all boundary temperature models that have been selected in
         * the input file (and are consequently currently active) and see if one
         * of them has the desired type specified by the template argument. If so,
         * return a pointer to it. If no boundary temperature model is active
         * that matches the given type, return a nullptr.
         */
        template <typename BoundaryTemperatureType>
        DEAL_II_DEPRECATED
        BoundaryTemperatureType *
        find_boundary_temperature_model () const;

        /**
         * Go through the list of all boundary temperature models that have been selected
         * in the input file (and are consequently currently active) and return
         * true if one of them has the desired type specified by the template
         * argument.
         */
        template <typename BoundaryTemperatureType>
        bool
        has_matching_boundary_temperature_model () const;

        /**
         * Go through the list of all boundary temperature models that have been selected
         * in the input file (and are consequently currently active) and see
         * if one of them has the type specified by the template
         * argument or can be casted to that type. If so, return a reference
         * to it. If no boundary temperature model is active that matches the given type,
         * throw an exception.
         */
        template <typename BoundaryTemperatureType>
        const BoundaryTemperatureType &
        get_matching_boundary_temperature_model () const;

        /*
         * Return a set of boundary indicators for which boundary
         * temperatures are prescribed.
         */
        const std::set<types::boundary_id> &
        get_fixed_temperature_boundary_indicators() const;

        /*
         * Return whether Dirichlet boundary conditions will be applied
         * on parts of the boundaries where material flows out.
         */
        bool
        allows_fixed_temperature_on_outflow_boundaries() const;

        /**
         * For the current plugin subsystem, write a connection graph of all of the
         * plugins we know about, in the format that the
         * programs dot and neato understand. This allows for a visualization of
         * how all of the plugins that ASPECT knows about are interconnected, and
         * connect to other parts of the ASPECT code.
         *
         * @param output_stream The stream to write the output to.
         */
        static
        void
        write_plugin_graph (std::ostream &output_stream);


        /**
         * Exception.
         */
        DeclException1 (ExcBoundaryTemperatureNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered boundary temperature objects.");
      private:
        /**
         * A list of boundary temperature objects that have been requested in the
         * parameter file.
         */
        std::vector<std::unique_ptr<Interface<dim>>> boundary_temperature_objects;

        /**
         * A list of names of boundary temperature objects that have been requested
         * in the parameter file.
         */
        std::vector<std::string> model_names;

        /**
         * A list of enums of boundary temperature operators that have been
         * requested in the parameter file. Each name is associated
         * with a model_name, and is used to modify the temperature
         * boundary with the values from the current plugin.
         */
        std::vector<aspect::Utilities::Operator> model_operators;

        /**
         * A set of boundary ids on which the boundary_temperature_objects
         * will be applied.
         */
        std::set<types::boundary_id> fixed_temperature_boundary_indicators;

        /**
         * Whether we allow the temperature to be fixed on parts of the boundary
         * where material flows out of the domain.
         */
        bool allow_fixed_temperature_on_outflow_boundaries;
    };



    template <int dim>
    template <typename BoundaryTemperatureType>
    inline
    BoundaryTemperatureType *
    Manager<dim>::find_boundary_temperature_model () const
    {
      for (const auto &p : boundary_temperature_objects)
        if (BoundaryTemperatureType *x = dynamic_cast<BoundaryTemperatureType *> ( p.get()) )
          return x;
      return nullptr;
    }


    template <int dim>
    template <typename BoundaryTemperatureType>
    inline
    bool
    Manager<dim>::has_matching_boundary_temperature_model () const
    {
      for (const auto &p : boundary_temperature_objects)
        if (Plugins::plugin_type_matches<BoundaryTemperatureType>(*p))
          return true;
      return false;
    }


    template <int dim>
    template <typename BoundaryTemperatureType>
    inline
    const BoundaryTemperatureType &
    Manager<dim>::get_matching_boundary_temperature_model () const
    {
      AssertThrow(has_matching_boundary_temperature_model<BoundaryTemperatureType> (),
                  ExcMessage("You asked BoundaryTemperature::Manager::get_boundary_temperature_model() for a "
                             "boundary temperature model of type <" + boost::core::demangle(typeid(BoundaryTemperatureType).name()) + "> "
                             "that could not be found in the current model. Activate this "
                             "boundary temperature model in the input file."));

      for (auto &p : boundary_temperature_objects)
        if (Plugins::plugin_type_matches<BoundaryTemperatureType>(*p))
          return Plugins::get_plugin_as_type<BoundaryTemperatureType>(*p);

      // We will never get here, because we had the Assert above. Just to avoid warnings.
      typename std::vector<std::unique_ptr<Interface<dim>>>::const_iterator boundary_temperature_model;
      return Plugins::get_plugin_as_type<BoundaryTemperatureType>(*(*boundary_temperature_model));
    }


    /**
     * Return a string that consists of the names of boundary temperature models that can
     * be selected. These names are separated by a vertical line '|' so
     * that the string can be an input to the deal.II classes
     * Patterns::Selection or Patterns::MultipleSelection.
     */
    template <int dim>
    std::string
    get_valid_model_names_pattern ();


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a boundary temperature model, register it with the functions that
     * can declare their parameters and create these objects.
     *
     * @ingroup BoundaryTemperatures
     */
#define ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(classname, name, description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryTemperature::Interface<2>,classname<2>> \
        dummy_ ## classname ## _2d (&aspect::BoundaryTemperature::Manager<2>::register_boundary_temperature, \
                                    name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::BoundaryTemperature::Interface<3>,classname<3>> \
        dummy_ ## classname ## _3d (&aspect::BoundaryTemperature::Manager<3>::register_boundary_temperature, \
                                    name, description); \
  }
  }
}


#endif
