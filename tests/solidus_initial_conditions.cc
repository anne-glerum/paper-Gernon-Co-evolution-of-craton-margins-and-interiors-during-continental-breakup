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


#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/boundary_temperature/interface.h>

#include <cmath>

namespace aspect
{
  namespace InitialTemperature
  {
    /**
     * Data class to handle the melting curve.
     */
    class MeltingCurve
    {
      public:

        /**
         * Read the data file into the class.
         */
        void read(const std::string &filename,
                  const MPI_Comm &comm);

        /**
         * Get the melting temperature.
         */
        double T(const double p, const double radius) const;

        /**
         * Is the melting curve dependent on radius. The melting curve can be
         * dependent on radius or pressure.
         */
        bool is_radius;

        /**
         * Number of data points in the melting curve data.
         */
        unsigned int n_points;
      private:
        /**
         * Data array for temperature.
         */
        std::vector<double> T_array;

        /**
         * Data array for pressure/radius.
         */
        std::vector<double> P_or_R_array;

        /**
         * Name of the data file.
         */
        std::string data_filename;
    };


    /**
     * A class that implements temperature initial conditions based on solidus
     * provided by a data file.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class Solidus : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Solidus ();

        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Returns spherical coordinates of a cartesian position.
         */
        const Tensor<1,dim>
        spherical_surface_coordinates(const Tensor<1,dim> &position) const;

        /**
         * Lithosphere thickness.
         */
        double       litho_thick;

        /**
         * Magnitude of temperature perturbation.
         */
        double       magnitude_T;

        /**
         * Magnitude of lithosphere thickness perturbation.
         */
        double       magnitude_lith;

        /**
         * Temperature difference from solidus, so the initial condition can
         * be super-solidus or sub-solidus.
         */
        double       deltaT;

        /**
         * The lateral wave number  of the harmonic perturbation in the first
         * dimension. This is the only lateral wave number in 2D and equals
         * the degree of the spherical harmonics in a 3D spherical shell.
         */
        int          lateral_wave_number_1;

        /**
         * The lateral wave number of the harmonic perturbation in the second
         * dimension. This is not used in 2D and equals the order of the
         * spherical harmonics in a 3D spherical shell.
         */
        int          lateral_wave_number_2;

        /**
         * The data file name for solidus data.
         */
        std::string  solidus_filename;

        /**
         * Data class for melting curve
         */
        MeltingCurve solidus_curve;
    };

    void MeltingCurve::read(const std::string &filename,
                            const MPI_Comm &comm)
    {
      data_filename=filename;
      // Read data from disk and distribute among processes
      std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

      std::string dummy;

      std::string T_Unit,P_Unit;
      n_points=0;
      getline (in, dummy);
      in>>T_Unit>>P_Unit;
      getline (in, dummy);
      while (!in.eof())
        {
          double T,p;
          in>>T>>p;
          if (!in.fail())
            {
              // Unit switching
              if (T_Unit=="C")
                {
                  T+=273.15;                      // Degree C to K
                }
              else if (T_Unit!="K")
                {
                  AssertThrow(false,ExcMessage ("Unit of the first column of melting curve data "
                                                "has to be one of the following: C/K."))
                }

              is_radius=false;                    // Second column is pressure
              if (P_Unit=="kbar")
                p*=1.e8;                        // kbar to Pa
              else if (P_Unit=="GPa")
                p*=1.e9;                        // GPa to Pa
              else
                {
                  is_radius=true;                 // Second column is radius instead of pressure
                  if (P_Unit=="km")
                    p*=1.e3;                    // km to meters
                  else if (P_Unit!="m")
                    AssertThrow(false,ExcMessage ("Unit of the second column of melting curve data "
                                                  "has to be one of the following: Pa/GPa/km/m."))
                  }
              T_array.push_back(T);
              P_or_R_array.push_back(p);
              n_points++;
            }
          getline (in, dummy);
        }
    }



    double MeltingCurve::T(const double p, const double radius) const
    {
      double T_value,P_or_R_value=is_radius?radius:p;
      if (T_array.size()==0)return (0);
      for (unsigned i=1; i<n_points; i++)
        {
          if (     (i==n_points-1) ||
                   (is_radius && P_or_R_value>P_or_R_array[i]) ||
                   (!is_radius && P_or_R_value<P_or_R_array[i]) )
            {
              T_value=T_array[i-1]+(T_array[i]-T_array[i-1])/(P_or_R_array[i]-P_or_R_array[i-1])*(P_or_R_value-P_or_R_array[i-1]);
              return (T_value);
            }
        }
      AssertThrow(false,ExcMessage(std::string("Something wrong with the melting curve data ")+ data_filename ));
      return (-1);
    }



    template <int dim>
    Solidus<dim>::Solidus ():
      solidus_curve()
    {}

    template <int dim>
    double
    Solidus<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      double T_min,T_litho;
      double T_solidus,T_perturbation;
      double litho_thick_theta;
      double lateral_perturbation;
      double Depth=this->get_geometry_model().depth(position);

      AssertThrow(solidus_curve.n_points!=0,ExcMessage("Error reading solidus file."));
      AssertThrow(solidus_curve.is_radius==true,ExcMessage("The solidus curve has to be radius dependent."));
      const GeometryModel::SphericalShell<dim> *spherical_geometry_model=
        dynamic_cast< const GeometryModel::SphericalShell<dim> *>(&this->get_geometry_model());

      AssertThrow(spherical_geometry_model!=0,
                  ExcMessage("This initial condition can only be used with spherical shell geometry model."));

      AssertThrow(this->has_boundary_temperature(),
                  ExcMessage("This initial condition can only be used with a prescribed boundary temperature."));

      T_min=(this->get_boundary_temperature_manager()).minimal_temperature();

      // In case of spherical shell calculate spherical coordinates
      const Tensor<1,dim> scoord = spherical_surface_coordinates(position);
      if (dim==2)
        {
          // Use a sine as lateral perturbation that is scaled to the opening angle of the geometry.
          // This way the perturbation is always 0 at the model boundaries.
          const double opening_angle = spherical_geometry_model->opening_angle()*numbers::PI/180.0;
          lateral_perturbation = std::sin(lateral_wave_number_1*scoord[1]*numbers::PI/opening_angle);
        }
      else if (dim==3)
        {
          // Spherical harmonics are only defined for order <= degree
          // and degree >= 0. Verify that it is indeed.
          Assert ( std::abs(lateral_wave_number_2) <= lateral_wave_number_1,
                   ExcMessage ("Spherical harmonics can only be computed for "
                               "order <= degree."));
          Assert ( lateral_wave_number_1 >= 0,
                   ExcMessage ("Spherical harmonics can only be computed for "
                               "degree >= 0."));
          // use a spherical harmonic function as lateral perturbation
          std::pair<double,double> sph_harm_vals = Utilities::real_spherical_harmonic( lateral_wave_number_1, lateral_wave_number_2, scoord[2], scoord[1] );
          lateral_perturbation = sph_harm_vals.first;
        }
      litho_thick_theta=litho_thick-magnitude_lith*lateral_perturbation;
      T_litho=solidus_curve.T(0,spherical_geometry_model->outer_radius()-litho_thick_theta)+deltaT;

      if (litho_thick_theta>0 && Depth<litho_thick_theta)
        T_solidus=T_min+(T_litho-T_min)*(Depth/litho_thick_theta);
      else
        T_solidus=solidus_curve.T(0,sqrt(position.square()))+deltaT;

      T_perturbation=Depth/( this->get_geometry_model().maximal_depth() )*magnitude_T*lateral_perturbation;
      return T_solidus+T_perturbation;
    }



    template <int dim>
    const Tensor<1,dim>
    Solidus<dim>::
    spherical_surface_coordinates(const Tensor<1,dim> &position) const
    {
      Tensor<1,dim> scoord;

      scoord[0] = std::sqrt(position.norm_square());                // R
      scoord[1] = std::atan2(position[1],position[0]);              // Phi
      if (scoord[1] < 0.0) scoord[1] = 2*numbers::PI + scoord[1];   // correct phi to [0,2*pi]
      if (dim==3)
        scoord[2] = std::acos(position[2]/std::sqrt(position.norm_square())); // Theta

      return scoord;
    }



    template <int dim>
    void
    Solidus<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Solidus");
        {
          prm.declare_entry ("Supersolidus","0e0",
                             Patterns::Double (),
                             "The difference from solidus, use this number to generate initial conditions "
                             "that close to solidus instead of exactly at solidus. Use small negative number"
                             " in this parameter to prevent large melting generation at the beginning. "
                             "  Units: \\si{\\kelvin} ");
          prm.declare_entry ("Lithosphere thickness","0",
                             Patterns::Double (0),
                             "The thickness of lithosphere thickness. Units: \\si{\\meter}");
          prm.enter_subsection("Perturbation");
          {
            prm.declare_entry ("Temperature amplitude", "0e0",
                               Patterns::Double (0),
                               "The amplitude of the initial spherical temperature perturbation in (K)");
            prm.declare_entry ("Lithosphere thickness amplitude", "0e0",
                               Patterns::Double (),
                               "The amplitude of the initial lithosphere thickness perturbation in (m)");
            prm.declare_entry ("Lateral wave number one","3",
                               Patterns::Integer(),
                               "Doubled first lateral wave number of the harmonic perturbation. "
                               "Equals the spherical harmonic degree in 3D spherical shells. "
                               "In all other cases one equals half of a sine period over "
                               "the model domain. This allows for single up-/downswings. "
                               "Negative numbers reverse the sign of the perturbation but are "
                               "not allowed for the spherical harmonic case.");
            prm.declare_entry ("Lateral wave number two", "2",
                               Patterns::Integer (),
                               "Doubled second lateral wave number of the harmonic perturbation. "
                               "Equals the spherical harmonic order in 3D spherical shells. "
                               "In all other cases one equals half of a sine period over "
                               "the model domain. This allows for single up-/downswings. "
                               "Negative numbers reverse the sign of the perturbation.");
          }
          prm.leave_subsection();
          prm.enter_subsection ("Data");
          {
            prm.declare_entry ("Solidus filename", "",
                               Patterns::Anything(),
                               "The solidus data filename. It is a function of radius or pressure "
                               "in the following format: \n"
                               "Line 1:  Header \n"
                               "Line 2:  Unit of temperature (C/K)        "
                               "Unit of pressure (GPa/kbar) or radius (km/m) \n"
                               "Line 3~: Column of solidus temperature    Column of radius/pressure \n"
                               "See data/initial-temperature/solidus.Mars as an example."
                               "\n\n"
                               "In order to facilitate placing input files in locations relative "
                               "to the ASPECT source directory, the file name may also include the special "
                               "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                               "in which the ASPECT source files were located when ASPECT was "
                               "compiled. This interpretation allows, for example, to reference "
                               "files located in the `data/' subdirectory of ASPECT.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Solidus<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Solidus");
        {
          deltaT=prm.get_double("Supersolidus");
          litho_thick=prm.get_double("Lithosphere thickness");
          prm.enter_subsection("Perturbation");
          {
            magnitude_T    = prm.get_double("Temperature amplitude");
            magnitude_lith = prm.get_double("Lithosphere thickness amplitude");
            lateral_wave_number_1 = prm.get_integer ("Lateral wave number one");
            lateral_wave_number_2 = prm.get_integer ("Lateral wave number two");
          }
          prm.leave_subsection();
          prm.enter_subsection("Data");
          solidus_filename = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Solidus filename"));
          solidus_curve.read(solidus_filename,this->get_mpi_communicator());
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }


}




// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(Solidus,
                                              "solidus",
                                              "This is a temperature initial condition that "
                                              "starts the model close to solidus, it also contains "
                                              "a user defined lithosphere thickness and with perturbations "
                                              " in both lithosphere thickness and temperature based on "
                                              "spherical harmonic functions. It was used as the initial condition "
                                              "of early Mars after the freezing of the magma ocean, "
                                              "using the solidus from Parmentier et al., "
                                              "Melt-solid segregation, Fractional magma ocean solidification, and implications for "
                                              "long term planetary evolution. Luna and Planetary Science, 2007.")
  }
}
