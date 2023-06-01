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


#ifndef _aspect_postprocess_stokes_residual_h
#define _aspect_postprocess_stokes_residual_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/data_out_base.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that outputs Stokes residuals to a file stokes_residuals.txt.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class StokesResidual : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        StokesResidual ();

        /**
         * This attaches to the Stokes solver signal.
         */
        void initialize() override;

        /**
         * Generate the output file.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Save the state of this object.
         */
        void save (std::map<std::string, std::string> &status_strings) const override;

        /**
         * Restore the state of the object.
         */
        void load (const std::map<std::string, std::string> &status_strings) override;

        /**
         * Serialize the contents of this class as far as they are not read
         * from input parameter files.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

      private:
        /**
         * A structure for a single time step record.
         */
        struct DataPoint
        {
          double time;
          unsigned int solve_index;
          std::vector<double> values;

          template <class Archive>
          void serialize (Archive &ar, const unsigned int version);
        };

        /**
         * Callback function to collect the data.
         */
        void stokes_solver_callback (const SolverControl &solver_control_cheap,
                                     const SolverControl &solver_control_expensive);

        /**
         * An array of all the past values
         */
        std::vector<DataPoint> entries;

    };
  }
}


#endif
