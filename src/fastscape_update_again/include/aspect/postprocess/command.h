/*
  Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_command_h
#define _aspect_postprocess_command_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <stdlib.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that runs a shell command at the end of each time step.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Command : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Declare new parameters
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Execute the command
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

      private:
        std::string command;
        bool terminate_on_failure;
        bool on_all_processes;
    };
  }
}


#endif
