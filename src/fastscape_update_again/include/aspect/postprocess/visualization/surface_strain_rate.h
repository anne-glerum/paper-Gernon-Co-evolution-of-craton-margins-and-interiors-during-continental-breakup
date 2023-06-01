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


#ifndef _aspect_postprocess_visualization_surface_strain_rate_h
#define _aspect_postprocess_visualization_surface_strain_rate_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessorScalar that computes a
       * variable that represents the second invariant of the strain
       * rate, and outputs it on all boundary faces.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class SurfaceStrainRate
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>,
          public SurfaceOnlyVisualization<dim>
      {
        public:
          SurfaceStrainRate();

          /**
           * Evaluate the strain rate for the current cell.
           *
           * @copydoc DataPostprocessorScalar<dim>::evaluate_vector_field()
           */
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double> > &computed_quantities) const override;
      };
    }
  }
}

#endif
