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


#ifndef _aspect_boundary_velocity_zero_velocity_h
#define _aspect_boundary_velocity_zero_velocity_h

#include <aspect/boundary_velocity/interface.h>

namespace aspect
{
  namespace BoundaryVelocity
  {
    using namespace dealii;

    /**
     * A class that implements zero velocity (stick) boundary conditions.
     *
     * @ingroup BoundaryVelocities
     */
    template <int dim>
    class ZeroVelocity : public Interface<dim>
    {
      public:
        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function obviously simply returns a zero
         * tensor.
         */
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
                           const Point<dim> &position) const override;

        // avoid -Woverloaded-virtual warning until the deprecated function
        // is removed from the interface:
        using Interface<dim>::boundary_velocity;
    };
  }
}


#endif
