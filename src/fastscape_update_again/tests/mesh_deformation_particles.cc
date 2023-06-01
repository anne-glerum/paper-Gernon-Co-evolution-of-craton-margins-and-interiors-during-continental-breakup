#include <aspect/simulator_signals.h>
#include <aspect/simulator_access.h>
#include <aspect/particle/world.h>
#include <deal.II/particles/particle_handler.h>

#include <iostream>

using namespace aspect;


template <int dim>
void post_mesh_deformation (const SimulatorAccess<dim> &simulator_access)
{
  // Get the reference location of the one particle,
  // and check that it has the correct value.
  const Point<dim> predicted_particle_reference_location =
    simulator_access.get_particle_world().get_particle_handler().begin()->get_reference_location();
  // At timestep 0, the mesh is not deformed
  // other than through the initial topography/
  // initial mesh deformation plugins, which
  // are not active in this test.
  // At timestep 1, the surface of the mesh is
  // upwardly displaced by dt*v_surf = 0.5*1.0 = 0.5.
  // The new surface lies therefore at 1.5.
  // Because the particle does not move (Stokes velocity
  // is set to zero), the vertical coordinate of the particle
  // reference location within the one cell of the mesh
  // should be updated from 0.5/1.0 to 0.5/1.5 = 1/3.
  Point<dim> correct_particle_reference_location;
  correct_particle_reference_location[0] = 0.5;
  if (simulator_access.get_timestep_number() == 0)
    correct_particle_reference_location[dim-1] = 0.5;
  else if (simulator_access.get_timestep_number() == 1)
    correct_particle_reference_location[dim-1] = 1./3.;
  else
    Assert(false, ExcNotImplemented());

  AssertThrow(predicted_particle_reference_location == correct_particle_reference_location,
              ExcMessage("The particle reference location is not correct."));
}


template <int dim>
void signal_connector (SimulatorSignals<dim> &signals)
{
  signals.post_mesh_deformation.connect (&post_mesh_deformation<dim>);
}


ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)
