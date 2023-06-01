/*
 Copyright (C) 2015 - 2020 by the authors of the ASPECT code.

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


#include <aspect/simulator.h>
#include <aspect/mesh_deformation/free_surface.h>
#include <aspect/utilities.h>
#include <aspect/postprocess/geoid.h>
#include <aspect/postprocess/dynamic_topography.h>
#include <aspect/postprocess/boundary_densities.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <aspect/citation_info.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::vector<double>,std::vector<double>>
                                                    Geoid<dim>::to_spherical_harmonic_coefficients(const std::vector<std::vector<double>> &spherical_function) const
    {
      std::vector<double> cosi(spherical_function.size(),0);
      std::vector<double> sini(spherical_function.size(),0);
      std::vector<double> coecos;
      std::vector<double> coesin;

      for (unsigned int ideg =  min_degree; ideg < max_degree+1; ++ideg)
        {
          for (unsigned int iord = 0; iord < ideg+1; ++iord)
            {
              // Do the spherical harmonic expansion.
              for (unsigned int ds_num = 0; ds_num < spherical_function.size(); ++ds_num)
                {
                  // Normalization after Dahlen and Tromp (1986) Appendix B.6.
                  const std::pair<double,double> sph_harm_vals = aspect::Utilities::real_spherical_harmonic(ideg,iord,spherical_function.at(ds_num).at(0),spherical_function.at(ds_num).at(1));
                  const double cos_component = sph_harm_vals.first; // real / cos part
                  const double sin_component = sph_harm_vals.second; // imaginary / sin part

                  cosi.at(ds_num) = (spherical_function.at(ds_num).at(3) * cos_component);
                  sini.at(ds_num) = (spherical_function.at(ds_num).at(3) * sin_component);
                }
              // Integrate the contribution of each spherical infinitesimal.
              double cosii = 0;
              double sinii = 0;
              for (unsigned int ds_num = 0; ds_num < spherical_function.size(); ++ds_num)
                {
                  cosii += cosi.at(ds_num) * spherical_function.at(ds_num).at(2);
                  sinii += sini.at(ds_num) * spherical_function.at(ds_num).at(2);
                }
              coecos.push_back(cosii);
              coesin.push_back(sinii);
            }
        }
      // Sum over each processor.
      dealii::Utilities::MPI::sum (coecos,this->get_mpi_communicator(),coecos);
      dealii::Utilities::MPI::sum (coesin,this->get_mpi_communicator(),coesin);

      return std::make_pair(coecos,coesin);
    }

    template <int dim>
    std::pair<std::vector<double>,std::vector<double>>
                                                    Geoid<dim>::density_contribution (const double &/*outer_radius*/) const
    {
      Assert(false, ExcNotImplemented());
      return std::make_pair(std::vector<double>(), std::vector<double>());

    }

    template <>
    std::pair<std::vector<double>,std::vector<double>>
                                                    Geoid<3>::density_contribution (const double &outer_radius) const
    {
      const unsigned int quadrature_degree = this->introspection().polynomial_degree.temperature;

      // Need to evaluate density contribution of each volume quadrature point.
      const QGauss<3> quadrature_formula(quadrature_degree);

      FEValues<3> fe_values (this->get_mapping(),
                             this->get_fe(),
                             quadrature_formula,
                             update_values |
                             update_quadrature_points |
                             update_JxW_values |
                             update_gradients);

      MaterialModel::MaterialModelInputs<3> in(fe_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<3> out(fe_values.n_quadrature_points, this->n_compositional_fields());

      std::vector<std::vector<double>> composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      // Directly do the global 3D integral over each quadrature point of every cell (different from traditional way to do layer integral).
      // This is necessary because of ASPECT's adaptive mesh refinement feature.
      std::vector<double> SH_density_coecos;
      std::vector<double> SH_density_coesin;
      for (unsigned int ideg =  min_degree; ideg < max_degree+1; ++ideg)
        {
          for (unsigned int iord = 0; iord < ideg+1; ++iord)
            {
              // Initialization of the density contribution integral per degree, order.
              double integrated_density_cos_component = 0;
              double integrated_density_sin_component = 0;

              // Loop over all of the cells.
              for (const auto &cell : this->get_dof_handler().active_cell_iterators())
                if (cell->is_locally_owned())
                  {
                    fe_values.reinit (cell);
                    // Set use_strain_rates to false since we don't need viscosity.
                    in.reinit(fe_values, cell, this->introspection(), this->get_solution(), false);

                    this->get_material_model().evaluate(in, out);

                    // Compute the integral of the density function
                    // over the cell, by looping over all quadrature points.
                    for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                      {
                        // Convert coordinates from [x,y,z] to [r, phi, theta].
                        const std::array<double,3> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(in.position[q]);

                        // Normalization after Dahlen and Tromp (1986) Appendix B.6.
                        const std::pair<double,double> sph_harm_vals = aspect::Utilities::real_spherical_harmonic(ideg,iord,scoord[2],scoord[1]);
                        const double cos_component = sph_harm_vals.first; // real / cos part
                        const double sin_component = sph_harm_vals.second; // imaginary / sin part

                        const double density = out.densities[q];
                        const double r_q = in.position[q].norm();

                        integrated_density_cos_component += density * (1./r_q) * std::pow(r_q/outer_radius,ideg+1) * cos_component * fe_values.JxW(q);
                        integrated_density_sin_component += density * (1./r_q) * std::pow(r_q/outer_radius,ideg+1) * sin_component * fe_values.JxW(q);
                      }
                  }
              SH_density_coecos.push_back(integrated_density_cos_component);
              SH_density_coesin.push_back(integrated_density_sin_component);
            }
        }
      // Sum over each processor.
      dealii::Utilities::MPI::sum (SH_density_coecos,this->get_mpi_communicator(),SH_density_coecos);
      dealii::Utilities::MPI::sum (SH_density_coesin,this->get_mpi_communicator(),SH_density_coesin);

      return std::make_pair(SH_density_coecos,SH_density_coesin);
    }

    template <int dim>
    std::pair<std::pair<double, std::pair<std::vector<double>,std::vector<double>>>, std::pair<double, std::pair<std::vector<double>,std::vector<double>>>>
    Geoid<dim>::topography_contribution(const double &/*outer_radius*/,
                                        const double &/*inner_radius*/) const
    {
      Assert(false, ExcNotImplemented());
      std::pair<double, std::pair<std::vector<double>,std::vector<double>>> temp;
      return std::make_pair(temp, temp);
    }

    template <>
    std::pair<std::pair<double, std::pair<std::vector<double>,std::vector<double>>>, std::pair<double, std::pair<std::vector<double>,std::vector<double>>>>
    Geoid<3>::topography_contribution(const double &outer_radius,
                                      const double &inner_radius) const
    {
      // Get a pointer to the boundary densities postprocessor.
      const Postprocess::BoundaryDensities<3> &boundary_densities =
        this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::BoundaryDensities<3>>();

      const double top_layer_average_density = boundary_densities.density_at_top();
      const double bottom_layer_average_density = boundary_densities.density_at_bottom();

      const types::boundary_id top_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      const types::boundary_id bottom_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("bottom");

      const unsigned int quadrature_degree = this->introspection().polynomial_degree.temperature;
      const QGauss<2> quadrature_formula_face(quadrature_degree);

      FEFaceValues<3> fe_face_values (this->get_mapping(),
                                      this->get_fe(),
                                      quadrature_formula_face,
                                      update_values |
                                      update_quadrature_points |
                                      update_JxW_values);

      // Vectors to store the location, infinitesimal area, and topography associated with each quadrature point of each surface and bottom cell respectively.
      std::vector<std::pair<Point<3>,std::pair<double,double>>> surface_stored_values;
      std::vector<std::pair<Point<3>,std::pair<double,double>>> CMB_stored_values;

      // Loop over all of the boundary cells and if one is at
      // surface or CMB, evaluate the topography vector there.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          {
            unsigned int face_idx = numbers::invalid_unsigned_int;
            bool at_upper_surface = false;
            {
              for (const unsigned int f : cell->face_indices())
                {
                  if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
                    {
                      // If the cell is at the top boundary, assign face_idx.
                      face_idx = f;
                      at_upper_surface = true;
                      break;
                    }
                  else if (cell->at_boundary(f) && cell->face(f)->boundary_id() == bottom_boundary_id)
                    {
                      // If the cell is at the bottom boundary, assign face_idx.
                      face_idx = f;
                      at_upper_surface = false;
                      break;
                    }
                }
              // If the cell is not at the boundary, jump to the next cell.
              if (face_idx == numbers::invalid_unsigned_int)
                continue;
            }

            // Focus on the boundary cell's upper face if on the top boundary and lower face if on the bottom boundary.
            fe_face_values.reinit(cell, face_idx);

            // Topography is evaluated at each quadrature
            // point on every top/bottom cell's boundary face. The
            // reason to do this -- as opposed to using a single
            // value per boundary face -- is that later in the
            // spherical harmonic expansion, we will calculate
            // sin(theta)*d_theta*d_phi by
            // infinitesimal_area/radius^2. The accuracy of this
            // transfer gets better as infinitesimal_area gets
            // closer to zero, so using every boundary quadrature
            // point's associated area (in the form of
            // FEFaceValues::JxW) will lead to better accuracy in
            // spherical harmonic expansion compared to using just
            // one average value per face, especially in the coarse
            // meshes.

            // If the cell is at the top boundary, add its contributions to the topography surface storage vector.
            if (at_upper_surface)
              {
                if (include_surface_topo_contribution == true)
                  {
                    if (use_free_surface_topography == true)
                      {
                        for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                          {
                            const Point<3> current_position = fe_face_values.quadrature_point(q);
                            const double topography = this->get_geometry_model().height_above_reference_surface(current_position);
                            surface_stored_values.emplace_back (current_position, std::make_pair(fe_face_values.JxW(q), topography));
                          }
                      }
                    else
                      {
                        // Get a reference to the dynamic topography postprocessor.
                        const Postprocess::DynamicTopography<3> &dynamic_topography =
                          this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::DynamicTopography<3>>();

                        // Get the already-computed dynamic topography solution.
                        const LinearAlgebra::BlockVector &topo_vector = dynamic_topography.topography_vector();

                        std::vector<double> topo_values(quadrature_formula_face.size());

                        fe_face_values[this->introspection().extractors.temperature].get_function_values(topo_vector, topo_values);

                        for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                          {
                            surface_stored_values.emplace_back (fe_face_values.quadrature_point(q), std::make_pair(fe_face_values.JxW(q), topo_values[q]));
                          }
                      }
                  }
                else
                  {
                    // If the cell is at the top boundary, but you do not include the surface topography in the geoid computation,
                    // add its zero contributions to the topography surface storage vector.
                    if (at_upper_surface)
                      {
                        for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                          {
                            surface_stored_values.emplace_back (fe_face_values.quadrature_point(q), std::make_pair(fe_face_values.JxW(q), 0.0));
                          }
                      }
                  }
              }

            // If the cell is at the bottom boundary, add its contributions to the topography CMB storage vector.
            if (at_upper_surface == false)
              {
                if (include_CMB_topo_contribution == true)
                  {
                    if (use_free_CMB_topography == true)
                      {
                        for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                          {
                            const Point<3> current_position = fe_face_values.quadrature_point(q);
                            const double topography = this->get_geometry_model().height_above_reference_surface(current_position) + (outer_radius - inner_radius);
                            CMB_stored_values.emplace_back (current_position, std::make_pair(fe_face_values.JxW(q), topography));
                          }
                      }
                    else
                      {
                        // Get a reference to the dynamic topography postprocessor.
                        const Postprocess::DynamicTopography<3> &dynamic_topography =
                          this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::DynamicTopography<3>>();

                        // Get the already-computed dynamic topography solution.
                        const LinearAlgebra::BlockVector &topo_vector = dynamic_topography.topography_vector();

                        std::vector<double> topo_values(quadrature_formula_face.size());

                        fe_face_values[this->introspection().extractors.temperature].get_function_values(topo_vector, topo_values);

                        for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                          {
                            CMB_stored_values.emplace_back (fe_face_values.quadrature_point(q), std::make_pair(fe_face_values.JxW(q), topo_values[q]));
                          }
                      }
                  }
                else
                  {
                    // If the cell is at the bottom boundary, but you do not include the CMB topography in the geoid computation,
                    // add its zero contributions to the topography CMB storage vector.
                    if (at_upper_surface == false)
                      {
                        for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                          {
                            CMB_stored_values.emplace_back (fe_face_values.quadrature_point(q), std::make_pair(fe_face_values.JxW(q), 0.0));
                          }
                      }
                  }
              }
          }

      std::vector<std::vector<double>> surface_topo_spherical_function;
      std::vector<std::vector<double>> CMB_topo_spherical_function;

      for (unsigned int i=0; i<surface_stored_values.size(); ++i)
        {
          const std::array<double,3> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(surface_stored_values[i].first);

          // Calculate spherical infinitesimal sin(theta)*d_theta*d_phi by infinitesimal_area/radius^2
          const double infinitesimal = surface_stored_values[i].second.first/(outer_radius*outer_radius);

          // Theta, phi, spherical infinitesimal, and surface topography
          surface_topo_spherical_function.emplace_back(std::vector<double> {scoord[2],
                                                                            scoord[1],
                                                                            infinitesimal,
                                                                            surface_stored_values[i].second.second
                                                                           });
        }

      for (unsigned int i=0; i<CMB_stored_values.size(); ++i)
        {
          const std::array<double,3> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(CMB_stored_values[i].first);

          // Calculate spherical infinitesimal sin(theta)*d_theta*d_phi by infinitesimal_area/radius^2
          const double infinitesimal = CMB_stored_values[i].second.first/(inner_radius*inner_radius);

          // Theta, phi, spherical infinitesimal, and CMB dynamic topography
          CMB_topo_spherical_function.emplace_back(std::vector<double> {scoord[2],
                                                                        scoord[1],
                                                                        infinitesimal,
                                                                        CMB_stored_values[i].second.second
                                                                       });
        }

      std::pair<double, std::pair<std::vector<double>,std::vector<double>>> SH_surface_topo_coes
        = std::make_pair(top_layer_average_density,to_spherical_harmonic_coefficients(surface_topo_spherical_function));
      std::pair<double, std::pair<std::vector<double>,std::vector<double>>> SH_CMB_topo_coes
        = std::make_pair(bottom_layer_average_density,to_spherical_harmonic_coefficients(CMB_topo_spherical_function));
      return std::make_pair(SH_surface_topo_coes,SH_CMB_topo_coes);
    }


    template <int dim>
    std::pair<std::string,std::string>
    Geoid<dim>::execute (TableHandler &)
    {
      // Current geoid code only works for spherical shell geometry.
      AssertThrow (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model())
                   &&
                   dim == 3,
                   ExcMessage("The geoid postprocessor is currently only implemented for the 3D spherical shell geometry model."));

      const GeometryModel::SphericalShell<dim> &geometry_model =
        Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model());

      // Get the value of the outer radius and inner radius.
      const double outer_radius = geometry_model.outer_radius();
      const double inner_radius = geometry_model.inner_radius();

      const types::boundary_id top_boundary_id = geometry_model.translate_symbolic_boundary_name_to_id("top");

      // Get the value of the surface gravity acceleration from the gravity model.
      Point<dim> surface_point;
      surface_point[0] = outer_radius;
      const double surface_gravity = this->get_gravity_model().gravity_vector(surface_point).norm();

      // Get the value of the universal gravitational constant.
      const double G = aspect::constants::big_g;

      // Get the spherical harmonic coefficients of the density contribution.
      std::pair<std::vector<double>,std::vector<double>> SH_density_coes = density_contribution(outer_radius);
      std::pair<double, std::pair<std::vector<double>,std::vector<double>>> SH_surface_topo_coes;
      std::pair<double, std::pair<std::vector<double>,std::vector<double>>> SH_CMB_topo_coes;

      // Initialize the surface and CMB density contrasts with NaNs because they may be unused in case of no topography contribution.
      double surface_delta_rho = numbers::signaling_nan<double>();
      double CMB_delta_rho = numbers::signaling_nan<double>();

      // Get the spherical harmonic coefficients of the surface and CMB topography.
      std::pair<std::pair<double, std::pair<std::vector<double>,std::vector<double>>>, std::pair<double, std::pair<std::vector<double>,std::vector<double>>>> SH_topo_coes;
      SH_topo_coes = topography_contribution(outer_radius,inner_radius);
      SH_surface_topo_coes = SH_topo_coes.first;
      SH_CMB_topo_coes = SH_topo_coes.second;

      // Get the density contrast at the surface and CMB to replace the initialized NaN values.
      // The surface and CMB density contrasts will be used later to calculate geoid,
      // and the spherical harmonic output of the surface and CMB topography contribution to geoid.
      surface_delta_rho =  SH_surface_topo_coes.first - density_above;
      CMB_delta_rho = density_below - SH_CMB_topo_coes.first;

      // Compute the spherical harmonic coefficients of geoid anomaly.
      std::vector<double> density_anomaly_contribution_coecos;
      std::vector<double> density_anomaly_contribution_coesin;
      std::vector<double> surface_topo_contribution_coecos;
      std::vector<double> surface_topo_contribution_coesin;
      std::vector<double> CMB_topo_contribution_coecos;
      std::vector<double> CMB_topo_contribution_coesin;
      geoid_coecos.clear();
      geoid_coesin.clear();

      // First compute the spherical harmonic contributions from density anomaly, surface topography and CMB topography.
      int ind = 0; // coefficients index
      for (unsigned int ideg =  min_degree; ideg < max_degree+1; ++ideg)
        {
          for (unsigned int iord = 0; iord < ideg+1; ++iord)
            {
              double coecos_density_anomaly = (4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1))) * SH_density_coes.first.at(ind);
              double coesin_density_anomaly = (4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1))) * SH_density_coes.second.at(ind);
              density_anomaly_contribution_coecos.push_back(coecos_density_anomaly);
              density_anomaly_contribution_coesin.push_back(coesin_density_anomaly);

              if (include_surface_topo_contribution == true || include_CMB_topo_contribution == true)
                {
                  const double coecos_surface_topo = (4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1)))
                                                     * surface_delta_rho*SH_surface_topo_coes.second.first.at(ind)*outer_radius;
                  const double coesin_surface_topo = (4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1)))
                                                     * surface_delta_rho*SH_surface_topo_coes.second.second.at(ind)*outer_radius;
                  surface_topo_contribution_coecos.push_back(coecos_surface_topo);
                  surface_topo_contribution_coesin.push_back(coesin_surface_topo);

                  const double coecos_CMB_topo = (4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1)))
                                                 * CMB_delta_rho*SH_CMB_topo_coes.second.first.at(ind)*inner_radius*std::pow(inner_radius/outer_radius,ideg+1);
                  const double coesin_CMB_topo = (4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1)))
                                                 * CMB_delta_rho*SH_CMB_topo_coes.second.second.at(ind)*inner_radius*std::pow(inner_radius/outer_radius,ideg+1);
                  CMB_topo_contribution_coecos.push_back(coecos_CMB_topo);
                  CMB_topo_contribution_coesin.push_back(coesin_CMB_topo);

                }

              ++ind;
            }
        }

      // Then sum the three contributions together to get the spherical harmonic coefficients of geoid anomaly.
      ind = 0; // coefficients index
      for (unsigned int ideg =  min_degree; ideg < max_degree+1; ++ideg)
        {
          for (unsigned int iord = 0; iord < ideg+1; ++iord)
            {
              if (include_surface_topo_contribution == true || include_CMB_topo_contribution == true)
                {
                  geoid_coecos.push_back(density_anomaly_contribution_coecos.at(ind)+surface_topo_contribution_coecos.at(ind)+CMB_topo_contribution_coecos.at(ind));
                  geoid_coesin.push_back(density_anomaly_contribution_coesin.at(ind)+surface_topo_contribution_coesin.at(ind)+CMB_topo_contribution_coesin.at(ind));
                }
              else
                {
                  geoid_coecos.push_back(density_anomaly_contribution_coecos.at(ind));
                  geoid_coesin.push_back(density_anomaly_contribution_coesin.at(ind));
                }

              ind += 1;
            }
        }

      const QMidpoint<dim-1> quadrature_formula_face_center;
      Assert(quadrature_formula_face_center.size() == 1, ExcInternalError());
      FEFaceValues<dim> fe_face_center_values (this->get_mapping(),
                                               this->get_fe(),
                                               quadrature_formula_face_center,
                                               update_values |
                                               update_quadrature_points|
                                               update_JxW_values);

      // Define a vector to store the location of the cells along the surface.
      std::vector<Point<dim>> surface_cell_locations;

      // Loop over all the cells to get the locations of the surface cells to prepare for the geoid computation.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          {
            // If the cell is at the top boundary, store the cell's upper face midpoint location.
            for (const unsigned int f : cell->face_indices())
              if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
                {
                  fe_face_center_values.reinit(cell,f);
                  const Point<dim> midpoint_at_top_face = fe_face_center_values.get_quadrature_points().at(0);
                  surface_cell_locations.push_back(midpoint_at_top_face);
                  break;
                }
          }

      // Transfer the geocentric coordinates of the surface cells to the surface spherical coordinates (theta,phi)
      std::vector<std::pair<double,double>> surface_cell_spherical_coordinates;
      for (unsigned int i=0; i<surface_cell_locations.size(); ++i)
        {
          const std::array<double,dim> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(surface_cell_locations.at(i));
          const double phi = scoord[1];
          const double theta = scoord[2];
          surface_cell_spherical_coordinates.emplace_back(theta,phi);
        }

      // Compute the grid geoid anomaly based on spherical harmonics.
      std::vector<double> geoid_anomaly;
      for (unsigned int i=0; i<surface_cell_spherical_coordinates.size(); ++i)
        {
          int ind = 0;
          double geoid_value = 0;
          for (unsigned int ideg =  min_degree; ideg < max_degree+1; ++ideg)
            {
              for (unsigned int iord = 0; iord < ideg+1; ++iord)
                {
                  // Normalization after Dahlen and Tromp (1986) Appendix B.6.
                  const std::pair<double,double> sph_harm_vals = aspect::Utilities::real_spherical_harmonic(ideg,iord,surface_cell_spherical_coordinates.at(i).first,surface_cell_spherical_coordinates.at(i).second);
                  const double cos_component = sph_harm_vals.first; // real / cos part
                  const double sin_component = sph_harm_vals.second; // imaginary / sin part

                  geoid_value += geoid_coecos.at(ind)*cos_component+geoid_coesin.at(ind)*sin_component;
                  ++ind;
                }
            }
          geoid_anomaly.push_back(geoid_value);
        }

      // The user can get the spherical harmonic coefficients of the density anomaly contribution if needed
      if (output_density_anomaly_contribution_SH_coes == true)
        {
          // Have a stream into which we write the SH coefficients data from density anomaly contribution.
          // The text stream is then later sent to processor 0.
          std::ostringstream output_density_anomaly_contribution_SH_coes;

          // Prepare the output SH coefficients data from density anomaly contribution.
          unsigned int SH_coes_ind = 0;
          for (unsigned int ideg =  min_degree; ideg < max_degree+1; ++ideg)
            {
              for (unsigned int iord = 0; iord < ideg+1; ++iord)
                {
                  output_density_anomaly_contribution_SH_coes << ideg
                                                              << ' '
                                                              << iord
                                                              << ' '
                                                              << density_anomaly_contribution_coecos.at(SH_coes_ind)
                                                              << ' '
                                                              << density_anomaly_contribution_coesin.at(SH_coes_ind)
                                                              << std::endl;
                  ++SH_coes_ind;
                }
            }

          const std::string density_anomaly_contribution_SH_coes_filename = this->get_output_directory() +
                                                                            "density_anomaly_contribution_SH_coefficients." +
                                                                            dealii::Utilities::int_to_string(this->get_timestep_number(), 5);

          // Because each processor already held all the SH coefficients from density anomaly contribution, we only need to stop by the processor 0 to get the data.
          // On processor 0, collect all the data and put them into the output density anomaly contribution SH coefficients file.
          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream density_anomaly_contribution_SH_coes_file (density_anomaly_contribution_SH_coes_filename.c_str());
              density_anomaly_contribution_SH_coes_file << "# "
                                                        << "degree order cosine_coefficient sine_coefficient"
                                                        << std::endl;

              // Write out the data on processor 0.
              density_anomaly_contribution_SH_coes_file << output_density_anomaly_contribution_SH_coes.str();
            }
        }

      // The user can get the spherical harmonic coefficients of the surface topography contribution if needed
      if (output_surface_topo_contribution_SH_coes == true)
        {
          // Have a stream into which we write the SH coefficients data from surface topography contribution.
          // The text stream is then later sent to processor 0.
          std::ostringstream output_surface_topo_contribution_SH_coes;

          // Prepare the output SH coefficients data from surface topography contribution.
          unsigned int SH_coes_ind = 0;
          for (unsigned int ideg =  min_degree; ideg < max_degree+1; ++ideg)
            {
              for (unsigned int iord = 0; iord < ideg+1; ++iord)
                {
                  output_surface_topo_contribution_SH_coes << ideg
                                                           << ' '
                                                           << iord
                                                           << ' '
                                                           << surface_topo_contribution_coecos.at(SH_coes_ind)
                                                           << ' '
                                                           << surface_topo_contribution_coesin.at(SH_coes_ind)
                                                           << std::endl;
                  ++SH_coes_ind;
                }
            }

          const std::string surface_topo_contribution_SH_coes_filename = this->get_output_directory() +
                                                                         "surface_topography_contribution_SH_coefficients." +
                                                                         dealii::Utilities::int_to_string(this->get_timestep_number(), 5);

          // Because each processor already held all the SH coefficients from surface topography contribution,
          // we only need to stop by the processor 0 to get the data. On processor 0, collect all the data
          // and put them into the output surface topography contribution SH coefficients file.
          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream surface_topo_contribution_SH_coes_file (surface_topo_contribution_SH_coes_filename.c_str());
              surface_topo_contribution_SH_coes_file << "# "
                                                     << "degree order cosine_coefficient sine_coefficient"
                                                     << std::endl;
              std::ostringstream output_surface_delta_rho;
              output_surface_delta_rho << surface_delta_rho;
              surface_topo_contribution_SH_coes_file << "surface density contrast(kg/m^3): "
                                                     << output_surface_delta_rho.str()
                                                     << std::endl;
              // Write out the data on processor 0
              surface_topo_contribution_SH_coes_file << output_surface_topo_contribution_SH_coes.str();
            }
        }

      // The user can get the spherical harmonic coefficients of the CMB topography contribution if needed.
      if (output_CMB_topo_contribution_SH_coes == true)
        {
          // Have a stream into which we write the SH coefficients data from CMB topography contribution.
          // The text stream is then later sent to processor 0.
          std::ostringstream output_CMB_topo_contribution_SH_coes;

          // Prepare the output SH coefficients data from CMB topography contribution.
          unsigned int SH_coes_ind = 0;
          for (unsigned int ideg =  min_degree; ideg < max_degree+1; ++ideg)
            {
              for (unsigned int iord = 0; iord < ideg+1; ++iord)
                {
                  output_CMB_topo_contribution_SH_coes << ideg
                                                       << ' '
                                                       << iord
                                                       << ' '
                                                       << CMB_topo_contribution_coecos.at(SH_coes_ind)
                                                       << ' '
                                                       << CMB_topo_contribution_coesin.at(SH_coes_ind)
                                                       << std::endl;
                  ++SH_coes_ind;
                }
            }

          const std::string CMB_topo_contribution_SH_coes_filename = this->get_output_directory() +
                                                                     "CMB_topography_contribution_SH_coefficients." +
                                                                     dealii::Utilities::int_to_string(this->get_timestep_number(), 5);

          // Because each processor already held all the SH coefficients from CMB topography contribution, we only need to stop by the processor 0
          // to get the data. On processor 0, collect all the data and put them into the output CMB topography contribution SH coefficients file.
          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream CMB_topo_contribution_SH_coes_file (CMB_topo_contribution_SH_coes_filename.c_str());
              CMB_topo_contribution_SH_coes_file << "# "
                                                 << "degree order cosine_coefficient sine_coefficient"
                                                 << std::endl;
              std::ostringstream output_CMB_delta_rho;
              output_CMB_delta_rho << CMB_delta_rho;
              CMB_topo_contribution_SH_coes_file << "CMB density contrast(kg/m^3): "
                                                 << output_CMB_delta_rho.str()
                                                 << std::endl;
              // Write out the data on processor 0
              CMB_topo_contribution_SH_coes_file << output_CMB_topo_contribution_SH_coes.str();
            }
        }

      // The user can get the spherical harmonic coefficients of the geoid anomaly if needed.
      if (output_geoid_anomaly_SH_coes == true)
        {
          // Have a stream into which we write the geoid anomaly SH coefficients data.
          // The text stream is then later sent to processor 0.
          std::ostringstream output_geoid_anomaly_SH_coes;

          // Prepare the output geoid anomaly SH coefficients data.
          unsigned int SH_coes_ind = 0;
          for (unsigned int ideg =  min_degree; ideg < max_degree+1; ++ideg)
            {
              for (unsigned int iord = 0; iord < ideg+1; ++iord)
                {
                  output_geoid_anomaly_SH_coes << ideg
                                               << ' '
                                               << iord
                                               << ' '
                                               << geoid_coecos.at(SH_coes_ind)
                                               << ' '
                                               << geoid_coesin.at(SH_coes_ind)
                                               << std::endl;
                  ++SH_coes_ind;
                }
            }

          const std::string geoid_anomaly_SH_coes_filename = this->get_output_directory() +
                                                             "geoid_anomaly_SH_coefficients." +
                                                             dealii::Utilities::int_to_string(this->get_timestep_number(), 5);

          // Because each processor already held all the geoid anomaly SH coefficients, we only need to stop by the processor 0 to get the data.
          // On processor 0, collect all the data and put them into the output geoid anomaly SH coefficients file.
          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream geoid_anomaly_SH_coes_file (geoid_anomaly_SH_coes_filename.c_str());
              geoid_anomaly_SH_coes_file << "# "
                                         << "degree order cosine_coefficient sine_coefficient"
                                         << std::endl;

              // Write out the data on processor 0.
              geoid_anomaly_SH_coes_file << output_geoid_anomaly_SH_coes.str();
            }
        }

      // Have a stream into which we write the geoid height data. the text stream is then
      // later sent to processor 0.
      std::ostringstream output;

      // Prepare the output data.
      if (output_in_lat_lon == true)
        {
          double lon, lat;
          for (unsigned int i=0; i<surface_cell_spherical_coordinates.size(); ++i)
            {
              // Transfer the spherical coordinates to geographical coordinates.
              lat = 90. - surface_cell_spherical_coordinates.at(i).first*(180./numbers::PI);
              lon = (surface_cell_spherical_coordinates.at(i).second <= numbers::PI
                     ?
                     surface_cell_spherical_coordinates.at(i).second*(180./numbers::PI)
                     :
                     surface_cell_spherical_coordinates.at(i).second*(180./numbers::PI) - 360.);

              // Write the solution to the stream output.
              output << lon
                     << ' '
                     << lat
                     << ' '
                     << geoid_anomaly.at(i)
                     << std::endl;
            }
        }
      else
        {
          for (unsigned int i=0; i<surface_cell_locations.size(); ++i)
            {
              // Write the solution to the stream output.
              output << surface_cell_locations.at(i)
                     << ' '
                     << geoid_anomaly.at(i)
                     << std::endl;
            }
        }

      const std::string filename = this->get_output_directory() +
                                   "geoid_anomaly." +
                                   dealii::Utilities::int_to_string(this->get_timestep_number(), 5);
      const std::vector<std::string> data = Utilities::MPI::gather(this->get_mpi_communicator(), output.str());

      // On processor 0, collect all of the data the individual processors sent
      // and concatenate them into one file:
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::ofstream file (filename.c_str());

          file << "# "
               << ((output_in_lat_lon == true)? "longitude latitude" : "x y z")
               << " geoid_anomaly" << std::endl;

          for (const auto &str : data)
            file << str;
        }

      // Prepare the free-air gravity anomaly output.
      if (output_gravity_anomaly == true)
        {
          // Have a stream into which we write the gravity anomaly data. the text stream is then
          // later sent to processor 0.
          std::ostringstream output_gravity_anomaly;

          // Compute the grid gravity anomaly based on spherical harmonics.
          std::vector<double> gravity_anomaly;
          gravity_anomaly.reserve(surface_cell_spherical_coordinates.size());

          for (unsigned int i=0; i<surface_cell_spherical_coordinates.size(); ++i)
            {
              int ind = 0;
              double gravity_value = 0;
              for (unsigned int ideg =  min_degree; ideg < max_degree+1; ++ideg)
                {
                  for (unsigned int iord = 0; iord < ideg+1; ++iord)
                    {
                      // Normalization after Dahlen and Tromp (1986) Appendix B.6.
                      const std::pair<double,double> sph_harm_vals = aspect::Utilities::real_spherical_harmonic(ideg,iord,surface_cell_spherical_coordinates.at(i).first,surface_cell_spherical_coordinates.at(i).second);
                      const double cos_component = sph_harm_vals.first; // real / cos part
                      const double sin_component = sph_harm_vals.second; // imaginary / sin part

                      // The conversion from geoid to gravity anomaly is given by gravity_anomaly = (l-1)*g/R_surface * geoid_anomaly
                      // based on Forte (2007) equation [97].
                      gravity_value += (geoid_coecos.at(ind)*cos_component+geoid_coesin.at(ind)*sin_component) * (ideg - 1) * surface_gravity / outer_radius;
                      ++ind;
                    }
                }
              gravity_anomaly.push_back(gravity_value);
            }

          // Prepare the output data.
          if (output_in_lat_lon == true)
            {
              double lon, lat;
              for (unsigned int i=0; i<surface_cell_spherical_coordinates.size(); ++i)
                {
                  // Transfer the spherical coordinates to geographical coordinates.
                  lat = 90. - surface_cell_spherical_coordinates.at(i).first*(180./numbers::PI);
                  lon = (surface_cell_spherical_coordinates.at(i).second <= numbers::PI
                         ?
                         surface_cell_spherical_coordinates.at(i).second*(180./numbers::PI)
                         :
                         surface_cell_spherical_coordinates.at(i).second*(180./numbers::PI) - 360.);

                  // Write the solution to the stream output.
                  output_gravity_anomaly << lon
                                         << ' '
                                         << lat
                                         << ' '
                                         << gravity_anomaly.at(i)
                                         << std::endl;
                }
            }
          else
            {
              for (unsigned int i=0; i<surface_cell_locations.size(); ++i)
                {
                  // Write the solution to the stream output.
                  output_gravity_anomaly << surface_cell_locations.at(i)
                                         << ' '
                                         << gravity_anomaly.at(i)
                                         << std::endl;
                }
            }

          const std::string filename = this->get_output_directory() +
                                       "gravity_anomaly." +
                                       dealii::Utilities::int_to_string(this->get_timestep_number(), 5);
          const std::vector<std::string> data = Utilities::MPI::gather(this->get_mpi_communicator(), output_gravity_anomaly.str());

          // On processor 0, collect all of the data the individual processors sent
          // and concatenate them into one file:
          if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream file (filename.c_str());

              file << "# "
                   << ((output_in_lat_lon == true)? "longitude latitude" : "x y z")
                   << " gravity_anomaly" << std::endl;

              for (const auto &str : data)
                file << str;
            }
        }

      return std::pair<std::string,std::string>("Writing geoid anomaly:",
                                                filename);
    }

    template <int dim>
    std::list<std::string>
    Geoid<dim>::required_other_postprocessors() const
    {
      std::list<std::string> deps;

      if ( (include_surface_topo_contribution == true && use_free_surface_topography == false) || (include_CMB_topo_contribution == true && use_free_CMB_topography == false) )
        deps.emplace_back("dynamic topography");

      deps.emplace_back("boundary densities");

      return deps;
    }

    template <int dim>
    double
    Geoid<dim>::evaluate (const Point<dim> &/*p*/) const
    {
      Assert(false, ExcNotImplemented());
      return 0;
    }

    template <>
    double
    Geoid<3>::evaluate (const Point<3> &p) const
    {
      const std::array<double,3> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(p);
      const double theta = scoord[2];
      const double phi = scoord[1];
      double value = 0.;

      for (unsigned int ideg=min_degree, k=0; ideg < max_degree+1; ++ideg)
        for (unsigned int iord = 0; iord < ideg+1; ++iord, ++k)
          {
            std::pair<double,double> val = aspect::Utilities::real_spherical_harmonic( ideg, iord, theta, phi );

            value += geoid_coecos[k] * val.first +
                     geoid_coesin[k] * val.second;

          }
      return value;
    }

    template <int dim>
    void
    Geoid<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Geoid");
        {
          prm.declare_entry("Include the contributon from dynamic topography", "true",
                            Patterns::Bool(),
                            "Option to include the contribution from dynamic topography on geoid. The default is true.");
          prm.declare_entry("Include surface topography contribution", "true",
                            Patterns::Bool(),
                            "Option to include the contribution from surface topography on geoid. The default is true.");
          prm.declare_entry("Include CMB topography contribution", "true",
                            Patterns::Bool(),
                            "Option to include the contribution from CMB topography on geoid. The default is true.");
          prm.declare_entry("Maximum degree","20",
                            Patterns::Integer (0),
                            "This parameter can be a random positive integer. However, the value normally should not exceed the maximum "
                            "degree of the initial perturbed temperature field. For example, if the initial temperature uses S40RTS, the "
                            "maximum degree should not be larger than 40.");
          prm.declare_entry("Minimum degree","2",
                            Patterns::Integer (0),
                            "This parameter normally is set to 2 since the perturbed gravitational potential at degree 1 always vanishes "
                            "in a reference frame with the planetary center of mass same as the center of figure.");
          prm.declare_entry("Output data in geographical coordinates", "false",
                            Patterns::Bool(),
                            "Option to output the geoid anomaly in geographical coordinates (latitude and longitude). "
                            "The default is false, so postprocess will output the data in geocentric coordinates (x,y,z) as normally.");
          prm.declare_entry("Density above","0.",
                            Patterns::Double (0.),
                            "The density value above the surface boundary.");
          prm.declare_entry("Density below","9900.",
                            Patterns::Double (0.),
                            "The density value below the CMB boundary.");
          prm.declare_entry("Output geoid anomaly coefficients", "false",
                            Patterns::Bool(),
                            "Option to output the spherical harmonic coefficients of the geoid anomaly up to the maximum degree. "
                            "The default is false, so postprocess will only output the geoid anomaly in grid format. ");
          prm.declare_entry("Output surface topography contribution coefficients", "false",
                            Patterns::Bool(),
                            "Option to output the spherical harmonic coefficients of the surface topography contribution "
                            "to the maximum degree. The default is false. ");
          prm.declare_entry("Output CMB topography contribution coefficients", "false",
                            Patterns::Bool(),
                            "Option to output the spherical harmonic coefficients of the CMB topography contribution "
                            "to the maximum degree. The default is false. ");
          prm.declare_entry("Output density anomaly contribution coefficients", "false",
                            Patterns::Bool(),
                            "Option to output the spherical harmonic coefficients of the density anomaly contribution to the "
                            "maximum degree. The default is false. ");
          prm.declare_entry("Output gravity anomaly", "false",
                            Patterns::Bool(),
                            "Option to output the free-air gravity anomaly up to the maximum degree. "
                            "The unit of the output is in SI, hence $m/s^2$ ($1mgal = 10^-5 m/s^2$). The default is false. ");

          prm.declare_alias("Output geoid anomaly coefficients","Also output the spherical harmonic coefficients of geoid anomaly");
          prm.declare_alias("Output surface topography contribution coefficients","Also output the spherical harmonic coefficients of surface dynamic topography contribution");
          prm.declare_alias("Output CMB topography contribution coefficients","Also output the spherical harmonic coefficients of CMB dynamic topography contribution");
          prm.declare_alias("Output density anomaly contribution coefficients","Also output the spherical harmonic coefficients of density anomaly contribution");
          prm.declare_alias("Output gravity anomaly","Also output the gravity anomaly");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void
    Geoid<dim>::parse_parameters (ParameterHandler &prm)
    {
      CitationInfo::add("geoid");

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Geoid");
        {
          const bool include_topo_contribution = prm.get_bool ("Include the contributon from dynamic topography");

          AssertThrow (include_topo_contribution == true,
                       ExcMessage("The parameter 'Include the contributon from dynamic topography' has been "
                                  " replaced by the two parameters 'Include surface topography contribution' and "
                                  "'Include CMB topography contribution'. Please use them instead."));

          include_surface_topo_contribution = prm.get_bool ("Include surface topography contribution");
          include_CMB_topo_contribution = prm.get_bool ("Include CMB topography contribution");
          max_degree = prm.get_integer ("Maximum degree");
          min_degree = prm.get_integer ("Minimum degree");
          output_in_lat_lon = prm.get_bool ("Output data in geographical coordinates");
          density_above = prm.get_double ("Density above");
          density_below = prm.get_double ("Density below");
          output_geoid_anomaly_SH_coes = prm.get_bool ("Output geoid anomaly coefficients");
          output_surface_topo_contribution_SH_coes = prm.get_bool ("Output surface topography contribution coefficients");
          output_CMB_topo_contribution_SH_coes = prm.get_bool ("Output CMB topography contribution coefficients");
          output_density_anomaly_contribution_SH_coes = prm.get_bool ("Output density anomaly contribution coefficients");
          output_gravity_anomaly = prm.get_bool ("Output gravity anomaly");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void
    Geoid<dim>::initialize ()
    {
      // Find if the included boundaries are active free surfaces
      if (include_surface_topo_contribution == true || include_CMB_topo_contribution == true)
        {
          if (this->get_parameters().mesh_deformation_enabled == true)
            {
              const  types::boundary_id surface_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
              const  types::boundary_id bottom_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("bottom");

              const std::set<types::boundary_id> mesh_deformation_boundaries = this->get_mesh_deformation_handler().get_active_mesh_deformation_boundary_indicators();

              use_free_surface_topography = mesh_deformation_boundaries.find(surface_id) != mesh_deformation_boundaries.end();
              use_free_CMB_topography = mesh_deformation_boundaries.find(bottom_id) != mesh_deformation_boundaries.end();
            }
        }
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(Geoid,
                                  "geoid",
                                  "A postprocessor that computes a representation of "
                                  "the geoid based on the density structure in the mantle, "
                                  "as well as the topography at the surface and "
                                  "core mantle boundary (CMB) if desired. The topography is based on the "
                                  "dynamic topography postprocessor in case of no free surface, "
                                  "and based on the real surface from the geometry model in case "
                                  "of a free surface. The geoid is computed "
                                  "from a spherical harmonic expansion, so the geometry "
                                  "of the domain must be a 3D spherical shell.")
  }
}
