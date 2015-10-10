/*
 * This file is part of nanolens, a free program to calculate microlensing 
 * magnification patterns.
 * Copyright (C) 2015  Aksel Alpay
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef FRAGMENTED_LENS_PLANE_HPP
#define	FRAGMENTED_LENS_PLANE_HPP

#include "lens_plane.hpp"
#include "spatial_grid_db.hpp"


namespace nanolens{

/// Implements a microlensing lens plane that is segmented along the x-Axis into
/// so called "fragments".
/// When the lens equation is evaluated, it loads the appropriate fragment and
/// all stars that are within a certain radius around the fragment center
/// into memory (if it is not already present). To correct for the jumps of 
/// the center of mass when a new fragment is loaded, a smooth mass contribution
/// is added with a hole where the stars are. This implementation only works 
/// for uniform star distributions. Note that due to fluctuations in the convergence
/// of stars a fragment, some visible inaccuracies may still exist. However, they
/// are usually so small that they can be convolved away.
template<class Deflection_engine_type>
class horizontally_fragmented_microlensing_lens_plane : public plane
{
public:
  typedef typename Deflection_engine_type::settings deflection_engine_settings_type;
  
  typedef microlensing_lens_plane<Deflection_engine_type> lens_plane_fragment;
  typedef std::shared_ptr<lens_plane_fragment> lens_plane_fragment_ptr;
  typedef spatial_grid_db<util::scalar, 2> star_db_type;
  
  /// Construct object
  /// \param star_db A spatial grid database object that contains the stars. Since
  /// this lens plane is built to calculate very long lightcurves, the amount of stars
  /// to include can exceed the available system memory. Hence, we need to load
  /// stars from the disk when we need them. For this, we use a \c spatial_grid_db
  /// as it can efficiently find all stars in the vicinity of a given point.
  /// \param settings The settings object for the deflection engine.
  /// \param shear The external shear
  /// \param sigma_star The convergence due to the stars. Since the stars are stored
  /// on the hard disk, we cannot efficiently calculate it here. Instead, we rely
  /// on the user having already calculated \c sigma_star.
  /// \param sigma_smooth The convergence due to smooth matter
  /// \param shear_rotation_angle The angle in degrees of the direction of the shear.
  /// 0 degrees means that the shear will be along the y-axis.
  /// \param fragment_size The size of a fragment. When a fragment is left, new
  /// new stars will be loaded from the disk, and a new deflection engine will be
  /// constructed. Hence, a too small value will decrease performance while a
  /// too large value will lead to decreased accuracy.
  /// \param fragment_star_distribution_radius The radius within which stars will
  /// be included when a new fragment is loaded into memory.
  explicit horizontally_fragmented_microlensing_lens_plane(const star_db_type& star_db,
                                              const deflection_engine_settings_type& settings,
                                              util::scalar y_fragment_center,
                                              util::scalar shear,
                                              util::scalar sigma_star,
                                              util::scalar sigma_smooth,
                                              util::scalar shear_rotation_angle,
                                              util::scalar fragment_size,
                                              util::scalar fragment_star_distribution_radius)
  : _star_db(star_db),
    _fragment_size(fragment_size),
    _fragment_star_distribution_radius(fragment_star_distribution_radius),
    _shear(shear),
    _sigma_smooth(sigma_smooth),
    _shear_rotation_angle(shear_rotation_angle),
    _y_fragment_center(y_fragment_center),
    _half_fragment_size(0.5 * fragment_size),
    _sigma_star(sigma_star),
    _settings(settings)
  {
    _current_fragment_center = {0.0, _y_fragment_center};
    
    _smooth_matter_fraction = _sigma_smooth / (_sigma_star + _sigma_smooth);
    
    load_new_fragment(_current_fragment_center);
  }
  
  /// Transforms from a position in the lens plane to a position in the source plane
  /// using the normalized lens equation.
  /// \param lens_plane_pos The position in the lens plane in Einstein radii.
  /// \return The resulting position in the source plane in Einstein radii.
  inline util::vector2 lensing_transformation(const util::vector2& lens_plane_pos) const
  {
    if(!is_within_fragment_range(lens_plane_pos))
      load_new_fragment(get_next_fragment_center(lens_plane_pos));
    
    util::vector2 result = _current_fragment->lensing_transformation(lens_plane_pos);
    // Correct result
    util::sub(result, get_deflection_correction(lens_plane_pos));
    
    
    return result;
  }
  
  /// Given a position in the source plane and a position in the lens plane,
  /// calculates the corresponding position in the image plane.
  /// \param lens_plane_pos The position in the plane in Einstein radii
  /// \param source_plane_pos The position in the source plane in Einstein radii
  /// \return The position in the image plane
  inline util::vector2 inverse_lensing_transformation(const util::vector2& lens_plane_pos,
                                                      const util::vector2& source_plane_pos) const
  {
    if(!is_within_fragment_range(lens_plane_pos))
      load_new_fragment(get_next_fragment_center(lens_plane_pos));
    
    util::vector2 result = _current_fragment->inverse_lensing_transformation(lens_plane_pos,
                                                                     source_plane_pos);
    // Correct result
    util::add(result, get_deflection_correction(lens_plane_pos));
    
    return result;
  }
  
  /// \return The external shear
  util::scalar get_shear() const
  {
    return _shear;
  }

  /// \return The fraction between smooth matter density and total density
  util::scalar get_smooth_matter_fraction() const
  {
    return _smooth_matter_fraction;
  }

  /// \return The convergence due to smooth matter.
  util::scalar get_sigma_smooth() const
  {
    return _sigma_smooth;
  }

  /// \return The convergence due to the stellar mass.
  util::scalar get_sigma_star() const
  {
    return _sigma_star;
  }

  /// \return The mean total convergence
  util::scalar get_mean_surface_density() const
  {
    return _sigma_smooth + _sigma_star;
  }
  
  /// Estimates the coordinates of a rectangle in the lens plane that is mapped
  /// to a given rectangle in the source plane
  /// \param source_plane_coordinates A matrix containing the coordinates of the
  /// source plane rectangle
  /// \param out A matrix into which the resulting coordinates in the lens plane
  /// will be written.
  void estimate_mapped_region_coordinates(const util::matrix_nxn<util::vector2, 2>& source_plane_coordinates,
                                          util::matrix_nxn<util::vector2, 2>& out) const
  {
    // Get shooting region using the global sigma*.
    this->_current_fragment->estimate_mapped_region_coordinates(source_plane_coordinates, out, _sigma_star);
    
  }
  
  /// Gathers information about the lens plane.
  /// \param out A map where the information will be stored. Each key in the map
  /// will specify the property name and the corresponding value will be the property
  /// value.
  void obtain_properties_set(std::map<std::string, util::scalar>& out) const
  {
    out.clear();

    out["N* (current fragment)"] = this->_current_fragment->get_num_stars();
    out["sigma* (estimated, current fragment)"] = this->get_sigma_star();
    out["sigma smooth"] = this->get_sigma_smooth();
    out["sigma total (current fragment)"] = this->get_mean_surface_density();
    out["shear"] = this->get_shear();
  }
  
private:
  
  /// Corrects the deflection for errors arising from the fragmented lens plane.
  /// When loading a new fragment, the center of the mass distribution of the
  /// stars will (in general) be shifted. If uncorrected, this leads to jumps
  /// in the deflection angle at the border between two fragments. We improve
  /// this error by adding a homogenous disk with constant convergence. Into this
  /// disk we "cut a hole" by adding second disk with negative convergence. All stars
  /// will then be put into this hole. This ensures that the center of mass always
  /// remains (more or less) constant at the coordinate origin.
  /// \param point The point for which the deflection shall be corrected
  /// \param A two dimensional vector that needs to be subtracted from a deflection
  /// angle to apply the correction.
  inline util::vector2 get_deflection_correction(const util::vector2& point) const
  { 
    util::vector2 r_hole = point;
    util::sub(r_hole, _center_of_mass);
    
    util::vector2 hole_contribution = r_hole;
 
    util::scalar R_squared = util::square(_fragment_star_distribution_radius);
    if(util::dot(r_hole, r_hole) <= R_squared)
    {
      util::scale(hole_contribution, -_sigma_star);
    }
    else
    {
      star hole_pseudo_star(_current_fragment_center, -R_squared * _sigma_star);
      hole_pseudo_star.calculate_deflection_angle(point, hole_contribution);
    }
    
    util::vector2 result = point;
    util::scale(result, _sigma_star);
    util::add(result, hole_contribution);

    return result;
  }
  
  /// Checks whether a point is within the current fragment. If it is not, a
  /// different fragment should be loaded.
  /// \param point The point that shall be checked
  /// \return whether the point is within the current fragment.
  inline bool is_within_fragment_range(const util::vector2& point) const
  {
    return std::abs(point[0] - _current_fragment_center[0]) <= _half_fragment_size;
  }
  
  /// Loads a new fragment.
  /// In more detail, loads all stars from the \c star_db that lie within a
  /// circle with its center at the given fragment center and the radius that
  /// has been specified in the constructor. These stars will then be used
  /// to create a new \c microlensing_lens_plane that serves as backend
  /// for calls to \c lensing_transformation and \c inverse_lensing_transformation.
  /// \param fragment_center The center of the fragment that shall be loaded
  void load_new_fragment(const util::vector2& fragment_center) const
  {
    _current_fragment_center = fragment_center;
    
    std::vector<star> stars;
    
    star_db_type::query_result_list_type query_result;
    this->_star_db.query_within_radius(_current_fragment_center, 
                                       _fragment_star_distribution_radius,
                                       query_result);

    stars.reserve(query_result.size());
    _total_mass = 0.0;
    _center_of_mass = {0.0, 0.0};
    for(const auto& element : query_result)
    {
      stars.push_back(star(element.first, element.second));
      _total_mass += element.second;
      
      util::vector2 center_of_mass_contribution = element.first;
      util::scale(center_of_mass_contribution, element.second);
      util::add(_center_of_mass, center_of_mass_contribution);
    }
    util::scale(_center_of_mass, 1.0 / _total_mass);
    
    _current_fragment = nullptr;
    
    _current_fragment = lens_plane_fragment_ptr(new lens_plane_fragment(stars, 
                                                                        _settings, 
                                                                        _shear, 
                                                                        _sigma_smooth,
                                                                        _shear_rotation_angle));
    
    
  }
  
  /// Calculates the center of the fragment closest to a given point
  /// \param position The point of which closest fragment center shall be found.
  /// \param The center of the closest fragment
  util::vector2 get_next_fragment_center(const util::vector2& position) const
  {
    util::scalar delta_x = position[0] - _current_fragment_center[0];
    
    long long int num_spanned_fragments = std::round(delta_x / _fragment_size);
    
    util::scalar x = _current_fragment_center[0] + num_spanned_fragments * _fragment_size;
    
    util::vector2 result = {x, _y_fragment_center};
    return result;
  }
  
  star_db_type _star_db;
  
  mutable util::vector2 _current_fragment_center;

  mutable lens_plane_fragment_ptr _current_fragment;
  
  util::scalar _shear;
  util::scalar _sigma_smooth;
  util::scalar _shear_rotation_angle;
  
  util::scalar _fragment_size;
  util::scalar _half_fragment_size;
  util::scalar _fragment_star_distribution_radius;
  mutable util::scalar _total_mass;
  
  mutable util::scalar _sigma_star;
  mutable util::scalar _smooth_matter_fraction;
  
  util::scalar _y_fragment_center;
  
  mutable util::vector2 _center_of_mass;
  
  deflection_engine_settings_type _settings;
};

}

#endif	/* FRAGMENTED_LENS_PLANE_HPP */

