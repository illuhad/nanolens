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


template<class Deflection_engine_type>
class horizontally_fragmented_microlensing_lens_plane : public plane
{
public:
  typedef typename Deflection_engine_type::settings deflection_engine_settings_type;
  
  typedef microlensing_lens_plane<Deflection_engine_type> lens_plane_fragment;
  typedef std::shared_ptr<lens_plane_fragment> lens_plane_fragment_ptr;
  typedef spatial_grid_db<util::scalar, 2> star_db_type;
  
  explicit horizontally_fragmented_microlensing_lens_plane(const star_db_type& star_db,
                                              const deflection_engine_settings_type& settings,
                                              util::scalar y_fragment_center,
                                              util::scalar shear,
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
    _settings(settings)
  {
    _current_fragment_center = {0.0, _y_fragment_center};
    
    load_new_fragment(_current_fragment_center);
  }
  
  inline util::vector2 lensing_transformation(const util::vector2& lens_plane_pos) const
  {
    if(!is_within_fragment_range(lens_plane_pos))
      load_new_fragment(get_next_fragment_center(lens_plane_pos));
    
    util::vector2 result = _current_fragment->lensing_transformation(lens_plane_pos);
    // Correct result
    util::sub(result, get_deflection_correction(lens_plane_pos));
    
    
    return result;
  }
  
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
  
  util::scalar get_shear() const
  {
    return _shear;
  }

  util::scalar get_smooth_matter_fraction() const
  {
    return _smooth_matter_fraction;
  }

  util::scalar get_sigma_smooth() const
  {
    return _sigma_smooth;
  }

  util::scalar get_sigma_star() const
  {
    return _sigma_star;
  }

  util::scalar get_mean_surface_density() const
  {
    return _sigma_smooth + _sigma_star;
  }
  
  /// Estimates the coordinates of a rectangle in the lens plane that is mapped
  /// to a given rectangle in the source plane
  /// @param source_plane_coordinates A matrix containing the coordinates of the
  /// source plane rectangle
  /// @param out A matrix into which the resulting coordinates in the lens plane
  /// will be written.
  void estimate_mapped_region_coordinates(const util::matrix_nxn<util::vector2, 2>& source_plane_coordinates,
                                          util::matrix_nxn<util::vector2, 2>& out) const
  {
    // Get shooting region, but ignore sigma_star.
    this->_current_fragment->estimate_mapped_region_coordinates(source_plane_coordinates, out, _sigma_star);
    
    /*
    util::matrix_nxn<util::vector2, 2> fragment_corners = source_plane_coordinates;
    fragment_corners[0][0][0] = _current_fragment_center[0] - 0.5 * _fragment_size;
    fragment_corners[0][1][0] = _current_fragment_center[0] - 0.5 * _fragment_size;
    fragment_corners[1][0][0] = _current_fragment_center[0] + 0.5 * _fragment_size;
    fragment_corners[1][1][0] = _current_fragment_center[0] + 0.5 * _fragment_size;
    for(std::size_t i = 0; i < 2; ++i)
    {
      for(std::size_t j = 0; j < 2; ++j)
      {
        // Add sigma star contribution based on the current fragment
        util::vector2 correction = fragment_corners[i][j];
        util::scale(correction, 1.0 / _sigma_star);
        util::add(out[i][j], correction);
      }
    }*/
    
  }
  
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

    return hole_contribution;
  }
  
  inline bool is_within_fragment_range(const util::vector2& point) const
  {
    return std::abs(point[0] - _current_fragment_center[0]) <= _half_fragment_size;
  }
  
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
    
    _sigma_star = _current_fragment->get_sigma_star();
    _smooth_matter_fraction = _current_fragment->get_smooth_matter_fraction();
  }
  
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

