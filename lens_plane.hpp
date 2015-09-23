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

#ifndef LENS_PLANE_HPP
#define	LENS_PLANE_HPP

#include <vector>
#include <cmath>
#include "util.hpp"
#include "plane.hpp"
#include "star.hpp"
#include "tree.hpp"
#include "input.hpp"

namespace nanolens
{

class tree_deflector
{
public:
  class settings 
  {
  public:
    
    explicit settings(const configuration& config)
    {
      _opening_angle = config.get_deflection_engine_property("opening_angle", 1.0);
    }
    
    settings(util::scalar opening_angle = 1.0)
    : _opening_angle(opening_angle) {}
    
    util::scalar get_opening_angle() const
    {
      return _opening_angle;
    }
  private:
    util::scalar _opening_angle;
  };
  
  tree_deflector(const std::vector<star>& stars, const settings& config)
  : _tree(stars, config.get_opening_angle())
  {
  }
  
  
  inline void operator()(const util::vector2& position, util::vector2& result) const
  {
    result = _tree.get_deflection(position);
  }
  
private:
  barnes_hut_tree _tree;
};

class exact_deflector
{
public:
  struct settings 
  {
    explicit settings(const configuration&){}
    settings(){}
  };
  
  exact_deflector(const std::vector<star>& stars, const settings&)
  : _stars(stars)
  {}
  
  
  inline void operator()(const util::vector2& position, util::vector2& result) const
  {
    result = {0.0, 0.0};
    
    util::vector2 current_deflection = {0.0, 0.0};
    for(std::size_t i = 0; i < _stars.size(); ++i)
    {
      _stars[i].calculate_deflection_angle(position, current_deflection);
      util::add(result, current_deflection);
    }
  }
  
private:
  std::vector<star> _stars;
};

template<class Deflector_type>
class microlensing_lens_plane : public plane
{
public:
  typedef typename Deflector_type::settings deflector_settings_type;

  typedef std::vector<star>::iterator star_iterator;
  typedef std::vector<star>::const_iterator const_star_iterator;

  explicit microlensing_lens_plane(const std::vector<star>& deflectors,
                      const deflector_settings_type& settings,
                      util::scalar shear = 0.,
                      util::scalar sigma_smooth = 0.,
                      util::scalar shear_rotation_angle = 0.)
  : _deflectors(deflectors), 
  _deflection_engine(deflectors, settings),
  _shear(shear),
  _sigma_smooth(sigma_smooth),
  _shear_rotation_angle(shear_rotation_angle)
  {
    // Get mean mass
    util::scalar mean_mass = 0.0;

    for(const star& s : deflectors)
      mean_mass += s.get_mass();
    mean_mass /= deflectors.size();

    // Estimate sigma*
    // First calculate distribution center
    util::vector2 distribution_center = {0.0,0.0};
    for(const star& s : deflectors)
      util::add(distribution_center, s.get_position());
    util::scale(distribution_center, 1.0 / deflectors.size());

    // Estimate average radius
    util::scalar avg_radius_squared = 0.0;
    for(const star& s : deflectors)
    {
      util::vector2 R = s.get_position();
      util::sub(R, distribution_center);
      avg_radius_squared += util::dot(R,R);
    }
    avg_radius_squared /= deflectors.size();

    util::scalar N_stars = 0.0;
    // Count stars within average radius
    for(const star& s: deflectors)
    {
      util::vector2 R = s.get_position();
      util::sub(R, distribution_center);
      if(util::dot(R,R) <= avg_radius_squared)
        ++N_stars;
    }
    //util::scalar area_radius_squared = 
    //  avg_radius_squared / util::square(calculate_einstein_radius(1.0,1.0,2.0));
    if(!deflectors.empty())
      _sigma_star = N_stars * mean_mass / avg_radius_squared;
    else _sigma_star = 0.0;


    this->_smooth_matter_fraction = _sigma_smooth / (_sigma_star + _sigma_smooth);

    // Set up shear matrix
    
    util::scalar shear_angle_rad = 
      _shear_rotation_angle * 2.0 * boost::math::constants::pi<util::scalar>() / 360.0;
    
    util::scalar cos_phi = std::cos(shear_angle_rad);
    util::scalar sin_phi = std::sin(shear_angle_rad);
    util::scalar cos2_phi = util::square(cos_phi);
    util::scalar sin2_phi = util::square(sin_phi);
    util::scalar A = 1.0 - _shear;
    util::scalar B = 1.0 + _shear;
    
    _shear_matrix = {{{A * cos2_phi + B * sin2_phi, (A - B)* cos_phi * sin_phi}, 
                      {(A -B )* cos_phi * sin_phi, A * sin2_phi + B * cos2_phi}}};

    _shear_matrix[0][0] -= _sigma_smooth - _sigma_star;
    _shear_matrix[1][1] -= _sigma_smooth - _sigma_star;

    //_shear_matrix[0][0] -= _sigma_smooth;
    //_shear_matrix[1][1] -= _sigma_smooth;

  }


  inline void get_deflection_angle(const util::vector2& position, util::vector2& result) const
  {
    this->_deflection_engine(position, result);
  }

  // Estimate for the radius of the mass distribution in the plane
  util::scalar get_radius_estimate() const
  {
    util::scalar result = 0.0;

    for(std::size_t i = 0; i < _deflectors.size(); ++i)
    {
      util::scalar distance = std::sqrt(util::dot(_deflectors[i].get_position(),
                                                  _deflectors[i].get_position()));

      if(distance > result)
        result = distance;
    }

    if(result == 0.0)
      return 1.0;

    return result;
  }

  // Transform from a position in the lens plane to a position in the source plane
  // at unit distance
  inline util::vector2 lensing_transformation(const util::vector2& lens_plane_pos) const
  {
    util::vector2 out;
    util::matrix_vector_mult(_shear_matrix, lens_plane_pos, out);

    util::vector2 deflection;
    get_deflection_angle(lens_plane_pos, deflection);
    util::add(out, deflection);

    return out;
  }

  inline util::vector2 inverse_lensing_transformation(const util::vector2& lens_plane_pos,
                                                  const util::vector2& source_plane_pos) const
  {
    util::vector2 out = source_plane_pos;

    util::vector2 lensing_transformation_result = lensing_transformation(lens_plane_pos);

    util::sub(out, lensing_transformation_result);

    return out;
  }

  inline star_iterator begin_stars()
  { return _deflectors.begin(); }

  inline const_star_iterator begin_stars() const
  { return _deflectors.begin(); }

  inline star_iterator end_stars()
  { return _deflectors.end(); }

  inline const_star_iterator end_stars() const
  { return _deflectors.end(); }

  inline std::size_t num_stars() const
  { return _deflectors.size(); }

  util::scalar find_squared_distance_to_nearest_star(const util::vector2& position) const
  {
    auto nearest_star = get_nearest_star(position);

    if(nearest_star == _deflectors.end())
      return std::numeric_limits<util::scalar>::max();

    util::vector2 diff_vector = position;
    util::sub(diff_vector, nearest_star->get_position());

    return util::dot(diff_vector, diff_vector);
  }

  util::vector2 get_nearest_star_position(const util::vector2& position) const
  {
    auto nearest_star = get_nearest_star(position);
    if(nearest_star != _deflectors.end())
      return nearest_star->get_position();
    else return {std::numeric_limits<util::scalar>::max(),
                 std::numeric_limits<util::scalar>::max()};

  }

  star& get_star_by_index(std::size_t idx)
  { return _deflectors[idx]; }

  const star& get_star_by_index(std::size_t idx) const
  { return _deflectors[idx]; }

  std::size_t get_num_stars() const
  { return _deflectors.size(); }

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

  void obtain_properties_set(std::map<std::string, util::scalar>& out) const
  {
    out.clear();

    out["N*"] = this->get_num_stars();
    out["sigma* (estimated)"] = this->get_sigma_star();
    out["sigma smooth"] = this->get_sigma_smooth();
    out["sigma total"] = this->get_mean_surface_density();
    out["shear"] = this->get_shear();
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
    estimate_mapped_region_coordinates(source_plane_coordinates, out, _sigma_star);
  }
  
  /// Estimates the coordinates of a rectangle in the lens plane that is mapped
  /// to a given rectangle in the source plane. This version uses a user supplied
  /// value for sigma star.
  /// @param source_plane_coordinates A matrix containing the coordinates of the
  /// source plane rectangle
  /// @param out A matrix into which the resulting coordinates in the lens plane
  /// will be written.
  /// @param sigma_star The stellar density that shall be used.
  void estimate_mapped_region_coordinates(const util::matrix_nxn<util::vector2, 2>& source_plane_coordinates,
                                        util::matrix_nxn<util::vector2, 2>& out,
                                        util::scalar sigma_star) const
  {
    // out = (shear_matrix - diag(sigma_total))^-1 * x
    // (The smooth matter term is already contained in the shear matrix)
    util::matrix_nxn<util::scalar, 2> lensing_transformation = _shear_matrix;
    lensing_transformation[0][0] -= sigma_star;
    lensing_transformation[1][1] -= sigma_star;
    
    util::matrix_nxn<util::scalar, 2> inverted_lensing_transformation;
    numeric::invert_matrix(lensing_transformation, inverted_lensing_transformation);
    
    for(std::size_t i = 0; i < 2; ++i)
      for(std::size_t j = 0; j < 2; ++j)
        util::matrix_vector_mult(inverted_lensing_transformation, source_plane_coordinates[i][j],
                                 out[i][j]);

  }
private:
  std::vector<star>::const_iterator get_nearest_star(const util::vector2& position) const
  {
    auto nearest_star = _deflectors.end();
    util::scalar min_distance = std::numeric_limits<util::scalar>::max();

    for(auto it = _deflectors.begin(); it != _deflectors.end(); ++it)
    {
      util::vector2 diff_vector = position;
      util::sub(diff_vector, it->get_position());

      util::scalar dist = util::dot(diff_vector, diff_vector);

      if(dist < min_distance && dist != 0.0)
      {
        min_distance = dist;
        nearest_star = it;
      }
    }

    return nearest_star;
  }

  std::size_t get_nearest_star_index(const util::vector2& position) const
  {
    std::size_t idx = _deflectors.size();
    util::scalar min_distance = std::numeric_limits<util::scalar>::max();

    for(std::size_t i = 0; i < _deflectors.size(); ++i)
    {
      util::vector2 diff_vector = position;
      util::sub(diff_vector, _deflectors[i].get_position());

      util::scalar dist = util::dot(diff_vector, diff_vector);

      if(dist < min_distance && dist != 0.0)
      {
        min_distance = dist;
        idx = i;
      }
    }

    return idx;
  }

  std::vector<star> _deflectors;

  util::scalar _shear;
  util::scalar _smooth_matter_fraction;

  util::scalar _sigma_star;
  util::scalar _sigma_smooth;

  util::scalar _shear_rotation_angle;
  
  util::matrix_nxn<util::scalar, 2> _shear_matrix;

  Deflector_type _deflection_engine;
};
}

#endif	/* LENS_PLANE_HPP */

