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

namespace nanolens
{

  //template<class Deflection_angle_method>
  class lens_plane : public plane
  {
  public:
    
    
    typedef std::vector<star>::iterator star_iterator;
    typedef std::vector<star>::const_iterator const_star_iterator;
    
    lens_plane()
    : plane(1.0)
    {}
    
    explicit lens_plane(const std::vector<star>& deflectors,
                        util::scalar distance_to_prev,
                        util::scalar shear = 0.,
                        util::scalar sigma_smooth = 0.)
    : _deflectors(deflectors), plane(distance_to_prev), 
    _tree(deflectors, 0.6),
    _shear(shear),
    _sigma_smooth(sigma_smooth)
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
      _sigma_star = N_stars * mean_mass / avg_radius_squared;
      
      
      this->_smooth_matter_fraction = _sigma_smooth / (_sigma_star + _sigma_smooth);
      
      _shear_and_smooth_matter_x = 1.0 - _shear - _sigma_smooth;
      _shear_and_smooth_matter_y = 1.0 + _shear - _sigma_smooth;
    }
    
    
    inline void get_deflection_angle(const util::vector2& position, util::vector2& result) const
    {
      result = {0.0, 0.0};
      
//      for(std::size_t i = 0; i < _deflectors.size(); ++i)
//      {
//        util::vector2 deflection;
//        _deflectors[i].calculate_deflection_angle(position, deflection);
//        util::add(result, deflection);
//      }
      result = _tree.get_deflection(position);
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
      util::vector2 out = lens_plane_pos;
      out[0] *= _shear_and_smooth_matter_x;
      out[1] *= _shear_and_smooth_matter_y;
      
      util::vector2 deflection;
      get_deflection_angle(lens_plane_pos, deflection);
      util::add(out, deflection);
      
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
    
    template<class RayBundleType>
    const_star_iterator find_star_in_bundle_area(const RayBundleType& bundle) const
    {
      for(auto it = _deflectors.begin(); it != _deflectors.end(); ++it)
        if(bundle.covered_area_contains_point_in_deflection_plane(it->get_position()))
          return it;
      return _deflectors.end();
    }
    
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
    
    util::scalar _shear_and_smooth_matter_x;
    util::scalar _shear_and_smooth_matter_y;
    //deflection_grid _grid;
    barnes_hut_tree _tree;
    //nested_interpolation_grid _grid;
  };
}

#endif	/* LENS_PLANE_HPP */

