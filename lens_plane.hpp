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

namespace nanolens
{
  
  class lens_plane : public plane
  {
  public:
    
    class star
    {
    public:
      star(const util::vector2& pos, util::scalar mass)
      : _mass(mass), _position(pos), _deflection_constant(0.)
      {
        _deflection_constant = 4 * util::G * _mass / util::square(util::c);
      }
      
      inline void calculate_deflection_angle(const util::vector2& position,
                                            util::vector2& result) const
      { 
        result = this->_position;
        util::sub(result, position);
        
        util::scalar squared_norm = util::dot(result, result);
        
        
        util::scalar factor = _deflection_constant / squared_norm;
        
        util::scale(result, factor);
      }
      
      const util::vector2& get_position() const
      { return _position; }
      
      util::scalar get_mass() const
      { return _mass; }
    private:
      util::vector2 _position;
      util::scalar _mass;
      util::scalar _deflection_constant;
      
    };
    
    typedef std::vector<star>::iterator star_iterator;
    typedef std::vector<star>::const_iterator const_star_iterator;
    
    lens_plane()
    : plane(1.0)
    {}
    
    explicit lens_plane(const std::vector<star>& deflectors,
                        util::scalar distance_to_prev)
    : _deflectors(deflectors), plane(distance_to_prev)
    {
      _distance_to_nearest_star.reserve(_deflectors.size());
      
      for(std::size_t i = 0; i < _deflectors.size(); ++i)
      {
        std::size_t nearest_star_idx = get_nearest_star_index(_deflectors[i].get_position());
        if(nearest_star_idx < _deflectors.size())
        {
          util::vector2 diff = _deflectors[i].get_position();
          util::sub(diff, _deflectors[nearest_star_idx].get_position());
          _distance_to_nearest_star.push_back(util::dot(diff, diff));
        }
        else
          _distance_to_nearest_star.push_back(std::numeric_limits<util::scalar>::max());
      }
    }
    
    void get_deflection_angle(const util::vector2& position, util::vector2& result) const
    {
      util::vector2 deflection;
      result = {0.0, 0.0};
      
      for(std::size_t i = 0; i < _deflectors.size(); ++i)
      {
        _deflectors[i].calculate_deflection_angle(position, deflection);
        util::add(result, deflection);
      }
    }
    
    // Einstein radius for object of mass 1
    static util::scalar calculate_einstein_radius(util::scalar d_ls,
                                           util::scalar d_l,
                                           util::scalar d_s)
    {
      return 2.0 / util::c * std::sqrt(util::G * d_ls / (d_l * d_s));
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
        return 2.0;
      
      return result;
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
    
    util::vector2 get_pseudo_deflection(const util::vector2& position,
                                               const util::scalar max_pull,
                                               const util::scalar cutoff) const
    {
      util::vector2 shift = {0.0,0.0};
      
      std::size_t nearest_star_idx = get_nearest_star_index(position);
      if(nearest_star_idx < _deflectors.size())
      {
        util::vector2 diff = _deflectors[nearest_star_idx].get_position();
        util::sub(diff, position);
        
        util::scalar next_nearest_star_distance =
                _distance_to_nearest_star[nearest_star_idx];
        
        util::scalar dist = util::dot(diff, diff);
        util::scalar decay = 10.0 * 0.5 * next_nearest_star_distance / util::square(cutoff);
        
        util::scale(diff,
                  std::min(max_pull, util::square(cutoff) / (decay * dist)));

        util::add(shift, diff);
      }
      return shift;
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
    std::vector<util::scalar> _distance_to_nearest_star;
    
  };
}

#endif	/* LENS_PLANE_HPP */

