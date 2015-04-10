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

#ifndef SYSTEM_HPP
#define	SYSTEM_HPP

#include <vector>
#include <memory>
#include <fstream>
#include <random>
#include <boost/mpi.hpp>
#include "util.hpp"
#include "ray.hpp"
#include "lens_plane.hpp"
#include "observer_plane.hpp"

namespace nanolens{

class system
{
public:
  static const std::size_t num_planes = 2;

  explicit system(const std::vector<star>& deflectors,
                  const std::array<util::scalar, num_planes>& plane_distances = {1.0, 1.0})
  : _deflector(new lens_plane(deflectors, plane_distances[1])),
    _observer(new observer_plane({0.0, 0.0}, plane_distances[0]))
  {
    init();
  }


  explicit system(const std::string& star_file_,
                  const std::array<util::scalar, num_planes>& plane_distances = {1.0, 1.0})
  : _observer(new observer_plane({0.0, 0.0}, plane_distances[1]))
  {
    init_einstein_radius(plane_distances[1],
                         plane_distances[0],
                         plane_distances[0] + plane_distances[1]);

    util::scalar lens_plane_einstein_radius = get_einstein_radius() 
            * _observer->distance_to_previous_plane();

    std::ifstream input_file;

    input_file.open(star_file_.c_str());

    std::vector<star> stars;

    if(input_file.is_open())
    {
      while(input_file.good())
      {
        util::vector2 position = {0.0, 0.0};
        util::scalar mass = 0.0;

        input_file >> position[0];
        input_file >> position[1];
        input_file >> mass;

        util::scale(position, lens_plane_einstein_radius);

        star new_star(position, mass);
        stars.push_back(new_star);
      }
    }

    _deflector = std::shared_ptr<lens_plane>(new lens_plane(stars, plane_distances[0]));
    
    init();


  }

  inline const observer_plane& get_observer() const
  { return *_observer; }

  inline const lens_plane& get_deflector() const
  { return *_deflector; }

  inline util::scalar get_distance_from_source_to_observer() const
  {
    return _deflector->distance_to_previous_plane()
            + _observer->distance_to_previous_plane();
  }

  template<typename RayType>
  inline void traverse(RayType& ray) const
  {
    ray.propagate(*_deflector);
    ray.propagate(*_observer);
  }

  inline util::scalar get_einstein_radius() const
  {
    return _einstein_radius;
  }

  system get_empty_system() const
  {
    return system(std::vector<star>(), 
      {_deflector->distance_to_previous_plane(),_observer->distance_to_previous_plane()});
  }
  
  /// Obtain the total deflection angle
  inline util::vector2 alpha(const util::vector2& lens_plane_position) const
  {
    util::vector2 out;
    _deflector->get_deflection_angle(lens_plane_position, out);
    return out;
  }
  
  inline util::vector2 ray_function(const util::vector2& lens_plane_pos) const
  {
    // ray equation:
    // 0 = (1+dL/dLS)ksi + alpha(ksi)*dL - pos * dL/dLS - pos_observer
    // ray function:
    // pos = (1.0 + dLS/dL)ksi + alpha(ksi) * d_LS - p_obs * d_LS/d_L
    
    util::vector2 out = lens_plane_pos;
    
    util::scale(out, 1.0 + this->_dLS_over_dL);
    
    util::vector2 deflection = alpha(lens_plane_pos);
    util::scale(deflection, this->get_deflector().distance_to_previous_plane());
    
    util::add(out, deflection);
    
    util::vector2 obs_position_term = _observer->get_observer_position();
    util::scale(obs_position_term, this->_dLS_over_dL);    
    util::sub(out, obs_position_term);
    
    return out;
  }
  
private:
  void init_einstein_radius(util::scalar d_ls, util::scalar d_l, util::scalar d_s)
  {
    _einstein_radius = lens_plane::calculate_einstein_radius(d_ls, d_l, d_s);
  }

  void init()
  {
    init_einstein_radius(_deflector->distance_to_previous_plane(),
            _observer->distance_to_previous_plane(),
            get_distance_from_source_to_observer());
    
        
    _dL_over_dLS = _observer->distance_to_previous_plane() /
                   _deflector->distance_to_previous_plane();
    
    _dLS_over_dL = _deflector->distance_to_previous_plane() /
                   _observer->distance_to_previous_plane();
  }
  
  util::scalar _dL_over_dLS;
  util::scalar _dLS_over_dL;
  std::shared_ptr<observer_plane> _observer;
  std::shared_ptr<lens_plane> _deflector;

  util::scalar _einstein_radius;
};

}

#endif	/* SYSTEM_HPP */

