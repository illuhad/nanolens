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
                  util::scalar shear,
                  util::scalar sigma_smooth,
                  const std::array<util::scalar, num_planes>& plane_distances = {1.0, 1.0})
  : _observer(new observer_plane(plane_distances[0]))
  {

    
    _deflector = std::shared_ptr<lens_plane>(new lens_plane(deflectors, 
                                                            plane_distances[1], 
                                                            shear, 
                                                            sigma_smooth));
    
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

  system get_empty_system() const
  {
    return system(std::vector<star>(), 0.0, 0.0,
      {_deflector->distance_to_previous_plane(),_observer->distance_to_previous_plane()});
  }
  
  inline util::vector2 lensing_transformation(const util::vector2& lens_plane_pos) const
  {
    util::vector2 source_plane_pos = _deflector->lensing_transformation(lens_plane_pos);
    util::scale(source_plane_pos, _deflector->distance_to_previous_plane());
    
    return source_plane_pos;
  }
  
private:

  void init()
  {
    
        
    _dL_over_dLS = _observer->distance_to_previous_plane() /
                   _deflector->distance_to_previous_plane();
    
    _dLS_over_dL = _deflector->distance_to_previous_plane() /
                   _observer->distance_to_previous_plane();
  }
  
  
  util::scalar _dL_over_dLS;
  util::scalar _dLS_over_dL;
  std::shared_ptr<observer_plane> _observer;
  std::shared_ptr<lens_plane> _deflector;

};

}

#endif	/* SYSTEM_HPP */

