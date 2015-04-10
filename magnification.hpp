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

#ifndef MAGNIFICATION_HPP
#define	MAGNIFICATION_HPP

#include "lensing_jacobian.hpp"
#include "util.hpp"
#include "geometry.hpp"
#include "ray.hpp"

namespace nanolens{


namespace magnification{

class by_lensing_jacobian
{
public:
  by_lensing_jacobian(util::scalar accuracy)
  : _accuracy(accuracy) {}
  

  template <class SystemType>
  util::scalar get_magnification(const SystemType& sys,
                                 const util::vector2& pixel_position,
                                 const util::vector2& image_position) const
  {
    lensing_jacobian<SystemType> jacobian(image_position, sys, _accuracy);
    
    return 1.0 / std::abs(jacobian.det());
  }  
  
private:
  util::scalar _accuracy;
  
};

template<std::size_t N_polygon_vertices>
class by_polygon_area
{
public:
  by_polygon_area(util::scalar polygon_radius)
  : _polygon_radius(polygon_radius), _template_shape(util::vector2({0.0,0.0}), polygon_radius)
  {}
  
  
  template <class SystemType>
  util::scalar get_magnification(const SystemType& sys, 
                                         const util::vector2& pixel_position,
                                         const util::vector2& image_position) const
  {
    util::vector2 angle = image_position;
    util::sub(angle, pixel_position);
    util::scale(angle, 1.0 / sys.get_deflector().distance_to_previous_plane());
      
    ray_bundle<N_polygon_vertices> bundle(_template_shape, pixel_position, angle);
    sys.traverse(bundle);
      
    return bundle.get_magnification();
  }  
    
private:
  util::scalar _polygon_radius;
  geometry::equilateral_polygon<N_polygon_vertices> _template_shape;
};

}

}

#endif	/* MAGNIFICATION_HPP */

