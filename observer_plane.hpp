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

#ifndef OBSERVER_PLANE_HPP
#define	OBSERVER_PLANE_HPP

#include "util.hpp"
#include "plane.hpp"

namespace nanolens
{
  class observer_plane : public plane
  {
  public:
    explicit observer_plane(const util::vector2& observer_position,
                            util::scalar distance_to_prev)
    : _pos(observer_position), plane(distance_to_prev)
    {}
    

    template<typename RayBundleType>
    bool is_hit(const RayBundleType& bundle) const
    {
      return bundle.covered_area_contains_point(_pos);
    }
    
    const util::vector2& get_observer_position() const
    {
      return _pos;
    }
  private:
    util::vector2 _pos;
  };
}

#endif	/* OBSERVER_PLANE_HPP */

