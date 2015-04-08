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

#ifndef PLANE_HPP
#define	PLANE_HPP

#include "util.hpp"

namespace nanolens
{
  class plane
  {
  public:
    explicit plane(util::scalar distance)
    : _distance_to_prev(distance)
    {}
    
    virtual ~plane(){}
    
    inline util::scalar distance_to_previous_plane() const
    {
      return _distance_to_prev;
    }
    
    virtual void get_deflection_angle(const util::vector2& position, 
                                      util::vector2& result) const
    {
      util::assign(result, 0.0);
    }
  protected:
    util::scalar _distance_to_prev;
  };
}

#endif	/* PLANE_HPP */

