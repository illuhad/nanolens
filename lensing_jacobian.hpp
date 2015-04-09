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

#ifndef LENSING_JACOBIAN_HPP
#define	LENSING_JACOBIAN_HPP

#include "numeric.hpp"

namespace nanolens{

template<class SystemType>
class lensing_jacobian : public numeric::jacobian<util::scalar, 2>
{
public:
  lensing_jacobian(const util::vector2& pos, const SystemType& sys, util::scalar accuracy)
  : numeric::jacobian<util::scalar, 2>(pos, accuracy,
    [&](const util::vector2& lens_plane_pos) -> util::vector2
    {
      util::vector2 source_plane_pos = sys.ray_function(lens_plane_pos);
      
      util::scalar dL_over_dS = sys.get_deflector().distance_to_previous_plane() 
        / sys.get_distance_from_source_to_observer();
      
      //convert to angle
      util::scale(source_plane_pos, dL_over_dS);
      
      return source_plane_pos;
    })
  {}
  
  inline util::scalar det() const
  {
    return (*this)[0][0] * (*this)[1][1] 
        - (*this)[0][1] * (*this)[1][0];
  }
};

}

#endif	/* LENSING_JACOBIAN_HPP */

