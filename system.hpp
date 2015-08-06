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
#include "lens_plane.hpp"
#include "observer_plane.hpp"

namespace nanolens{

template<class Lens_plane_type>
class system
{
public:
  static const std::size_t num_planes = 2;

  explicit system(const std::shared_ptr<Lens_plane_type>& deflector)
  : _observer(new observer_plane()), _deflector(deflector)
  {
  }


  inline const observer_plane& get_observer() const
  { return *_observer; }

  inline const Lens_plane_type& get_deflector() const
  { return *_deflector; }
  

  inline util::vector2 lensing_transformation(const util::vector2& lens_plane_pos) const
  { 
    return _deflector->lensing_transformation(lens_plane_pos);
  }
  
  inline util::vector2 inverse_lensing_transformation(const util::vector2& lens_plane_pos,
                                                      const util::vector2& source_plane_pos) const
  {
    return _deflector->inverse_lensing_transformation(lens_plane_pos, source_plane_pos);
  }
  
private:

  std::shared_ptr<observer_plane> _observer;
  std::shared_ptr<Lens_plane_type> _deflector;

};

}

#endif	/* SYSTEM_HPP */

