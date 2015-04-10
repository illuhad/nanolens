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

#ifndef STAR_HPP
#define	STAR_HPP

namespace nanolens {

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

}

#endif	/* STAR_HPP */

