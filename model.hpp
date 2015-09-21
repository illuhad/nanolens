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

#ifndef MODEL_HPP
#define	MODEL_HPP

#include "util.hpp"
#include "input.hpp"

namespace nanolens{
namespace models{

template<class Return_type, class Vector_type>
class basic_model
{
public:
  virtual ~basic_model(){}
  
  virtual Return_type operator()(const Vector_type& x) const = 0;
  virtual util::scalar get_evaluation_diameter() const = 0;
};

class shakura_sunyaev : public basic_model<util::scalar, util::vector2>
{
public:
  shakura_sunyaev(util::scalar radius)
  : _radius(radius) {}
  
  shakura_sunyaev(const configuration& config,
                  const configuration::config_node& node)
  {
    _radius = config.get_config_node_property(node, "radius", 1.0);
    _shift = _radius / 100.0;
  }
  
  virtual ~shakura_sunyaev(){}
  
  virtual util::scalar operator()(const util::vector2& x) const
  {
    // I(R) ~ 1/[exp((R/r)^3/4 -1]
    util::vector2 x_shifted = x;
    x_shifted[0] += _shift;
    x_shifted[1] += _shift;
    
    util::vector2 zero = {static_cast<util::scalar>(0.0), 
                          static_cast<util::scalar>(0.0)};
    if(x == zero)
      return std::numeric_limits<util::scalar>::max();
    
    util::scalar r = std::sqrt(util::dot(x_shifted,x_shifted));
    return 1.0 / (std::exp(std::pow(r / _radius, 0.75)) - 1);
  }
  
  virtual util::scalar get_evaluation_diameter() const
  {
    return 10 * _radius;
  }

  
private:
  util::scalar _radius;
  util::scalar _shift;
};

} // models
} // nanolens


#endif	/* MODEL_HPP */

