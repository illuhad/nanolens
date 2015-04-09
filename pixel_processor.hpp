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

#ifndef PIXEL_PROCESSOR_HPP
#define	PIXEL_PROCESSOR_HPP

#include <array>
#include <vector>
#include "util.hpp"
#include "image_finder.hpp"
#include "ray.hpp"
#include "numeric.hpp"
#include "magnification.hpp"

namespace nanolens{

template<class MagnificationCalculatorType>
class pixel_processor
{
public:
  pixel_processor(util::scalar accuracy)
  : _accuracy(accuracy),
    _magnification_calculator(accuracy)
  {}
  
  template<class SystemType>
  util::scalar get_pixel_magnification(const util::vector2& pixel_position,
                                       SystemType& sys,
                                      image_finder<SystemType>& img_finder)
  {
    std::vector<util::vector2> image_positions;
    img_finder.get_images(pixel_position, image_positions);
    
    
    util::scalar magnification = 0.0;
    
    for(const util::vector2& img : image_positions)
      magnification += _magnification_calculator.get_magnification(sys, pixel_position, img);
    
    return magnification;
  }
  
private:
  MagnificationCalculatorType _magnification_calculator;
  util::scalar _accuracy;
};

}

#endif	/* PIXEL_PROCESSOR_HPP */

