/* 
 * File:   screen.hpp
 * Author: aksel
 *
 * Created on 9. April 2015, 23:40
 */

#ifndef SCREEN_HPP
#define	SCREEN_HPP

#include "util.hpp"
#include "status.hpp"

namespace nanolens {

class screen_descriptor
{
public:
  screen_descriptor(const std::array<std::size_t, 2>& num_pixels,
                    const util::vector2& screen_size,
                    const util::vector2& screen_position)
  : _num_pixels(num_pixels), _screen_position(screen_position), _physical_size(screen_size)
  {
    util::vector2 half_size = _physical_size;
    util::scale(half_size, 0.5);
    
    _pixel_sizes = {_physical_size[0] / static_cast<util::scalar>(num_pixels[0]),
                    _physical_size[1] / static_cast<util::scalar>(num_pixels[1])};
    
    _min_corner = _screen_position;
    util::sub(_min_corner, half_size);
    
    _max_corner = _min_corner;
    util::add(_max_corner, _physical_size);
    
    util::vector2 half_pixel_sizes =  _pixel_sizes;
    util::scale(half_pixel_sizes, 0.5);
    
    _pixel_origin = _min_corner;
    util::add(_pixel_origin, half_pixel_sizes);
  }
  
  const std::array<std::size_t, 2>& get_num_pixels() const
  { return _num_pixels; }
  
  const util::vector2& get_physical_size() const
  { return _physical_size; }
  
  const util::vector2& get_screen_position() const
  { return _screen_position; }
  
  const util::vector2& get_corner_of_min_extent() const
  { return _min_corner; }
  
  const util::vector2& get_corner_of_max_extent() const
  { return _max_corner; }
  
  inline util::vector2 get_pixel_coordinates(const std::array<std::size_t, 2>& pixel_index) const
  {
    assert(pixel_index[0] < _num_pixels[0]);
    assert(pixel_index[1] < _num_pixels[1]);
    
    util::vector2 result = _pixel_origin;
    result[0] += _pixel_sizes[0] * static_cast<util::scalar>(pixel_index[0]);
    result[1] += _pixel_sizes[1] * static_cast<util::scalar>(pixel_index[1]);
    
    return result;
  }
  
private:
  util::vector2 _pixel_origin;
  
  util::vector2 _min_corner;
  util::vector2 _max_corner;
  
  util::vector2 _pixel_sizes;
  std::array<std::size_t, 2> _num_pixels;
  util::vector2 _physical_size;
  util::vector2 _screen_position;
};

}


#endif	/* SCREEN_HPP */

