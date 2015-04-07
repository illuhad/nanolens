/* 
 * File:   pixel_processor.hpp
 * Author: aksel
 *
 * Created on 7. April 2015, 05:12
 */

#ifndef PIXEL_PROCESSOR_HPP
#define	PIXEL_PROCESSOR_HPP

#include <array>
#include <vector>
#include "util.hpp"
#include "image_finder.hpp"
#include "ray.hpp"

namespace nanolens{

template<std::size_t N_rays>
class pixel_processor
{
public:
  pixel_processor(util::scalar polygon_radius)
  : _template_shape({0.0, 0.0}, polygon_radius)
  {}
  
  template<class SystemType>
  util::scalar get_pixel_magnification(const util::vector2& pixel_position,
                                       SystemType& sys
                                      image_finder<SystemType>& img_finder) const
  {
    std::vector<util::vector2> image_positions;
    img_finder.get_images(_pos, image_positions);
    
    util::scalar magnification = 0.0;
    for(const util::vector2& img : image_positions)
    {
      // Convert positions in the lens plane to angles
      util::vector2 angle = img;
      util::sub(angle, pixel_position);
      util::scale(angle, 1.0 / sys.get_deflector().distance_to_previous_plane());
      
      ray_bundle<N_rays> bundle(_template_shape, pixel_position, angle);
      sys.traverse(bundle);
      
      magnification += bundle.get_magnification();
    }
    return magnification;
  }
  
private:
  geometry::equilateral_polygon<N_rays> _template_shape;
};

}

#endif	/* PIXEL_PROCESSOR_HPP */

