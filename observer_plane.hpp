/* 
 * File:   observer_plane.hpp
 * Author: aksel
 *
 * Created on 8. Dezember 2014, 00:31
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

