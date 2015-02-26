/* 
 * File:   plane.hpp
 * Author: aksel
 *
 * Created on 9. Dezember 2014, 17:41
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

