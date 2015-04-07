/* 
 * File:   ray.hpp
 * Author: aksel
 *
 * Created on 7. Dezember 2014, 14:55
 */

#ifndef RAY_HPP
#define	RAY_HPP

#include <vector>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "util.hpp"
#include "lens_plane.hpp"
#include "geometry.hpp"

namespace nanolens
{
  
  class ray
  {
  public:
    ray()
    : _ready(false)
    {}
    
    ray(const util::vector2& initial_position, const util::vector2& angle)
    : _current_angle(angle), _current_position(initial_position),
      _initial_position(initial_position), _initial_angle(angle),
      _current_reference_position(initial_position),
      _deflection_angle({0.0, 0.0}), _deflection_plane_position({0.0,0.0}),
      _ready(true)
    {}
    
    inline void deflect(const util::vector2& deflection_angle)
    {
      util::add(_current_angle, deflection_angle);
      util::add(_deflection_angle, deflection_angle);
    }
    
    template<typename PlaneType>
    inline void propagate(const PlaneType& destination_plane)
    {
      perform_propagation(destination_plane);
    }
    
    inline void propagate(const lens_plane& destination_plane)
    {
      perform_propagation(destination_plane);
      
      _deflection_plane_position = _current_position;
    }
    
    void shift_position(const util::vector2& shift)
    {
      util::add(_current_position, shift);
    }
    
    const util::vector2& get_impact_position() const
    {
      return _current_position;
    }
    
    const util::vector2& get_impact_reference_position() const
    {
      return _current_reference_position;
    }
    
    const util::vector2& get_initial_position() const
    {
      return _initial_position;
    }
    
    const util::vector2& get_initial_angle() const
    {
      return _initial_angle;
    }
    
    const util::vector2& get_deflection() const
    {
      return _deflection_angle;
    }
    
    const util::vector2& get_deflection_plane_impact_position() const
    {
      return _deflection_plane_position;
    }
    
    const util::vector2& get_current_angle() const
    {
      return _current_angle;
    }
    
    void set_position(const util::vector2& pos)
    {
      _current_position = pos;
    }
    
    void set_angle(const util::vector2& angle)
    {
      _current_angle = angle;
    }
    
    void reset()
    {
      _current_position = _initial_position;
      _current_angle = _initial_angle;
      _current_reference_position = _initial_position;
    }
    
    bool is_ready() const
    {
      return _ready;
    }
  private:
    template<class PlaneType>
    inline void perform_propagation(const PlaneType& destination_plane)
    {
      util::scalar distance = destination_plane.distance_to_previous_plane();
      
      util::scale_add(_current_position, _current_angle, distance);
      util::scale_add(_current_reference_position, _initial_angle, distance);
      
      util::vector2 deflection_angle;
      destination_plane.get_deflection_angle(_current_position, deflection_angle);
      
      deflect(deflection_angle);
    }
    
    bool _ready;
    util::vector2 _current_position;
    util::vector2 _current_angle;
    util::vector2 _initial_position;
    util::vector2 _initial_angle;
    util::vector2 _current_reference_position;
    util::vector2 _deflection_angle;
    util::vector2 _deflection_plane_position;
  };
  

  template<std::size_t N_rays>
  class ray_bundle
  {
  public:
    ray_bundle() = default;
    
    /// Create new ray bundle
    /// @param template_polygon A model for the shape of the polygon. Its center
    /// must be at (0, 0) and its radius will equal the radius of the ray_bundle
    ray_bundle(const geometry::polygon<N_rays>& template_polygon,
              const util::vector2& position,
              const util::vector2& angle)
    : _shape(template_polygon), _initial_area(0.0), _travelled_distance(0.0)
    {
      _shape.shift_coordinates(angle);
     
      for(std::size_t i = 0; i < _rays.size(); ++i)
        _rays[i] = ray(position, _shape[i]);
      
      _initial_area = _shape.area();
      
    }
    
    util::scalar get_magnification() const
    {
      util::scalar reference_area = _travelled_distance * _initial_area;
      
      geometry::polygon<N_rays> final_shape;
      for(std::size_t i = 0; i < final_shape.num_vertices(); ++i)
        final_shape[i] = _rays[i].get_impact_position();
      
      util::scalar current_area = final_shape.area();
      
      return reference_area / current_area;
    }
    
    const ray& get_ray(std::size_t ray_id) const
    {
      return *(_rays[ray_id]);
    }
    
    template<typename PlaneType>
    inline void propagate(const PlaneType& destination_plane)
    {
      for(std::size_t i = 0; i < _rays.size(); ++i)
        _rays[i].propagate(destination_plane);
      
      _travelled_distance += destination_plane.distance_to_previous_plane();
    }
  private:
 
    std::array<ray, geometry::polygon<N_rays>::num_vertices()> _rays;
    util::scalar _initial_area;
    util::scalar _travelled_distance;
    
    geometry::polygon<N_rays>  _shape;
  };
  
}

#endif	/* RAY_HPP */

