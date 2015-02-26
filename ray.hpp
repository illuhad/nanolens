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
#include "mesh.hpp"

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
  
  template<class AdaptiveMeshPolicy>
  class ray_bundle
  {
  public:
    static constexpr std::size_t num_rays = 3;
    static constexpr std::size_t num_triangles = 1;
    
    ray_bundle() = default;
    
    ray_bundle(ray* r0, ray* r1, ray* r2)
    {
      assert(r0 != nullptr && r1 != nullptr && r2 != nullptr);
      
      _rays[0] = r0;
      _rays[1] = r1;
      _rays[2] = r2;
    }
    
    util::scalar get_magnification() const
    {
      return retrieve_magnification();
    }
    
    // after one call, following calls must be with the same point as argument
    // since the result of the calculation is cached.
    bool covered_area_contains_point(const util::vector2& point) const
    {
      return _is_hit.retrieve([this, point]() -> util::scalar
      {
        return ray_triangle_contains_point(*_rays[0],
                                     *_rays[1],
                                     *_rays[2],
                                     point);
      });
      
    }
    
    
    bool covered_area_contains_point_in_deflection_plane(const util::vector2& point) const
    {
      const util::vector2& a = _rays[0]->get_deflection_plane_impact_position();
      const util::vector2& b = _rays[1]->get_deflection_plane_impact_position();
      const util::vector2& c = _rays[2]->get_deflection_plane_impact_position();

      return triangle_contains_point(a, b, c, point);
    }

    
    template<typename SystemType>
    util::scalar refinement_rating(const SystemType& sys) const
    {
      util::scalar result = 0.0;

      if(AdaptiveMeshPolicy::refine_if_strong_deflection)
        result += AdaptiveMeshPolicy::strong_deflection_weight * strong_deflection_rating();

      if(AdaptiveMeshPolicy::refine_if_strong_distortion)
        result += AdaptiveMeshPolicy::strong_distortion_weight * strong_distortion_rating();

      if(AdaptiveMeshPolicy::refine_if_strong_magnification)
        result += AdaptiveMeshPolicy::strong_magnification_weight * strong_magnification_rating();

      if(AdaptiveMeshPolicy::refine_if_strong_attenuation)
        result += AdaptiveMeshPolicy::strong_attenuation_weight * strong_attenuation_rating();

      if(AdaptiveMeshPolicy::refine_if_triangle_twisted)
        result += AdaptiveMeshPolicy::triangle_twisted_weight * triangle_twisted_rating();

      if(AdaptiveMeshPolicy::refine_if_star_nearby)
        result += AdaptiveMeshPolicy::star_nearby_weight * star_nearby_rating(sys);

      if(AdaptiveMeshPolicy::refine_if_close_to_observer)
        result += AdaptiveMeshPolicy::close_to_observer_weight * close_to_observer_rating(sys);

      return result;
    }
    
    
    const ray& get_ray(std::size_t ray_id) const
    {
      return *(_rays[ray_id]);
    }
     
  private:
    
    util::scalar strong_distortion_rating() const
    {
      util::scalar limit = AdaptiveMeshPolicy::strong_distortion_maximum;
      
      std::array<util::vector2, 3> angle_differences;
      angle_differences[0] = _rays[0]->get_current_angle();
      angle_differences[1] = _rays[1]->get_current_angle();
      angle_differences[2] = _rays[2]->get_current_angle();
      
      util::sub(angle_differences[0], _rays[1]->get_current_angle());
      util::sub(angle_differences[1], _rays[2]->get_current_angle());
      util::sub(angle_differences[2], _rays[0]->get_current_angle());
      
      util::scalar sum_of_differences = 0.0;
      for(std::size_t j = 0; j < angle_differences.size(); ++j)
        sum_of_differences += util::dot(angle_differences[j], angle_differences[j]);
      
      //std::cout << sum_of_differences << std::endl;
      
      if(sum_of_differences > util::square(limit))
        return 1.0;

      return sum_of_differences / util::square(limit);
    }
    
    template<typename SystemType>
    util::scalar star_nearby_rating(const SystemType& sys) const
    {
      util::vector2 central_deflection_plane_position;
      util::average(_rays[0]->get_deflection_plane_impact_position(),
                    _rays[1]->get_deflection_plane_impact_position(),
                    _rays[2]->get_deflection_plane_impact_position(),
                    central_deflection_plane_position);
      
      util::scalar distance 
              = sys.get_deflector().find_squared_distance_to_nearest_star(
                            central_deflection_plane_position);
      
      util::scalar threshold = AdaptiveMeshPolicy::star_nearby_threshold;
      if(distance < threshold)
      {
        if(distance == 0.0)
          return 1.0;
        
        util::scalar cutoff = util::square(AdaptiveMeshPolicy::star_nearby_close_cutoff);
        
        if(distance < cutoff)
          return 1.0;
        
        return 1.0 - 1.0 / (threshold - cutoff) * (distance - cutoff);
      }
      
      return 0.0;
    }
    
    util::scalar triangle_twisted_rating() const
    {
      // TODO
      return 0.0;
    }
    
    bool strong_deflection_rating() const
    {
      util::scalar max_deflection = AdaptiveMeshPolicy::strong_deflection_maximum;
      
      util::scalar deflection_squared = 0.0;
      
      for(std::size_t i = 0; i < _rays.size(); ++i)
      {
        deflection_squared += 
                util::dot(_rays[i]->get_deflection(), 
                          _rays[i]->get_deflection());
        
      }
      
      return deflection_squared / max_deflection;
    }
    
    util::scalar strong_magnification_rating() const
    {
      util::scalar mag = retrieve_magnification();
      util::scalar max_magnification = AdaptiveMeshPolicy::strong_magnification_maximum;
      
      if(mag > max_magnification)
        return 1.0;
      
      return mag / max_magnification;
    }
    
    util::scalar strong_attenuation_rating() const
    {
      util::scalar mag = retrieve_magnification();
      
      util::scalar cutoff = AdaptiveMeshPolicy::strong_attenuation_cutoff;
      
      if(mag <= 0.0 || mag <= cutoff)
        return 1.0;
      
      return cutoff / mag;
    }
    
    template<class SystemType>
    util::scalar close_to_observer_rating(const SystemType& sys) const
    {
      util::scalar distance_squared = std::numeric_limits<util::scalar>::max();
      util::vector2 difference;
      for(std::size_t i = 0; i < _rays.size(); ++i)
      {
        difference = _rays[i]->get_impact_position();
        util::sub(difference, sys.get_observer().get_observer_position());
        
        util::scalar ray_distance = util::dot(difference, difference);
        
        if(ray_distance < distance_squared)
          distance_squared = ray_distance;
      }
      
      if(distance_squared > util::square(AdaptiveMeshPolicy::close_to_observer_threshold))
        return 0.0;

      if(distance_squared < util::square(AdaptiveMeshPolicy::close_to_observer_cutoff)
        || distance_squared == 0.0)
        return 1.0;
      
      return util::square(AdaptiveMeshPolicy::close_to_observer_cutoff) / distance_squared;
    }
    
    util::scalar area_of_ray_triangle(const ray& r0, const ray& r1, const ray& r2) const
    {
      const util::vector2& p0 = r0.get_impact_position();
      const util::vector2& p1 = r1.get_impact_position();
      const util::vector2& p2 = r2.get_impact_position();
      
      return area_of_triangle(p0, p1, p2);
    }
    
    inline util::scalar area_of_triangle(const util::vector2& p0,
                                  const util::vector2& p1,
                                  const util::vector2& p2) const
    {
      return std::abs(signed_area_of_triangle(p0, p1, p2));
    }
    
    inline util::scalar signed_area_of_triangle(const util::vector2& p0,
                                  const util::vector2& p1,
                                  const util::vector2& p2) const
    {
      return 0.5 * ((p1[0] - p0[0]) * (p2[1] - p0[1]) 
                 - (p2[0] - p0[0]) * (p1[1] - p0[1]));
    }
    
    bool ray_triangle_contains_point(const ray& r0, 
                                     const ray& r1, 
                                     const ray& r2, 
                                     const util::vector2& point) const
    {
      const util::vector2& a = r0.get_impact_position();
      const util::vector2& b = r1.get_impact_position();
      const util::vector2& c = r2.get_impact_position();
     
      return triangle_contains_point(a, b, c, point);
    }
    
    inline bool triangle_contains_point(const util::vector2& a,
                                        const util::vector2& b,
                                        const util::vector2& c,
                                        const util::vector2& point) const
    {
      util::vector2 v0 = c;
      util::vector2 v1 = b;
      util::vector2 v2 = point;
      util::sub(v0,a);
      util::sub(v1,a);
      util::sub(v2,a);
      
      util::scalar dot00, dot01, dot02, dot11, dot12;
      dot00 = util::dot(v0, v0);
      dot01 = util::dot(v0, v1);
      dot02 = util::dot(v0, v2);
      dot11 = util::dot(v1, v1);
      dot12 = util::dot(v1, v2);
      
      util::scalar factor = 1. / (dot00 * dot11 - dot01 * dot01);
 
      // transform to u,v coordinates
      util::scalar u = factor * (dot11 * dot02 - dot01 * dot12);
      util::scalar v = factor * (dot00 * dot12 - dot01 * dot02);

      return u > 0. && v > 0. && ((u + v) < 1.);
    }
    
    inline util::scalar retrieve_magnification() const
    {  
      return _magnification.retrieve([this]() -> util::scalar
      {
        return calculate_magnification();
      });
    }
    
    util::scalar calculate_magnification() const
    {
      util::scalar area = area_of_ray_triangle(*_rays[0],
                                         *_rays[1],
                                         *_rays[2]);
          
      util::scalar reference_area
              = area_of_triangle(_rays[0]->get_impact_reference_position(),
                                 _rays[1]->get_impact_reference_position(),
                                 _rays[2]->get_impact_reference_position());
          
      return reference_area / area;
    }
    
    std::array<ray*, num_rays> _rays;

    mutable util::cached_value<util::scalar> _magnification;
    mutable util::cached_value<util::scalar> _is_hit;
  };
  
  template<class AdaptiveMeshPolicy>
  class ray_dome
  {
  public:
    
    typedef mesh_ptr<ray, ray_bundle<AdaptiveMeshPolicy>> ray_mesh_ptr;
    typedef mesh<ray, ray_bundle<AdaptiveMeshPolicy>> ray_mesh;
    
    static constexpr std::size_t min_rays_per_dim = 2;
    
    ray_dome(const util::vector2& central_position,
             const util::vector2& central_angle,
             util::scalar z_axis_offset,
             const util::scalar& angular_radius,
             unsigned num_rays_per_dim,
             util::scalar refinement_level)
    : _num_rays(0),
      _central_position(central_position),
      _central_angle(central_angle),
      _angular_radius(angular_radius),
      _z_axis_offset(z_axis_offset),
      // Refinement will bring in 2 new refinement levels
      _max_refinement_level(2 * refinement_level)
    {
      if(num_rays_per_dim < min_rays_per_dim)
      {
        _num_rays_per_dim = min_rays_per_dim;
      }
      else _num_rays_per_dim = num_rays_per_dim;
    }
    
    
    /// Sends this dome through the lensing system and calculates the magnification
    template<typename SystemType>
    inline void traverse(const SystemType& sys, bool dump=false)
    {
      _solid_angle_mesh = ray_mesh_ptr(new ray_mesh());
      _solid_angle_mesh->create_square_mesh(_central_angle,
                                            _angular_radius,
                                            _num_rays_per_dim,
                                            0.0);
      
      _magnification = 0.;
      
      // Set up rays and ray bundles
      for(auto current_vertex = _solid_angle_mesh->get_meshdb().begin_vertices();
              current_vertex != _solid_angle_mesh->get_meshdb().end_vertices();
              ++current_vertex)
      {
        util::vector2 ray_initial_position = _central_position;
        util::vector2 ray_initial_angle = (*current_vertex)->get_position();
        
        util::scale_add(ray_initial_position, ray_initial_angle, _z_axis_offset);
        
        util::vector2 lens_plane_impact_position = ray_initial_position;
        util::scale_add(lens_plane_impact_position,
                        ray_initial_angle,
                        sys.get_deflector().distance_to_previous_plane());
        
        // Pull vertices to stars 
        /*
        util::vector2 pull_vector = sys.get_deflector().get_pseudo_deflection(
                lens_plane_impact_position,
                AdaptiveMeshPolicy::vertex_pull_max,
                AdaptiveMeshPolicy::vertex_pull_cutoff);
        
        
        util::vector2 new_vertex_position = (*current_vertex)->get_position();
        util::add(new_vertex_position, pull_vector);
        (*current_vertex)->set_position(new_vertex_position);
        ray_initial_angle = new_vertex_position;
        */
        // Attach a ray to the vertex
        (*current_vertex)->set_payload(ray(ray_initial_position, ray_initial_angle));
        
        // Send ray through the system
        sys.traverse((*current_vertex)->get_payload());
        
      }
      
      _solid_angle_mesh->get_meshdb().for_each_cell([] (const typename ray_mesh::cell_ptr& cell)
      {
        assert(cell);
        ray* r0 = &(cell->get_vertex(0)->get_payload());
        ray* r1 = &(cell->get_vertex(1)->get_payload());
        ray* r2 = &(cell->get_vertex(2)->get_payload());

        assert(r0);
        assert(r1);
        assert(r2);
        
        cell->set_payload(ray_bundle<AdaptiveMeshPolicy>(r0, r1, r2));
      });
      
      // refine cells if neccessary
      for(std::size_t refinement_level = 0;
              refinement_level < _solid_angle_mesh->get_meshdb().get_num_contained_levels(); 
              ++refinement_level)
      {
        if(_solid_angle_mesh->get_meshdb().contains_refinement_level(refinement_level))
        {
              
          // Iterate over all cells of the current refinement level
          _solid_angle_mesh->get_meshdb().for_each_cell(refinement_level,
          [&] (const typename ray_mesh::cell_ptr& cell)
          {
            assert(cell);
            assert(cell->get_refinement_level() == refinement_level);
            
            ray_bundle<AdaptiveMeshPolicy>& bundle = cell->get_payload();
            bool is_hit = sys.get_observer().is_hit(bundle);

            bool refine = false;
            // Find out if we need to refine this cell
            if(refinement_level + 2 <= _max_refinement_level)
            {
              if(is_hit)
                refine = true;
              else
              {
                // evaluate refinement rating
                util::scalar desired_level = bundle.refinement_rating(sys) * _max_refinement_level;
                
                if(std::round(desired_level) > refinement_level)
                  refine = true;
              }
            }
            
            if(refine)
            {
              typename ray_mesh::cell::refined_cells_container new_cells;
              typename ray_mesh::cell::new_vertices_container new_vertices;

              cell->refine(new_cells, new_vertices);

              for(auto&& new_vertex : new_vertices)
              {
                assert(new_vertex);
                // Setup ray and send it through the system
                util::vector2 ray_initial_position = _central_position;
                util::vector2 ray_initial_angle = new_vertex->get_position();
        
                util::scale_add(ray_initial_position, ray_initial_angle, _z_axis_offset);
        
                new_vertex->set_payload(ray(ray_initial_position, ray_initial_angle));
        
                // Send ray through the system
                sys.traverse(new_vertex->get_payload());
                
                //std::cout << new_vertex->get_payload().get_impact_position()[0]
                //<< new_vertex->get_payload().get_impact_position()[0] << std::endl;
              }
              
              for(auto && new_cell : new_cells)
              {
                assert(new_cell);
                ray* r0 = &(new_cell->get_vertex(0)->get_payload());
                ray* r1 = &(new_cell->get_vertex(1)->get_payload());
                ray* r2 = &(new_cell->get_vertex(2)->get_payload());
                
                new_cell->set_payload(ray_bundle<AdaptiveMeshPolicy>(r0, r1, r2));
              }
            }

          });
          
          
        }
      }
      
      std::size_t num_hits = 0;
      // gather magnifications
      _solid_angle_mesh->get_meshdb().for_each_cell(
      [&] (const typename ray_mesh::cell_ptr& cell)
      {
        assert(cell);
        ray_bundle<AdaptiveMeshPolicy>& bundle = cell->get_payload();

        if(dump)
        {
          const ray* r0 = &(bundle.get_ray(0));
          const ray* r1 = &(bundle.get_ray(1));
          const ray* r2 = &(bundle.get_ray(2));
          
          std::cerr << r0->get_initial_angle()[0] << " " << r0->get_initial_angle()[1] << std::endl;
          std::cerr << r1->get_initial_angle()[0] << " " << r1->get_initial_angle()[1] << std::endl;
          std::cerr << r2->get_initial_angle()[0] << " " << r2->get_initial_angle()[1] << std::endl;
          std::cerr << r0->get_initial_angle()[0] << " " << r0->get_initial_angle()[1] << std::endl;
          std::cerr << std::endl << std::endl;

          /*
          std::cerr << r0->get_impact_position()[0] << " " << r0->get_impact_position()[1] << std::endl;
          std::cerr << r1->get_impact_position()[0] << " " << r1->get_impact_position()[1] << std::endl;
          std::cerr << r2->get_impact_position()[0] << " " << r2->get_impact_position()[1] << std::endl;
          std::cerr << r0->get_impact_position()[0] << " " << r0->get_impact_position()[1] << std::endl;
          std::cerr << std::endl << std::endl;*/
          
        }

        bool is_hit = sys.get_observer().is_hit(bundle);

        //if(is_hit)
        //            std::cout << cell->get_refinement_level() << std::endl;
        
        if(is_hit)
        {
          _magnification += bundle.get_magnification();
          ++num_hits;
        }

      });
      
      //std::cout << num_hits << " " << _magnification <<  std::endl;
    }
    
    util::scalar get_magnification() const
    {
      return _magnification;
    }
    
    unsigned get_num_rays() const
    {
      return _solid_angle_mesh->get_meshdb().get_num_vertices();
    }
  private:

    util::scalar _magnification;
    util::scalar _angular_radius;
    
    std::size_t _num_rays_per_dim;
    
    ray_mesh_ptr _solid_angle_mesh;
    
    unsigned _num_rays;
    unsigned _max_refinement_level;
    
    util::vector2 _central_position;
    util::vector2 _central_angle;
    util::scalar& _z_axis_offset;
  };
}

#endif	/* RAY_HPP */

