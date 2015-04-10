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

#ifndef IMAGE_FINDER_HPP
#define	IMAGE_FINDER_HPP

#include <complex>
#include <functional>
#include "numeric.hpp"
#include "util.hpp"
#include "status.hpp"
#include "geometry.hpp"
#include "screen.hpp"
#include "magnification.hpp"

namespace nanolens{

template<class SystemType>
class image_finder
{
public:
  typedef std::function<void(const status_info&)> status_handler_type;
  
  image_finder(const SystemType* sys, util::scalar accuracy, const status_handler_type& handler)
  : _system(sys), _accuracy(accuracy), _handler(handler)
  {
    assert(sys != nullptr);
  }
  
  virtual ~image_finder()
  {}
  
  
  virtual void get_images(const util::vector2& source_plane_pos,
                          std::vector<util::vector2>& out) = 0;
protected:
  
  const SystemType* _system;
  const util::scalar _accuracy;
  status_handler_type _handler;
};

namespace image_finders {

template<class SystemType>
class inversion_table : public image_finder<SystemType>
{
public:
  typedef numeric::function_inverter<util::scalar, 2, util::scalar, 2> inverter_type;
  typedef std::shared_ptr<inverter_type> inverter_ptr_type;
  
  static const std::size_t benchmark_size = 10000;
  
  inversion_table(const SystemType* sys,
                               const util::vector2& sampling_center,
                               const util::scalar sampling_radius,
                               const util::vector2& physical_source_plane_size,
                               const util::vector2& screen_position,
                               const std::array<std::size_t, 2>& num_pixels,
                               std::size_t num_samples_per_dim,
                               util::scalar accuracy,
                               const boost::mpi::communicator& comm,
                               const typename image_finder<SystemType>::status_handler_type& handler)
  : image_finder<SystemType>(sys, accuracy, handler),
    _differential_delta(0.25 * accuracy),
    _newton_tolerance(accuracy),
    _max_iterations(100),
    _roots_epsilon(4 * accuracy)
  {
    util::vector2 pixel_sizes = {physical_source_plane_size[0] / static_cast<util::scalar>(num_pixels[0]),
                                 physical_source_plane_size[1] / static_cast<util::scalar>(num_pixels[1])};
    
    util::vector2 screen_start = screen_position;
    screen_start[0] -= 0.5 * physical_source_plane_size[0];
    screen_start[1] -= 0.5 * physical_source_plane_size[1];
    
    inverter_type::codomain_vector pixel_min_coordinates = screen_start;
    inverter_type::codomain_vector pixel_max_coordinates = screen_start;
    pixel_max_coordinates[0] += physical_source_plane_size[0];
    pixel_max_coordinates[1] += physical_source_plane_size[1];
    
    inverter_type::domain_vector sampling_start_vector = sampling_center;
    inverter_type::domain_vector sampling_end_vector = sampling_center;
    for(std::size_t i = 0; i < sampling_start_vector.size(); ++i)
    {
      sampling_start_vector[i] -= sampling_radius;
      sampling_end_vector[i] += sampling_radius;
    }
    
    std::size_t avg_pixel_number = 0.5 * (num_pixels[0] + num_pixels[1]);
    
    this->_inverter = inverter_ptr_type(new inverter_type(pixel_min_coordinates,
                                                          pixel_max_coordinates,
                                                          sampling_start_vector,
                                                          sampling_end_vector,
                                                          avg_pixel_number,
                                                          comm));
    
    
    auto function_evaluator = [this](const inverter_type::domain_vector& domain_vec) -> inverter_type::codomain_vector
    {
      return this->_system->ray_function(domain_vec);
    }; 
    
    scheduler schedule = _inverter->create_schedule(
            function_evaluator,
            num_samples_per_dim,
            benchmark_size);
            
    // run inverter
            
    _inverter->run(function_evaluator, num_samples_per_dim, schedule);
    
    _inverter->allcombine();
  }
  
  
  
  virtual void get_images(const util::vector2& source_plane_pos,
                          std::vector<util::vector2>& out)
  {
    out.clear();
    std::vector<util::vector2> start_positions = _inverter->inverse(source_plane_pos);
    
    // Define the function of which we want determine the root using Newton's method
    auto ray_equation = [&](const util::vector2& lens_plane_pos) -> util::vector2
    {
      // left hand side of the ray equation
      util::vector2 lhs = this->_system->ray_function(lens_plane_pos);
      //subtract right hand side from left hand side
      util::sub(lhs, source_plane_pos);
      
      return lhs;
    };
    // run newton
    for(const util::vector2& pos : start_positions)
    {
      numeric::newton<util::scalar, 2> newton2d(pos, _differential_delta, ray_equation);
      
      newton2d.run(_newton_tolerance, _max_iterations);
      
      // Save root if Newton was successful and we haven't saved this root already
      if(newton2d.was_successful())
        if(is_new_root(out, newton2d.get_position()))
          out.push_back(newton2d.get_position());
    }
  }
  
private:
  inline bool is_new_root(const std::vector<util::vector2>& root_list,
                          const util::vector2& root) const
  {
    for(const util::vector2& r : root_list)
    {
      if((std::abs(r[0] - root[0]) + std::abs(r[1] - root[1])) < _roots_epsilon)
        return false;
    }
    return true;
  }
  
  util::scalar _roots_epsilon;
  util::scalar _differential_delta;
  util::scalar _newton_tolerance;
  std::size_t _max_iterations;
  inverter_ptr_type _inverter;
};


template<class SystemType>
class root_tracing : public image_finder<SystemType>
{
public:
  root_tracing(const SystemType* sys,
               util::scalar accuracy,
               const  typename image_finder<SystemType>::status_handler_type& handler,
               image_finder<SystemType>* reliable_image_finder,
               const screen_descriptor& screen)
  : image_finder<SystemType>(sys, accuracy, handler),
    _backend_finder(reliable_image_finder),
    _screen(screen),
    _roots_at_pixel(screen.get_corner_of_min_extent(), 
                    screen.get_corner_of_max_extent(),
                    screen.get_num_pixels()),
    _mag_calculator(accuracy),
    _differential_delta(0.25 * accuracy)
  {
    util::vector2 pixel_coordinates = {0.0, 0.0};
    
    auto ray_equation = [&](const util::vector2& lens_plane_pos) -> util::vector2
    {
      // left hand side of the ray equation
      util::vector2 lhs = this->_system->ray_function(lens_plane_pos);
      //subtract right hand side from left hand side
      util::sub(lhs, pixel_coordinates);
      
      return lhs;
    };
    
    
    for(std::size_t px_x = 0; px_x < _screen.get_num_pixels()[0]; ++px_x)
    {
      std::vector<util::vector2> current_root_list;
      
      for(std::size_t px_y = 0; px_y < _screen.get_num_pixels()[1]; ++px_y)
      {
        pixel_coordinates = _screen.get_pixel_coordinates({px_x, px_y});
        
        if(current_root_list.empty())
          _backend_finder->get_images(pixel_coordinates, current_root_list);
        else
        {
          std::vector<util::vector2>& new_root_list = _roots_at_pixel[pixel_coordinates];
          new_root_list.clear();
          
          // Propagate roots
          for(const util::vector2& root : current_root_list)
          {
            numeric::newton<util::scalar, 2> newton2d(root, 
                                                      _differential_delta, 
                                                      ray_equation);
            
            newton2d.run(this->_accuracy, 50);
            
            if(newton2d.was_successful())
              new_root_list.push_back(newton2d.get_position());
          }
          
          if(caustic_crossing(current_root_list, new_root_list))
            _backend_finder->get_images(pixel_coordinates, new_root_list);

          current_root_list = new_root_list;
        }
        
      }
    }
  }
  
  virtual ~root_tracing(){}
  
  virtual void get_images(const util::vector2& source_plane_pos,
                          std::vector<util::vector2>& out)
  {
    out = _roots_at_pixel[source_plane_pos];
  }
  
private:
  inline bool caustic_crossing(const std::vector<util::vector2>& root_list0,
                        const std::vector<util::vector2>& root_list1) const
  {
    util::scalar mag0 = get_magnification(root_list0);
    util::scalar mag1 = get_magnification(root_list1);

    return std::abs(mag0 - mag1) > 25.0;
  }
  
  inline util::scalar get_magnification(const std::vector<util::vector2>& root_list) const
  {
    util::scalar mag = 0.0;
    
    // magnification by lensing jacobian does not care about the pixel position,
    // hence we just assume 0.0, 0.0
    util::vector2 pixel_position = {0.0, 0.0};
    
    for(const util::vector2& root : root_list)
      mag += _mag_calculator.get_magnification(*(this->_system), pixel_position, root);
      
    return mag;
  }
  
  util::scalar _differential_delta;
  
  magnification::by_lensing_jacobian _mag_calculator;
  
  screen_descriptor _screen;
  util::grid2d<util::scalar, std::vector<util::vector2>> _roots_at_pixel;
  image_finder<SystemType>* _backend_finder;
};


// TODO
template<class SystemType>
class complex_polynomial : public image_finder<SystemType>
{
  struct root
  {
    util::vector2 position;
    unsigned multiplicity;
  };
  
public:
  typedef std::complex<util::scalar> complex_type;
  
  complex_polynomial(const SystemType* sys, util::scalar accuracy,
                     const typename image_finder<SystemType>::status_handler_type& handler)
  : image_finder<SystemType>(sys, accuracy, handler),
    _newton_tolerance(accuracy),
    _differential_delta(0.25 * accuracy)
  {
  }
  
  virtual void get_images(const util::vector2& source_plane_pos,
                          std::vector<util::vector2>& out) const
  {
    out.clear();
    
    std::vector<root> found_roots;
    
    while(!all_roots_found(found_roots))
    {
      // Create random starting point
    }
    
    out.reserve(found_roots.size());
    for(const root& r : found_roots)
      out.push_back(r.position);
  }
  
private:
  
  inline std::size_t get_multiplicity(const util::vector2& root) const
  {
    // TODO
    return 0;
  }

  inline bool all_roots_found(const std::vector<root>& root_list) const
  {
    std::size_t num_weighted_roots = 0;
    for(const root& r : root_list)
      num_weighted_roots += r.multiplicity;
    
    return num_weighted_roots >= this->_system->get_deflector().num_stars();
  }
  
  util::scalar _newton_tolerance;
  util::scalar _differential_delta;
};

template<class SystemType, std::size_t N_start_points>
class newton_crown : public image_finder<SystemType>
{
public:
  static_assert(N_start_points > 0, "newton_crown: At least 1 start point is required");
  
  typedef util::grid2d<util::scalar, bool> grid_type;
  
  newton_crown(const SystemType* sys,
               util::scalar accuracy, 
               const  typename image_finder<SystemType>::status_handler_type& handler)
  : image_finder<SystemType>(sys, accuracy, handler),
    _shape_template(util::vector2({0.0, 0.0}), 1.0),
    _differential_delta(0.25 * accuracy),
    _newton_tolerance(accuracy),
    _max_iterations(50),
    _roots_epsilon(4 * accuracy)
  {
  }
  
  virtual void get_images(const util::vector2& source_plane_pos,
                          std::vector<util::vector2>& out)
  {
    out.clear();
    
    // Define the function of which we want determine the root using Newton's method
    auto ray_equation = [&](const util::vector2& lens_plane_pos) -> util::vector2
    {
      // left hand side of the ray equation
      util::vector2 lhs = this->_system->ray_function(lens_plane_pos);
      //subtract right hand side from left hand side
      util::sub(lhs, source_plane_pos);
      
      return lhs;
    };
    
    for(std::size_t star_idx = 0; star_idx < this->_system->get_deflector().get_num_stars();
      ++star_idx)
    {
      util::vector2 star_position 
        = this->_system->get_deflector().get_star_by_index(star_idx).get_position();
      
      geometry::equilateral_polygon<N_start_points> newton_start_points = 
        _shape_template;
      
      //Scale radius
      util::scalar radius 
        = 0.2 * this->_system->get_deflector().get_distance_to_nearest_star_for_star(star_idx);
      
      for(std::size_t vertex = 0; vertex < newton_start_points.num_hull_vertices(); ++vertex)
        util::scale(newton_start_points[vertex], radius);
      
      newton_start_points.shift_coordinates(star_position);
      
      for(std::size_t vertex_index = 0; 
        vertex_index < newton_start_points.num_hull_vertices();
        ++vertex_index)
      {
        numeric::newton<util::scalar, 2> newton2d(newton_start_points[vertex_index], 
                                                  _differential_delta, ray_equation);
      
        newton2d.run(_newton_tolerance, _max_iterations);

        if(newton2d.was_successful())
        {
          if(is_new_root(out, newton2d.get_position()))
          {
            out.push_back(newton2d.get_position());
          }
        }
      }
    }
  }
private:
  inline bool is_new_root(const std::vector<util::vector2>& root_list, 
                             const util::vector2& root) const
  {
    for(const util::vector2& r : root_list)
    {
      if((std::abs(r[0] - root[0]) + std::abs(r[1] - root[1])) < _roots_epsilon)
        return false;
    }
    return true;
  }
  
  geometry::equilateral_polygon<N_start_points> _shape_template;
  
  util::scalar _differential_delta;
  util::scalar _newton_tolerance;
  util::scalar _max_iterations;
  util::scalar _roots_epsilon;
};

}

}

#endif	/* IMAGE_FINDER_HPP */

