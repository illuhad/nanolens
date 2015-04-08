/* 
 * File:   image_finder.hpp
 * Author: aksel
 *
 * Created on 24. März 2015, 03:19
 */

#ifndef IMAGE_FINDER_HPP
#define	IMAGE_FINDER_HPP

#include <complex>
#include "numeric.hpp"
#include "util.hpp"

namespace nanolens{

template<class SystemType>
class image_finder
{
public:
  image_finder(const SystemType* sys, util::scalar accuracy)
  : _system(sys), _accuracy(accuracy)
  {
    assert(sys != nullptr);
    
    _dL_over_dLS = _system->get_observer().distance_to_previous_plane() /
                   _system->get_deflector().distance_to_previous_plane();
    
    _dLS_over_dL = _system->get_deflector().distance_to_previous_plane() /
                   _system->get_observer().distance_to_previous_plane();
  }
  
  virtual ~image_finder()
  {}
  
  
  virtual void get_images(const util::vector2& source_plane_pos,
                          std::vector<util::vector2>& out) = 0;
protected:
  /// Obtain the deflection angle
  inline util::vector2 alpha(const util::vector2& lens_plane_position) const
  {
    util::vector2 out;
    _system->get_deflector().get_deflection_angle(lens_plane_position, out);
    return out;
  }
  
  inline util::vector2 ray_function(const util::vector2& lens_plane_pos) const
  {
    // ray equation:
    // 0 = (1+dL/dLS)ksi + alpha(ksi)*dL - pos * dL/dLS - pos_observer
    // ray function:
    // pos = (1.0 + dLS/dL)ksi + alpha(ksi) * d_LS - p_obs * d_LS/d_L
    
    util::vector2 out = lens_plane_pos;
    
    util::scale(out, 1.0 + this->_dLS_over_dL);
    
    util::vector2 deflection = alpha(lens_plane_pos);
    util::scale(deflection, _system->get_deflector().distance_to_previous_plane());
    
    util::add(out, deflection);
    
    util::vector2 obs_position_term = _system->get_observer().get_observer_position();
    util::scale(obs_position_term, this->_dLS_over_dL);    
    util::sub(out, obs_position_term);
    
    return out;
  }
  
  util::scalar _dL_over_dLS;
  util::scalar _dLS_over_dL;
  
  const SystemType* _system;
  util::scalar _accuracy;
};

template<class SystemType>
class inversion_table_image_finder : public image_finder<SystemType>
{
public:
  typedef numeric::function_inverter<util::scalar, 2, util::scalar, 2> inverter_type;
  typedef std::shared_ptr<inverter_type> inverter_ptr_type;
  
  static const std::size_t benchmark_size = 10000;
  
  inversion_table_image_finder(const SystemType* sys,
                               const util::vector2& sampling_center,
                               const util::scalar sampling_radius,
                               const util::vector2& physical_source_plane_size,
                               const util::vector2& screen_position,
                               const std::array<std::size_t, 2>& num_pixels,
                               std::size_t num_samples_per_dim,
                               util::scalar accuracy)
  : image_finder<SystemType>(sys, accuracy),
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
    
    this->_inverter = inverter_ptr_type(new inverter_type(pixel_min_coordinates,
                                                          pixel_max_coordinates,
                                                          sampling_start_vector,
                                                          sampling_end_vector));
    
    
    auto function_evaluator = [this](const inverter_type::domain_vector& domain_vec) -> inverter_type::codomain_vector
    {
      return this->ray_function(domain_vec);
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
                          std::vector<util::vector2>& out) const
  {
    out.clear();
    std::vector<util::vector2> start_positions = _inverter.inverse(source_plane_pos);
    
    // Define the function of which we want determine the root using Newton's method
    auto ray_equation = [&](const util::vector2& lens_plane_pos) -> util::vector2
    {
      // left hand side of the ray equation
      util::vector2 lhs = this->ray_function(lens_plane_pos);
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
  numeric::function_inverter<util::scalar, 2, util::scalar, 2> _inverter;
};

//TODO
template<class SystemType>
class root_tracing_image_finder : public image_finder<SystemType>
{
public:
  root_tracing_image_finder(const SystemType* sys)
  : image_finder(sys)
  {
    
  }
  
  virtual void get_images(const util::vector2& source_plane_pos,
                          std::vector<util::vector2>& out) const
  {
    
  }
};


template<class SystemType>
class complex_polynomial_image_finder : public image_finder<SystemType>
{
  struct root
  {
    util::vector2 position;
    unsigned multiplicity;
  };
  
public:
  typedef std::complex<util::scalar> complex_type;
  
  complex_polynomial_image_finder(const SystemType* sys, util::scalar accuracy)
  : image_finder(sys, accuracy),
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

}

#endif	/* IMAGE_FINDER_HPP */

