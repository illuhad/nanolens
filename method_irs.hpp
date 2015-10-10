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

/// @file This file contains the implementation of the widely used inverse
/// ray shooting method. I.e. rays are shot backwards from the observer
/// to the source and counted in the pixels of a screen in the source plane.
/// The number of rays in one pixel is then proportional to the pixel's magnification.

#ifndef INVERSE_RAY_SHOOTING_HPP
#define	INVERSE_RAY_SHOOTING_HPP

#include "util.hpp"
#include "method.hpp"
#include "input.hpp"
#include "scheduler.hpp"
#include "grid.hpp"

namespace nanolens{
namespace magnification_pattern_generation{

template<class System_type>
class inverse_ray_shooting_region
{
public:
  static constexpr util::scalar overshooting_area_size = 6.0;
  
  inverse_ray_shooting_region(const System_type& sys,
                              const screen<util::scalar>& s,
                              util::scalar n_rays_per_pixel)
  {
    _center = s.get_properties().get_screen_position();
    
    util::vector2 screen_size = s.get_properties().get_physical_size();

    util::vector2 screen_min_corner = s.get_properties().get_corner_of_min_extent();
    
    util::matrix_nxn<util::vector2, 2> corner_coordinates;
    for(std::size_t i = 0; i < 2; ++i)
      for(std::size_t j = 0; j < 2; ++j)
      {
        corner_coordinates[i][j] = screen_min_corner;
        util::vector2 offset = screen_size;
        offset[0] *= i;
        offset[1] *= j;
        
        util::add(corner_coordinates[i][j], offset);
      }
    
    util::matrix_nxn<util::vector2, 2> shooting_region_corners;
    sys.get_deflector().estimate_mapped_region_coordinates(corner_coordinates, 
                                                           shooting_region_corners);
    
    _shooting_region_size = 
      {std::abs(shooting_region_corners[1][0][0] - shooting_region_corners[0][0][0]) + overshooting_area_size,
       std::abs(shooting_region_corners[0][1][1] - shooting_region_corners[0][0][1]) + overshooting_area_size};
    
    
    util::vector2 half_shooting_region_size = _shooting_region_size;
    util::scale(half_shooting_region_size, 0.5);
    
    _min_corner = _center;
    util::sub(_min_corner, half_shooting_region_size);
    
    util::vector2 pixel_size = {_shooting_region_size[0] / s.get_properties().get_num_pixels()[0],
                                _shooting_region_size[1] / s.get_properties().get_num_pixels()[1]};

    _ray_density 
      =  n_rays_per_pixel / (pixel_size[0] * pixel_size[1]);
    
    util::scalar rays_per_px_y = std::sqrt(n_rays_per_pixel * pixel_size[1] / pixel_size[0]);
    util::scalar rays_per_px_x = pixel_size[0] / pixel_size[1] * rays_per_px_y;
    
    _ray_density_per_dim = {rays_per_px_x / pixel_size[0], rays_per_px_y / pixel_size[1]};
    
    _num_rays = {static_cast<std::size_t>(_ray_density_per_dim[0] * _shooting_region_size[0]),
                 static_cast<std::size_t>(_ray_density_per_dim[1] * _shooting_region_size[1])};
    
    _ray_distances = {_shooting_region_size[0] / _num_rays[0],
                      _shooting_region_size[1] / _num_rays[1]};
    
  }
  
  const util::vector2& get_ray_separation() const
  { return _ray_distances; }
  
  const std::array<std::size_t, 2>& get_num_rays() const
  { return _num_rays; }
  
  util::scalar get_ray_density() const
  { return _ray_density; }
  
  util::scalar get_ray_density_x() const
  { return _ray_density_per_dim[0]; }
  
  util::scalar get_ray_density_y() const
  { return _ray_density_per_dim[1]; }
  
  const util::vector2& get_region_center() const
  { return _center; }
  
  const util::vector2& get_region_size() const
  { return _shooting_region_size; }
  
  const util::vector2& get_region_min_corner() const
  { return _min_corner; }
private:
  
  util::vector2 _ray_density_per_dim;
  util::vector2 _ray_distances;
  std::array<std::size_t, 2> _num_rays;
  util::scalar _ray_density;
  util::vector2 _center;
  util::vector2 _shooting_region_size;
  util::vector2 _min_corner;
};

class inverse_ray_shooting_settings
{
  util::scalar _n_rays_per_px;
public:

  explicit inverse_ray_shooting_settings(const configuration& config)
  {
    _n_rays_per_px = config.get_method_property("num_rays_per_pixel", 1.0);
  }

  explicit inverse_ray_shooting_settings(util::scalar n_rays_per_px)
  : _n_rays_per_px(n_rays_per_px) {}

  util::scalar get_n_rays_per_px() const
  { return _n_rays_per_px; }
  
};

template<class System_type>
class inverse_ray_shooting : public method<System_type, 
                                           inverse_ray_shooting_settings, 
                                           screen<util::scalar>
                                          >
{  
public:
  typedef inverse_ray_shooting_settings settings;
  

  inverse_ray_shooting(const boost::mpi::communicator& comm,
                       const settings& config,
                       screen<util::scalar>* s)
  : method<System_type, 
          inverse_ray_shooting<System_type>::settings, 
          screen<util::scalar>>(comm, config, s)
  {}
  
  virtual void run(const System_type& sys, status_handler_type handler)
  {
    // Determine shooting region
    inverse_ray_shooting_region<System_type> shooting_region(sys, *(this->_screen), this->_config.get_n_rays_per_px());
    
    if(this->_comm.rank() == 0)
      handler(status_info("", nullptr, std::string("Shooting region size: ") 
             + std::to_string(shooting_region.get_region_size()[0])
             + "x"
             + std::to_string(shooting_region.get_region_size()[1])));
    
    std::size_t n_rays_x = shooting_region.get_num_rays()[0];
    std::size_t n_rays_y = shooting_region.get_num_rays()[1];
    
    util::vector2 step_sizes = shooting_region.get_ray_separation();
    
    util::vector2 min_corner = shooting_region.get_region_min_corner();
    util::vector2 ray_position = min_corner;
    
    util::scalar gamma = sys.get_deflector().get_shear();
    util::scalar sigma = sys.get_deflector().get_mean_surface_density();
    
    util::scalar magnification_per_ray 
      = 1.0
      / (this->_config.get_n_rays_per_px() * std::abs((util::square(1.0 - sigma) - util::square(gamma))));
    
    scheduler schedule(this->_comm);
    
    std::fill(this->_screen->get_pixels().begin(), this->_screen->get_pixels().end(), 0.0);
    
    handler(status_info("Scheduling inverse ray shooting"));
    
    util::grid_translator<util::scalar, 2, false> pixel_grid_translator(
      this->_screen->get_properties().get_corner_of_min_extent(),
      this->_screen->get_properties().get_corner_of_max_extent(),
      this->_screen->get_properties().get_num_pixels());
    
    
    // the dummy counter prevents the compiler from optimizing stuff away
    std::size_t dummy_counter = 0;
    schedule.autosized_run(n_rays_x, [&](std::size_t test_ray_index, std::size_t num_tests)
    {
      util::vector2 ray_position = shooting_region.get_region_center();
      ray_position[1] = shooting_region.get_region_min_corner()[1];
      
      util::scalar step_size = shooting_region.get_region_size()[1] / static_cast<util::scalar>(num_tests);
      
      ray_position[1] += step_size * test_ray_index;
      
      util::vector2 impact_position = sys.lensing_transformation(ray_position);
      if(pixel_grid_translator.contains_point(impact_position))
        ++dummy_counter;
      
    }, 2.0);
      
    
    handler(status_info("Scheduling complete", &schedule, ""));
    
    // Use double precision for the position of the ray to avoid rounding errors
    // for large shooting regions
    std::array<double, 2> double_min_corner = {static_cast<double>(min_corner[0]),
                                               static_cast<double>(min_corner[1])};
    std::array<double, 2> double_stepwidth = {static_cast<double>(step_sizes[0]),
                                              static_cast<double>(step_sizes[1])};
    
    std::size_t num_jobs_completed = 0;
    for(std::size_t x_idx = schedule.get_ownership_range_begin();
            x_idx <= schedule.get_ownership_range_end(); ++x_idx)
    {
      ray_position[0] = double_min_corner[0]
                      + x_idx * double_stepwidth[0];
      
      util::scalar progress = static_cast<util::scalar>(num_jobs_completed) /
        static_cast<util::scalar>(schedule.get_num_assigned_jobs());
      
      handler(status_info("Evaluating ray function", nullptr, "", progress));
      
      for(std::size_t y_idx = 0; y_idx < n_rays_y; ++y_idx)
      {
        ray_position[1] = double_min_corner[1]
                        + y_idx * double_stepwidth[1];
        
        util::vector2 cast_ray_position = {static_cast<util::scalar>(ray_position[0]),
                                           static_cast<util::scalar>(ray_position[1])};
        
        util::vector2 impact_position = sys.lensing_transformation(cast_ray_position);

        std::array<std::size_t, 2> pixel_idx = pixel_grid_translator(impact_position);

        if(this->_screen->get_pixels().is_within_bounds(pixel_idx.data()))
          this->_screen->get_pixels()[pixel_idx.data()] += magnification_per_ray;

      }
      
      ++num_jobs_completed;
    }
    
    handler(status_info("", nullptr, "Finished computation"));
    
    handler(status_info("Waiting for processes...\n"));
    
    auto add =  [](util::scalar& a, util::scalar b){ a += b;};
    this->_screen->get_pixels().reduce_parallel_array(this->_comm, 0, add);
  }
  
private:
};

}
}


#endif	/* INVERSE_RAY_SHOOTING_HPP */

