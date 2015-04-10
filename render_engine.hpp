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

#ifndef RAYTRACER_HPP
#define	RAYTRACER_HPP

#include <vector>
#include <memory>
#include <fstream>
#include <random>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/math/constants/constants.hpp>
#include "util.hpp"
#include "ray.hpp"
#include "lens_plane.hpp"
#include "observer_plane.hpp"
#include "system.hpp"
#include "image_finder.hpp"
#include "pixel_processor.hpp"
#include "status.hpp"
#include "screen.hpp"

namespace nanolens{

class render_engine
{
public:
  typedef std::size_t count_type;

  static constexpr int master_rank = 0;
  
  render_engine(const std::array<std::size_t, 2>& num_pixels,
              const util::vector2& physical_size,
              const util::vector2& screen_position,
              util::scalar numerical_accuracy = 1.e-7)
  : _num_pixels(num_pixels),
    _size(physical_size),
    _position(screen_position),
    _accuracy(numerical_accuracy),
    _scaled_screen_size(physical_size)
  {}

  template<class Function>
  void run(system& sys, const boost::mpi::communicator& comm, Function status_handler)
  {

    init_screen(sys);

    std::size_t npixels_x = _num_pixels[0];
    std::size_t npixels_y = _num_pixels[1];

    // Find images
    status_handler(status_info("Initializing image finding algorithm"));
    
/*    
    std::shared_ptr<image_finder<system>> backend_img_finder(
          new image_finders::newton_crown<system, 256>(&sys, _accuracy, status_handler));
    
    std::shared_ptr<image_finder<system>> img_finder(
          new image_finders::root_tracing<system>(&sys, 
                                                  _accuracy,
                                                   status_handler,
                                                   backend_img_finder.get(),
                                                   *_screen));
    
 */
    std::shared_ptr<image_finder<system>> img_finder(
      new image_finders::inverse_ray_shooting<system>(&sys, status_handler, *_screen, {0.0,0.0}, {8.0, 8.0}, {10000, 10000}, comm));
      
    pixel_processor<magnification::by_compact_image_count> pixel_evaluator(2.0 * _accuracy);
    
    scheduler pixel_schedule(comm);
    
    status_handler(status_info("Scheduling pixel processing"));
    // Calculate 50 pixel to get performance estimates for the scheduler
    std::size_t benchmark_size = 1000;
    auto scheduler_test_function = [&]()
    {
      std::size_t px_x = std::max(npixels_x / 2, static_cast<std::size_t>(1));
      
      for(std::size_t px_y = 0; px_y < benchmark_size; ++px_y)
      {
        util::vector2 pixel_position = _screen->get_pixel_coordinates({px_x, px_y});

        pixel_evaluator.get_pixel_magnification(pixel_position,
                                                 sys,
                                                 *img_finder);
      }
    };
    
    pixel_schedule.run(npixels_x, scheduler_test_function);
    
    status_handler(status_info("Pixel processing complete", &pixel_schedule, ""));
    
    std::vector<util::scalar> send_data;
    // Reserve a generous amount of space so that the vector won't have
    // to grow often
    send_data.reserve(0.7 * _num_pixels[0] * _num_pixels[1]);
    
    // Compute pixel values in parallel
    for(std::size_t px_x = 0; px_x < npixels_x; ++px_x)
    {
      double progress = pixel_schedule.get_progress_at_job(px_x);
      
      status_handler(status_info("Processing pixels", nullptr, "", progress));
      if(pixel_schedule.is_scheduled_to_this_process(px_x))
      {
        for(std::size_t px_y = 0; px_y < npixels_y; ++px_y)
        { 
          util::vector2 pixel_position = _screen->get_pixel_coordinates({px_x, px_y});

          util::scalar magnification = pixel_evaluator.get_pixel_magnification(pixel_position,
                                                                               sys,
                                                                               *img_finder);

          send_data.push_back(magnification);

        }
      }
    }
    
    std::vector<std::vector<util::scalar>> recv_data;
    boost::mpi::gather(comm, send_data, recv_data, master_rank);
    
    if(comm.rank() == master_rank)
    {
      _pixels = util::multi_array<util::scalar>(_num_pixels[0], _num_pixels[1]);
      
      std::fill(_pixels.begin(), _pixels.end(), 0.0);
      
      // Sort received pixel values from processes and compose image
      std::vector<std::size_t> rank_pixel_counters(comm.size(), 0);

      for(std::size_t px_x = 0; px_x < npixels_x; ++px_x)
      {
        int assigned_process_rank = pixel_schedule.get_assigned_process_rank(px_x);
        assert(assigned_process_rank >= 0);

        for(std::size_t px_y = 0; px_y < npixels_y; ++px_y)
        {
          std::size_t pixel_index [] = {px_x, px_y};

          std::size_t current_pixel_count_for_rank = rank_pixel_counters[assigned_process_rank];

          _pixels[pixel_index] = recv_data[assigned_process_rank][current_pixel_count_for_rank];

            ++rank_pixel_counters[assigned_process_rank];
        }
      }
    }
  }


  /// works only on the master process
  void save_pixels(const std::string& filename)
  {
    
    std::ofstream file(filename.c_str(), 
                         std::ofstream::out|std::ofstream::binary|std::ofstream::trunc);
    
    if(_pixels.size() != 0)
    {

      if(file.is_open())
      {
        file.write(reinterpret_cast<char*>(_pixels.begin()), 
                   sizeof(util::scalar) * _pixels.size());
      }
    }
  }

private:
  
  void init_screen(const system& sys)
  {
    util::scalar einstein_radius = sys.get_einstein_radius();

    _scaled_screen_size = _size;
    util::scale(_scaled_screen_size, einstein_radius * sys.get_distance_from_source_to_observer());
    
    _screen = std::shared_ptr<screen_descriptor>(new screen_descriptor(_num_pixels, 
                                                                       _scaled_screen_size, 
                                                                       _position));
    
  }
  

  std::shared_ptr<screen_descriptor> _screen;
  std::array<std::size_t, 2> _num_pixels;
  util::scalar _accuracy;
  util::multi_array<util::scalar> _pixels;
  util::vector2 _size;
  util::vector2 _pixel_size;
  util::vector2 _scaled_screen_size;
  util::vector2 _position;
};

}

#endif	/* SOURCE_PLANE_HPP */

