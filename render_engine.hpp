/* 
 * File:   source_plane.hpp
 * Author: aksel
 *
 * Created on 8. Dezember 2014, 00:40
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

namespace nanolens{

  class render_engine
  {
  public:
    typedef std::size_t count_type;
    
    render_engine(const util::vector2& num_pixels,
                const util::vector2& physical_size,
                const util::vector2& screen_position)
    : _pixels(num_pixels[0], num_pixels[1]),
      _size(physical_size),
      _position(screen_position),
    {}
    
    template<class Function>
    void run(system& sys, const boost::mpi::communicator& comm, Function status_handler)
    {
      this->_num_rays_total = 0;
      std::fill(_pixels.begin(), _pixels.end(), 0.0);
      
      
      util::scalar einstein_radius = sys.get_einstein_radius();
      
      
      util::vector2 scaled_size = _size;
      util::scale(scaled_size, einstein_radius * sys.get_distance_from_source_to_observer());
      
      const util::scalar pi = boost::math::constants::pi<util::scalar>();
      
      util::vector2 position = {0.0, 0.0};
      
      std::size_t npixels_x = _pixels.get_extent_of_dimension(0);
      std::size_t npixels_y = _pixels.get_extent_of_dimension(1);
      //util::scalar pixel_size = (scaled_size[0] + scaled_size[1]) / (npixels_x + npixels_y);
      
      util::scalar distance_to_observer = sys.get_deflector().distance_to_previous_plane()
              + sys.get_observer().distance_to_previous_plane();
      
      util::scalar angular_radius = std::max(3.5 * sys.get_deflector().get_radius_estimate() / 
                sys.get_deflector().distance_to_previous_plane(),
              0.5 * (scaled_size[0] + scaled_size[1]) /
                sys.get_deflector().distance_to_previous_plane());
      
      
      util::scalar z_offset = 0.0;
      
      
      unsigned num_rays_per_dim = static_cast<unsigned>(std::sqrt(_num_rays));

      std::default_random_engine rng_engine(_rd());
      std::uniform_real_distribution<util::scalar> random_angle(0.0, 2 * pi);
      
      int rank = comm.rank();
      int num_procs = comm.size();   
      
      for(std::size_t x = 0; x < npixels_x; ++x)
      {
        
        if(x % num_procs == static_cast<unsigned>(rank))
        {
          status_handler(static_cast<double>(x)/npixels_x);

          util::scalar prev_magnification = 1.0;
          for(std::size_t y = 0; y < npixels_y; ++y)
          {
            std::size_t pixel_idx[] = {x, y};

            position[0] = _position[0] - 0.5 * scaled_size[0]  
                    + (static_cast<util::scalar>(x) + 0.5) / npixels_x * scaled_size[0];

            position[1] = _position[1] - 0.5 * scaled_size[1]  
                    + (static_cast<util::scalar>(y) + 0.5) / npixels_y * scaled_size[1];


            std::size_t num_rays = std::min(prev_magnification / 3.0, 
                                            static_cast<util::scalar>(3)) * num_rays_per_dim;
            if(num_rays < num_rays_per_dim)
              num_rays = num_rays_per_dim;
            
            util::vector2 central_angle = sys.get_observer().get_observer_position();

            util::sub(central_angle, position);
            util::scale(central_angle, 1.0 / distance_to_observer);
            //central_angle={0.0,0.0};
            ray_dome<AdaptiveMeshPolicy> dome(position, 
                                              central_angle, 
                                              z_offset, 
                                              angular_radius, 
                                              num_rays,
                                              _refinement_level);
            bool dump = false;
            if(x == 0 && y == 0)
              dump = false;
            
            dome.traverse(sys, dump);
            _num_rays_total += dome.get_num_rays();
            ++_num_domes;

            util::scalar mag = dome.get_magnification();
            _pixels[pixel_idx] = mag;
            
            prev_magnification = 1.0;
          }
        }

      }
      
      // collect data
      std::vector<util::scalar> data;
      data.reserve(1.2 * npixels_x / num_procs * npixels_y);
      
      for(std::size_t x = 0; x < npixels_x; ++x)
      {
        for(std::size_t y = 0; y < npixels_y; ++y)
        {
          std::size_t idx [] = {x, y};
          if(x % num_procs == static_cast<unsigned>(rank))
            data.push_back(_pixels[idx]);
        }
      }
      
      std::vector<std::vector<util::scalar> > received_data;
      if(rank == 0)
        received_data.reserve(num_procs);
      
      boost::mpi::gather(comm, data, received_data, 0);
      
      if(rank == 0)
      {
        std::vector<std::size_t> read_positions(num_procs, 0);
        for(std::size_t x = 0; x < npixels_x; ++x)
        {
          for(std::size_t y = 0; y < npixels_y; ++y)
          {
            int process_id = x % num_procs;
            std::size_t idx [] = {x, y};
            _pixels[idx] = received_data[process_id][read_positions[process_id]];
            ++read_positions[process_id];
          }
        }
      }

    }
    
    count_type get_num_traced_rays() const
    {
      return _num_rays_total;
    }
    
    count_type get_num_domes() const
    {
      return _num_domes;
    }
    
    void save_pixels(const std::string& filename)
    {
      std::ofstream file(filename.c_str(), 
                         std::ofstream::out|std::ofstream::binary|std::ofstream::trunc);
      
      if(file.is_open())
      {
        file.write(reinterpret_cast<char*>(_pixels.begin()), 
                   sizeof(util::scalar) * _pixels.size());
      }
    }
    
  private:

    util::multi_array<util::scalar> _pixels;
    util::vector2 _size;
    util::vector2 _position;
  };
}

#endif	/* SOURCE_PLANE_HPP */

