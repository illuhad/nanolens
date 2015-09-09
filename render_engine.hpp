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

#ifndef RENDERENGINE_HPP
#define	RENDERENGINE_HPP

#include <vector>
#include <memory>
#include <fstream>
#include <random>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/math/constants/constants.hpp>

#include "util.hpp"
#include "lens_plane.hpp"
#include "observer_plane.hpp"
#include "system.hpp"
#include "image_finder.hpp"
#include "pixel_processor.hpp"
#include "status.hpp"
#include "fits.hpp"
#include "screen.hpp"

namespace nanolens{
namespace render_engines{

template<class System_type, class Method_type, class Screen_type>
class renderer
{
public:
  typedef typename Method_type::settings method_settings_type;
  static constexpr int master_rank = 0;
  
  renderer(const boost::mpi::communicator& comm,
           const method_settings_type& config,
           const Screen_type& s)
  : _config(config), _screen(std::move(s)), _comm(comm), _my_rank(comm.rank())
  {
    _method = std::make_shared<Method_type>(_comm, _config, &_screen);
  }
  
  virtual ~renderer(){}
  
  virtual void run(System_type& sys, status_handler_type handler) = 0;
  
  virtual void save_as_fits(const std::string& filename)
  {
    if(_my_rank == master_rank)
    {
      util::fits<typename Screen_type::pixel_type> fits_file(filename);
      fits_file.save(_screen.get_pixels());
    }
  }
  
  void save_as_raw(const std::string& filename)
  {
    if(_my_rank == master_rank)
    {
      std::ofstream file(filename.c_str(), 
                           std::ofstream::out|std::ofstream::binary|std::ofstream::trunc);

      if(_screen.get_pixels().size() != 0)
      {

        if(file.is_open())
        {
          file.write(reinterpret_cast<const char*>(_screen.get_pixels().begin()), 
                     sizeof(typename Screen_type::pixel_type) * _screen.get_pixels().size());
        }
      }
    }
  }
protected:
  method_settings_type _config;
  Screen_type _screen;
  boost::mpi::communicator _comm;
  std::shared_ptr<Method_type> _method;
  int _my_rank;
};

template<class System_type, class Method_type, class Pixel_type>
class standard_renderer : public renderer<System_type, Method_type, screen<Pixel_type>>
{
public:
  using typename renderer<System_type, Method_type, screen<Pixel_type>>::method_settings_type;
  
  /// Construct object.
  /// @param comm The mpi communicator
  /// @param method_config The settings for the chosen method
  /// @param num_pixels The number of pixels in x and y direction
  /// @param physical_size The physical size in Einstein radii of the system
  /// @param screen_position The x and y coordinates of the center of the screen
  standard_renderer(const boost::mpi::communicator& comm,
                    const method_settings_type& method_config,
              const std::array<std::size_t, 2>& num_pixels,
              const util::vector2& physical_size,
              const util::vector2& screen_position)
  : renderer<System_type, Method_type, screen<Pixel_type>>(comm, 
                                       method_config, 
                                       std::move(screen<Pixel_type>(num_pixels, physical_size, screen_position)))
  {}

  void run(System_type& sys, status_handler_type handler) 
  {
    this->_method->run(sys, handler);
  }
  
  virtual ~standard_renderer(){}
  
private:
  util::scalar _accuracy;
};

template<class System_type, class Method_type, class Screen_type, class Pixel_type>
class system_modifying_renderer : public renderer<System_type, Method_type, Screen_type>
{
public:
  /// Functor that will be called upon every step to modify the system. If it returns
  /// false, computation terminates.
  typedef std::function<bool (System_type&)> system_modification_function;

  using typename renderer<System_type, Method_type, screen<Pixel_type>>::method_settings_type;
  
  
  /// Construct object.
  /// @param comm The mpi communicator
  /// @param method_config The settings for the chosen method
  /// @param num_pixels The number of pixels in x and y direction
  /// @param physical_size The physical size in Einstein radii of the system
  /// @param screen_position The x and y coordinates of the center of the screen
  /// @param mod_fun A function that will be applied on the system after each
  /// run of the method. It must have the signature bool (System_type&). If it returns
  /// false, the computation will finish.
  /// @param shared_screen Whether the individual runs of the method will share one
  /// pixel screen. If false, all methods will operate on different screens
  /// and the result will be a stack of screens.
  system_modifying_renderer(const boost::mpi::communicator& comm,
                    const method_settings_type& method_config,
                    const std::array<std::size_t, 2>& num_pixels,
                    const util::vector2& physical_size,
                    const util::vector2& screen_position,
                    system_modification_function mod_fun,
                    bool shared_screen)
  : renderer<System_type, Method_type, screen<Pixel_type>>(comm, 
                                       method_config, 
                                       std::move(screen<Pixel_type>(num_pixels, physical_size, screen_position))),
    _mod_fun(mod_fun),
    _num_pixels(num_pixels),
    _physical_screen_size(physical_size)
  {
  }
  
  virtual ~system_modifying_renderer(){}
  
  void run(System_type& sys, status_handler_type handler)
  {
    do
    {
      if(!_shared_screen)
        this->_screen = screen<Pixel_type>(_num_pixels, _physical_screen_size, _screen_position);
      
      this->_method->run(sys, handler);
      
      if(!_shared_screen)
        _results.push_back(this->_screen);
    }
    while(_mod_fun(sys));
    
    // Make sure the screen stack contains the single screen in case we have
    // worked with a shared screen
    if(_shared_screen)
      _results.push_back(this->_screen);
  }
  
  const std::vector<screen<Pixel_type>>& get_screen_stack() const
  {
    return _results;
  }
  
private:
  system_modification_function _mod_fun;
  bool _shared_screen;
  std::vector<screen<Pixel_type>> _results;
  
  std::array<std::size_t, 2> _num_pixels;
  util::vector2 _physical_screen_size;
  util::vector2 _screen_position;
};

}
}

#endif	/* SOURCE_PLANE_HPP */

