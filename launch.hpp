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

#ifndef LAUNCH_HPP
#define	LAUNCH_HPP

#include "input.hpp"
#include "method_irs.hpp"

namespace nanolens {

template<class System_type>
class standard_launcher
{
public:
  
  template<class Method_type>
  void run_method(const boost::mpi::communicator& comm,
                           const configuration& config,
                           System_type& sys) const
  {
    typename Method_type::settings method_settings(config);

    nanolens::render_engines::standard_renderer
    <
      System_type, 
      Method_type,
      util::scalar
    > engine(comm,
             method_settings,
             config.get_resolution(),
             config.get_physical_screen_size(),
             config.get_screen_pos());

    nanolens::standard_terminal_output status_handler(comm, 0);

    engine.run(sys, status_handler);
    
    if(!config.get_fits_output().empty())
      engine.save_as_fits(config.get_fits_output());
  
    if(!config.get_raw_output().empty())
      engine.save_as_raw(config.get_raw_output());

  }


  void run_configured_method(const boost::mpi::communicator& comm,
                               const configuration& config,
                               System_type& sys) const
  {
    switch(config.get_method_type())
    {
    case nanolens::configuration::INVERSE_RAY_SHOOTING:
      run_method
      <
        magnification_pattern_generation::inverse_ray_shooting<System_type>
      >(comm, config, sys);
      break;
    }
  }
  
};


}

#endif	/* LAUNCH_HPP */

