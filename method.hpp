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

#ifndef METHOD_HPP
#define	METHOD_HPP

#include "screen.hpp"
#include "status.hpp"
#include "input.hpp"

namespace nanolens {

template<class System_type, class Settings_type, class Screen_type>
class method
{
public:
  
  method(const boost::mpi::communicator& comm,
         const Settings_type config,
         Screen_type* s,
         int master_rank = 0)
  : _screen(s), _comm(comm), _config(config)
  {
    assert(s != nullptr);
    
    boost::mpi::broadcast(_comm, _config, master_rank);
  }
  
  virtual void run(const System_type& sys, status_handler_type handler) = 0;
  
  virtual ~method()
  {}
protected:
  Screen_type* _screen;
  boost::mpi::communicator _comm;
  Settings_type _config;
};




}

#endif	/* METHOD_HPP */

