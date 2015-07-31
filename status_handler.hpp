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

#ifndef STATUS_HANDLER_HPP
#define	STATUS_HANDLER_HPP

#include "status.hpp"
#include "util.hpp"
#include <boost/mpi.hpp>

namespace nanolens{

class status_handler
{
public:
  virtual ~status_handler(){}
  virtual void operator() (const status_info& status) = 0;
};


class standard_terminal_output : public status_handler
{
public:
  standard_terminal_output(const boost::mpi::communicator& comm, int master_rank)
  : _comm(comm), _master_cout(std::cout, master_rank), _max_line_length(0)
  {
    std::cout.width(6);
    std::cout.precision(5);
  }
  
  virtual void operator() (const status_info& status)
  {
    if(status.get_current_notification() != "")
    {
      std::cout << "\nMessage from process " << _comm.rank() << ": " << status.get_current_notification() << std::endl;
    }
    
    if(_comm.rank() == 0)
    {
      if(status.get_schedule() != nullptr)
      {
        _master_cout << "\n\nPerformance metric of current task schedule:\n";

        std::vector<double> performances = status.get_schedule()->get_relative_performances();

        for(std::size_t i = 0; i < performances.size(); ++i)
          _master_cout << "Process " << i << ": " << performances[i]/performances[0] << "x\n";
        
        _master_cout << "(relative to process 0)\n\n";
      }

      std::string clear_string = "\r";
      for(std::size_t i = 0; i < _max_line_length; ++i)
        clear_string += " ";

      _master_cout << clear_string;
      std::stringstream sstr;
      sstr << "\r" << "Status: " << status.get_task();
      if(status.progress_known())
        sstr << " | " << status.get_progress() * 100 << "%";
        
      std::string msg =  sstr.str();

      if(msg.length() > _max_line_length)
        _max_line_length = msg.length();

      _master_cout << "\r" << msg;
      std::cout.flush();
    }
  }
  
private:
  const boost::mpi::communicator _comm;
  util::master_ostream _master_cout;
  unsigned _max_line_length;
};

}


#endif	/* STATUS_HANDLER_HPP */

