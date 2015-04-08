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

#ifndef STATUS_HPP
#define	STATUS_HPP

#include <string>

#include "scheduler.hpp"

namespace nanolens{

class status_info
{
public:
  status_info(const std::string& current_task,
         const scheduler* current_schedule,
         const std::string& notification,
         double progress)
  : _current_task(current_task),
    _current_notification(notification),
    _current_schedule(current_schedule),
    _progress(progress),
    _with_progress(true){}
  
 status_info(const std::string& current_task,
         const scheduler* current_schedule,
         const std::string& notification)
  : _current_task(current_task),
    _current_notification(notification),
    _current_schedule(current_schedule),
    _progress(0.0),
    _with_progress(false){}
  
  status_info(const std::string& current_task)
  : _current_task(current_task),
    _current_notification(""),
    _current_schedule(nullptr),
    _progress(0.0),
    _with_progress(false){}
  
  const std::string& get_task() const
  { return _current_task; }
  
  const std::string& get_current_notification() const
  { return _current_notification; }
  
  const scheduler* get_schedule() const
  { return _current_schedule; }
  
  double get_progress() const
  {
    return _progress;
  }
  
  bool progress_known() const
  {
    return _with_progress;
  }
  
private:
  bool _with_progress;
  std::string _current_task;
  std::string _current_notification;
  const scheduler* _current_schedule;
  double _progress;
};

}

#endif	/* STATUS_HPP */

