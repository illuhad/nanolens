/* 
 * File:   status.hpp
 * Author: aksel
 *
 * Created on 8. April 2015, 04:00
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

