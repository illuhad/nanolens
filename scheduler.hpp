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

#ifndef SCHEDULER_HPP
#define	SCHEDULER_HPP

#include <vector>
#include "boost/mpi.hpp"
#include "timer.hpp"

namespace nanolens{

class scheduler
{
public:
  scheduler(const boost::mpi::communicator& comm)
  : _comm(comm) {}
  
  template<class Function>
  void run(size_t num_job_entities, Function test_job)
  { 
    util::timer t;
    
    t.start();
    test_job();
    double time = t.stop();
    double speed = 1.0 / time;
    
    boost::mpi::all_gather(_comm, speed, this->_performance_statistic);
    
    // normalize speeds
    double sum = std::accumulate(_performance_statistic.begin(), _performance_statistic.end(), 0.0);
    
    for(double& element : _performance_statistic)
      element /= sum; 
    
    _num_jobs = num_job_entities;
    
    calculate_ownership_range();
    
    boost::mpi::all_gather(_comm, _own_beg, this->_all_begins);
    boost::mpi::all_gather(_comm, _own_end, this->_all_ends);
    
    assert(_all_begins.size() == _all_ends.size());
  }
  
  const std::vector<double>& get_relative_performances() const
  {
    return _performance_statistic;
  }
  
  bool is_scheduled_to_this_process(size_t job_entity_index) const
  {
    assert(job_entity_index < _num_jobs);
    
    return _own_beg <= job_entity_index && job_entity_index <= _own_end;
  }
  
  bool is_scheduled_to_process(std::size_t job_entity_index, int process_rank) const
  {
    assert(static_cast<std::size_t>(process_rank) < _all_begins.size());
    
    return _all_begins[process_rank] <= job_entity_index && 
          job_entity_index <= _all_ends[process_rank];
  }
  
  int get_assigned_process_rank(std::size_t job_entity_index) const
  {
    for(std::size_t i = 0; i < _all_begins.size(); ++i)
      if(is_scheduled_to_process(job_entity_index, static_cast<int>(i)))
        return static_cast<int>(i);
    return -1;
  }
  
  std::size_t get_num_assigned_jobs() const
  {
    return _own_end - _own_beg + 1;
  }
  
  std::size_t get_num_assigned_jobs(int rank) const
  {
    assert(static_cast<std::size_t>(rank) < _all_begins.size());
    return _all_ends[rank] - _all_begins[rank] + 1;
  }
  
  double get_progress_at_job(std::size_t job_entity_index) const
  {
    return static_cast<double>(job_entity_index - _own_beg) /
      static_cast<double>(get_num_assigned_jobs());
  }
private:
  void calculate_ownership_range()
  {
    int rank = _comm.rank();
    
    size_t current_beg = 0;
    for(int i = 0; i < rank; ++i)
      current_beg += static_cast<size_t>(_performance_statistic[i] * static_cast<double>(_num_jobs)) + 1;
    
    _own_end = current_beg + static_cast<size_t>(_performance_statistic[rank] * static_cast<double>(_num_jobs));
    if(_own_end >= _num_jobs || rank == (_comm.size() - 1))
      _own_end = _num_jobs - 1;

    _own_beg = current_beg;
  }
  
  boost::mpi::communicator _comm;
  std::size_t _num_jobs;
  
  std::size_t _own_beg;
  std::size_t _own_end;
  
  std::vector<std::size_t> _all_begins;
  std::vector<std::size_t> _all_ends;
  
  std::vector<double> _performance_statistic;
};

}

#endif	/* SCHEDULER_HPP */

