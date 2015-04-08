/* 
 * File:   scheduler.hpp
 * Author: aksel
 *
 * Created on 21. März 2015, 05:07
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
    assert(process_rank < _all_begins.size());
    
    return _all_begins[process_rank] <= job_entity_index && 
          job_entity_index <= _all_ends[process_rank];
  }
  
  int get_assigned_process_rank(std::size_t job_entity_index) const
  {
    for(int i = 0; i < _all_begins.size(); ++i)
      if(is_scheduled_to_process(job_entity_index, i))
        return i;
    return -1;
  }
  
  std::size_t get_num_assigned_jobs() const
  {
    return _own_end - _own_beg + 1;
  }
  
  std::size_t get_num_assigned_jobs(int rank) const
  {
    assert(rank < _all_begins.size());
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

