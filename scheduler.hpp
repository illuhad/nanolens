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
  size_t _num_jobs;
  
  size_t _own_beg;
  size_t _own_end;
  
  std::vector<double> _performance_statistic;
};

}

#endif	/* SCHEDULER_HPP */

