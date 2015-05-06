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

#ifndef NUMERIC_HPP
#define	NUMERIC_HPP

#include<array>
#include<functional>
#include<cmath>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

#include "grid.hpp"
#include "scheduler.hpp"

namespace nanolens{
namespace numeric{

/// Inverts a function f: Domain -> Codomain via an inversion of tabulated values
template<typename DomainFieldType, unsigned DomainDimension,
         typename CodomainFieldType, unsigned CodomainDimension>
class function_inverter
{
  static_assert(DomainDimension > 0, "function_inverter: dimension of domain must be > 0");
  static_assert(CodomainDimension > 0, "function_inverter: dimension of codomain must be > 0");
  
public:
  typedef std::array<DomainFieldType, DomainDimension> domain_vector;
  typedef std::array<CodomainFieldType, CodomainDimension> codomain_vector;
  typedef std::function<codomain_vector (const domain_vector&)> function_type;
  
  typedef util::grid<CodomainFieldType, std::vector<domain_vector>, CodomainDimension> table_type;
  typedef std::shared_ptr<table_type> table_ptr_type;
  
  typedef util::grid<CodomainFieldType, std::vector<domain_vector>, CodomainDimension> grid_type;
  
  function_inverter(const codomain_vector& codomain_min_extent,
                    const codomain_vector& codomain_max_extent,
                    const domain_vector& domain_min_extent,
                    const domain_vector& domain_max_extent,
                    std::size_t num_codomain_bins_per_dim,
                    const boost::mpi::communicator& comm)
  : _lookup_table(nullptr),
    _min_extent(codomain_min_extent),
    _max_extent(codomain_max_extent),
    _domain_min_extent(domain_min_extent),
    _domain_max_extent(domain_max_extent),
    _num_codomain_bins(num_codomain_bins_per_dim),
    _comm(comm)
  {
    for(std::size_t i = 0; i < DomainDimension; ++i)
      assert(_domain_min_extent[i] <= _domain_max_extent[i]);
    
    _rank = _comm.rank();
  }
  
  scheduler create_schedule(const function_type& f,
                            std::size_t num_steps_per_dim,
                            std::size_t benchmark_size) const
  {
    table_type test_table(_min_extent, _max_extent, _num_codomain_bins);
    
    scheduler job_scheduler;
    
    std::size_t num_jobs = num_steps_per_dim;
    
    job_scheduler.run(num_jobs, 
                      [&]()
    {
                        
      domain_vector sampling_center;
      for(std::size_t i = 0; i < DomainDimension; ++i)
        sampling_center[i] = static_cast<DomainFieldType>(0.5 * (_domain_min_extent[i] 
                                                               + _domain_max_extent[i]));
      
      for(std::size_t i = 0; i < benchmark_size; ++i)
      {
        domain_vector current_position = sampling_center;
        
        DomainFieldType dim0_range = _domain_max_extent[0] - _domain_min_extent[0];
        
        DomainFieldType relative_position = static_cast<DomainFieldType>(i) / 
                                            static_cast<DomainFieldType>(benchmark_size);
        
        current_position[0] = _domain_min_extent[0] + relative_position * dim0_range;
        
        codomain_vector function_evaluation_result = f(current_position);
        
        if(is_codomain_vector_within_range(function_evaluation_result)) 
          test_table[function_evaluation_result].push_back(current_position);
      }
    });
    
    return job_scheduler;
  }
  
  void run(const function_type& f,
           std::size_t num_steps_per_dim,
           std::size_t scheduler_benchmark_size)
  {
    scheduler job_scheduler = create_schedule(f, num_steps_per_dim, 
                                              scheduler_benchmark_size);
    
    run(f, num_steps_per_dim, job_scheduler);
  }
  
  void run(const function_type& f,
           std::size_t num_steps_per_dim,
           const scheduler& process_schedule)
  {
    _lookup_table = table_ptr_type(new table_type(_min_extent, _max_extent, _num_codomain_bins));
    
    std::size_t num_jobs = num_steps_per_dim;
    
    domain_vector step_sizes;
    for(std::size_t i = 0; i < DomainDimension; ++i)
    {
      step_sizes[i] = (_domain_max_extent[i] - _domain_min_extent[i]) / 
              (static_cast<DomainFieldType>(num_steps_per_dim));
    }
    
    scheduler job_scheduler(_comm);
    
    
    std::size_t steps_per_job = 1;

    for(std::size_t dimension = 1; dimension < DomainDimension; ++dimension)
      steps_per_job *= num_steps_per_dim;   

    // job id corresponds to steps in x direction
    for(std::size_t job_id = 0; job_id < num_jobs; ++job_id)
    {
      domain_vector current_position = _domain_min_extent;
      current_position[0] = static_cast<DomainFieldType>(job_id * step_sizes[0] + _domain_min_extent[0]);     
      
      std::array<std::size_t, DomainDimension> step_indices;
      for(std::size_t i = 0; i < DomainDimension; ++i)
        step_indices[i] = 0;
      
      if(job_scheduler.is_scheduled_to_this_process(job_id))
      {
        for(std::size_t total_step_id = 0; total_step_id < steps_per_job; ++total_step_id)
        {
          
          for(std::size_t i = 1; i < DomainDimension; ++i)
            current_position[i] = step_indices[i] * step_sizes[i] + _domain_min_extent[i];
          
          codomain_vector function_evaluation_result = f(current_position);
        
          if(is_codomain_vector_within_range(function_evaluation_result)) 
            _lookup_table[function_evaluation_result].push_back(current_position);          
          
          // Increment step id
          increment_step_indices(step_indices, num_steps_per_dim);
        } 
      }
    }
  }
  
  void combine_on_master(int master_rank)
  {
    assert(_lookup_table != nullptr);
    
    auto concatenate = [](typename grid_type::value_type& a, const typename grid_type::value_type& b)
    {
      for(auto element : b)
        a.push_back(element);
    };
    
    _lookup_table.combine_parallel_grid(_comm, master_rank,concatenate);
  }
  
  void allcombine()
  {
    assert(_lookup_table != nullptr);
    
    auto concatenate = [](typename grid_type::value_type& a, const typename grid_type::value_type& b)
    {
      for(auto element : b)
        a.push_back(element);
    };
    
    _lookup_table.allcombine_parallel_grid(_comm, concatenate);
  }
  
  const std::vector<domain_vector>& inverse(const codomain_vector& v) const
  {
    return _lookup_table[v];
  }
private:
  inline void increment_step_indices(std::array<std::size_t, DomainDimension>& idx,
                                     std::size_t steps_per_dim) const
  {
    for(int i = DomainDimension - 1; i >= 0; --i)
    {
      if((idx[i] + 1) < steps_per_dim)
      {
        ++(idx[i]);
        return;
      }
      else
        idx[i] = 0;
    }
  }
  
  inline bool is_codomain_vector_within_range(const codomain_vector& v) const
  {
    for(std::size_t i = 0; i < CodomainDimension; ++i)
      if(v[i] > _max_extent[i] || v[i] < _min_extent[i])
        return false;
    
    return true;
  }
  
  grid_type _lookup_table;
  codomain_vector _min_extent;
  codomain_vector _max_extent;
  domain_vector _domain_min_extent;
  domain_vector _domain_max_extent;
  std::size_t _num_codomain_bins;
  
  boost::mpi::communicator _comm;
  int _rank;
};

template<class DomainScalarType, std::size_t Domain_dim,
  class CodomainScalarType, std::size_t Codomain_dim>
class derivative
{
public:
  typedef util::vector<DomainScalarType, Domain_dim> domain_vector;
  typedef util::vector<CodomainScalarType, Codomain_dim> codomain_vector;
  typedef std::function<codomain_vector (const domain_vector&)> function_type;
  
  derivative(const function_type& f, DomainScalarType accuracy)
  : _f(f), _accuracy(accuracy), _half_accuracy(0.5 * accuracy)
  {}
  
  // @param direction A normalized vector indicating the direction of the derivative
  codomain_vector operator()(const domain_vector& direction, const domain_vector& pos) const
  {
    domain_vector p0 = pos;
    domain_vector p1 = pos;
    
    for(std::size_t i = 0; i < Domain_dim; ++i)
    {
      p0[i] -= _half_accuracy * direction[i];
      p1[i] += _half_accuracy * direction[i];
    }
    
    codomain_vector v0 = _f(p0);
    codomain_vector v1 = _f(p1);
    
    for(std::size_t i = 0; i < Codomain_dim; ++i)
    {
      v1[i] -= v0[i];
      v1[i] /= _accuracy;
    }
    
    return v1;
  }
  
protected:
  function_type _f;
  const DomainScalarType _accuracy;
  const DomainScalarType _half_accuracy;
};

// derivative in one coordinate direction
template<class DomainScalarType, std::size_t Domain_dim,
  class CodomainScalarType, std::size_t Codomain_dim, std::size_t I>
class d_dx_i : public derivative<DomainScalarType, Domain_dim, CodomainScalarType, Codomain_dim>
{
  static_assert(I < Domain_dim, "Dimension too small");
public:
  typedef util::vector<DomainScalarType, Domain_dim> domain_vector;
  typedef util::vector<CodomainScalarType, Codomain_dim> codomain_vector;
  typedef std::function<codomain_vector (const domain_vector&)> function_type;
  
  d_dx_i(const function_type& f, DomainScalarType accuracy)
  : derivative<DomainScalarType, Domain_dim, CodomainScalarType, Codomain_dim>(f, accuracy)
  {
    for(std::size_t i = 0; i < Domain_dim; ++i)
      _direction[i] = 0.;
    _direction[I] = 1.;
  }
  
  codomain_vector operator()(const domain_vector& pos) const
  {
    return derivative<DomainScalarType, 
                      Domain_dim, 
                      CodomainScalarType,
                      Codomain_dim>::operator ()(_direction, pos);
  }
  
private:
  domain_vector _direction;
};

template<class DomainScalarType, std::size_t Domain_dim,
  class CodomainScalarType, std::size_t Codomain_dim>
using d_dx = d_dx_i<DomainScalarType, Domain_dim, CodomainScalarType, Codomain_dim, 0>;

template<class DomainScalarType, std::size_t Domain_dim,
  class CodomainScalarType, std::size_t Codomain_dim>
using d_dy = d_dx_i<DomainScalarType, Domain_dim, CodomainScalarType, Codomain_dim, 1>;

template<class DomainScalarType, std::size_t Domain_dim,
  class CodomainScalarType, std::size_t Codomain_dim>
using d_dz = d_dx_i<DomainScalarType, Domain_dim, CodomainScalarType, Codomain_dim, 2>;

template<typename T>
constexpr T kronecker_delta(T a, T b)
{
  return (a == b) ? 1. : 0.;
}

template<class ScalarType, std::size_t N>
class jacobian
{
public:
  typedef util::vector<ScalarType, N> Rn_vector;
  typedef std::function<Rn_vector (const Rn_vector&)> function_type;
  
  jacobian() = default;
  
  jacobian(const Rn_vector& pos, ScalarType differential_delta,
           const function_type& f)
  : _f(f)
  {
    calculate_jacobian(pos, differential_delta);
  }
  
  std::array<ScalarType, N>& operator[](std::size_t i)
  {
    assert(i < N);
    return _j[i];
  }
  
  const std::array<ScalarType, N>& operator[](std::size_t i) const
  {
    assert(i < N);
    return _j[i];
  }
  
  util::matrix_nxn<ScalarType, N>& matrix()
  {
    return _j;
  }
  
  const util::matrix_nxn<ScalarType, N>& matrix() const
  {
    return _j;
  }
  
private:
  void calculate_jacobian(const Rn_vector& pos, ScalarType differential_delta)
  { 
    for(std::size_t i = 0; i < N; ++i)
    {
      Rn_vector point1 = pos;
      Rn_vector point2 = pos;
      
      point1[i] += differential_delta;
      point2[i] -= differential_delta;
      
      Rn_vector point1_value = _f(point1);
      Rn_vector point2_value = _f(point2);
      
      for(std::size_t j = 0; j < N; ++j)
        // df_j / dx_i
        _j[j][i] = 0.5 * (point1_value[j] - point2_value[j]) / differential_delta;
      
    }
  }
  
  function_type _f;
  util::matrix_nxn<ScalarType, N> _j;
};

template<class ScalarType, std::size_t N>
void invert_matrix(const util::matrix_nxn<ScalarType, N>& mat,
                   util::matrix_nxn<ScalarType, N>& out)
{
  static_assert(N == 2, "Matrix inversion currently only supported for 2x2 matrices.");
  
  ScalarType a = mat[0][0];
  ScalarType b = mat[0][1];
  ScalarType c = mat[1][0];
  ScalarType d = mat[1][1];
  
  ScalarType det = a * d - b * c;
  
  out[0][0] = d / det;
  out[0][1] = -b / det;
  out[1][0] = -c / det;
  out[1][1] = a / det;
}

// Newton method for a function f: R^n -> R^n
template<class ScalarType, std::size_t Dimension>
class newton
{
public:
  typedef std::array<ScalarType, Dimension> Rn_vector;
  typedef std::function<Rn_vector (const Rn_vector&)> function_type;
  
  
  newton(const Rn_vector& starting_point,
          ScalarType differential_delta,
          const function_type& f)
  : _delta(differential_delta),
    _current_position(starting_point),
    _resid(0.0),
    _f(f),
    _requested_tolerance(0.0),
    _last_step_size(0.0),
    _num_iter(0)
  {}
  
  
  void single_iteration()
  {
    this->_jacobian = jacobian<ScalarType,Dimension>(this->_current_position,
                                     this->_delta, _f);
    
    util::matrix_nxn<ScalarType, Dimension> J_inverse;
    
    invert_matrix(_jacobian.matrix(), J_inverse);
    
    Rn_vector function_value = _f(this->_current_position);
    
    // correction = -J^(-1)*function_value
    Rn_vector correction;
    for(std::size_t i = 0; i < Dimension; ++i)
    {
      correction[i] = 0.0;
      for(std::size_t j = 0; j < Dimension; ++j)
        correction[i] += (J_inverse[i][j] * function_value[j]);
    }
    
    _last_step_size = 0.0;
    
    for(std::size_t i = 0; i < Dimension; ++i)
    {
      this->_current_position[i] -= correction[i];
      this->_last_step_size += util::square(correction[i]);
    }
    
        // Update residual
    Rn_vector residual = _f(_current_position);
    this->_resid = 0.0;

    for(std::size_t i = 0; i < Dimension; ++i)
      this->_resid += std::abs(residual[i]);
      
  } 
  
  
  /// @param tolerance the squared distance to the root that the result is allowed to have
  void run(ScalarType tolerance, std::size_t max_iterations)
  {
    _requested_tolerance = tolerance;
    _last_step_size = std::numeric_limits<util::scalar>::max();
    _resid = std::numeric_limits<util::scalar>::max();
    
    for(_num_iter = 0; _num_iter < max_iterations && _resid > tolerance;
      ++_num_iter)
    {
      single_iteration();
    }
  }
  
  /// Once a root has been approximated, returns f(approximation)
  ScalarType get_residual() const
  {
    return _resid;
  }
  
  ScalarType get_squared_distance_error() const
  {
    return _last_step_size;
  }
  
  inline bool was_successful() const
  {
    if(_requested_tolerance == 0.0)
      return false;
    
    return _resid <= _requested_tolerance;
  }
  
  const util::vector2& get_position() const
  {
    return _current_position;
  }
  
  std::size_t get_num_iterations() const
  {
    return _num_iter;
  }
  
  const jacobian<ScalarType, Dimension>& get_current_jacobian() const
  {
    return _jacobian;
  }
private:
  jacobian<ScalarType, Dimension> _jacobian;
  ScalarType _requested_tolerance;
  ScalarType _delta;
  ScalarType _resid;
  ScalarType _last_step_size;
  std::size_t _num_iter;
  Rn_vector _current_position;
  function_type _f;
};

}
}

#endif	/* NUMERIC_HPP */

