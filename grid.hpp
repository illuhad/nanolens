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

#ifndef GRID_HPP
#define	GRID_HPP

#include <array>
#include "util.hpp"


namespace nanolens{
namespace util{


template<typename FieldType, class ValueType, std::size_t Dimension, bool EdgeLayer = true>
class grid
{
public:
  typedef std::array<FieldType, Dimension> scalar_array_type;
  typedef std::array<std::size_t, Dimension> index_type;
  
  typedef typename util::multi_array<ValueType>::iterator iterator;
  typedef typename util::multi_array<ValueType>::const_iterator const_iterator;
  
  typedef ValueType value_type;
  
  grid(){}
  
  grid(const scalar_array_type& min_extent,
       const scalar_array_type& max_extent,
       const std::array<std::size_t, Dimension> num_buckets)
  : _min(min_extent), _max(max_extent),
    _num_buckets(num_buckets)
  {
    std::vector<std::size_t> extended_num_buckets;
    if(EdgeLayer)
      for(std::size_t i = 0; i < Dimension; ++i)
        extended_num_buckets.push_back(num_buckets[i] + 2);
    else
      for(std::size_t i = 0; i < Dimension; ++i)
        extended_num_buckets.push_back(num_buckets[i]);
    
    _data = util::multi_array<ValueType>(extended_num_buckets);
    for(std::size_t i = 0; i < Dimension; ++i)
    {
      assert(_max[i] > _min[i]);
      _stepwidths[i] = (_max[i] - _min[i]) / static_cast<FieldType>(_num_buckets[i]);
      _half_stepwidth[i] = 0.5 * _stepwidths[i];
    }
  }
  
  ValueType& operator[](const scalar_array_type& position)
  {
    return _data[get_index(position).data()];
  }
  
  const ValueType& operator[](const scalar_array_type& position) const
  {
    return _data[get_index(position).data()];
  }
  
  ValueType& operator[](const index_type& position)
  {    
    return _data[position.data()];
  }
  
  const ValueType& operator[](const index_type& position) const
  {
    return _data[position.data()];
  }
  
  iterator begin()
  {
    return _data.begin();
  }
  
  const_iterator begin() const
  {
    return _data.begin();
  }

  iterator end()
  {
    return _data.end();
  }

  const_iterator end() const
  {
    return _data.end();
  }

  size_t get_num_entries_per_dim() const
  {
    return _num_buckets;
  }
  
  size_t get_num_entries() const
  {
    return _data.get_num_elements();
  }
  
  util::multi_array<ValueType>& data() 
  {
    return _data;
  }
  
  const util::multi_array<ValueType>& data() const
  {
    return _data;
  }
  
  template<class Function>
  void combine_parallel_grid(const boost::mpi::communicator& comm,
                             int master_rank,
                             Function combination_method)
  {
    int rank = comm.rank();
    int nprocs = comm.size();
    
    if(rank != master_rank)
    {
      // Use synchronous communication to save memory (the data tables can be quite large)
      comm.send(master_rank, 0, _data);
    }
    else
    {
      for(int proc = 0; proc < nprocs; ++proc)
      {
        if(proc != master_rank)
        {
          util::multi_array<ValueType> recv_data;
          comm.recv(proc, 0, recv_data);
          
          assert(recv_data.get_dimension() == _data.get_dimension());
          for(std::size_t dim = 0; dim < recv_data.get_dimension(); ++dim)
            assert(recv_data.get_extent_of_dimension(dim) == _data.get_extent_of_dimension(dim));
          
          auto recv_element = recv_data.begin();
          auto data_element = _data.begin();
          
          for(; recv_element != recv_data.end(); ++recv_element, ++data_element)
            combination_method(*data_element, *recv_element);
          
        }
      }
    }   
  }
  
  template<class Function>
  void allcombine_parallel_grid(const boost::mpi::communicator& comm,
                             Function combination_method)
  {
    int master_rank = 0;
    
    combine_parallel_grid(comm, master_rank, combination_method);
    
    boost::mpi::broadcast(comm, this->_data, master_rank);
  }
  
  inline bool contains_point(const scalar_array_type& point) const
  {
    for(std::size_t i = 0; i < Dimension; ++i)
      if(point[i] < _min[i] || point[i] >= _max[i])
        return false;
    
    return true;
  }
  
  inline std::size_t get_total_num_buckets(std::size_t dimension) const
  {
    return _data.get_extent_of_dimension(dimension);
  }
  
    
  inline index_type get_index(const scalar_array_type& point) const
  {
    index_type result;
    
    if(EdgeLayer)
    {
      for(size_t i = 0; i < Dimension; ++i)
      {
        if(point[i] < _min[i])
          result[i] = 0;
        else if(point[i] >= _max[i])
          result[i] = _data.get_extent_of_dimension(i) - 1;
        else
          result[i] = static_cast<std::size_t>((point[i] - _min[i]) / _stepwidths[i]) + 1;
      }
    }
    else
    {
      for(size_t i = 0; i < Dimension; ++i)
      {
        result[i] = static_cast<std::size_t>((point[i] - _min[i]) / _stepwidths[i]);
      }
    }
    
    return result;
  }
  
  index_type get_interior_start_bucket() const
  {
    index_type result;
    std::fill(result.begin(), result.end(), 0);
    
    if(EdgeLayer)
      for(std::size_t i = 0; i < Dimension; ++i)
        result[i] = 1;
    
    return result;
  }
  
  // one bucket beyond the final interior bucket
  index_type get_interior_end_bucket() const
  {
    index_type result;
    
    for(std::size_t i = 0; i < Dimension; ++i)
      if(EdgeLayer)
        result[i] = _data.get_extent_of_dimension(i) - 1;
      else
        result[i] = _data.get_extent_of_dimension(i);
    
    return result;
  }
  
  scalar_array_type get_min_position_of_bucket(const index_type& p) const
  {
    std::array<int, Dimension> exterior_layer_corrected_position;
    for(std::size_t i = 0; i < Dimension; ++i)
      exterior_layer_corrected_position[i] = static_cast<int>(p[i]);
    
    if(EdgeLayer)
      for(std::size_t i = 0; i < Dimension; ++i)
        --exterior_layer_corrected_position[i];
    
    scalar_array_type position;
    
    for(std::size_t i = 0; i < Dimension; ++i)
    {
      position[i] = _min[i] + exterior_layer_corrected_position[i] * _stepwidths[i];
    }
    
    return position;
  }
  
  scalar_array_type get_central_position_of_bucket(const index_type& p) const
  {
    scalar_array_type pos = get_min_position_of_bucket(p);
    util::add(pos, _half_stepwidth);
    
    return pos;
  }
  
  scalar_array_type get_bucket_size() const
  {
    return _stepwidths;
  }
private:
  
  util::multi_array<ValueType> _data;
  std::array<std::size_t, Dimension> _num_buckets;
  
  scalar_array_type _min;
  scalar_array_type _max;
  scalar_array_type _stepwidths;
  scalar_array_type _half_stepwidth;
  
};

template<typename KeyType, class ValueType>
using grid2d = grid<KeyType, ValueType, 2>;

template<typename KeyType, class ValueType>
using grid3d = grid<KeyType, ValueType, 3>;

}
}

#endif	/* GRID_HPP */

