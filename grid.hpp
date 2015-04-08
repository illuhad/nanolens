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

template<typename FieldType, class ValueType, size_t Dimension>
class grid
{
public:
  typedef std::array<FieldType, Dimension> scalar_array_type;
  typedef std::array<std::size_t, Dimension> index_type;
  
  typedef typename util::multi_array<ValueType>::iterator iterator;
  typedef typename util::multi_array<ValueType>::const_iterator const_iterator;
  
  typedef ValueType value_type;
  
  grid(const scalar_array_type& min_extent,
       const scalar_array_type& max_extent,
       size_t num_buckets_per_dim)
  : _data(std::vector<size_t>(Dimension, num_buckets_per_dim + 2)),
    _min(min_extent), _max(max_extent),
    _num_buckets(num_buckets_per_dim)
  {
    for(std::size_t i = 0; i < Dimension; ++i)
    {
      assert(_max[i] > _min[i]);
      _stepwidths[i] = (_max[i] - _min[i]) / num_buckets_per_dim;
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
private:
  
  
  index_type get_index(const scalar_array_type& point) const
  {
    index_type result;
    
    for(size_t i = 0; i < Dimension; ++i)
    {
      if(point[i] < _min[i])
        result[i] = 0;
      else if(point[i] >= _max[i])
        result[i] = _data.get_extent_of_dimension(i) - 1;
      else
        result[i] = static_cast<size_t>((point[i] - _min[i]) / _stepwidths[i]) + 1;
    }
    
    return result;
  }
  
  util::multi_array<ValueType> _data;
  FieldType _num_buckets;
  
  scalar_array_type _min;
  scalar_array_type _max;
  scalar_array_type _stepwidths;
  
};

template<typename KeyType, class ValueType>
using grid2d = grid<KeyType, ValueType, 2>;

template<typename KeyType, class ValueType>
using grid3d = grid<KeyType, ValueType, 3>;

}
}

#endif	/* GRID_HPP */

