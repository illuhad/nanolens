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
#include <boost/circular_buffer.hpp>
#include "util.hpp"


namespace nanolens{
namespace util{

namespace grid_impl_{

template<class Vector_type, class Index_type, std::size_t Dimension, bool With_edge_layer>
struct translator_impl
{
};

template<class Vector_type, class Index_type, std::size_t Dimension>
struct translator_impl<Vector_type, Index_type, Dimension, true>
{
  static inline Index_type to_grid_index(const Vector_type& position,
                                         const Vector_type& min_extent,
                                         const Vector_type& max_extent,
                                         const Vector_type& stepwidths,
                                         const Index_type& num_buckets)
  {
    Index_type result;
    for(size_t i = 0; i < Dimension; ++i)
    {
      if(position[i] <= min_extent[i])
        result[i] = 0;
      else if(position[i] >= max_extent[i])
        result[i] = num_buckets[i] - 1;
      else
        result[i] = static_cast<std::size_t>((position[i] - min_extent[i]) / stepwidths[i]) + 1;
    }
    return result;
  }
  
  static inline Vector_type to_bucket_min_position(const Index_type& idx,
                                                   const Vector_type& min_extent,
                                                   const Vector_type& stepwidths)
  {
    Vector_type position;
    
    for(std::size_t i = 0; i < Dimension; ++i)
      position[i] = min_extent[i] + (static_cast<long long int>(idx[i]) - 1) * stepwidths[i];
    
    return position;
  }
};

template<class Vector_type, class Index_type, std::size_t Dimension>
struct translator_impl<Vector_type, Index_type, Dimension, false>
{
  static inline Index_type to_grid_index(const Vector_type& position,
                                         const Vector_type& min_extent,
                                         const Vector_type& max_extent,
                                         const Vector_type& stepwidths,
                                         const Index_type& num_buckets)
  {
    Index_type result;
    
    for(size_t i = 0; i < Dimension; ++i)
    {
      result[i] = static_cast<std::size_t>((position[i] - min_extent[i]) / stepwidths[i]);
    }
    
    return result;
  }
  
  static inline Vector_type to_bucket_min_position(const Index_type& idx,
                                                   const Vector_type& min_extent,
                                                   const Vector_type& stepwidths)
  { 
    Vector_type position;
    
    for(std::size_t i = 0; i < Dimension; ++i)
      position[i] = min_extent[i] + idx[i] * stepwidths[i];
    
    return position;
  }
};

}
// Translates a point into a grid index
template<class FieldType, std::size_t Dimension, bool Edge_layer>
class grid_translator
{
public:
  typedef std::array<FieldType, Dimension> scalar_array_type;
  typedef std::array<std::size_t, Dimension> index_type;
  
  grid_translator() = default;
  
  grid_translator(const scalar_array_type& min_extent,
                  const scalar_array_type& max_extent,
                  const index_type& num_buckets)
  : _min_extent(min_extent), _max_extent(max_extent),
    _num_buckets(num_buckets)
  {
    for(std::size_t i = 0; i < Dimension; ++i)
      assert(_max_extent[i] > _min_extent[i]);
    
    _interior_num_buckets = _num_buckets;
    
    if(Edge_layer)
      for(std::size_t i = 0; i < _num_buckets.size(); ++i)
        _num_buckets[i] += 2;
    
    
    for(std::size_t dim = 0; dim < Dimension; ++dim)
    {
      _stepwidths[dim] = (_max_extent[dim] - _min_extent[dim]) / static_cast<FieldType>(_interior_num_buckets[dim]);
      _half_stepwidth[dim] = 0.5 * _stepwidths[dim];
    }

  }
  
  index_type operator()(const scalar_array_type& position) const
  {
    return grid_impl_::translator_impl
    <
      scalar_array_type, 
      index_type, 
      Dimension, 
      Edge_layer
    >::to_grid_index(position, _min_extent, _max_extent, _stepwidths, _num_buckets);
  }
  
  const scalar_array_type& get_min_extent() const
  { return _min_extent; }
  
  const scalar_array_type& get_max_extent() const
  { return _max_extent; }
  
  const scalar_array_type& get_bucket_size() const
  { return _stepwidths; }
  
  const scalar_array_type& get_half_bucket_size() const
  { return _half_stepwidth; }
  
  const index_type& get_num_buckets_per_dim() const
  { return _num_buckets; }
  
  const index_type& get_interior_num_buckets_per_dim() const
  { return _interior_num_buckets; }
  
  index_type get_interior_start_bucket() const
  {
    index_type result;
    
    if(Edge_layer)
      std::fill(result.begin(), result.end(), 1);
    else
      std::fill(result.begin(), result.end(), 0);
    
    return result;
  }
  
    // one bucket beyond the final interior bucket
  index_type get_interior_end_bucket() const
  {
    index_type result;
    
    for(std::size_t i = 0; i < Dimension; ++i)
      if(Edge_layer)
        result[i] = _num_buckets[i] - 1;
      else
        result[i] = _num_buckets[i];
    
    return result;
  }
  
  scalar_array_type get_min_position_of_bucket(const index_type& p) const
  {
    return grid_impl_::translator_impl
    <
      scalar_array_type, 
      index_type, 
      Dimension, 
      Edge_layer
    >::to_bucket_min_position(p, this->_min_extent, this->_stepwidths);
  }
  
    
  scalar_array_type get_central_position_of_bucket(const index_type& p) const
  {
    scalar_array_type pos = get_min_position_of_bucket(p);
    for(std::size_t i = 0; i < Dimension; ++i)
      pos[i] += _half_stepwidth[i];
    
    return pos;
  }
  
  inline bool contains_point(const scalar_array_type& point) const
  {
    for(std::size_t i = 0; i < Dimension; ++i)
      if(point[i] < _min_extent[i] || point[i] >= _max_extent[i])
        return false;
    
    return true;
  }
private:
  scalar_array_type _min_extent;
  scalar_array_type _max_extent;
  scalar_array_type _stepwidths;
  scalar_array_type _half_stepwidth;
  index_type _num_buckets;
  index_type _interior_num_buckets;
};

template<typename Field_type, class Value_type, std::size_t Dimension, bool Edge_layer = true>
class grid
{
public:
  typedef typename grid_translator<Field_type, Dimension, Edge_layer>::scalar_array_type scalar_array_type;
  typedef typename grid_translator<Field_type, Dimension, Edge_layer>::index_type index_type;

  typedef typename util::multi_array<Value_type>::iterator iterator;
  typedef typename util::multi_array<Value_type>::const_iterator const_iterator;
  
  typedef Value_type value_type;
  
  grid(){}
  
  grid(const scalar_array_type& min_extent,
       const scalar_array_type& max_extent,
       const std::array<std::size_t, Dimension> num_buckets)
  {
    
    _translator = grid_translator<Field_type,Dimension,Edge_layer>(min_extent, 
                                                                 max_extent, 
                                                                 num_buckets);
   
    index_type extended_num_buckets = num_buckets;
    
    // Make sure we have enough space for the EdgeLayer
    if(Edge_layer)
      for(std::size_t i = 0; i < Dimension; ++i)
        extended_num_buckets[i] += 2;
    
    _data = util::multi_array<Value_type>(
      std::move(util::array_to_vector(extended_num_buckets)));
    
  }
  
  Value_type& operator[](const scalar_array_type& position)
  {
    return _data[get_index(position).data()];
  }
  
  const Value_type& operator[](const scalar_array_type& position) const
  {
    return _data[get_index(position).data()];
  }
  
  Value_type& operator[](const index_type& position)
  {    
    return _data[position.data()];
  }
  
  const Value_type& operator[](const index_type& position) const
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
  
  size_t get_num_entries() const
  {
    return _data.get_num_elements();
  }
  
  util::multi_array<Value_type>& data() 
  {
    return _data;
  }
  
  const util::multi_array<Value_type>& data() const
  {
    return _data;
  }
  
  template<class Function>
  void combine_parallel_grid(const boost::mpi::communicator& comm,
                             int master_rank,
                             Function combination_method)
  {
    _data.reduce_parallel_array(comm, master_rank, combination_method);  
  }
  
  template<class Function>
  void allcombine_parallel_grid(const boost::mpi::communicator& comm,
                             Function combination_method)
  {
    _data.allreduce_parallel_array(comm, combination_method);
  }
  
  inline bool contains_point(const scalar_array_type& point) const
  {
    return _translator.contains_point(point);
  }
  
  inline std::size_t get_total_num_buckets(std::size_t dimension) const
  {
    return _data.get_extent_of_dimension(dimension);
  }
  
  inline index_type get_index(const scalar_array_type& point) const
  {
    return _translator(point);
  }
  
  index_type get_interior_start_bucket() const
  {
    return _translator.get_interior_start_bucket();
  }
  
  // one bucket beyond the final interior bucket
  index_type get_interior_end_bucket() const
  {
    return _translator.get_interior_end_bucket();
  }
  
  scalar_array_type get_min_position_of_bucket(const index_type& p) const
  {
    return _translator.get_min_position_of_bucket(p);
  }

  scalar_array_type get_central_position_of_bucket(const index_type& p) const
  {
    return _translator.get_central_position_of_bucket(p);
  }
  
  scalar_array_type get_bucket_size() const
  {
    return _translator.get_bucket_size();
  }
private:
  grid_translator<Field_type,Dimension,Edge_layer> _translator;
  
  util::multi_array<Value_type> _data;

  scalar_array_type _stepwidths;
  scalar_array_type _half_stepwidth;
  
};

template<typename KeyType, class ValueType>
using grid2d = grid<KeyType, ValueType, 2>;

template<typename KeyType, class ValueType>
using grid3d = grid<KeyType, ValueType, 3>;

template<class Field_type, std::size_t Dimension, class T>
class buffered_grid_cache
{
public:
    
  static_assert(Dimension > 0, "Dimension of cache grid must be at least 1");
  
  class buffer_entry;
  
  typedef boost::circular_buffer<buffer_entry> buffer_type;
  typedef buffer_entry* data_ptr_type;
  
  typedef grid<Field_type, data_ptr_type, Dimension, true> grid_type;
  
  
  
  typedef typename grid_type::scalar_array_type scalar_array_type;
  typedef typename grid_type::index_type index_type;
  
  class buffer_entry
  {
    index_type _grid_position;
    T _data;
    
  public:
    buffer_entry(const index_type& grid_position,
                 const T& data)
    : _grid_position(grid_position),
      _data(data)
    {}
    
    buffer_entry()
    : _grid_position({-1, -1})
    {}
    
    const index_type& get_grid_position() const
    { return _grid_position; }
    
    T& data()
    { return _data; }
    
    const T& data() const
    { return _data; }
  };
  
  typedef std::function<void (T&, const scalar_array_type&)> entry_initialization_function;
  
  buffered_grid_cache() = default;
  buffered_grid_cache(const buffered_grid_cache&) = delete;
  buffered_grid_cache& operator=(const buffered_grid_cache&) = delete;
  
  
  buffered_grid_cache(const scalar_array_type& min_extent,
                      const scalar_array_type& max_extent,
                      const index_type& num_buckets,
                      entry_initialization_function f,
                      std::size_t cache_size)
  : _grid(min_extent, max_extent, num_buckets),
    _init_func(f), _buffer(cache_size)
  {
    std::fill(_grid.begin(), _grid.end(), nullptr);
  }
  
  const scalar_array_type get_cell_size() const
  { return _grid.get_bucket_size(); }
  
    
  T& operator[](const scalar_array_type& pos)
  {
    index_type grid_index = _grid.get_index(pos);
    data_ptr_type data_iterator = _grid[grid_index];
    
    if(data_iterator == nullptr)
    {
      // Create new cell
      return create_new_cell(grid_index);
    }
    else
    {
      if(data_iterator->get_grid_position() == grid_index)
        return data_iterator->data();
      else
      {
        // Create new cell
        return create_new_cell(grid_index);
      }
    }
    
  }
private:
  T& create_new_cell(const index_type& grid_index)
  {
    scalar_array_type center = _grid.get_central_position_of_bucket(grid_index);
    
    T data;
    _init_func(data, center);
    
    _buffer.push_back(buffer_entry(grid_index, data));
    _grid[grid_index] = &(_buffer[_buffer.size() - 1]);
    
    return _buffer.rbegin()->data();
  }
  
  entry_initialization_function _init_func;
  grid_type _grid;
  buffer_type _buffer;
  
};

}
}

#endif	/* GRID_HPP */

