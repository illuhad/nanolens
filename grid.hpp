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

/// Implementation details of the grid - should not be used by the user.
namespace grid_impl_{

/// Implementation of a grid tranlator core, i.e. functionality to transform
/// between \f$ R^n \f$ coordinates and the index of a grid bin.
/// There are two template specializations: With enabled and with disabled edge layer.
/// The edge layer is a layer of additional grid cells around the actual grid.
/// If enabled, all coordinates that lie outside the actual grid will be mapped
/// to the closest edge layer.
/// \tparam Vector_type the data type to use for vectors
/// \tparam Index_type the data type used for grid indices (must be some sort of
/// array to support multidimensional grids)
/// \tparam Dimension The dimension of the grid
/// \tparam With_edge_layer Whether to use an edge layer cell
template<class Vector_type, class Index_type, std::size_t Dimension, bool With_edge_layer>
struct translator_impl
{
};

/// Specialization for an enabled edge layer
template<class Vector_type, class Index_type, std::size_t Dimension>
struct translator_impl<Vector_type, Index_type, Dimension, true>
{
  /// Translates a coordinate vector into a grid index vector.
  /// \param position the coordinate vector
  /// \param min_extent The coordinates of the corner of the grid with the
  /// mimimum coordinate values
  /// \param max_extent The coordinates of the corner of the grid with the
  /// maximum coordinate values
  /// \param stepwidths A vector containing the step width between grid cells
  /// for each dimension
  /// \param num_buckets A vector containing the number of grid cells for each dimension
  /// \return The index of the grid cell that corresponds to the given coordinates
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
  
  
  /// Calculates for a given grid index the coordinates of the corner of the grid
  /// cell with the minimum coordinate values
  /// \param idx the grid cell index
  /// \param min_extent The coordinates of the corner of the grid with the
  /// mimimum coordinate values
  /// \param stepwidths A vector containing the step width between grid cells
  /// for each dimension
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

/// Specialization for a disabled edge layer
template<class Vector_type, class Index_type, std::size_t Dimension>
struct translator_impl<Vector_type, Index_type, Dimension, false>
{
  
  /// Translates a coordinate vector into a grid index vector.
  /// \param position the coordinate vector
  /// \param min_extent The coordinates of the corner of the grid with the
  /// mimimum coordinate values
  /// \param max_extent The coordinates of the corner of the grid with the
  /// maximum coordinate values
  /// \param stepwidths A vector containing the step width between grid cells
  /// for each dimension
  /// \param num_buckets A vector containing the number of grid cells for each dimension
  /// \return The index of the grid cell that corresponds to the given coordinates
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
  
  /// Calculates for a given grid index the coordinates of the corner of the grid
  /// cell with the minimum coordinate values
  /// \param idx the grid cell index
  /// \param min_extent The coordinates of the corner of the grid with the
  /// mimimum coordinate values
  /// \param stepwidths A vector containing the step width between grid cells
  /// for each dimension
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


/// Implementation of a grid tranlator, i.e. functionality to transform
/// between \f$ R^n \f$ coordinates and the index of a grid bin.
/// \tparam FieldType the scalar data type to use for vectors
/// \tparam Dimension The dimension of the grid
/// \tparam With_edge_layer Whether to use an edge layer.
/// The edge layer is a layer of additional grid cells around the actual grid.
/// If enabled, all coordinates that lie outside the actual grid will be mapped
/// to the closest edge layer cell. The difference between enabled and
/// disabled edge layer can visualized like this:
///  _ _               _|_|_|_
/// |_|_|              _|_|_|_ 
/// |_|_|              _|_|_|_
///                     | | |
/// A 2x2 grid without  A 2x2 grid with
/// edge layer          edge layer
///
/// Enabling an edge layer guarantees that all coordinates can be mapped to
/// a grid cell, whereas without an edge layer, coordinates outside the bounds
/// of the grid have no corresponding grid cell.
template<class FieldType, std::size_t Dimension, bool Edge_layer>
class grid_translator
{
public:
  typedef std::array<FieldType, Dimension> scalar_array_type;
  typedef std::array<std::size_t, Dimension> index_type;
  
  /// Construct object
  grid_translator() = default;
  
  /// Construct object
  /// \param min_extent The coordinates of the corner of the grid with the
  /// mimimum coordinate values
  /// \param max_extent The coordinates of the corner of the grid with the
  /// maximum coordinate values
  /// \param num_buckets A vector containing the number of grid cells for each dimension
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
  
  /// Translates a coordinate vector to its corresponding grid cell index
  /// \param position The coordinate vector that shall be translated
  /// \return The grid index vector corresponding to the given position
  /// If the edge layer has been disabled, it is the user's responsibility
  /// to check that the resulting grid index is valid and lies within
  /// the bounds of the grid.
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
  
  /// \return The coordinates of the corner of the grid with the
  /// mimimum coordinate values
  const scalar_array_type& get_min_extent() const
  { return _min_extent; }
  
  /// \return The coordinates of the corner of the grid with the
  /// maximum coordinate values
  const scalar_array_type& get_max_extent() const
  { return _max_extent; }
  
  /// \return A vector containing the size of a grid cell in each dimension
  const scalar_array_type& get_bucket_size() const
  { return _stepwidths; }
  
  /// \return A vector containing half the size of a grid cell in each dimension
  const scalar_array_type& get_half_bucket_size() const
  { return _half_stepwidth; }
  
  /// \return A vector containing the total number of grid cells for each dimension,
  /// including edge layer cells, if enabled.
  const index_type& get_num_buckets_per_dim() const
  { return _num_buckets; }
  
  /// \return A vector containing the number of non-edge layer cells for each dimension
  const index_type& get_interior_num_buckets_per_dim() const
  { return _interior_num_buckets; }
  
  /// \return A vector containing the indices of the first non-edge layer cell
  index_type get_interior_start_bucket() const
  {
    index_type result;
    
    if(Edge_layer)
      std::fill(result.begin(), result.end(), 1);
    else
      std::fill(result.begin(), result.end(), 0);
    
    return result;
  }
  
  /// \return A vector containing the indices of the first grid cell after the last
  /// non-edge layer cell
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
  
  /// Calculates the coordinates of the corner of a grid cell with the minimum
  /// coordinate values
  /// \param p The index of the grid cell
  /// \return The coordinates of the cell corner with the minimum coordinate values
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
  
  /// Calculates the coordinates of the center of a given grid cell
  /// \param p The index of the grid cell
  /// \return The coordinates of the center of the given cell
  scalar_array_type get_central_position_of_bucket(const index_type& p) const
  {
    scalar_array_type pos = get_min_position_of_bucket(p);
    for(std::size_t i = 0; i < Dimension; ++i)
      pos[i] += _half_stepwidth[i];
    
    return pos;
  }
  
  /// \return Whether a given point lies within the bounds of the grid
  /// \param point The point that shall be checked
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

/// Implements a grid using a \c grid_translator<> object to translate coordinates
/// to grid cells. In constrast to the \c grid_translator, the \c grid class
/// actually allocates memory to store the data in the grid cells.
/// \tparam Field_type The scalar type of the coordinate vectors
/// \tparam Value_type The data type stored in the grid cells
/// \tparam Dimension The dimensions of the grid
/// \tparam Edge_layer Whether to use an edge layer. See the documentation of
/// the \c grid_translator for more details.
template<typename Field_type, class Value_type, std::size_t Dimension, bool Edge_layer = true>
class grid
{
public:
  typedef typename grid_translator<Field_type, Dimension, Edge_layer>::scalar_array_type scalar_array_type;
  typedef typename grid_translator<Field_type, Dimension, Edge_layer>::index_type index_type;

  typedef typename util::multi_array<Value_type>::iterator iterator;
  typedef typename util::multi_array<Value_type>::const_iterator const_iterator;
  
  typedef Value_type value_type;
  
  /// Default construct object. The resulting grid should not be used!
  grid(){}
  
  /// Construct object
  /// \param min_extent The coordinates of the corner of the grid with the minimum
  /// coordinate values
  /// \param max_extent The coordinates of the corner of the grid with the maximum
  /// coordinate values
  /// \param num_vectors An array containing the number of grid cells in each
  /// dimension
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
  
  /// Access the content of a grid cell.
  /// \param position The coordinates of which the corresponding cell shall be
  /// retrieved. If the edge layer has been enabled, this succeeds for arbitry
  /// coordinates. Without an edge layer, it is the responsibility of the user
  /// to make sure that \c position lies within the bounds of the grid. If the
  /// coordinates are outside the grid and no edge layer is used, the operation
  /// results in undefined behavior.
  /// \return A reference to the content of the grid cell in which \c position
  /// lies.
  Value_type& operator[](const scalar_array_type& position)
  {
    return _data[get_index(position).data()];
  }
  
  /// Access the content of a grid cell.
  /// \param position The coordinates of which the corresponding cell shall be
  /// retrieved. If the edge layer has been enabled, this succeeds for arbitry
  /// coordinates. Without an edge layer, it is the responsibility of the user
  /// to make sure that \c position lies within the bounds of the grid. If the
  /// coordinates are outside the grid and no edge layer is used, the operation
  /// results in undefined behavior.
  /// \return A reference to the content of the grid cell in which \c position
  /// lies.
  const Value_type& operator[](const scalar_array_type& position) const
  {
    return _data[get_index(position).data()];
  }
  
  /// Access the content of a grid cell by the index of the cell.
  /// \param position The index of the cell that shall be accessed. It is the 
  /// user's responsibility to ensure that the index is within the bounds.
  /// \return A reference to the content of the grid cell
  Value_type& operator[](const index_type& position)
  {    
    return _data[position.data()];
  }
  
  /// Access the content of a grid cell by the index of the cell.
  /// \param position The index of the cell that shall be accessed. It is the 
  /// user's responsibility to ensure that the index is within the bounds.
  /// \return A reference to the content of the grid cell
  const Value_type& operator[](const index_type& position) const
  {
    return _data[position.data()];
  }
  
  /// \return An iterator to the first cell
  iterator begin()
  {
    return _data.begin();
  }
  
  /// \return An iterator to the first cell
  const_iterator begin() const
  {
    return _data.begin();
  }

  /// \return An iterator to one element beyond the last cell
  iterator end()
  {
    return _data.end();
  }

  /// \return An iterator to one element beyond the last cell
  const_iterator end() const
  {
    return _data.end();
  }
  
  /// \return The total number of grid cells
  std::size_t get_num_entries() const
  {
    return _data.get_num_elements();
  }
  
  /// \return A direct reference to the array containing the contents of the
  /// grid cells
  util::multi_array<Value_type>& data() 
  {
    return _data;
  }
  
  /// \return A direct reference to the array containing the contents of the
  /// grid cells
  const util::multi_array<Value_type>& data() const
  {
    return _data;
  }
  
  /// Combines each element within all parallel instances of the grid on the
  /// master process by applying a user specified combination function.
  /// \param comm The mpi communicator
  /// \param master_rank The rank of the master process. Only this process
  /// will have access to the combined grid.
  /// \param combination_method The function that will be applied to combine
  /// all elements at the same position in the parallel array instances.
  /// It must have the signature \c void(T&, const T&) and be an in-place operation
  template<class Function>
  void combine_parallel_grid(const boost::mpi::communicator& comm,
                             int master_rank,
                             Function combination_method)
  {
    _data.reduce_parallel_array(comm, master_rank, combination_method);  
  }
  
  /// Like \c combine_parallel_grid, but distributes the result to all processes
  /// \param comm The mpi communicator
  /// \param combination_method The function that will be applied to combine
  /// all elements at the same position in the parallel array instances.
  /// It must have the signature \c void(T&, const T&) and be an in-place operation
  template<class Function>
  void allcombine_parallel_grid(const boost::mpi::communicator& comm,
                             Function combination_method)
  {
    _data.allreduce_parallel_array(comm, combination_method);
  }
  
  /// \return Whether a given point lies within the bounds of the grid
  /// \param point The point that shall be checked
  inline bool contains_point(const scalar_array_type& point) const
  {
    return _translator.contains_point(point);
  }
  
  /// \return The total number of grid cells in one dimension
  /// \param dimension The index of the dimension
  inline std::size_t get_total_num_buckets(std::size_t dimension) const
  {
    return _data.get_extent_of_dimension(dimension);
  }
  
  /// \return The index of the grid cell that corresponds to a given point in space
  /// \param point The point for which the corresponding cell index shall
  /// be calculated
  inline index_type get_index(const scalar_array_type& point) const
  {
    return _translator(point);
  }
  
  /// \return A vector containing the indices of the first non-edge layer cell
  index_type get_interior_start_bucket() const
  {
    return _translator.get_interior_start_bucket();
  }
  
  /// \return A vector containing the indices of the first grid cell after the last
  /// non-edge layer cell
  index_type get_interior_end_bucket() const
  {
    return _translator.get_interior_end_bucket();
  }
    
  /// Calculates the coordinates of the corner of a grid cell with the minimum
  /// coordinate values
  /// \param p The index of the grid cell
  /// \return The coordinates of the cell corner with the minimum coordinate values
  scalar_array_type get_min_position_of_bucket(const index_type& p) const
  {
    return _translator.get_min_position_of_bucket(p);
  }

  /// Calculates the coordinates of the center of a given grid cell
  /// \param p The index of the grid cell
  /// \return The coordinates of the center of the given cell
  scalar_array_type get_central_position_of_bucket(const index_type& p) const
  {
    return _translator.get_central_position_of_bucket(p);
  }
  
  /// \return A vector containing the size of a grid cell in each dimension
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

/// A two dimensional grid
/// \tparam KeyType The scalar type of coordinate vectors
/// \tparam ValueType The data type stored in the grid cells
/// \tparam Edge_layer Wheter to use an edge layer (see \c grid_translator for
/// more details)
template<typename KeyType, class ValueType, bool Edge_layer=true>
using grid2d = grid<KeyType, ValueType, 2, Edge_layer>;

/// A three dimensional grid
/// \tparam KeyType The scalar type of coordinate vectors
/// \tparam ValueType The data type stored in the grid cells
/// \tparam Edge_layer Wheter to use an edge layer (see \c grid_translator for
/// more details)
template<typename KeyType, class ValueType, bool Edge_layer=true>
using grid3d = grid<KeyType, ValueType, 3, Edge_layer>;

/// Implements a low-latency FIFO cache for grid cells. 
/// The core of the implemtation is a circular array, that overwrites elements 
/// starting from the beginning if it gets full. 
/// This array (called "the buffer") is the actual cache and stores
/// the cell data. The important property of this array that we exploit is that
/// pointers to its elements always remain valid - the only thing that can happen
/// is that the element that the pointer is pointing to has been overwritten by
/// data from a different cell. This allows for a very efficient, low-latency
/// implementation of a cell cache: We create a grid of equal size and extent
/// as the grid whose cells we want to cache. But instead of storing the full
/// cell data in each cell, we only store a pointer to a position in the buffer
/// where the cell has been last stored. Hence, by following the pointer
/// we immediately have found the cell in the cache, or we know that the
/// cell is not in the cache if the cell at the pointer turns out to be
/// not the cell we are looking for. In this case, we need to load cell from
/// the actual grid. This is done by calling a user-supplied 
/// \c entry_initialization_function that may for example look up the element
/// from storage, or calculate its content on the fly.
/// This class is non-copyable.
/// \tparam Field_type the scalar type of coordinate vectors
/// \tparam Dimension the dimension of the grid
/// \tparam T The type of the data to be stored in the grid cells.
template<class Field_type, std::size_t Dimension, class T>
class buffered_grid_cache
{
public:
    
  static_assert(Dimension > 0, "Dimension of cache grid must be at least 1");
  
  class buffer_entry;
  
  /// This is the type used for the central cache buffer. We use the
  /// \c boost::circular_buffer class as provides everything we need.
  typedef boost::circular_buffer<buffer_entry> buffer_type;
  
  typedef buffer_entry* data_ptr_type;
  
  typedef grid<Field_type, data_ptr_type, Dimension, true> grid_type;
  
  typedef typename grid_type::scalar_array_type scalar_array_type;
  typedef typename grid_type::index_type index_type;
  
  /// This class represents entries in the cache buffer. In addition to
  /// the user data of type \c T, we also need to store the index of the
  /// grid cell to which this entry belongs, so that we can check if
  /// the entry is the one we were looking for.
  class buffer_entry
  {
    index_type _grid_position;
    T _data;
    
  public:
    /// Construct object
    /// \param grid_position The index of the grid cell to which this cache
    /// entry belongs
    /// \param data The user data that is stored in this grid cell
    buffer_entry(const index_type& grid_position,
                 const T& data)
    : _grid_position(grid_position),
      _data(data)
    {}
    
    /// Default construct object.
    buffer_entry()
    : _grid_position({-1, -1})
    {}
    
    /// \return The index of the grid cell to which this cache entry belongs. If
    /// this object has been constructed without a valid index as parameter, 
    /// will return \c {std::size_t(-1),std::size_t(-1)}
    const index_type& get_grid_position() const
    { return _grid_position; }
    
    /// \return The user data of the grid cell corresponding to this cache entry
    T& data()
    { return _data; }
    
    /// \return The user data of the grid cell corresponding to this cache entry
    const T& data() const
    { return _data; }
  };
  
  /// The type of the function that will be used to load new cells into the cache
  typedef std::function<void (T&, const scalar_array_type&)> entry_initialization_function;
  
  /// Default construct object
  buffered_grid_cache() = default;
  
  buffered_grid_cache(const buffered_grid_cache&) = delete;
  buffered_grid_cache& operator=(const buffered_grid_cache&) = delete;
  
  /// Construct object
  /// \param min_extent the coordinates of the corner of the grid with
  /// the minimum coordinate values
  /// \param max_extent the coordinates of the corner of the grid with
  /// the maximum coordinate values
  /// \param num_buckets the number of grid cells in each dimension
  /// \param f the function that will be called to load grid cells
  /// that are not yet in the cache. It must be of the signature
  /// \c void(T&,const scalar_array_type&). The first parameter
  /// will be used to pass a reference to the cell data that shall be initialized,
  /// while the second parameter will contain the coordinates of the center
  /// of the grid cell that shall be loaded into the cache.
  /// \param cache_size the maximum number of cached elements
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
  
  /// \return the size of a grid cell in each dimension
  const scalar_array_type get_cell_size() const
  { return _grid.get_bucket_size(); }
  
  /// Given a position, accesses the grid cell in which the supplied position lies.
  /// If the grid cell cannot be found in the cache, the entry initialization
  /// function that that has been supplied in the constructor will be called to load
  /// the requested grid cell into the cache.
  /// \param pos The position of which the corresponding grid cell shall be
  /// accessed
  /// \return The user data that has been stored in the requested grid cell
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
  /// Loads a new grid cell using the entry initialization funtion and inserts
  /// it into the cache.
  /// \param grid_index The index of the grid cell that shall be inserted into
  /// the cache
  /// \return A reference to the newly inserted element in the cache
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

