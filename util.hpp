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

#ifndef UTIL_HPP
#define	UTIL_HPP

#include <array>
#include <vector>
#include <cassert>
#include <algorithm>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>

namespace nanolens{
namespace util{

typedef double scalar;

template<typename ScalarType, std::size_t N>
using vector = std::array<ScalarType, N>;

typedef vector<scalar, 3> vector3;
typedef vector<scalar, 2> vector2;

template<typename ScalarType, std::size_t M, std::size_t N>
using matrix = std::array<std::array<ScalarType, N>, M>;

template<typename ScalarType, std::size_t N>
using matrix_nxn = matrix<ScalarType, N, N>;

const scalar G = 1.;
const scalar c = 1.;

constexpr scalar square(scalar x)
{ return x*x; }

inline void add(vector2& a, const vector2& b)
{
  a[0] += b[0];
  a[1] += b[1];
}

inline void sub(vector2& a, const vector2& b)
{
  a[0] -= b[0];
  a[1] -= b[1];
}

inline void scale_add(vector2& a, const vector2& b, scalar s)
{
  a[0] += s * b[0];
  a[1] += s * b[1];
}

inline void scale(vector2& a, scalar s)
{
  a[0] *= s;
  a[1] *= s;
}

inline scalar dot(const vector2& a, const vector2& b)
{
  return a[0] * b[0] + a[1] * b[1];
}

inline void assign(vector2& a, scalar s)
{
  a[0] = s;
  a[1] = s;
}

// Determinant of a matrix defined by two column vectors
inline scalar det(const vector2& a, const vector2& b)
{
  return a[0] * b[1] - a[1] * b[0];
}

inline void average(const vector2& a, const vector2& b, vector2& out)
{
  out[0] = 0.5 * (a[0] + b[0]);
  out[1] = 0.5 * (a[1] + b[1]);
}

inline void average(const vector2& a, const vector2& b, const vector2& c, vector2& out)
{
  out[0] = 1.0 / 3.0 * (a[0] + b[0] + c[0]);
  out[1] = 1.0 / 3.0 * (a[1] + b[1] + c[1]);
}

inline void normalize(vector2& a)
{
  scalar length = std::sqrt(dot(a, a));
  scale(a, 1.0 / length);
}

template<typename ScalarType, std::size_t M, std::size_t N>
inline void matrix_vector_mult(const matrix<ScalarType, M, N>& mat,
                        const vector<ScalarType, N>& x,
                        vector<ScalarType,M>& out)
{
  for(std::size_t i = 0; i < M; ++i)
  {
    out[i] = 0.0;
    for(std::size_t j = 0; j < N; ++j)
    {
      out[i] += mat[i][j] * x[j];
    }
  }
}


template<typename T>
class multi_array
{
public:
  typedef std::size_t size_type;
  typedef std::size_t index_type;

  typedef T* iterator;
  typedef const T* const_iterator;

  /// Construct empty array with no dimensions

  multi_array()
  : data_(nullptr), buffer_size_(0)
  {
  }

  /// Construct multi dimensional array with the dimensions given
  /// as a \c std::vector.
  /// @param sizes Specifies the dimensions. \c sizes.size() is the number
  /// of dimensions, and \c sizes[i] the extent in the i-th dimension.
  /// E.g, to construct a 2x3 array, \c sizes has to contain the elements {2, 3}.

  explicit multi_array(const std::vector<size_type>& sizes)
  : sizes_(sizes), data_(nullptr)
  {
    init();
  }

  /// Construct multi dimensional array with the dimensions given
  /// as a simple stack-based C-style array.
  /// @param sizes Specifies the dimensions. \c sizes.size() is the number
  /// of dimensions, and \c sizes[i] the extent in the i-th dimension.
  /// E.g, to construct a 2x3 array, \c sizes has to contain the elements {2, 3}.

  template<size_type N>
  explicit multi_array(const size_type(&sizes) [N])
  : data_(nullptr), sizes_(sizes + 0, sizes + N)
  {
    init();
  }

  /// Construct multi dimensional array with the dimensions given
  /// as a simple C-style array.
  /// @param sizes Specifies the dimensions. \c sizes.size() is the number
  /// of dimensions, and \c sizes[i] the extent in the i-th dimension.
  /// E.g, to construct a 2x3 array, \c sizes has to contain the elements {2, 3}.
  /// @param num_dimensions The number of elements of \c sizes and thus
  /// the number of dimensions of the \c multi_array

  multi_array(const size_type* sizes, size_type num_dimensions)
  : data_(nullptr), sizes_(sizes, sizes + num_dimensions)
  {
    init();
  }

  /// Construct two dimensional array
  /// @param size_x The extent of the array in dimension 0
  /// @param size_y The extent of the array in dimension 1

  multi_array(size_type size_x, size_type size_y)
  : data_(nullptr)
  {
    sizes_.reserve(2);
    sizes_.push_back(size_x);
    sizes_.push_back(size_y);
    init();
  }

  /// Construct three dimensional array
  /// @param size_x The extent of the array in dimension 0
  /// @param size_y The extent of the array in dimension 1
  /// @param size_z The extent of the array in dimension 2

  multi_array(size_type size_x, size_type size_y, size_type size_z)
  : data_(nullptr)
  {
    sizes_.reserve(3);
    sizes_.push_back(size_x);
    sizes_.push_back(size_y);
    sizes_.push_back(size_z);
    init();
  }

  /// Copy Constructor. May Throw.

  multi_array(const multi_array<T>& other)
  : sizes_(other.sizes_), buffer_size_(other.buffer_size_), data_(nullptr),
  position_increments_(other.position_increments_)
  {
    if (sizes_.size() != 0)
    {
      init();
      std::copy(other.begin(), other.end(), data_);
    }
  }

  /// Assignment operator. Provides strong exception guarantee.

  multi_array<T>& operator=(multi_array<T> other)
  {
    // using the copy and swap idiom (note the by-value function
    // parameter!) grants us a strong exception guarantee
    swap(*this, other);
    return *this;
  }

  ~multi_array()
  {
    if (data_)
      delete [] data_;
  }

  /// Swap two multi arrays. Their sizes do not have to equal.
  /// @param a The first array
  /// @param b The second array

  friend void swap(multi_array<T>& a, multi_array<T>& b)
  {
    using std::swap;

    swap(a.sizes_, b.sizes_);
    swap(a.buffer_size_, b.buffer_size_);
    swap(a.data_, b.data_);
    swap(a.position_increments_, b.position_increments_);
  }

  /// @return The extent of a dimension
  /// @param dim The index of the dimension

  size_type get_extent_of_dimension(std::size_t dim) const
  {
    assert(dim < sizes_.size());
    return sizes_[dim];
  }

  /// @return The dimension of the array

  size_type get_dimension() const
  {
    return sizes_.size();
  }

  /// @return An iterator type to the beginning of the array

  iterator begin()
  {
    return data_;
  }

  /// @return An iterator type to the beginning of the array

  const_iterator begin() const
  {
    return data_;
  }

  /// @return An iterator type pointing to one element beyond the
  /// last element of the array

  iterator end()
  {
    return data_ + buffer_size_;
  }

  /// @return An iterator type pointing to one element beyond the
  /// last element of the array

  const_iterator end() const
  {
    return data_ + buffer_size_;
  }

  /// @return The total number of elements in the array

  size_type get_num_elements() const
  {
    return buffer_size_;
  }

  /// @return The total number of elements in the array

  size_type size() const
  {
    return get_num_elements();
  }

  /// Grants access to the raw data buffer
  /// @return the raw data buffer
  T* data()
  {
    return data_;
  }

  /// Grants access to the raw data buffer
  /// @return the raw data buffer
  const T* data() const
  {
    return data_;
  }

  /// Checks if an (multi-dimensional) index is within the bounds
  // of the multi_array
  /// @return whether the index is within the bounds
  /// @param position the index
  bool is_within_bounds(const index_type* position) const
  {
    for(size_type i = 0; i < get_dimension(); ++i)
      if(position[i] >= sizes_[i])
        return false;
    return true;
  }

  /// Checks if an (multi-dimensional) index is within the bounds
  // of the multi_array
  /// @return whether the index is within the bounds
  /// @param position the index
  bool is_within_bounds(const std::vector<index_type>& position) const
  {
    for(size_type i = 0; i < get_dimension(); ++i)
      if(position[i] >= sizes_[i])
        return false;
    return true;
  }



  /// Access an element of the array
  /// @param position Contains the indices of the element to look up
  /// @return A reference to the specified element

  T& operator[](const std::vector<index_type>& position)
  {
    assert(position.size() == get_dimension());
    assert(data_ != nullptr);

    size_type pos = calculate_position(position);
    assert(pos < buffer_size_);
    return data_[pos];
  }

  /// Access an element of the array
  /// @param position Contains the indices of the element to look up
  /// @return A reference to the specified element

  const T& operator[](const std::vector<index_type>& position) const
  {
    assert(position.size() == get_dimension());
    assert(data_ != nullptr);

    size_type pos = calculate_position(position);
    assert(pos < buffer_size_);
    return data_[pos];
  }

  /// Access an element of the array
  /// @param position Contains the indices of the element to look up
  /// @return A reference to the specified element

  T& operator[](const std::vector<int>& position)
  {
    assert(position.size() == get_dimension());
    assert(data_ != nullptr);

    size_type pos = calculate_position(position);
    assert(pos < buffer_size_);
    return data_[pos];
  }

  /// Access an element of the array
  /// @param position Contains the indices of the element to look up
  /// @return A reference to the specified element

  const T& operator[](const std::vector<int>& position) const
  {
    assert(position.size() == get_dimension());
    assert(data_ != nullptr);

    size_type pos = calculate_position(position);
    assert(pos < buffer_size_);
    return data_[pos];
  }

  /// Access an element of the array
  /// @param position A simple C-array containing the indices of the
  /// element to look up
  /// @return A reference to the specified element

  template<size_type N>
  T& operator[](const index_type(&position) [N])
  {
    assert(N == get_dimension());
    assert(data_ != nullptr);

    size_type pos = calculate_position(position);
    assert(pos < buffer_size_);
    return data_[pos];
  }

  /// Access an element of the array
  /// @param position A simple C-array containing the indices of the
  /// element to look up
  /// @return A reference to the specified element

  template<size_type N>
  const T& operator[](const index_type(&position) [N]) const
  {
    assert(N == get_dimension());
    assert(data_ != nullptr);

    size_type pos = calculate_position(position);
    assert(pos < buffer_size_);
    return data_[pos];
  }

  /// Access an element of the array
  /// @param position A pointer to a simple C-array containing the
  /// indices of the element to look up
  /// @return A reference to the specified element

  T& operator[](const index_type* position)
  {
    assert(data_ != nullptr);

    size_type pos = calculate_position(position);
    if(pos >= buffer_size_)
    {
      std::cout << pos << " " << buffer_size_ << " " << position[0] << " " << position[1] << std::endl;
      std::cout << this->sizes_[0] << " " << this->sizes_[1] << std::endl;
    }

    assert(pos < buffer_size_);
    return data_[pos];
  }

  /// Access an element of the array
  /// @param position A pointer to asimple C-array containing the
  /// indices of the element to look up
  /// @return A reference to the specified element

  const T& operator[](const index_type* position) const
  {
    assert(data_ != nullptr);

    size_type pos = calculate_position(position);
    assert(pos < buffer_size_);
    return data_[pos];
  }

#ifndef WITHOUT_MPI
  /// Reduces each element within all parallel instances of the array on the
  /// master process by applying a user specified combination function.
  /// @param comm The mpi communicator
  /// @param master_rank The rank of the master process. Only this process
  /// will have access to the reduces array.
  /// @param combination_method The function that will be applied to combine
  /// all elements at the same position in the parallel array instances.
  /// It must have the signature T (const T&, const T&)

  template<class Function>
  void reduce_parallel_array(const boost::mpi::communicator& comm,
                             int master_rank,
                             Function combination_method)
  {
    int rank = comm.rank();
    int nprocs = comm.size();

    if(rank != master_rank)
    {
      // Use synchronous communication to save memory (the data tables can be quite large)
      comm.send(master_rank, 0, *this);
    }
    else
    {
      for(int proc = 0; proc < nprocs; ++proc)
      {
        if(proc != master_rank)
        {
          util::multi_array<T> recv_data;
          comm.recv(proc, 0, recv_data);

          assert(recv_data.get_dimension() == this->get_dimension());
          for(std::size_t dim = 0; dim < recv_data.get_dimension(); ++dim)
            assert(recv_data.get_extent_of_dimension(dim) == this->get_extent_of_dimension(dim));

          auto recv_element = recv_data.begin();
          auto data_element = this->begin();

          for(; recv_element != recv_data.end(); ++recv_element, ++data_element)
            combination_method(*data_element, *recv_element);

        }
      }
    }   
  }

  /// Like \c reduce_parallel_array, but makes the result available to each
  /// process.
  /// @param comm The mpi communicator
  /// @param combination_method The function that will be applied to combine
  /// all elements at the same position in the parallel array instances.
  /// It must have the signature T (const T&, const T&)
  template<class Function>
  void allreduce_parallel_array(const boost::mpi::communicator& comm,
                             Function combination_method)
  {
    int master_rank = 0;

    reduce_parallel_array(comm, master_rank, combination_method);

    boost::mpi::broadcast(comm, *this, master_rank);
  }

#endif

private:
  /// Based on the indices of an element, calculates the position
  /// of the element in the flat, one-dimensional data array.
  /// @param position An object of a type offering operator[], that
  /// stores the indices of the element
  /// @return the position in the one dimensional data array.

  template<typename Container>
  inline size_type calculate_position(const Container& position) const
  {
    size_type pos = 0;
    for (size_type i = 0; i < get_dimension(); ++i)
      pos += position[i] * position_increments_[i];

    return pos;
  }

  /// Initializes the data array and the position increments of each dimension

  void init()
  {
    assert(sizes_.size() != 0);
    for (std::size_t i = 0; i < sizes_.size(); ++i)
      assert(sizes_[i] != 0);

    if(data_ != nullptr)
      delete [] data_;

    buffer_size_ = std::accumulate(sizes_.begin(), sizes_.end(), 1,
                                   std::multiplies<size_type>());

    data_ = new T [buffer_size_];

    position_increments_.resize(get_dimension());
    position_increments_[0] = 1;

    for (std::size_t i = 1; i < sizes_.size(); ++i)
    {
      position_increments_[i] = position_increments_[i - 1] * sizes_[i - 1];
    }
  }

  std::vector<size_type> sizes_;
  std::vector<index_type> position_increments_;

  size_type buffer_size_;

  T* data_;

  friend class boost::serialization::access;

  template<class Archive>
  void save(Archive& ar, const unsigned int version) const
  {
    ar & sizes_;
    for(size_t i = 0; i < buffer_size_; ++i)
      ar & data_[i];
  }

  template<class Archive>
  void load(Archive& ar, const unsigned int version)
  {
    ar & sizes_;
    init();
    for(size_t i = 0; i < buffer_size_; ++i)
      ar & data_[i];
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
};


/// This class implements an ostream that will only write to the stream
/// on the master level. Calls from other processes will be ignored.
/// E.g. the following code will print "Hello world" only the process 
/// of rank 0 in MPI_COMM_WORLD:
///
/// \code{.cpp}
/// master_ostream stream(std::cout, 0);
/// stream << "Hello World!\n";
/// \endcode

class master_ostream
{
public:
  /// Initializes the object. Collective on MPI_COMM_WORLD.
  /// @param ostr The output stream that shall be used. The reference
  /// must remain valid throughout the existence of the master_ostream
  /// object.
  /// @param master_rank The rank of the process (in MPI_COMM_WORLD) 
  /// on which data shall be written.

  master_ostream(std::ostream& ostr, int master_rank)
  : ostr_(ostr), master_rank_(master_rank)
  {
    rank_ = _comm.rank();
  }

  /// Conversion to std::ostream&

  operator std::ostream& ()
  {
    return ostr_;
  }

  operator const std::ostream& () const
  {
    return ostr_;
  }

  /// Return the master rank.

  int get_master_rank() const
  {
    return master_rank_;
  }

  /// Get the process rank

  int get_rank() const
  {
    return rank_;
  }

  /// @return whether the calling process is the master process

  inline bool is_master_process() const
  {
    return rank_ == master_rank_;
  }

  typedef std::ostream& (*io_manip_type)(std::ostream&);
private:
  boost::mpi::communicator _comm;
  std::ostream& ostr_;
  int master_rank_;
  int rank_;
};

/// Implements \c operator<<.
/// Only calls from the process specified as \c master_process during
/// the construction of the \c master_ostream object will lead to output
/// being written to the output stream. Calls from other processes will
/// be ignored. This call is not collective, but obviously at least
/// the specified master process has to call the function if any
/// output is to be written at all.

template<class T>
master_ostream& operator<<(master_ostream& ostr, const T& msg)
{
  if (ostr.get_rank() == ostr.get_master_rank())
    (std::ostream&)ostr << msg;

  return ostr;
}

/// This version enables the use of io manips

master_ostream& operator<<(master_ostream& ostr,
        master_ostream::io_manip_type io_manip)
{
  if (ostr.get_rank() == ostr.get_master_rank())
    io_manip((std::ostream&)ostr);

  return ostr;
}

template<class T>
class cached_value
{
  bool _value_present;
  T _value;
public:
  cached_value()
  : _value_present(false) {}

  template<class GeneratorFunction>
  inline const T& retrieve(GeneratorFunction f)
  {
    if(_value_present)
      return _value;

    _value = f();
    _value_present = true;
    return _value;
  }

};  


template<typename T, std::size_t N>
std::vector<T> array_to_vector(const std::array<T,N>& data)
{
  return std::vector<T>(data.begin(), data.end());
}


}
}

#endif	/* UTIL_HPP */

