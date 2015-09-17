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

#ifndef INTERPOLATOR_HPP
#define	INTERPOLATOR_HPP

#include "util.hpp"


namespace nanolens{

template<typename T, typename Coordinate_type>
class bicubic_interpolator
{
public:
  typedef std::array<Coordinate_type, 2> position_vector_type;
  typedef std::array<T, 16> eval_values_storage_type;
  typedef util::matrix_nxn<T, 4> eval_values_matrix_type;
  
  static constexpr unsigned num_samples_x = 4;
  static constexpr unsigned num_samples_y = 4;
  
  bicubic_interpolator() = default;
  bicubic_interpolator(position_vector_type min_extent, position_vector_type max_extent)
  : _min_extent(min_extent)
  {
    _area_size = max_extent;
    for(std::size_t i = 0; i < 2; ++i)
      _area_size[i] -= min_extent[i];
    
    for(std::size_t i = 0; i < 4; ++i)
      for(std::size_t j = 0; j < 4; ++j)
        _bicubic_interpolation_coefficients[i][j] = T();
  }

  T interpolate(const position_vector_type& v) const
  {
    position_vector_type relative_pos = v;
    
    for(std::size_t i=0; i < 2; ++i)
    {
      relative_pos[i] -= _min_extent[i];
      relative_pos[i] /= _area_size[i];
    }
    
    return interpolate_at_relative(relative_pos);
  }
  
  template<class Function>
  void foreach_evaluation_point(Function fn) const
  {
    unsigned eval_point_id = 0;
    for(unsigned i = 0; i < 4; ++i)
      for(unsigned j = 0; j < 4; ++j)
      {
        position_vector_type position = _min_extent;
        //position[0] += (static_cast<int>(i) - 1) * _area_size[0];
        //position[1] += (static_cast<int>(j) - 1) * _area_size[1];
        
        position[0] += static_cast<int>(i) * _area_size[0] / 3.0;
        position[1] += static_cast<int>(j) * _area_size[1] / 3.0;
        fn(position, eval_point_id, i, j);
        
        ++eval_point_id;
      }
  }
  
  
  template<class Eval_function>
  void init(Eval_function fn)
  {
    eval_values_matrix_type p;
     
    for(std::size_t i = 0; i < 4; ++i)
      for(std::size_t j = 0; j < 4; ++j)
      {
        position_vector_type position = _min_extent;

        position[0] += static_cast<int>(i) * _area_size[0] / 3.0;
        position[1] += static_cast<int>(j) * _area_size[1] / 3.0;
        p[i][j] = fn(position);
      }

    init(p);
  }
  
  inline void init(const eval_values_storage_type& eval)
  {
    eval_values_matrix_type eval_matrix;
    
    unsigned array_pos = 0;
    for(unsigned i=0; i < 4; ++i)
      for(unsigned j=0; j < 4; ++j)
      {
        eval_matrix[i][j] = eval[array_pos];
        ++array_pos;
      }
    
    init(eval_matrix);
  }
  
  
  /// Initializes the interpolator object by calculating the bicubic
  /// interpolation coefficients.
  /// @param p A 4x4 matrix containing the function values at the evaluation
  /// positions.
  inline void init(const eval_values_matrix_type& p)
  {
    util::matrix_nxn<T, 4>& a = this->_bicubic_interpolation_coefficients;
    
    a[0][0] =        p[1][1];
    a[0][1] = -0.5 * p[1][0] + 0.5 * p[1][2];
    a[0][2] =        p[1][0] - 2.5 * p[1][1] + 2   * p[1][2] - 0.5 * p[1][3];
    a[0][3] = -0.5 * p[1][0] + 1.5 * p[1][1] - 1.5 * p[1][2] + 0.5 * p[1][3];
    a[1][0] = -0.5 * p[0][1] + 0.5 * p[2][1];
    a[1][1] = 0.25 * p[0][0] - 0.25* p[0][2] - 0.25* p[2][0] + 0.25* p[2][2];
    a[1][2] = -0.5 * p[0][0] + 1.25* p[0][1] -       p[0][2] + 0.25* p[0][3] 
              +0.5 * p[2][0] - 1.25* p[2][1] +       p[2][2] - 0.25* p[2][3];
    a[1][3] = 0.25 * p[0][0] - 0.75* p[0][1] + 0.75* p[0][2] - 0.25* p[0][3] 
             -0.25 * p[2][0] + 0.75* p[2][1] - 0.75* p[2][2] + 0.25* p[2][3];
    a[2][0] =        p[0][1] - 2.5 * p[1][1] + 2   * p[2][1] - 0.5 * p[3][1];
    a[2][1] = -0.5 * p[0][0] + 0.5 * p[0][2] + 1.25* p[1][0] - 1.25* p[1][2]
                   - p[2][0] +       p[2][2] + 0.25* p[3][0] - 0.25* p[3][2];
    a[2][2] =        p[0][0] - 2.5 * p[0][1] + 2   * p[0][2] - 0.5 * p[0][3] 
              -2.5 * p[1][0] + 6.25* p[1][1] - 5   * p[1][2] + 1.25* p[1][3] 
              +2   * p[2][0] - 5   * p[2][1] + 4   * p[2][2] -       p[2][3] 
              -0.5 * p[3][0] + 1.25* p[3][1] -       p[3][2] + 0.25* p[3][3];
    a[2][3] = -0.5 * p[0][0] + 1.5 * p[0][1] - 1.5 * p[0][2] + 0.5 * p[0][3] 
              +1.25* p[1][0] - 3.75* p[1][1] + 3.75* p[1][2] - 1.25* p[1][3] 
              -      p[2][0] + 3   * p[2][1] - 3   * p[2][2] +       p[2][3] 
              +0.25* p[3][0] - 0.75* p[3][1] + 0.75* p[3][2] - 0.25* p[3][3];
    a[3][0] = -0.5 * p[0][1] + 1.5 * p[1][1] - 1.5 * p[2][1] + 0.5 * p[3][1];
    a[3][1] =  0.25* p[0][0] - 0.25* p[0][2] - 0.75* p[1][0] + 0.75* p[1][2] 
              +0.75* p[2][0] - 0.75* p[2][2] - 0.25* p[3][0] + 0.25* p[3][2];
    a[3][2] = -0.5 * p[0][0] + 1.25* p[0][1] -       p[0][2] + 0.25* p[0][3] 
              +1.5 * p[1][0] - 3.75* p[1][1] + 3   * p[1][2] - 0.75* p[1][3] 
              -1.5 * p[2][0] + 3.75* p[2][1] - 3   * p[2][2] + 0.75* p[2][3] 
              +0.5 * p[3][0] - 1.25* p[3][1] +       p[3][2] - 0.25* p[3][3];
    a[3][3] =  0.25* p[0][0] - 0.75* p[0][1] + 0.75* p[0][2] - 0.25* p[0][3] 
              -0.75* p[1][0] + 2.25* p[1][1] - 2.25* p[1][2] + 0.75* p[1][3] 
              +0.75* p[2][0] - 2.25* p[2][1] + 2.25* p[2][2] - 0.75* p[2][3] 
              -0.25* p[3][0] + 0.75* p[3][1] - 0.75* p[3][2] + 0.25* p[3][3];

  }
private:
  
  inline T interpolate_at_relative(const position_vector_type& relative_position) const
  {
    position_vector_type transformed_rel_pos;
    
    //Transform coordinates from relative position on the entire sampling area
    //to relative position on the central square of evaluation points
    //(neccessary due to the calculation formulas of the interpolation matrix entries)
    
    for(std::size_t i = 0; i < transformed_rel_pos.size(); ++i)
      transformed_rel_pos[i] = relative_position[i] * 3.0 - 1.0;
    
    Coordinate_type x = transformed_rel_pos[0];
    Coordinate_type y = transformed_rel_pos[1];

    std::array<Coordinate_type, 4> x_powers = {1., x, 0., 0.};
    std::array<Coordinate_type, 4> y_powers = {1., y, 0., 0.};

    for(std::size_t i = 2; i < 4; ++i)
    {
      x_powers[i] = x * x_powers[i - 1];
      y_powers[i] = y * y_powers[i - 1];
    }

    T result = T();

    for(std::size_t i = 0; i < 4; ++i)
      for(std::size_t j = 0; j < 4; ++j)
        result += x_powers[i] * y_powers[j] * _bicubic_interpolation_coefficients[i][j];
    
    return result;
  }
  
  position_vector_type _area_size;
  position_vector_type _min_extent;
  util::matrix_nxn<T, 4> _bicubic_interpolation_coefficients;
};

template<typename T,
         typename Coordinate_type,
         unsigned Dimension,
         class Vector_type,
         template<class, class> class Interpolator_type>
class vector_interpolator
{
public:
  
  static_assert(Dimension > 0, "Dimension must be greater than 0.");
  
  typedef Interpolator_type<T, Coordinate_type> interpolator_type;
  
  typedef typename interpolator_type::position_vector_type position_vector_type;
  
  static constexpr unsigned num_samples_x = interpolator_type::num_samples_x;
  static constexpr unsigned num_samples_y = interpolator_type::num_samples_y;
  
  vector_interpolator() = default;
  vector_interpolator(const position_vector_type& min_extent,
                      const position_vector_type& max_extent)
  {
    for(unsigned i = 0; i < Dimension; ++i)
      _interpolators[i] = interpolator_type(min_extent, max_extent);
  }
  
  typedef typename interpolator_type::eval_values_matrix_type backend_eval_values_matrix_type;
  typedef std::array
  <
    backend_eval_values_matrix_type, 
    Dimension
  > eval_values_matrix_type;
  
  template<class Eval_function>
  void init(Eval_function fn)
  {
    eval_values_matrix_type eval;
    
    _interpolators[0].foreach_evaluation_point(
    [&eval, &fn](const position_vector_type& pos, unsigned eval_point_id, unsigned i, unsigned j)
    {
      auto function_value = fn(pos);
      
      for(std::size_t k = 0; k < Dimension; ++k)
        eval[k][i][j] = function_value[k];
    });
    
    init(eval);
  }
  

  void init(const eval_values_matrix_type& eval)
  {
    for(std::size_t i = 0; i < Dimension; ++i)
    {
      _interpolators[i].init(eval[i]);
    }
  }
  
  static inline void set_eval_value(eval_values_matrix_type& m,
                                    unsigned i, 
                                    unsigned j, 
                                    const Vector_type& x)
  {
    for(unsigned k = 0; k < Dimension; ++k)
      m[k][i][j] = x[k];
  }
  
  static inline Vector_type get_eval_value(const eval_values_matrix_type& m,
                                           unsigned i,
                                           unsigned j)
  {
    Vector_type result;
    for(unsigned k = 0; k < Dimension; ++k)
      result[k] = m[k][i][j];
    return result;
  }
  
    
  template<class Function>
  void foreach_evaluation_point(Function fn) const
  {
    _interpolators[0].foreach_evaluation_point(fn);
  }

  Vector_type interpolate(const position_vector_type& v) const
  {
    Vector_type result;
    
    for(unsigned i = 0; i < Dimension; ++i)
      result[i] = _interpolators[i].interpolate(v);
    
    return result;
  }
  
  std::size_t get_dimension() const
  { return Dimension; }
  
  interpolator_type& operator[](std::size_t i)
  { return _interpolators[i]; }
  
  const interpolator_type& operator[](std::size_t i) const
  { return _interpolators[i]; }
  
private:
  std::array<interpolator_type, Dimension> _interpolators;
  
};

typedef vector_interpolator<util::scalar, util::scalar, 2, util::vector2, bicubic_interpolator>
standard_vector2_interpolator;

// Uses a taylor series to approximate the deflection angle. Expects Vector_type
// to be two-dimensional.
template<typename T, class Vector_type>
class deflection_angle_taylor_interpolator
{
public:
  typedef Vector_type position_vector_type;
  typedef std::array<Vector_type, 4> eval_values_storage_type;
  typedef util::matrix_nxn<Vector_type, 4> eval_values_matrix_type;
  
  static constexpr unsigned num_samples_x = 2;
  static constexpr unsigned num_samples_y = 2;
  
  deflection_angle_taylor_interpolator() = default;
  deflection_angle_taylor_interpolator(const position_vector_type& min_extent,
                                       const position_vector_type& max_extent)
  : _min_corner(min_extent)
  {
    // Currently, only square cells are supported (the assertion had to be 
    // commented out because it may trigger due to rounding errors, especially
    // when using float)
    // assert((max_extent[0] - min_extent[0]) == (max_extent[1] - min_extent[0]));

    _delta = max_extent[0] - min_extent[0];
    _center = _min_corner;
    for(std::size_t i = 0; i < _center.size(); ++i)
      _center[i] += 0.5 * _delta;
  }
  
  deflection_angle_taylor_interpolator(const position_vector_type& center,
                                       T cell_size)
  : _center(center), _delta(cell_size)
  {
    for(std::size_t i = 0; i < _center.size(); ++i)
      _min_corner[i] = _center[i] - 0.5 * _delta;
  }
  
  template<class Eval_function>
  void init(Eval_function fn)
  {
    eval_values_storage_type corner_values;
    corner_values[0] = fn({_min_corner[0] + _delta, _min_corner[1] + _delta});
    corner_values[1] = fn({_min_corner[0],          _min_corner[1] + _delta});
    corner_values[2] = fn(_min_corner);
    corner_values[3] = fn({_min_corner[0] + _delta, _min_corner[1]});
    
    init(corner_values);
  }
  
  void init(const eval_values_matrix_type& corner_values)
  {
    eval_values_storage_type values;
    values[2] = corner_values[0][0];
    values[3] = corner_values[1][0];
    values[1] = corner_values[0][1];
    values[0] = corner_values[1][1];
    init(values);
  }
  
  template<class Function>
  void foreach_evaluation_point(Function fn) const
  {
    util::matrix_nxn<unsigned, 2> index_map;
    index_map[0][0] = 2;
    index_map[1][0] = 3;
    index_map[0][1] = 1;
    index_map[1][1] = 0;
      
    for(unsigned i = 0; i < 2; ++i)
      for(unsigned j = 0; j < 2; ++j)
      {
        position_vector_type position = _min_corner;
        
        position[0] += i * _delta;
        position[1] += j * _delta;
        fn(position, index_map[i][j], i, j);
      }
  }
  
  void init(const eval_values_storage_type& corner_values)
  {
    
    T delta_inverse = 1.0 / _delta;
    T delta_inverse2 = util::square(delta_inverse);
    T delta_inverse3 = delta_inverse2 * delta_inverse;
    
    // Average of the x components
    _coefficients[0] = 0.25 * (corner_values[0][0] + corner_values[1][0] 
                             + corner_values[2][0] + corner_values[3][0]);
   // Average of the y components
    _coefficients[1] = 0.25 * (corner_values[0][1] + corner_values[1][1] 
                             + corner_values[2][1] + corner_values[3][1]);
    
    // d/dx(alpha_x)
    _coefficients[2] = 0.125 * delta_inverse * (corner_values[0][0] - corner_values[1][0]
                                              + corner_values[3][0] - corner_values[2][0]
                                              + corner_values[2][1] - corner_values[1][1]
                                              + corner_values[3][1] - corner_values[0][1]);
    
    // d/dx (alpha_y) = d/dy (alpha_x)
    _coefficients[3] = 0.125 * delta_inverse * (corner_values[1][0] - corner_values[2][0]
                                              + corner_values[0][0] - corner_values[3][0]
                                              + corner_values[0][1] - corner_values[1][1]
                                              + corner_values[3][1] - corner_values[2][1]);
    
    // d/dy(d/dx(alpha_y))
    _coefficients[4] = 0.25 * delta_inverse2 * (corner_values[3][1] - corner_values[2][1]
                                             +  corner_values[1][1] - corner_values[0][1]);
    
    // d/dy(d/dx(alpha_x))
    _coefficients[5] = 0.25 * delta_inverse2 * (corner_values[0][0] - corner_values[1][0]
                                             +  corner_values[2][0] - corner_values[3][0]);
    
    _coefficients[6] = 0.375 * delta_inverse3 * (corner_values[1][0] - corner_values[0][0]
                                               + corner_values[2][0] - corner_values[3][0]
                                               + corner_values[3][1] - corner_values[0][1]
                                               + corner_values[2][1] - corner_values[1][1]);
    
    _coefficients[7] = 0.375 * delta_inverse3 * (corner_values[0][0] - corner_values[3][0]
                                               + corner_values[1][0] - corner_values[2][0]
                                               + corner_values[2][1] - corner_values[3][1]
                                               + corner_values[1][1] - corner_values[0][1]);
  }
  
  Vector_type interpolate(const position_vector_type& v) const
  {
    Vector_type result;
    
    position_vector_type delta = v;
    util::sub(delta, _center);
    
    position_vector_type delta2 = {util::square(delta[0]), util::square(delta[1])};
    
    std::array<T, 4> eval_polynomials;
    eval_polynomials[0] = 0.5 * (delta2[0] - delta2[1]);
    eval_polynomials[1] = delta[0] * delta[1];
    eval_polynomials[2] = 0.5 * delta[0] * (delta2[0] / 3.0 - delta2[1]);
    eval_polynomials[3] = 0.5 * delta[1] * (delta2[0] - delta2[1] / 3.0);
    
    result[0] = _coefficients[0] + delta[0] * _coefficients[2] + delta[1] * _coefficients[3]
      + eval_polynomials[0] * _coefficients[4] + eval_polynomials[1] * _coefficients[5]
      + eval_polynomials[2] * _coefficients[6] + eval_polynomials[3] * _coefficients[7];
    
    result[1] = _coefficients[1] - delta[1] * _coefficients[2] + delta[0] * _coefficients[3]
      - eval_polynomials[1] * _coefficients[4] + eval_polynomials[0] * _coefficients[5]
      - eval_polynomials[3] * _coefficients[6] + eval_polynomials[2] * _coefficients[7];
   
    
    return result;
  }
  
  static inline void set_eval_value(eval_values_matrix_type& m,
                                    unsigned i, 
                                    unsigned j, 
                                    const Vector_type& x)
  {
    m[i][j] = x;
  }
  
  static inline Vector_type get_eval_value(const eval_values_matrix_type& m,
                                           unsigned i,
                                           unsigned j)
  {
    return m[i][j];
  }
  
  
private:
  T _delta;
  std::array<T, 8> _coefficients;
  position_vector_type _min_corner;
  position_vector_type _center;
};


template<class Interpolator_type, std::size_t N>
struct interpolator_grid_bulk_initializer
{
  typedef typename Interpolator_type::position_vector_type position_vector_type;
  
  template<class Eval_function>
  inline void shared_border_init(
                            std::array<std::array<Interpolator_type,N>,N>& interpolators,
                            Eval_function fn)
  {
    std::array<std::array<typename Interpolator_type::eval_values_matrix_type, N>,N> cell_eval_values;
    
    unsigned num_samples_x = Interpolator_type::num_samples_x;
    unsigned num_samples_y = Interpolator_type::num_samples_y;
    
    assert(N > 0);

    for(std::size_t i = 0; i < N; ++i)
    {
      for(std::size_t j = 0; j < N; ++j)
      {
        
        interpolators[i][j].foreach_evaluation_point([&](const position_vector_type& pos,
                                                     unsigned eval_id, 
                                                     unsigned local_i, 
                                                     unsigned local_j)
        {
          if(local_i == 0 && i != 0)
          {
            auto border_value = Interpolator_type::get_eval_value(cell_eval_values[i-1][j],
                                                                  num_samples_x - 1, 
                                                                  local_j);
            
            Interpolator_type::set_eval_value(cell_eval_values[i][j],
                                              local_i, 
                                              local_j, 
                                              border_value);
          }
          else if(local_j == 0 && j != 0)
          {
            
            auto border_value = Interpolator_type::get_eval_value(cell_eval_values[i][j-1],
                                                      local_i, 
                                                      num_samples_y - 1);
            
            Interpolator_type::set_eval_value(cell_eval_values[i][j], 
                                              local_i, 
                                              local_j, 
                                              border_value);
          }
          else
            Interpolator_type::set_eval_value(cell_eval_values[i][j], 
                                              local_i, 
                                              local_j, 
                                              fn(pos));
          
        });

        interpolators[i][j].init(cell_eval_values[i][j]);
      }
    }
    
    
  }
};
  

}

#endif	/* INTERPOLATOR_HPP */

