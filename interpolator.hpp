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
  typedef std::array<Coordinate_type, 2> vector_type;
  typedef std::array<T, 16> eval_values_storage_type;
  
  bicubic_interpolator() = default;
  bicubic_interpolator(vector_type min_extent, vector_type max_extent)
  : _min_extent(min_extent)
  {
    _area_size = max_extent;
    for(std::size_t i = 0; i < 2; ++i)
      _area_size[i] -= min_extent[i];
    
    for(std::size_t i = 0; i < 4; ++i)
      for(std::size_t j = 0; j < 4; ++j)
        _bicubic_interpolation_coefficients[i][j] = T();
  }

  T interpolate(const vector_type& v) const
  {
    vector_type relative_pos = v;
    
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
    for(std::size_t i = 0; i < 4; ++i)
      for(std::size_t j = 0; j < 4; ++j)
      {
        vector_type position = _min_extent;
        //position[0] += (static_cast<int>(i) - 1) * _area_size[0];
        //position[1] += (static_cast<int>(j) - 1) * _area_size[1];
        
        position[0] += static_cast<int>(i) * _area_size[0] / 3.0;
        position[1] += static_cast<int>(j) * _area_size[1] / 3.0;
        fn(position, eval_point_id);
        
        ++eval_point_id;
      }
  }
  
  
  template<class Eval_function>
  void init(Eval_function fn)
  {
    util::matrix_nxn<T, 4> p;
     
    for(std::size_t i = 0; i < 4; ++i)
      for(std::size_t j = 0; j < 4; ++j)
      {
        vector_type position = _min_extent;
        //position[0] += (static_cast<int>(i) - 1) * _area_size[0];
        //position[1] += (static_cast<int>(j) - 1) * _area_size[1];
        
        position[0] += static_cast<int>(i) * _area_size[0] / 3.0;
        position[1] += static_cast<int>(j) * _area_size[1] / 3.0;
        p[i][j] = fn(position);
      }

    init(p);
  }
  
  inline void init(const eval_values_storage_type& eval)
  {
    util::matrix_nxn<T, 4> eval_matrix;
    
    unsigned array_pos = 0;
    for(unsigned i=0; i < 4; ++i)
      for(unsigned j=0; j < 4; ++j)
      {
        eval_matrix[i][j] = eval[array_pos];
        ++array_pos;
      }
    
    init(eval_matrix);
  }
  
  inline void init(const util::matrix_nxn<T, 4>& p)
  {
    util::matrix_nxn<T, 4>& a = this->_bicubic_interpolation_coefficients;
    
    a[0][0] = p[1][1];
    a[0][1]= -.5*p[1][0] + .5*p[1][2];
    a[0][2] = p[1][0] - 2.5*p[1][1] + 2*p[1][2] - .5*p[1][3];
    a[0][3] = -.5*p[1][0] + 1.5*p[1][1] - 1.5*p[1][2] + .5*p[1][3];
    a[1][0] = -.5*p[0][1] + .5*p[2][1];
    a[1][1] = .25*p[0][0] - .25*p[0][2] - .25*p[2][0] + .25*p[2][2];
    a[1][2] = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + .5*p[2][0] - 1.25*p[2][1] + p[2][2] - .25*p[2][3];
    a[1][3] = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .25*p[2][0] + .75*p[2][1] - .75*p[2][2] + .25*p[2][3];
    a[2][0] = p[0][1] - 2.5*p[1][1] + 2*p[2][1] - .5*p[3][1];
    a[2][1] = -.5*p[0][0] + .5*p[0][2] + 1.25*p[1][0] - 1.25*p[1][2] - p[2][0] + p[2][2] + .25*p[3][0] - .25*p[3][2];
    a[2][2] = p[0][0] - 2.5*p[0][1] + 2*p[0][2] - .5*p[0][3] - 2.5*p[1][0] + 6.25*p[1][1] - 5*p[1][2] + 1.25*p[1][3] + 2*p[2][0] - 5*p[2][1] + 4*p[2][2] - p[2][3] - .5*p[3][0] + 1.25*p[3][1] - p[3][2] + .25*p[3][3];
    a[2][3] = -.5*p[0][0] + 1.5*p[0][1] - 1.5*p[0][2] + .5*p[0][3] + 1.25*p[1][0] - 3.75*p[1][1] + 3.75*p[1][2] - 1.25*p[1][3] - p[2][0] + 3*p[2][1] - 3*p[2][2] + p[2][3] + .25*p[3][0] - .75*p[3][1] + .75*p[3][2] - .25*p[3][3];
    a[3][0] = -.5*p[0][1] + 1.5*p[1][1] - 1.5*p[2][1] + .5*p[3][1];
    a[3][1] = .25*p[0][0] - .25*p[0][2] - .75*p[1][0] + .75*p[1][2] + .75*p[2][0] - .75*p[2][2] - .25*p[3][0] + .25*p[3][2];
    a[3][2] = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + 1.5*p[1][0] - 3.75*p[1][1] + 3*p[1][2] - .75*p[1][3] - 1.5*p[2][0] + 3.75*p[2][1] - 3*p[2][2] + .75*p[2][3] + .5*p[3][0] - 1.25*p[3][1] + p[3][2] - .25*p[3][3];
    a[3][3] = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .75*p[1][0] + 2.25*p[1][1] - 2.25*p[1][2] + .75*p[1][3] + .75*p[2][0] - 2.25*p[2][1] + 2.25*p[2][2] - .75*p[2][3] - .25*p[3][0] + .75*p[3][1] - .75*p[3][2] + .25*p[3][3];

  }
private:
  
  inline T interpolate_at_relative(const vector_type& relative_position) const
  {
    vector_type transformed_rel_pos;
    
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
  
  vector_type _area_size;
  vector_type _min_extent;
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
  
  typedef typename interpolator_type::vector_type position_vector_type;
  
  vector_interpolator() = default;
  vector_interpolator(const position_vector_type& min_extent,
                      const position_vector_type& max_extent)
  {
    for(unsigned i = 0; i < Dimension; ++i)
      _interpolators[i] = interpolator_type(min_extent, max_extent);
  }
  
  template<class Eval_function>
  void init(Eval_function fn)
  {
    std::array<typename interpolator_type::eval_values_storage_type, Dimension> eval;
    
    _interpolators[0].foreach_evaluation_point(
    [&eval, &fn](const position_vector_type& pos, unsigned eval_point_id)
    {
      auto function_value = fn(pos);
      
      for(std::size_t i = 0; i < Dimension; ++i)
        eval[i][eval_point_id] = function_value[i];
    });
    
    for(std::size_t i = 0; i < Dimension; ++i)
    {
      _interpolators[i].init(eval[i]);
    }
  }
  
  Vector_type interpolate(const position_vector_type& v) const
  {
    Vector_type result;
    
    for(unsigned i = 0; i < Dimension; ++i)
      result[i] = _interpolators[i].interpolate(v);
    
    return result;
  }
  
private:
  std::array<interpolator_type, Dimension> _interpolators;
  
};

typedef vector_interpolator<util::scalar, util::scalar, 2, util::vector2, bicubic_interpolator>
standard_vector2_interpolator;

}

#endif	/* INTERPOLATOR_HPP */

