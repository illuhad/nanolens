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

#ifndef GEOMETRY_HPP
#define	GEOMETRY_HPP

#include <cmath>
#include "util.hpp"
#include "mesh.hpp"

namespace nanolens{
namespace geometry{



class triangle
{
public:
  triangle(const util::vector2& point_a,
           const util::vector2& point_b,
           const util::vector2& point_c)
  : _points{point_a, point_b, point_c}
  {}
  
  util::scalar area() const
  {
    return std::abs(signed_area());
  }
  
  inline util::scalar signed_area() const
  {
    return 0.5 * ((_points[1][0] - _points[0][0]) * (_points[2][1] - _points[0][1]) 
                 - (_points[2][0] - _points[0][0]) * (_points[1][1] - _points[0][1]));   
  }
  
  bool contains_point(const util::vector2& point) const
  {
    util::vector2 v0 = _points[2];
    util::vector2 v1 = _points[1];
    util::vector2 v2 = point;
    util::sub(v0, _points[0]);
    util::sub(v1, _points[0]);
    util::sub(v2, _points[0]);

    util::scalar dot00, dot01, dot02, dot11, dot12;
    dot00 = util::dot(v0, v0);
    dot01 = util::dot(v0, v1);
    dot02 = util::dot(v0, v2);
    dot11 = util::dot(v1, v1);
    dot12 = util::dot(v1, v2);

    util::scalar factor = 1. / (dot00 * dot11 - dot01 * dot01);

    // transform to u,v coordinates
    util::scalar u = factor * (dot11 * dot02 - dot01 * dot12);
    util::scalar v = factor * (dot00 * dot12 - dot01 * dot02);

    return u > 0. && v > 0. && ((u + v) < 1.);
  }
  
  virtual std::size_t num_vertices() const
  {
    return 3;
  }
  
  
  const util::vector2& operator[](std::size_t index) const
  {
    assert(index < 3);
    return _points[index];
  }

  util::vector2& operator[](std::size_t index)
  {
    assert(index < 3);
    return _points[index];
  }
  

  
private:
  std::array<util::vector2, 3> _points;
  util::vector2 _a;  
  util::vector2 _b;
  util::vector2 _c;
};

template<std::size_t N>
class polygon
{
public:
  static constexpr std::size_t num_hull_vertices()
  { return N; }
  
  static constexpr std::size_t num_vertices()
  { return N + 1; }
  
  virtual ~polygon(){}
                 
  util::scalar area() const
  {
    util::scalar result = 0.0;
    
    for(std::size_t i = 0; i < _triangles.size(); ++i)
    {
      triangle tri(_vertices[_triangles[i][0]],
                   _vertices[_triangles[i][1]],
                   _vertices[_triangles[i][2]]);
      
      result += tri.area();
    }
  }
  
  bool contains_point(const util::vector2& point) const
  {
    for(std::size_t i = 0; i < _triangles.size(); ++i)
    {
      triangle tri(_vertices[_triangles[i][0]],
                   _vertices[_triangles[i][1]],
                   _vertices[_triangles[i][2]]);
      
      if(tri.contains_point(point))
        return true;
    }
    
    return false;
  }
  
  std::size_t num_vertices() const
  {
    return _vertices.size();
  }
  
  std::size_t num_triangles() const
  {
    return _triangles.size();
  }
  
  const std::array<std::size_t, 3>& get_triangle_vertex_indices(std::size_t triangle_index) const
  {
    assert(triangle_index < _triangles.size());
    return _triangles[triangle_index];
  }
  
  const util::vector2& operator[](std::size_t index) const
  {
    assert(index < _vertices.size());
    return _points[index];
  }

  util::vector2& operator[](std::size_t index)
  {
    assert(index < _vertices.size());
    return _points[index];
  }

  
  void shift_coordinates(const util::vector2& offset)
  {
    for(util::vector2& vertex : _vertices)
      util::add(vertex, offset);
  }
  
  const util::vector2& get_center() const
  {
    return _vertices[N];
  }
  
protected:
  // There are N vertices on the sides of the polygon and one at the center
  std::array<util::vector2, N + 1> _vertices;
  
  std::array<std::array<std::size_t, 3>, N> _triangles;
};

template<std::size_t N>
class equilateral_polygon : public polygon<N>
{
public:
  
  equilateral_polygon(const util::vector2& center,
                      util::scalar radius)
  {
    _vertices[N] = center;
    
    util::scalar pi = boost::math::constants::pi<util::scalar>();
    
    for(std::size_t i = 0; i < N; ++i)
    {
      util::scalar current_angle = i * 2.0 * pi / static_cast<util::scalar>(N);
      
      _vertices[i] = {radius * std::cos(current_angle), radius * std::sin(current_angle)};
    }
    
    for(std::size_t i = 0; i < N - 1; ++i)
    {
      // define triangle between the current vertex on the hull i, the
      // adjecent vertex at i + 1 and the central vertex at N
      _triangles[i] = {i, i + 1, N};
    }
    // special treatment for the last triangle
    _triangles[N - 1] = {N - 1, 0, N};
  }
  
  virtual ~equilateral_polygon(){}
private:
  // There are N vertices on the sides of the polygon and one at the center
  std::array<util::vector2, N + 1> _vertices;
  
  std::array<std::array<std::size_t, 3>, N> _triangles;
};

}
}

#endif	/* GEOMETRY_HPP */

