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


/// Implements a triangle
class triangle
{
public:
  /// Construct object
  /// \param point_a The first vertex
  /// \param point_b The second vertex
  /// \param point_c The third vertex
  triangle(const util::vector2& point_a,
           const util::vector2& point_b,
           const util::vector2& point_c)
  : _points{point_a, point_b, point_c}
  {}
  
  /// \return The area of the triangle
  util::scalar area() const
  {
    return std::abs(signed_area());
  }
  
  /// \return The area of the triangle. The sign of the area will be positive
  /// the vertices are oriented counter-clockwise and negative if they
  /// are oriented clockwise
  inline util::scalar signed_area() const
  {
    return 0.5 * ((_points[1][0] - _points[0][0]) * (_points[2][1] - _points[0][1]) 
                 - (_points[2][0] - _points[0][0]) * (_points[1][1] - _points[0][1]));   
  }
  
  /// \return whether a given point is inside the triangle
  /// \param point The point that shall be checked
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
  
  /// \return The number of vertices (always 3 for a triangle)
  virtual std::size_t num_vertices() const
  {
    return 3;
  }
  
  /// Access a vertex by its index. Valid indices are in \c [0, num_vertices-1]
  /// \param index The vertex index
  /// \return The position of the vertex with the given index
  const util::vector2& operator[](std::size_t index) const
  {
    assert(index < 3);
    return _points[index];
  }

  /// Access a vertex by its index. Valid indices are in \c [0, num_vertices-1]
  /// \param index The vertex index
  /// \return The position of the vertex with the given index
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

/// Implements a polygon by decomposing it into triangles. Decomposition is done
/// by adding a vertex at the center of the polygon and connecting every two adjacent
/// vertices on the hull of the polygon with the central vertex.
/// \tparam N The number of (hull) vertices of the polygon
template<std::size_t N>
class polygon
{
public:
  /// \return The number of vertices on the hull of the polygon
  static constexpr std::size_t num_hull_vertices()
  { return N; }
  
  /// \return The total number of vertices used, i.e. the number of hull vertices
  /// plus one central vertex.
  static constexpr std::size_t num_vertices()
  { return N + 1; }
  
  /// Construct object
  polygon()
  {
    // Initialize default triangulation
        
    for(std::size_t i = 0; i < N - 1; ++i)
    {
      // define triangle between the current vertex on the hull i, the
      // adjecent vertex at i + 1 and the central vertex at N
      this->_triangles[i] = {i, i + 1, N};
    }
    // special treatment for the last triangle
    this->_triangles[N - 1] = {N - 1, 0, N};
  }
  
  virtual ~polygon(){}
                 
  /// Calculates the area of the polygon by adding the areas of all its triangles
  /// \return The polygon area
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
    
    return result;
  }
  
  /// \return whether a given point is inside the polygon
  /// \param point The point that shall be checked
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
  
  /// \return The number of triangles that make up the polygon
  std::size_t num_triangles() const
  {
    return _triangles.size();
  }
  
  /// Gets the indices of the vertices of one of the triangles that make up
  /// the polygon
  /// \param triangle_index The index of the triangle. Valid indices are in
  /// \c [0,num_triangles()-1]
  /// \return An array with the indices of the three vertices that define
  /// the specified triangle.
  const std::array<std::size_t, 3>& get_triangle_vertex_indices(std::size_t triangle_index) const
  {
    assert(triangle_index < _triangles.size());
    return _triangles[triangle_index];
  }
  
  /// Accesses the coordinates of one of the vertices
  /// \param index The vertex index
  /// \return The coordinates of the specified vertex
  const util::vector2& operator[](std::size_t index) const
  {
    assert(index < _vertices.size());
    return _vertices[index];
  }

  /// Accesses the coordinates of one of the vertices
  /// \param index The vertex index
  /// \return The coordinates of the specified vertex
  util::vector2& operator[](std::size_t index)
  {
    assert(index < _vertices.size());
    return _vertices[index];
  }

  /// Shifts the coordinates of all vertices by a given offset vector
  /// \param offset A vector that will be added to all vertex coordinates
  void shift_coordinates(const util::vector2& offset)
  {
    for(util::vector2& vertex : _vertices)
      util::add(vertex, offset);
  }
  
  /// \return The coordinates of the central vertex
  const util::vector2& get_center() const
  {
    return _vertices[N];
  }
  
protected:
  // There are N vertices on the sides of the polygon and one at the center
  std::array<util::vector2, N + 1> _vertices;
  
  std::array<std::array<std::size_t, 3>, N> _triangles;
};

/// Implements an equilateral polygon
/// \tparam N The number of vertices of the polygon
template<std::size_t N>
class equilateral_polygon : public polygon<N>
{
public:
  
  /// Construct object
  /// \param center The coordinates of the center of the equilateral polygon
  /// \param radius The radius of the polygon, i.e. the distance of all its
  /// vertices from the center of the polygon.
  equilateral_polygon(const util::vector2& center,
                      util::scalar radius)
  {
    this->_vertices[N] = center;
    
    util::scalar pi = boost::math::constants::pi<util::scalar>();
    
    for(std::size_t i = 0; i < N; ++i)
    {
      util::scalar current_angle = i * 2.0 * pi / static_cast<util::scalar>(N);
      
      this->_vertices[i] = {radius * std::cos(current_angle), radius * std::sin(current_angle)};
    }

  }
  
  virtual ~equilateral_polygon(){}

};

}
}

#endif	/* GEOMETRY_HPP */

