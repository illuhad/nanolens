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

#ifndef TREE_HPP
#define	TREE_HPP

#include "util.hpp"
#include "grid.hpp"
#include "star.hpp"
#include "numeric.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/index/rtree.hpp>

namespace nanolens{

//class nested_interpolation_grid 
//{
//
//  class grid_layer
//  {
//  public:
//    typedef std::shared_ptr<grid_layer> layer_ptr_type;
//    
//    struct grid_entry
//    {
//      // The external deflections at the four corners are stored in a matrix
//      // where the 0,0 entry corresponds to the contributions from the upper left
//      // corner, 0,1 to the contributions from the upper right corner and so on
//      util::matrix_nxn<util::vector2, 4> external_deflections;
//      util::matrix_nxn<util::vector2, 2> corner_positions;
//      
//      layer_ptr_type sub_grid;
//      
//      // Interpolation coefficients are vectors because we need to interpolate
//      // both x and y components
//      util::matrix_nxn<util::vector2, 4> bicubic_interpolation_coefficients;
//      
//      std::vector<star> interior_stars;
//      // holds a list of the exterior stars. This vector is no longer needed
//      // once the interpolation coefficients have been determined and will
//      // then be emptied.
//      std::vector<star> exterior_stars;
//      
//      grid_entry()
//      {
//        for(std::size_t i = 0; i < external_deflections.size(); ++i)
//          for(std::size_t j = 0; j < external_deflections[i].size(); ++j)
//            external_deflections[i][j] = {0., 0.};
//      }
//    };
//    
//    typedef util::grid2d<util::scalar, grid_entry> grid_type;
//    
//    grid_layer() = default;
//    
//    grid_layer(std::size_t layer_depth, 
//               std::size_t max_layer_depth,
//               std::size_t num_allowed_stars_per_bucket,
//               const util::vector2& min_extent,
//               const util::vector2& max_extent,
//               const std::array<std::size_t, 2>& num_cells,
//               const std::array<std::size_t, 2>& num_subgrid_cells,
//               util::scalar accuracy = 1.e-8)
//    : _layer_idx(layer_depth), _max_layer_idx(max_layer_depth - 1),
//      _num_cells(num_cells), _num_subgrid_cells(num_subgrid_cells),
//      _num_allowed_stars(num_allowed_stars_per_bucket),
//      _accuracy(accuracy)
//    {
//      
//      _grid = grid_type(min_extent, max_extent, num_cells);
//      
//      grid_type::scalar_array_type bucket_size = _grid.get_bucket_size();
//      
//      for(std::size_t grid_i = 0; grid_i < _grid.get_total_num_buckets(0); ++grid_i)
//      {
//        for(std::size_t grid_j = 0; grid_j < _grid.get_total_num_buckets(1); ++grid_j)
//        {
//          grid_type::index_type bucket_index = {grid_i, grid_j};
//          
//          grid_entry& grid_element = _grid[bucket_index];
//          
//          util::vector2 min_corner_position = _grid.get_min_position_of_bucket(bucket_index);
//          for(std::size_t i = 0; i < grid_element.corner_positions.size(); ++i)
//          {
//            for(std::size_t j = 0; j < grid_element.corner_positions[i].size(); ++j)
//            { 
//              grid_element.corner_positions[i][j] = min_corner_position;
//              grid_element.corner_positions[i][j][0] += i * bucket_size[0];
//              grid_element.corner_positions[i][j][1] += j * bucket_size[1];
//            }
//          }
//          
//        }
//      }
//    }
//    
//    util::vector2 get_deflection(const util::vector2& position) const
//    {
//      util::vector2 result = {0.0, 0.0};
//         
//      grid_type::index_type bucket = _grid.get_index(position);
//      
//      //Check if we are a leaf in the tree structure
//      if(_grid[bucket].sub_grid == nullptr)
//        result = exact_deflection(bucket, position);
//      else       
//        result = _grid[bucket].sub_grid->get_deflection(position);
//      
//      
//      util::vector2 interpolation_contribution = interpolate(bucket, position);
//      util::add(result, interpolation_contribution);
//      
//      //std::cout << _layer_idx << "  " << result[0] << " " << result[1] << std::endl;
//      
//      return result;
//    }
//    
//    void add_star(const star& s)
//    {
//      grid_type::index_type bucket = _grid.get_index(s.get_position());
//      grid_entry& entry = _grid[bucket];
//      
//      entry.interior_stars.push_back(s);
//      
//      // Adapt stored exterior deflections 
//      grid_type::scalar_array_type bucket_size = _grid.get_bucket_size();
//      grid_type::scalar_array_type step_size = bucket_size;
//      step_size[0] /= static_cast<util::scalar>(entry.external_deflections.size()-1);
//      step_size[1] /= static_cast<util::scalar>(entry.external_deflections.size()-1);
//      
//      for(std::size_t i = 0; i < _grid.get_total_num_buckets(0); ++i)
//      {
//        for(std::size_t j = 0; j < _grid.get_total_num_buckets(1); ++j)
//        {
//          grid_type::index_type idx = {i, j};
//          
//          if(idx != bucket)
//          {
//            grid_entry& current_entry = _grid[idx];
//          
//            // Add the star's deflection to the interpolation values
//            for(std::size_t interpol_i = 0;
//              interpol_i < current_entry.external_deflections.size();
//              ++interpol_i)
//            {
//              for(std::size_t interpol_j = 0;
//                interpol_j < current_entry.external_deflections[interpol_i].size();
//                ++interpol_j)
//              {
//                util::vector2 evaluation_position 
//                  = current_entry.corner_positions[0][0];
//                util::vector2 correction = {interpol_i * step_size[0],
//                                            interpol_j * step_size[1]};
//                util::add(evaluation_position, correction);
//                
//                util::vector2 deflection;
//                s.calculate_deflection_angle(evaluation_position, deflection);
//                
//                util::add(current_entry.external_deflections[interpol_i][interpol_j], deflection);
//              }
//            }
//            
//            current_entry.exterior_stars.push_back(s);
//          }
//        }
//      }
//      
//      if(entry.interior_stars.size() > _num_allowed_stars)
//      {  
//        // Attempt to refine grid
//        if(_layer_idx < _max_layer_idx)
//        {
//          // If we have just exceeded the maximum number of stars,
//          // all stars need to be redistributed into the refined grid
//          if(entry.interior_stars.size() == _num_allowed_stars + 1)
//            for(const star& interior_star : entry.interior_stars)
//              insert_star_into_bucket(bucket, interior_star);
//          else
//            // Otherwise we only need to take care of the new star
//            insert_star_into_bucket(bucket, s);
//        }
//      }
//    }
//    
//    void init_interpolation()
//    {
//      for(std::size_t i = 0; i < _grid.get_total_num_buckets(0); ++i)
//      {
//        for(std::size_t j = 0; j < _grid.get_total_num_buckets(1); ++j)
//        {
//          grid_type::index_type bucket_index = {i, j};
//          grid_entry& bucket = _grid[bucket_index];
//          
//          update_bicubic_interpolation_coefficients(bucket);
//          
//          if(bucket.sub_grid != nullptr)
//            bucket.sub_grid->init_interpolation();
//        }
//      }
//    }
//  private:
//    inline void insert_star_into_bucket(const grid_type::index_type& bucket_for_star,
//                                        const star& s)
//    {
//      grid_type::scalar_array_type bucket_size = _grid.get_bucket_size();
//      
//      util::vector2 bucket_min_corner = _grid.get_min_position_of_bucket(bucket_for_star);
//      util::vector2 bucket_max_corner = bucket_min_corner;
//      util::add(bucket_max_corner, bucket_size);
//
//      //Create new subgrid if necessary
//      if(!_grid[bucket_for_star].sub_grid)
//        _grid[bucket_for_star].sub_grid = layer_ptr_type(new grid_layer(_layer_idx + 1, 
//                                                   _max_layer_idx + 1, 
//                                                   _num_allowed_stars, 
//                                                   bucket_min_corner,
//                                                   bucket_max_corner,
//                                                   _num_subgrid_cells,
//                                                   _num_subgrid_cells));
//
//      _grid[bucket_for_star].sub_grid->add_star(s);
// 
//    }
//    
//    inline util::vector2 exact_deflection(const std::vector<star>& star_list,
//                                          const util::vector2& position) const
//    {
//      util::vector2 result = {0.0, 0.0};
//      
//      for(const auto& star : star_list)
//      {
//        util::vector2 contribution;
//        star.calculate_deflection_angle(position, contribution);
//        util::add(result, contribution);
//      }
//      
//      return result;
//    }
//    
//    /// @return exact deflection of interior stars
//    inline util::vector2 exact_deflection(const grid_type::index_type& bucket,
//                                          const util::vector2& position) const
//    {
//      return exact_deflection(_grid[bucket].interior_stars, position);
//    }
//    
//    inline util::vector2 bicubic_interpolation(const grid_entry& bucket,
//                                               const util::vector2& relative_position) const
//    {
//      util::scalar x = relative_position[0];
//      util::scalar y = relative_position[1];
//    
//      std::array<util::scalar, 4> x_powers = {1., x, x, x};
//      std::array<util::scalar, 4> y_powers = {1., y, y, y};
//      
//      for(std::size_t i = 2; i < 4; ++i)
//      {
//        x_powers[i] *= x_powers[i - 1];
//        y_powers[i] *= y_powers[i - 1];
//      }
//     
//      util::vector2 result = {0.0, 0.0};
//      
//      for(std::size_t i = 0; i < 4; ++i)
//      {
//        for(std::size_t j = 0; j < 4; ++j)
//        {
//          util::vector2 contribution = bucket.bicubic_interpolation_coefficients[i][j];
//          util::scale(contribution, x_powers[i] * y_powers[j]);
//          util::add(result, contribution);
//        }
//      }
//      return result;
//    }
//    
//    void update_bicubic_interpolation_coefficients(grid_entry& bucket)
//    {
//      for(std::size_t component = 0; component < 2; ++component)
//      {
//        util::matrix_nxn<util::vector2, 4>& a = bucket.bicubic_interpolation_coefficients;
//        util::matrix_nxn<util::scalar, 4> p;
//        for(std::size_t i = 0; i < 4; ++i)
//          for(std::size_t j = 0; j < 4; ++j)
//            p[i][j] = bucket.external_deflections[i][j][component];
//        
//        a[0][0][component] = p[1][1];
//        a[0][1][component] = -.5*p[1][0] + .5*p[1][2];
//        a[0][2][component] = p[1][0] - 2.5*p[1][1] + 2*p[1][2] - .5*p[1][3];
//        a[0][3][component] = -.5*p[1][0] + 1.5*p[1][1] - 1.5*p[1][2] + .5*p[1][3];
//        a[1][0][component] = -.5*p[0][1] + .5*p[2][1];
//        a[1][1][component] = .25*p[0][0] - .25*p[0][2] - .25*p[2][0] + .25*p[2][2];
//        a[1][2][component] = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + .5*p[2][0] - 1.25*p[2][1] + p[2][2] - .25*p[2][3];
//        a[1][3][component] = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .25*p[2][0] + .75*p[2][1] - .75*p[2][2] + .25*p[2][3];
//        a[2][0][component] = p[0][1] - 2.5*p[1][1] + 2*p[2][1] - .5*p[3][1];
//        a[2][1][component] = -.5*p[0][0] + .5*p[0][2] + 1.25*p[1][0] - 1.25*p[1][2] - p[2][0] + p[2][2] + .25*p[3][0] - .25*p[3][2];
//        a[2][2][component] = p[0][0] - 2.5*p[0][1] + 2*p[0][2] - .5*p[0][3] - 2.5*p[1][0] + 6.25*p[1][1] - 5*p[1][2] + 1.25*p[1][3] + 2*p[2][0] - 5*p[2][1] + 4*p[2][2] - p[2][3] - .5*p[3][0] + 1.25*p[3][1] - p[3][2] + .25*p[3][3];
//        a[2][3][component] = -.5*p[0][0] + 1.5*p[0][1] - 1.5*p[0][2] + .5*p[0][3] + 1.25*p[1][0] - 3.75*p[1][1] + 3.75*p[1][2] - 1.25*p[1][3] - p[2][0] + 3*p[2][1] - 3*p[2][2] + p[2][3] + .25*p[3][0] - .75*p[3][1] + .75*p[3][2] - .25*p[3][3];
//        a[3][0][component] = -.5*p[0][1] + 1.5*p[1][1] - 1.5*p[2][1] + .5*p[3][1];
//        a[3][1][component] = .25*p[0][0] - .25*p[0][2] - .75*p[1][0] + .75*p[1][2] + .75*p[2][0] - .75*p[2][2] - .25*p[3][0] + .25*p[3][2];
//        a[3][2][component] = -.5*p[0][0] + 1.25*p[0][1] - p[0][2] + .25*p[0][3] + 1.5*p[1][0] - 3.75*p[1][1] + 3*p[1][2] - .75*p[1][3] - 1.5*p[2][0] + 3.75*p[2][1] - 3*p[2][2] + .75*p[2][3] + .5*p[3][0] - 1.25*p[3][1] + p[3][2] - .25*p[3][3];
//        a[3][3][component] = .25*p[0][0] - .75*p[0][1] + .75*p[0][2] - .25*p[0][3] - .75*p[1][0] + 2.25*p[1][1] - 2.25*p[1][2] + .75*p[1][3] + .75*p[2][0] - 2.25*p[2][1] + 2.25*p[2][2] - .75*p[2][3] - .25*p[3][0] + .75*p[3][1] - .75*p[3][2] + .25*p[3][3];
//      }
//      
//    }
//    
//    inline util::vector2 bilinear_interpolation(const grid_entry& bucket,
//                                                const util::vector2& relative_position) const
//    { 
//      util::scalar x = relative_position[0];
//      util::scalar y = relative_position[1];
//      
//      util::vector2 result = bucket.external_deflections[0][0];
//      util::scale(result, (1.-x)*(1.-y));
//      
//      util::scale_add(result, bucket.external_deflections[1][0], x * (1. - y));
//      util::scale_add(result, bucket.external_deflections[0][1], y * (1. - x));
//      util::scale_add(result, bucket.external_deflections[1][1], x * y);
//      
//      return result;
//    }
//    
//    inline util::vector2 interpolate(const grid_type::index_type& bucket_index,
//                                     const util::vector2& position) const
//    {
//      const grid_entry& bucket = _grid[bucket_index];
//      
//      util::vector2 relative_position = position;
//      util::sub(relative_position, _grid.get_min_position_of_bucket(bucket_index));
//      
//      grid_type::scalar_array_type bucket_size = _grid.get_bucket_size();
//      
//      for(std::size_t i = 0; i < relative_position.size(); ++i)
//        relative_position[i] /= bucket_size[i];
//      
//      return bilinear_interpolation(bucket, relative_position);
//    }
//    
//    std::array<std::size_t, 2> _num_cells;
//    std::array<std::size_t, 2> _num_subgrid_cells;
//    
//    grid_type _grid;
//    
//    
//    std::size_t _layer_idx;
//    std::size_t _max_layer_idx;
//    std::size_t _num_allowed_stars;
//    
//    util::scalar _accuracy;
//  };
//public:
//  
//  nested_interpolation_grid() = default;
//  
//  nested_interpolation_grid(const std::vector<star>& stars,
//                            std::size_t max_num_layers,
//                            const util::vector2& min_extent,
//                            const util::vector2& max_extent,
//                            const std::array<std::size_t, 2>& num_toplevel_cells,
//                            const std::array<std::size_t, 2>& num_subgrid_cells)
//  : _top_layer(0, max_num_layers, 1, min_extent, max_extent, num_toplevel_cells, num_subgrid_cells)
//  {
//    assert(max_num_layers > 0);
//    
//    for(const star& s : stars)
//      _top_layer.add_star(s);
//    
//    _top_layer.init_interpolation();
//    
//    std::ofstream ostr("defl.dat", std::ios::trunc);
//    for(util::scalar x = min_extent[0]; x < max_extent[0]; x+= 0.1)
//      for(util::scalar y = min_extent[1]; y < max_extent[1]; y+=0.1)
//      {
//        util::vector2 defl = get_deflection({x,y});
//        ostr << x << " " << y << " " << (util::square(defl[0]) + util::square(defl[1])) << "\n";
//      }
//    ostr.close();
//    
//  }
//  
//  util::vector2 get_deflection(const util::vector2& pos) const
//  {
//    //std::cout << "------------------\n";
//    return _top_layer.get_deflection(pos);
//  }
//  
//  
//  
//private:
//  
//  grid_layer _top_layer;
//};  
//  
//class deflection_grid
//{
//public:
//  class indexed_star
//  {
//  public:
//    indexed_star() = default;
//    
//    indexed_star(const star& s, std::size_t id)
//    : _star(s), _id(id)
//    {}
//    
//    const star& get_star() const
//    { return _star; }
//    
//    std::size_t get_id() const
//    { return _id; }
//  private:
//    star _star;
//    std::size_t _id;
//  };
//
//  typedef std::vector<pseudo_star> grid_entry_type;
//  
//  typedef util::grid2d<util::scalar, grid_entry_type> grid_type;
//  
//  typedef boost::geometry::model::point<util::scalar, 2, boost::geometry::cs::cartesian> point;
//  typedef std::pair<point, indexed_star> tree_value_type;
//  typedef boost::geometry::index::rtree<tree_value_type, boost::geometry::index::rstar<16> > tree_type;
//  
//  deflection_grid() = default;
//  
//  
//  deflection_grid(const std::vector<star>& stars,
//                  const util::vector2& min_extent,
//                  const util::vector2& max_extent,
//                  const std::array<std::size_t, 2>& num_cells,
//                  util::scalar max_deflection_error)
//  : _grid(min_extent, max_extent, num_cells),
//    _deflection_tolerance(max_deflection_error)
//  {
//    // Create grid stars
//    std::vector<indexed_star> indexed_stars;
//    for(std::size_t i = 0; i < stars.size(); ++i)
//      indexed_stars.push_back(indexed_star(stars[i], i));
//    
//    std::vector<std::pair<point, indexed_star>> tree_data;
//    for(auto& s : indexed_stars)
//    {
//      point pt(s.get_star().get_position()[0],
//               s.get_star().get_position()[1]);
//      
//      tree_data.push_back(std::make_pair(pt, s));
//    }
//    
//    tree_type tree(tree_data);
//    
//    for(std::size_t i = 0; i < _grid.get_total_num_buckets(0); ++i)
//      for(std::size_t j = 0; j < _grid.get_total_num_buckets(1); ++j)
//      {
//        grid_type::index_type cell = {i, j};
//        init_cell(tree, cell, indexed_stars);
//      }
//    
//    std::size_t n_pseudo_stars = 0;
//    util::scalar ncells = 0.;
//    for(std::size_t i = 0; i < _grid.get_total_num_buckets(0); ++i)
//      for(std::size_t j = 0; j < _grid.get_total_num_buckets(1); ++j)
//      {
//        grid_type::index_type cell = {i, j};
//        n_pseudo_stars += _grid[cell].size();
//        ncells += 1.;
//      }
//    std::cout << "avg stars per cell "<< (double)n_pseudo_stars/ncells << std::endl;
//    
//  }
//  
//  util::vector2 get_deflection(const util::vector2& pos) const
//  {
//    const grid_entry_type& cell_elements = _grid[pos];
//    
//    util::vector2 result = {0.0, 0.0};
//    
//    for(auto& ps : cell_elements)
//    {
//      util::vector2 deflection;
//      ps.calculate_deflection_angle(pos, deflection);
//      util::add(result, deflection);
//    }
//    
//    return result;
//  }
//    
//private:
//  void init_cell(const tree_type& tree,
//                 const grid_type::index_type& cell_index,
//                 const std::vector<indexed_star>& stars)
//  {
//    std::vector<bool> is_star_processed(stars.size(), false);
//    
//    grid_entry_type& cell = _grid[cell_index];
//    util::vector2 cell_position = _grid.get_min_position_of_bucket(cell_index);
//    util::vector2 cell_size = _grid.get_bucket_size();
//    util::scale_add(cell_position, cell_size, 0.5);
//
//    point pt(cell_position[0], cell_position[1]);
//    
//    cell.clear();
//
//    std::vector<indexed_star> sorted_stars;
//    for ( tree_type::const_query_iterator it 
//      = tree.qbegin(boost::geometry::index::nearest(pt, stars.size())) ;
//          it != tree.qend() ; ++it )
//    {
//      sorted_stars.push_back(it->second);
//    }
//    
//    // Start from the farthest star
//    for(int i = sorted_stars.size() - 1; i >= 0; --i)
//    {
//      if(!is_star_processed[sorted_stars[i].get_id()])
//      {
//        pseudo_star combined_star;
//        group_stars(sorted_stars[i], cell_position, tree, is_star_processed, combined_star);
//        cell.push_back(combined_star);
//      }
//    }
//
//  }
//  
//  void group_stars(const indexed_star& grouping_center,
//                   const util::vector2& reference_position,
//                   const tree_type& tree,
//                   std::vector<bool>& processed_stars,
//                   pseudo_star& out) const
//  {
//    point pt(grouping_center.get_star().get_position()[0], 
//             grouping_center.get_star().get_position()[1]);
//    
//    util::vector2 center_of_mass_sum = grouping_center.get_star().get_position();
//    
//    util::scalar masses = grouping_center.get_star().get_mass();
//    util::scale(center_of_mass_sum, masses);
//    
//    util::vector2 exact_deflection;
//    grouping_center.get_star().calculate_deflection_angle(reference_position,
//                                                          exact_deflection);
//
//    std::vector<star> grouped_stars;
//    grouped_stars.push_back(grouping_center.get_star());
//    
//    for ( tree_type::const_query_iterator it 
//      = tree.qbegin(boost::geometry::index::nearest(pt, processed_stars.size())) ;
//          it != tree.qend() ; ++it )
//    {
//      // Don't count the central star twice
//      if(it->second.get_id() != grouping_center.get_id() && !processed_stars[it->second.get_id()])
//      {
//        util::vector2 new_cm_sum = center_of_mass_sum;
//        util::scale_add(new_cm_sum, it->second.get_star().get_position(),
//                        it->second.get_star().get_mass());
//        util::scalar new_total_mass = masses + it->second.get_star().get_mass();
//        
//        util::vector2 new_cm = new_cm_sum;
//        util::scale(new_cm, 1.0 / new_total_mass);
//        
//        star new_pseudo_star(new_cm, new_total_mass);
//        
//        util::vector2 new_deflection;
//        new_pseudo_star.calculate_deflection_angle(reference_position, new_deflection);
//        
//        util::vector2 individual_deflection;
//        it->second.get_star().calculate_deflection_angle(reference_position, individual_deflection);
//        util::vector2 new_exact_deflection = exact_deflection;
//        util::add(new_exact_deflection, individual_deflection);
//        
//        if(std::abs(new_exact_deflection[0] - new_deflection[0])
//          +std::abs(new_exact_deflection[1] - new_deflection[1]) < _deflection_tolerance)
//        {
//          center_of_mass_sum = new_cm_sum;
//          masses = new_total_mass;
//          exact_deflection = new_exact_deflection;
//          processed_stars[it->second.get_id()] = true;
//          grouped_stars.push_back(it->second.get_star());
//        }
//        else
//          break;
//        
//      }
//    }
//    
//    util::vector2 center_of_mass = center_of_mass_sum;
//    util::scale(center_of_mass, 1.0 / masses);
//    out = pseudo_star(grouped_stars);
//  }
//
//
//  util::scalar _deflection_tolerance;
//  grid_type _grid;
//};

namespace impl_{

template<unsigned Multipole_order>
class barnes_hut_grid_cell
{
public:
  static const std::size_t max_allowed_stars_per_cell = 1;
  
  typedef util::grid<util::scalar, barnes_hut_grid_cell<Multipole_order>, 2 , false> grid_layer_type;
  
  barnes_hut_grid_cell()
  : _min_extent({0.0,0.0}), _max_extent({0.0,0.0}) {}

  barnes_hut_grid_cell(const util::vector2& min_extent,
            const util::vector2& max_extent)
  {
    set_cell_size(min_extent, max_extent);
  }

  inline void set_cell_size(const util::vector2& min_extent,
                            const util::vector2& max_extent)
  {
    _min_extent = min_extent;
    _max_extent = max_extent;
    
    _center = min_extent;
    util::add(_center, max_extent);
    util::scale(_center, 0.5);
    
    _diameter = 0.0;
    for(std::size_t i = 0; i < max_extent.size(); ++i)
      _diameter += max_extent[i] - min_extent[i];
    
    _diameter /= static_cast<util::scalar>(max_extent.size());
  }

  inline const util::vector2& get_cell_center() const
  { return _center; }
  
  inline util::scalar get_cell_diameter() const
  { return _diameter; }
  
  inline bool is_leaf() const
  { return _sub_grid == nullptr; }
  
  inline const grid_layer_type& get_sub_cells() const
  {
    assert(_sub_grid != nullptr);
    return *_sub_grid;
  }
  
  const std::vector<star>& get_interior_stars() const
  {
    return _interior_stars;
  }
  
  const pseudo_star<Multipole_order>& get_pseudo_star() const
  {
    return _approx_pseudo_star;
  }
  
  const util::vector2& get_min_extent() const
  {
    return _min_extent;
  }
  
  const util::vector2& get_max_extent() const
  {
    return _max_extent;
  }
  
  inline bool cell_contains_point(const util::vector2& pos)
  {
    for(std::size_t i = 0; i < Dimension; ++i)
      if(!(_min_extent[i] <= pos[i] && pos[i] < _max_extent[i]))
        return false;
    return true;
  }

  void add_star(const star& s)
  {
    assert(cell_contains_point(s.get_position()));
    
    _interior_stars.push_back(s);

    if(_interior_stars.size() > max_allowed_stars_per_cell)
    {
      // Subdivide cell
      if(_sub_grid == nullptr)
      {
        std::array<std::size_t, 2> num_cells = {2, 2};

        _sub_grid = std::shared_ptr<grid_layer_type>(new grid_layer_type(_min_extent,
                                                                         _max_extent,
                                                                         num_cells));
        for(std::size_t x = 0; x < Dimension; ++x)
          for(std::size_t y = 0; y < Dimension; ++y)
          {
            // Set cell positions and sizes
            typename grid_layer_type::index_type grid_idx = {x,y};

            util::vector2 min_extent = _sub_grid->get_min_position_of_bucket(grid_idx);
            util::vector2 bucket_size = _sub_grid->get_bucket_size();

            util::vector2 max_extent = min_extent;
            util::add(max_extent, bucket_size);

            (*_sub_grid)[grid_idx].set_cell_size(min_extent, max_extent);
          }

        // Redistribute the one existing star
        for(std::size_t i  = 0; i < max_allowed_stars_per_cell; ++i)
          insert_star_into_subcell(_interior_stars[i]);

      }

      insert_star_into_subcell(s);
    }
  }

  void init_pseudo_stars()
  {
    if(!_interior_stars.empty())
    {
      _approx_pseudo_star = pseudo_star<Multipole_order>(this->_interior_stars);
      _center = _approx_pseudo_star.center_of_mass();
    }

    if(_sub_grid != nullptr)
    {
      for(auto sub_cell = _sub_grid->begin();
          sub_cell != _sub_grid->end();
          ++sub_cell)
      {
        sub_cell->init_pseudo_stars();
      }
    }
  }
  
private:
  util::vector2 _min_extent;
  util::vector2 _max_extent;
  util::scalar _diameter;
  util::vector2 _center;
  
  static const std::size_t Dimension = 2;

  std::shared_ptr<grid_layer_type> _sub_grid;
  std::vector<star> _interior_stars;
  pseudo_star<Multipole_order> _approx_pseudo_star;


  inline void insert_star_into_subcell(const star& s)
  {
    assert(_sub_grid != nullptr);

    const util::vector2& pos = s.get_position();      

    for(auto sub_cell = _sub_grid->begin();
        sub_cell != _sub_grid->end();
        ++sub_cell)
    {
      if(sub_cell->cell_contains_point(pos))
      {
        sub_cell->add_star(s);
        return;
      }
    }
  }

};
}

class barnes_hut_tree
{
public:
  typedef impl_::barnes_hut_grid_cell<3> grid_cell;
  
  barnes_hut_tree() = default;
  
  barnes_hut_tree(const std::vector<star>& stars,
             util::scalar accuracy)
  : _accuracy_squared(util::square(accuracy))
  {
    if(!stars.empty())
    {
      util::scalar scalar_max = std::numeric_limits<util::scalar>::max();
      util::scalar scalar_min = std::numeric_limits<util::scalar>::min();

      util::vector2 min_extent = {scalar_max, scalar_max};
      util::vector2 max_extent = {scalar_min, scalar_min};
      for(const star& s : stars)
      {
        const util::vector2& pos = s.get_position();
        for(std::size_t i = 0; i < 2; ++i)
        {
          if(pos[i] < min_extent[0])
          {
            min_extent[0] = pos[i]-0.01;
            min_extent[1] = pos[i]-0.01;
          }
          if(pos[i] > max_extent[0])
          {
            max_extent[0] = pos[i]+0.01;
            max_extent[1] = pos[i]+0.01;
          }
        }
      }
      
      _tree_root.set_cell_size(min_extent, max_extent);
      
      for(const star& s : stars)
        _tree_root.add_star(s);
      
      _tree_root.init_pseudo_stars();
    }
  }
  
  util::vector2 get_deflection(const util::vector2& pos) const
  {
    //std::cerr << "----------------------\n";
    
    util::vector2 defl  = get_deflection_of_cell(this->_tree_root, pos);
    
    //std::cerr << defl[0] << " " << defl[1] << std::endl;
    return defl;
  }
  
private:
  util::vector2 get_deflection_of_cell(const grid_cell& cell, const util::vector2& pos) const
  {
   // std::cerr << "cell_min = " << cell.get_min_extent()[0] << " "<<cell.get_min_extent()[1] <<
   //              " cell_max = " << cell.get_max_extent()[0] << " " << cell.get_max_extent()[1] <<
   //              " nstars = "  << cell.get_interior_stars().size() << std::endl;
    
    util::vector2 result = {0.0, 0.0};
    
    if(cell.is_leaf())
    {
      for(const star& s : cell.get_interior_stars())
      {
        assert(cell.get_interior_stars().size() <= grid_cell::max_allowed_stars_per_cell);
        
        util::vector2 contribution;
        s.calculate_deflection_angle(pos, contribution);
        util::add(result, contribution);
      }
    }
    else
    {
      // Check if we can use the approximation
      util::vector2 diff = pos;
      util::sub(diff, cell.get_cell_center());
      
      util::scalar distance_squared = util::dot(diff, diff);
      
      if(util::square(cell.get_cell_diameter()) / distance_squared < _accuracy_squared)
      {
        // Get approximated deflection
        cell.get_pseudo_star().calculate_deflection_angle(pos, result);
      }
      else
      {
        // Try the children
        auto& sub_cells = cell.get_sub_cells();
        
        std::size_t n = 0;
        for(auto current_cell = sub_cells.begin(); 
          current_cell != sub_cells.end(); 
          ++current_cell, ++n)
        {
          //std::cerr << "subcell " << n << ": ";
          util::vector2 contribution = get_deflection_of_cell(*current_cell, pos);
          util::add(result, contribution);
        }
        //std::cerr << "\n";
      }
    }
    
    return result;
  }
  
  grid_cell _tree_root;
  util::scalar _accuracy_squared;
};




}

#endif

