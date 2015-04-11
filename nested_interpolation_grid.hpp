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

#ifndef NESTED_INTERPOLATION_GRID_HPP
#define	NESTED_INTERPOLATION_GRID_HPP

#include "util.hpp"
#include "grid.hpp"
#include "star.hpp"

namespace nanolens{

class nested_interpolation_grid 
{

  class grid_layer
  {
  public:
    typedef std::shared_ptr<grid_layer> layer_ptr_type;
    
    struct grid_entry
    {
      // The external deflections at the four corners are stored in a matrix
      // where the 0,0 entry corresponds to the contributions from the upper left
      // corner, 0,1 to the contributions from the upper right corner and so on
      util::matrix_nxn<util::vector2, 2> external_deflections;
      layer_ptr_type sub_grid;
      
      std::vector<star> interior_stars;
    };
    
    typedef util::grid2d<util::scalar, grid_entry> grid_type;
    
    grid_layer() = default;
    
    grid_layer(std::size_t layer_depth, 
               std::size_t max_layer_depth,
               std::size_t num_allowed_stars_per_bucket,
               const util::vector2& min_extent,
               const util::vector2& max_extent,
               const std::array<std::size_t, 2>& num_cells,
               const std::array<std::size_t, 2>& num_subgrid_cells)
    : _layer_idx(layer_depth), _max_layer_idx(max_layer_depth - 1),
      _num_cells(num_cells), _num_subgrid_cells(num_subgrid_cells),
      _num_allowed_stars(num_allowed_stars_per_bucket)
    {
      
      _grid = grid_type(min_extent, max_extent, num_cells);
      
      for(auto grid_element = _grid.begin(); grid_element != _grid.end(); ++grid_element)
      {
        for(std::size_t i = 0; i < grid_element->external_deflections.size(); ++i)
          for(std::size_t j = 0; j < grid_element->external_deflections[i].size(); ++j)
            grid_element->external_deflections[i][j] = {0.0, 0.0};
      }
    }
    
    util::vector2 get_deflection(const util::vector2& position) const
    {
      util::vector2 result = {0.0, 0.0};
         
      grid_type::index_type bucket = _grid.get_index(position);
      
      //Check if we are a leaf in the tree structure
      if(_grid[bucket].sub_grid == nullptr)
        result = exact_deflection(bucket, position);
      else       
        result = _grid[bucket].sub_grid->get_deflection(position);
      
      
      util::vector2 interpolation_contribution = interpolate(bucket, position);
      util::add(result, interpolation_contribution);
      
      return result;
    }
    
    void add_star(const star& s)
    {
      grid_type::index_type bucket = _grid.get_index(s.get_position());
      grid_entry& entry = _grid[bucket];
      
      entry.interior_stars.push_back(s);
      
      // Adapt stored exterior deflections 
      grid_type::scalar_array_type bucket_size = _grid.get_bucket_size();
      
      util::matrix_nxn<util::vector2, 2> bucket_corner_positions;

      for(std::size_t i = 0; i < _grid.get_total_num_buckets(0); ++i)
      {
        for(std::size_t j = 0; j < _grid.get_total_num_buckets(1); ++j)
        {
          grid_type::index_type idx = {i, j};
          
          if(idx != bucket)
          {
            util::vector2 bucket_min_position = _grid.get_min_position_of_bucket(idx); 

            for(std::size_t matrix_i = 0; matrix_i < bucket_corner_positions.size(); ++ matrix_i)
              for(std::size_t matrix_j = 0; matrix_j < bucket_corner_positions[matrix_i].size(); ++ matrix_j)
              {
                bucket_corner_positions[matrix_i][matrix_j] = bucket_min_position;
                bucket_corner_positions[matrix_i][matrix_j][0] += matrix_i * bucket_size[0];
                bucket_corner_positions[matrix_i][matrix_j][1] += matrix_j * bucket_size[1];
              }
            
            // Add the star's deflection to the corner interpolation values
            for(std::size_t corner_i = 0; corner_i < bucket_corner_positions.size(); ++corner_i)
              for(std::size_t corner_j = 0; corner_j < bucket_corner_positions.size(); ++corner_j)
              {
                util::vector2 deflection;
                s.calculate_deflection_angle(bucket_corner_positions[corner_i][corner_j], deflection);
                util::add(_grid[idx].external_deflections[corner_i][corner_j], deflection);
              }
          }
        }
      }
          
      
      if(entry.interior_stars.size() > _num_allowed_stars)
      {  
        // Attempt to refine grid
        if(_layer_idx != _max_layer_idx)
        {
          // If we have just exceeded the maximum number of stars,
          // all stars need to be redistributed into the refined grid
          if(entry.interior_stars.size() == _num_allowed_stars + 1)
            for(const star& interior_star : entry.interior_stars)
              insert_star_into_bucket(bucket, interior_star);
          else
            // Otherwise we only need to take care of the new star
            insert_star_into_bucket(bucket, s);
        }
      }
    }
  private:
    inline void insert_star_into_bucket(const grid_type::index_type& bucket_for_star,
                                        const star& s)
    {
      grid_type::scalar_array_type bucket_size = _grid.get_bucket_size();
      
      util::vector2 bucket_min_corner = _grid.get_min_position_of_bucket(bucket_for_star);
      util::vector2 bucket_max_corner = bucket_min_corner;
      util::add(bucket_max_corner, bucket_size);

      //Create new subgrid if necessary
      if(!_grid[bucket_for_star].sub_grid)
        _grid[bucket_for_star].sub_grid = layer_ptr_type(new grid_layer(_layer_idx + 1, 
                                                   _max_layer_idx + 1, 
                                                   _num_allowed_stars, 
                                                   bucket_min_corner,
                                                   bucket_max_corner,
                                                   _num_subgrid_cells,
                                                   _num_subgrid_cells));

      _grid[bucket_for_star].sub_grid->add_star(s);
 
    }
    
    /// @return exact deflection of interior stars
    inline util::vector2 exact_deflection(const grid_type::index_type& bucket,
                                          const util::vector2& position) const
    {
      util::vector2 result = {0.0, 0.0};
      
      for(const auto& star : _grid[bucket].interior_stars)
      {
        util::vector2 contribution;
        star.calculate_deflection_angle(position, contribution);
        util::add(result, contribution);
      }
      
      return result;
    }
    
    inline util::vector2 interpolate(const grid_type::index_type bucket_index,
                                     const util::vector2& position) const
    {
      const grid_entry& bucket = _grid[bucket_index];
      
      util::vector2 relative_position = position;
      util::sub(relative_position, _grid.get_min_position_of_bucket(bucket_index));
      
      grid_type::scalar_array_type bucket_size = _grid.get_bucket_size();
      
      for(std::size_t i = 0; i < relative_position.size(); ++i)
        relative_position[i] /= bucket_size[i];
      
      util::scalar x = relative_position[0];
      util::scalar y = relative_position[1];
      
      util::vector2 result = bucket.external_deflections[0][0];
      util::scale(result, (1.-x)*(1.-y));
      
      util::scale_add(result, bucket.external_deflections[1][0], x * (1. - y));
      util::scale_add(result, bucket.external_deflections[0][1], y * (1. - x));
      util::scale_add(result, bucket.external_deflections[1][1], x * y);
      
      return result;
    }
    
    std::array<std::size_t, 2> _num_cells;
    std::array<std::size_t, 2> _num_subgrid_cells;
    
    grid_type _grid;
    
    
    std::size_t _layer_idx;
    std::size_t _max_layer_idx;
    std::size_t _num_allowed_stars;
  };
public:
  
  nested_interpolation_grid() = default;
  
  nested_interpolation_grid(const std::vector<star>& stars,
                            std::size_t max_num_layers,
                            const util::vector2& min_extent,
                            const util::vector2& max_extent,
                            const std::array<std::size_t, 2>& num_toplevel_cells,
                            const std::array<std::size_t, 2>& num_subgrid_cells)
  : _top_layer(0, max_num_layers, 1, min_extent, max_extent, num_toplevel_cells, num_subgrid_cells)
  {
    assert(max_num_layers > 0);
    
    for(const star& s : stars)
      _top_layer.add_star(s);
  }
  
  util::vector2 get_deflection(const util::vector2& pos) const
  {
    return _top_layer.get_deflection(pos);
  }
  
  
  
private:
  
  grid_layer _top_layer;
};  
  
}

#endif	/* NESTED_INTERPOLATION_GRID_HPP */

